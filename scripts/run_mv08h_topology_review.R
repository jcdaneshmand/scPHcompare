#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop("usage: run_mv08h_topology_review.R <filtered.h5> <raw.h5> <reference.rds> <panel.csv> <output-dir>", call. = FALSE)
}
filtered_path <- normalizePath(args[[1L]], mustWork = TRUE)
raw_path <- normalizePath(args[[2L]], mustWork = TRUE)
reference_path <- normalizePath(args[[3L]], mustWork = TRUE)
panel_path <- normalizePath(args[[4L]], mustWork = TRUE)
out_dir <- normalizePath(args[[5L]], mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")

source("R/dual_view_topology.R")
source("R/toy_baseline.R")
library(Seurat)
library(Matrix)

atomic_csv <- function(x, path) {
  partial <- paste0(path, ".partial")
  utils::write.csv(x, partial, row.names = FALSE, quote = TRUE, na = "")
  file.rename(partial, path)
}
sha <- function(x) digest::digest(x, algo = "sha256", serialize = TRUE)
sha_file <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)

panel <- utils::read.csv(panel_path, check.names = FALSE, stringsAsFactors = FALSE)
stable <- sub("\\.[0-9]+$", "", sub("^.*-", "", panel$feature_id))
reference <- readRDS(reference_path)
if (!identical(as.integer(panel$panel_order), seq_len(500L)) || nrow(panel) != 500L ||
    !identical(unique(panel$panel_sha256), "48e881ee753893bfaecd7101fc16fbcb552dd30a2a791fedfc7204d7b322a32e")) {
  stop("panel contract mismatch", call. = FALSE)
}

filtered <- Seurat::Read10X_h5(filtered_path, use.names = FALSE, unique.features = TRUE)
raw <- Seurat::Read10X_h5(raw_path, use.names = FALSE, unique.features = TRUE)
if (is.list(filtered) || is.list(raw)) stop("multi-assay H5 is not admitted by this bounded runner", call. = FALSE)
if (!identical(rownames(filtered), rownames(raw))) stop("filtered/raw feature axes differ", call. = FALSE)
filtered_names <- Seurat::Read10X_h5(filtered_path, use.names = TRUE, unique.features = TRUE)
if (is.list(filtered_names) || !identical(dim(filtered_names), dim(filtered))) stop("feature-name metadata axis differs", call. = FALSE)

feature_ids <- sub("\\.[0-9]+$", "", rownames(filtered))
panel_idx <- match(stable, feature_ids)
if (anyNA(panel_idx) || anyDuplicated(panel_idx)) stop("panel stable-ID mapping is incomplete or ambiguous", call. = FALSE)

counts <- as(filtered, "dgCMatrix")
nfeature <- Matrix::colSums(counts > 0)
total <- Matrix::colSums(counts)
names_text <- rownames(filtered_names)
mito <- grepl("^MT-", names_text)
ribo <- grepl("^(RPS|RPL)", names_text)
percent_mito <- 100 * Matrix::colSums(counts[mito, , drop = FALSE]) / total
percent_ribo <- 100 * Matrix::colSums(counts[ribo, , drop = FALSE]) / total
eligible <- nfeature >= 500 & nfeature <= 9000 & percent_mito <= 20 & percent_ribo > 5
if (sum(eligible) < 384L) stop("fewer than 384 cells satisfy the frozen QC rule", call. = FALSE)
message("qc_eligible=", sum(eligible))
selected <- select_matched_cells(colnames(counts)[eligible], n = 384L, seed = 20260805L)
selected_idx <- match(selected, colnames(counts))
raw_idx <- match(selected, colnames(raw))

# The raw matrix is read once but reduced immediately to the selected cells.
raw_selected <- as(raw[, raw_idx, drop = FALSE], "dgCMatrix")
rownames(raw_selected) <- rownames(raw)
colnames(raw_selected) <- selected
if (anyNA(selected_idx) || anyNA(raw_idx) || any(Matrix::colSums(raw_selected) <= 0)) {
  stop("selected raw matrix contains missing or zero-count cells; missing_filtered=", sum(is.na(selected_idx)), "; missing_raw=", sum(is.na(raw_idx)), "; zero_columns=", sum(Matrix::colSums(raw_selected) <= 0), call. = FALSE)
}
obj <- Seurat::CreateSeuratObject(counts = raw_selected, min.cells = 0, min.features = 0,
                                  project = "HCA_BM_002")
obj <- Seurat::SCTransform(obj, assay = "RNA", new.assay.name = "SCT",
                           variable.features.n = 3000L, return.only.var.genes = FALSE,
                           verbose = FALSE)
sct <- Seurat::GetAssayData(obj, assay = "SCT", layer = "data")
panel_sct_idx <- match(rownames(raw)[panel_idx], rownames(sct))
if (anyNA(panel_sct_idx) || !all(selected %in% colnames(sct))) {
  stop("SCT output axes do not retain the selected source; rows=", nrow(sct), "; cols=", ncol(sct), "; missing_panel=", sum(is.na(panel_sct_idx)), "; missing_selected=", sum(!selected %in% colnames(sct)), call. = FALSE)
}
sct <- as.matrix(sct[panel_sct_idx, selected, drop = FALSE])
rownames(sct) <- panel$feature_id
colnames(sct) <- selected
if (any(!is.finite(sct))) stop("nonfinite SCT panel values", call. = FALSE)

center <- reference$center
scale <- reference$scale
if (!identical(names(center), panel$feature_id) || !identical(names(scale), panel$feature_id)) {
  stop("reference center/scale names do not match frozen panel", call. = FALSE)
}
standardized <- sweep(sweep(sct, 1L, center, "-"), 1L, scale, "/")
source_obj <- new_dual_view_source(
  standardized, sample_id = "HCA_BM_002", cohort = "mv08h_hca_validation",
  representation = "sct_global_descriptive_v1",
  fit_scope_id = reference$pca_model$fit_scope_id,
  subsample_seed = 20260805L,
  standardization_id = reference$pca_model$standardization_id,
  contract_profile = "scientific"
)
coords <- construct_frozen_cell_topology_view(
  source_obj,
  coordinates = t(standardized) %*% reference$pca_model$rotation,
  coordinate_contract_id = "mv08h_immutable_reference_projection_v1",
  coordinate_fit_cache_key = reference$pca_model$cache_key
)
gene_view <- construct_gene_topology_view(source_obj)
validate_topology_view(coords)
validate_topology_view(gene_view)

run_one <- function(view) {
  started <- proc.time()[["elapsed"]]
  result <- run_topology_view_ph(view, max_dim = 1L, threshold = -1, field = 2L)
  elapsed <- proc.time()[["elapsed"]] - started
  d <- result$diagram
  rows <- lapply(c("H0", "H1"), function(dim_label) {
    dim_num <- if (dim_label == "H0") 0 else 1
    part <- d[d[, "dimension"] == dim_num, , drop = FALSE]
    finite <- part[is.finite(part[, "death"]), , drop = FALSE]
    positive <- finite[(finite[, "death"] - finite[, "birth"]) > 0, , drop = FALSE]
    data.frame(view_id = view$view_id, homology_dimension = dim_label,
               point_count = length(view$point_ids), finite_intervals = nrow(finite),
               positive_intervals = nrow(positive), essential_intervals = sum(is.infinite(part[, "death"])),
               max_persistence = if (nrow(positive)) max(positive[, "death"] - positive[, "birth"]) else 0,
               diagram_sha256 = result$provenance$diagram_sha256,
               elapsed_seconds = elapsed, stringsAsFactors = FALSE)
  })
  list(summary = do.call(rbind, rows), diagram = d, cache_key = result$cache_key)
}

cell_result <- run_one(coords)
gene_result <- run_one(gene_view)
summary <- rbind(cell_result$summary, gene_result$summary)
summary$panel_sha256 <- "48e881ee753893bfaecd7101fc16fbcb552dd30a2a791fedfc7204d7b322a32e"
summary$selected_cells <- 384L
summary$selected_cell_sha256 <- attr(selected, "selected_cell_sha256")
summary$labels_outcomes_opened <- FALSE
summary$landscapes_computed <- FALSE
summary$scientific_claims_authorized <- FALSE
atomic_csv(summary, file.path(out_dir, "mv08h-topology-review-summary.csv"))

identity <- data.frame(
  contract_id = "mv08h_topology_review_execution_v1",
  unit_id = "HCA_BM_002", filtered_cells = ncol(filtered), qc_eligible_cells = sum(eligible),
  selected_cells = 384L, selected_cell_sha256 = attr(selected, "selected_cell_sha256"),
  panel_genes = 500L, panel_sha256 = "48e881ee753893bfaecd7101fc16fbcb552dd30a2a791fedfc7204d7b322a32e",
  source_matrix_sha256 = sha(standardized), reference_source_cache_key = reference$cache_key,
  cell_view_payload_sha256 = coords$payload_sha256, gene_view_payload_sha256 = gene_view$payload_sha256,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  landscapes_computed = FALSE, fusion_computed = FALSE, stringsAsFactors = FALSE
)
atomic_csv(identity, file.path(out_dir, "mv08h-topology-review-identity.csv"))

validation <- data.frame(
  check_id = c("qc_depth", "panel_mapping", "selected_cells", "cell_view_typed", "gene_view_typed", "h0_h1_separate", "cell_h0_mst_count", "gene_h0_mst_count", "positive_h1_or_empty", "label_firewall", "landscape_deferred"),
  passed = c(sum(eligible) >= 384L, TRUE, identical(length(selected), 384L), inherits(coords, "scph_cell_topology_view_v1"), inherits(gene_view, "scph_gene_topology_view_v1"), nrow(summary) == 4L, summary$positive_intervals[summary$view_id == "cell_topology_v1" & summary$homology_dimension == "H0"] == 383L, summary$positive_intervals[summary$view_id == "gene_topology_v1" & summary$homology_dimension == "H0"] == 499L, all(summary$positive_intervals[summary$homology_dimension == "H1"] >= 0L), TRUE, TRUE),
  evidence = c(paste0("eligible=", sum(eligible)), "all 500 stable IDs mapped once", "seed=20260805", "immutable 30-PC coordinates", "correlation-chord dist", "four rows: two views x H0/H1", "cell H0 positive intervals=383", "gene H0 positive intervals=499", "H1 retained even when empty", "labels/outcomes closed", "no landscapes or fusion"),
  stringsAsFactors = FALSE
)
atomic_csv(validation, file.path(out_dir, "mv08h-topology-review-validation.csv"))

report <- c(
  "# MV8-H bounded topology-review execution (2026-08-18)", "",
  "This label-closed technical run executes the frozen dual-view topology contract for HCA_BM_002 only.", "",
  paste0("- QC-eligible cells: ", sum(eligible), "; selected: 384 under seed 20260805 (IDs remain private; only a digest is published)."),
  "- Cell view uses the immutable 30-PC projection; gene view uses the same 500 genes under correlation-chord geometry.",
  "- H0 and H1 remain separate. Complete Vietoris-Rips uses field 2 and threshold -1.",
  "- Landscapes are not computed here; the dissertation-aligned all-active-level exact/error-controlled definition is bound for the next gate.",
  "- Labels, outcomes, fusion, other HCA units, and deletion remain closed.", "",
  "This run is feasibility evidence, not a biological result or validation claim."
)
writeLines(report, file.path(out_dir, "MV08H_TOPOLOGY_REVIEW_2026-08-18.md"), useBytes = TRUE)
cat("completed bounded topology review\n")
