#!/usr/bin/env Rscript

# Label-closed exact-500 transform comparison for the recovered HCA raw-read
# sentinel. This script deliberately stops before persistence or landscapes.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("usage: run_mv08k_exact500_transform_contract.R <filtered.h5> <reference.rds> <panel.csv> <output-dir>", call. = FALSE)
}

filtered_path <- normalizePath(args[[1L]], mustWork = TRUE)
reference_path <- normalizePath(args[[2L]], mustWork = TRUE)
panel_path <- normalizePath(args[[3L]], mustWork = TRUE)
output_dir <- normalizePath(args[[4L]], mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
source("R/toy_baseline.R")
source("R/dual_view_topology.R")
library(Seurat)
library(Matrix)

atomic_csv <- function(x, path) {
  partial <- paste0(path, ".partial")
  utils::write.csv(x, partial, row.names = FALSE, quote = TRUE, na = "")
  if (!file.rename(partial, path)) stop("failed to publish ", basename(path), call. = FALSE)
}

sha_file <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
sha_object <- function(x) digest::digest(x, algo = "sha256", serialize = TRUE)

panel <- utils::read.csv(panel_path, check.names = FALSE, stringsAsFactors = FALSE)
if (nrow(panel) != 500L || !all(c("panel_order", "feature_id", "panel_sha256") %in% names(panel))) {
  stop("expected the ordered exact-500 panel", call. = FALSE)
}
panel <- panel[order(panel$panel_order), , drop = FALSE]
panel_sha <- unique(tolower(panel$panel_sha256))
if (length(panel_sha) != 1L || !identical(as.integer(panel$panel_order), seq_len(500L))) {
  stop("exact-500 panel order or SHA-256 is invalid", call. = FALSE)
}
stable_ids <- sub("\\.[0-9]+$", "", sub("^.*-", "", panel$feature_id))

reference <- readRDS(reference_path)
if (!identical(length(reference$panel$feature_id), 500L) ||
    !identical(tolower(reference$identity$panel_sha256), panel_sha) ||
    !identical(reference$panel$feature_id, panel$feature_id) ||
    nrow(reference$pca_model$rotation) != 500L || ncol(reference$pca_model$rotation) != 30L ||
    !identical(rownames(reference$pca_model$rotation), panel$feature_id)) {
  stop("frozen exact-500 reference contract is incompatible", call. = FALSE)
}

filtered <- Seurat::Read10X_h5(filtered_path, use.names = FALSE, unique.features = TRUE)
filtered_names <- Seurat::Read10X_h5(filtered_path, use.names = TRUE, unique.features = TRUE)
if (is.list(filtered) || is.list(filtered_names) || !identical(dim(filtered), dim(filtered_names))) {
  stop("input must be one feature-by-cell Gene Expression matrix", call. = FALSE)
}
counts <- as(filtered, "dgCMatrix")
feature_ids <- sub("\\.[0-9]+$", "", rownames(counts))
panel_index <- match(stable_ids, feature_ids)
if (anyNA(panel_index) || anyDuplicated(panel_index)) {
  stop("exact-500 stable-ID mapping is incomplete or ambiguous", call. = FALSE)
}

nfeature <- Matrix::colSums(counts > 0)
total <- Matrix::colSums(counts)
names_text <- rownames(filtered_names)
mito <- grepl("^MT-", names_text)
ribo <- grepl("^(RPS|RPL)", names_text)
percent_mito <- 100 * Matrix::colSums(counts[mito, , drop = FALSE]) / total
percent_ribo <- 100 * Matrix::colSums(counts[ribo, , drop = FALSE]) / total
eligible <- nfeature >= 500 & nfeature <= 9000 & percent_mito <= 20 & percent_ribo > 5 & total > 0
if (sum(eligible) < 384L) stop("fewer than 384 cells satisfy the frozen QC rule", call. = FALSE)
selected <- select_matched_cells(colnames(counts)[eligible], n = 384L, seed = 20260805L)
selected_index <- match(selected, colnames(counts))
selected_counts <- as(counts[, selected_index, drop = FALSE], "dgCMatrix")
rownames(selected_counts) <- rownames(counts)
colnames(selected_counts) <- selected
if (anyNA(selected_index) || any(Matrix::colSums(selected_counts) <= 0)) {
  stop("selected cells are missing or zero-count", call. = FALSE)
}

run_transform <- function(config_id, min_cells) {
  started <- proc.time()[["elapsed"]]
  object <- Seurat::CreateSeuratObject(counts = selected_counts, min.cells = 0, min.features = 0,
                                        project = "HCA_BM_002")
  object <- Seurat::SCTransform(
    object, assay = "RNA", new.assay.name = "SCT", variable.features.n = 3000L,
    return.only.var.genes = FALSE, verbose = FALSE, min_cells = min_cells
  )
  sct <- Seurat::GetAssayData(object, assay = "SCT", layer = "data")
  panel_source_ids <- rownames(selected_counts)[panel_index]
  sct_index <- match(panel_source_ids, rownames(sct))
  retained <- !is.na(sct_index)
  elapsed <- proc.time()[["elapsed"]] - started
  result <- list(
    config_id = config_id, min_cells = min_cells, elapsed_seconds = elapsed,
    sct_feature_count = nrow(sct), panel_retained = sum(retained),
    missing_panel_count = sum(!retained), exact500_retained = all(retained),
    standardized_finite = FALSE, pca_compatible = FALSE, cell_view_valid = FALSE,
    gene_view_valid = FALSE, source_matrix_sha256 = NA_character_,
    cell_view_payload_sha256 = NA_character_, gene_view_payload_sha256 = NA_character_,
    zero_variance_gene_count = NA_integer_, minimum_gene_sd = NA_real_,
    view_error = NA_character_
  )
  if (!all(retained)) return(result)
  sct_panel <- as.matrix(sct[sct_index, selected, drop = FALSE])
  rownames(sct_panel) <- panel$feature_id
  colnames(sct_panel) <- selected
  standardized <- sweep(sweep(sct_panel, 1L, reference$center, "-"), 1L, reference$scale, "/")
  result$standardized_finite <- all(is.finite(standardized))
  if (!result$standardized_finite) return(result)
  coordinates <- t(standardized) %*% reference$pca_model$rotation
  result$pca_compatible <- identical(dim(coordinates), c(384L, 30L)) && all(is.finite(coordinates))
  if (!result$pca_compatible) return(result)
  # Retain a digest even when the downstream geometry gate declines admission;
  # it is the deterministic-repeat evidence for the immutable source matrix.
  result$source_matrix_sha256 <- sha_object(standardized)
  gene_sd <- apply(standardized, 1L, stats::sd)
  result$zero_variance_gene_count <- sum(!is.finite(gene_sd) | gene_sd <= sqrt(.Machine$double.eps))
  result$minimum_gene_sd <- suppressWarnings(min(gene_sd, na.rm = TRUE))
  # This is a scientific admission gate: gene correlation-chord geometry is
  # undefined for an effectively constant feature. Record the result rather
  # than weakening new_dual_view_source() or silently dropping the feature.
  if (result$zero_variance_gene_count > 0L) return(result)
  source_obj <- new_dual_view_source(
    standardized, sample_id = "HCA_BM_002", cohort = "mv08k_hca_exact500",
    representation = "sct_global_descriptive_v1", fit_scope_id = reference$pca_model$fit_scope_id,
    subsample_seed = 20260805L, standardization_id = reference$pca_model$standardization_id,
    contract_profile = "scientific", expected_genes = 500L
  )
  views <- tryCatch({
    cell_view <- construct_frozen_cell_topology_view(
      source_obj, coordinates = coordinates,
      coordinate_contract_id = "mv08k_immutable_exact500_reference_projection_v1",
      coordinate_fit_cache_key = reference$pca_model$cache_key
    )
    gene_view <- construct_gene_topology_view(source_obj)
    validate_topology_view(cell_view)
    validate_topology_view(gene_view)
    list(cell = cell_view, gene = gene_view)
  }, error = function(error) {
    result$view_error <<- conditionMessage(error)
    NULL
  })
  if (is.null(views)) return(result)
  cell_view <- views$cell
  gene_view <- views$gene
  result$cell_view_valid <- inherits(cell_view, "scph_cell_topology_view_v1")
  result$gene_view_valid <- inherits(gene_view, "scph_gene_topology_view_v1")
  result$cell_view_payload_sha256 <- cell_view$payload_sha256
  result$gene_view_payload_sha256 <- gene_view$payload_sha256
  result
}

results <- lapply(
  list(list("sct_default_min_cells5", 5L), list("sct_lowcount_min_cells1", 1L)),
  function(config) run_transform(config[[1L]], config[[2L]])
)
summary <- do.call(rbind, lapply(results, function(x) as.data.frame(x, stringsAsFactors = FALSE)))
summary$selected_cells <- 384L
summary$qc_eligible_cells <- sum(eligible)
summary$selected_cell_sha256 <- attr(selected, "selected_cell_sha256")
summary$panel_sha256 <- panel_sha
summary$reference_pca_cache_key <- reference$pca_model$cache_key
summary$labels_outcomes_opened <- FALSE
summary$persistence_computed <- FALSE
summary$landscapes_computed <- FALSE
summary$fusion_computed <- FALSE
atomic_csv(summary, file.path(output_dir, "mv08k-exact500-transform-summary.csv"))

identity <- data.frame(
  contract_id = "mv08k_exact500_transform_contract_v1",
  unit_id = "HCA_BM_002", filtered_h5_sha256 = sha_file(filtered_path),
  exact500_panel_sha256 = panel_sha, reference_record_sha256 = sha_file(reference_path),
  reference_pca_cache_key = reference$pca_model$cache_key,
  selected_cells = 384L, selected_cell_sha256 = attr(selected, "selected_cell_sha256"),
  labels_outcomes_opened = FALSE, persistence_computed = FALSE,
  landscapes_computed = FALSE, fusion_computed = FALSE,
  stringsAsFactors = FALSE
)
atomic_csv(identity, file.path(output_dir, "mv08k-exact500-transform-identity.csv"))

default_row <- summary[summary$config_id == "sct_default_min_cells5", , drop = FALSE]
low_row <- summary[summary$config_id == "sct_lowcount_min_cells1", , drop = FALSE]
validation <- data.frame(
  check_id = c("input_panel_mapping", "frozen_reference_binding", "frozen_qc_selection",
               "standard_sct_insufficient_exact500", "lowcount_sct_retains_exact500",
               "lowcount_standardized_finite", "lowcount_frozen_pca_compatible",
               "lowcount_no_effective_zero_variance_genes", "lowcount_separate_views_valid",
               "label_firewall", "persistence_deferred"),
  passed = c(TRUE, TRUE, identical(length(selected), 384L),
             default_row$panel_retained[[1L]] < 500L,
             low_row$exact500_retained[[1L]], low_row$standardized_finite[[1L]],
             low_row$pca_compatible[[1L]], low_row$zero_variance_gene_count[[1L]] == 0L,
             low_row$cell_view_valid[[1L]] && low_row$gene_view_valid[[1L]],
             TRUE, TRUE),
  evidence = c("all 500 stable IDs map exactly once", "panel, center/scale, and 30-PC model are frozen",
               paste0("eligible=", sum(eligible), "; selected=384; seed=20260805"),
               paste0("retained=", default_row$panel_retained[[1L]], "/500"),
               paste0("retained=", low_row$panel_retained[[1L]], "/500"),
               paste0("finite=", low_row$standardized_finite[[1L]]),
               paste0("coordinates=384x30; finite=", low_row$pca_compatible[[1L]]),
               paste0("count=", low_row$zero_variance_gene_count[[1L]], "; minimum_sd=", low_row$minimum_gene_sd[[1L]]),
               "cell Euclidean and gene correlation-chord views typed without PH",
               "labels/outcomes remain closed", "persistence and landscapes remain deferred"),
  stringsAsFactors = FALSE
)
atomic_csv(validation, file.path(output_dir, "mv08k-exact500-transform-validation.csv"))

report <- c(
  paste0("# MV8-K exact-500 transform-contract comparison (", format(Sys.Date()), ")"), "",
  "This label-closed comparison determines whether the recovered exact-500 raw-read panel can retain its full axis under the frozen HCA topology transform.", "",
  "- Both configurations use the same frozen 384-cell selection, exact-500 panel, center/scale, and immutable 30-PC reference model.",
  "- Standard SCT uses its default five-cell inclusion rule; low-count SCT changes only the inclusion rule to one cell.",
  "- No persistence diagrams, landscapes, clustering, labels, outcomes, fusion, or manuscript claims are computed.",
  "- A scale-up decision requires every validation check to pass and a separate deterministic repeat.",
  "- Retention alone is not enough: every standardized feature must have effective nonzero variance before gene correlation-chord topology is admissible."
)
writeLines(report, file.path(output_dir, paste0("MV08K_EXACT500_TRANSFORM_CONTRACT_", format(Sys.Date()), ".md")), useBytes = TRUE)

passed <- all(validation$passed)
cat("MV8-K exact500 transform contract checks=", sum(validation$passed), "/", nrow(validation), "\n", sep = "")
quit(status = if (passed) 0L else 2L)
