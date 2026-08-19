#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop("usage: run_mv08i_hca_common475_validation.R <file-manifest.csv> <h5-cache-dir> <reference.rds> <panel.csv> <output-dir>", call. = FALSE)
}

manifest_path <- normalizePath(args[[1L]], mustWork = TRUE)
cache_dir <- normalizePath(args[[2L]], mustWork = TRUE)
reference_path <- normalizePath(args[[3L]], mustWork = TRUE)
panel_path <- normalizePath(args[[4L]], mustWork = TRUE)
out_dir <- normalizePath(args[[5L]], mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "private"), recursive = TRUE, showWarnings = FALSE)

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
source("R/dual_view_topology.R")
source("R/toy_baseline.R")

atomic_csv <- function(x, path) {
  partial <- paste0(path, ".partial")
  utils::write.csv(x, partial, row.names = FALSE, quote = TRUE, na = "")
  file.rename(partial, path)
}
sha_file <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)

manifest <- utils::read.csv(manifest_path, check.names = FALSE, stringsAsFactors = FALSE)
required <- c("file_order", "unit_id", "file_name", "expected_bytes", "expected_sha256")
if (!all(required %in% names(manifest)) || nrow(manifest) != 8L) {
  stop("MV8-I requires the frozen eight-row HCA file manifest", call. = FALSE)
}
manifest <- manifest[order(manifest$file_order), , drop = FALSE]
if (!identical(as.integer(manifest$file_order), seq_len(8L)) ||
    anyDuplicated(manifest$unit_id) || anyDuplicated(manifest$file_name)) {
  stop("HCA file manifest order or unit identity is invalid", call. = FALSE)
}

panel <- utils::read.csv(panel_path, check.names = FALSE, stringsAsFactors = FALSE)
if (nrow(panel) != 475L || !all(c("panel_order_475", "feature_id", "common475_axis_sha256") %in% names(panel))) {
  stop("common475 panel contract is invalid", call. = FALSE)
}
panel <- panel[order(panel$panel_order_475), , drop = FALSE]
panel_digest <- unique(panel$common475_axis_sha256)
if (length(panel_digest) != 1L || !identical(as.integer(panel$panel_order_475), seq_len(475L))) {
  stop("common475 panel order or digest is invalid", call. = FALSE)
}
stable <- sub("\\.[0-9]+$", "", sub("^.*-", "", panel$feature_id))
reference <- readRDS(reference_path)
if (!identical(length(reference$panel$feature_id), 475L) ||
    !identical(reference$identity$common475_axis_sha256, panel_digest)) {
  stop("common475 reference identity does not match the panel", call. = FALSE)
}

run_one <- function(row) {
  unit_id <- row$unit_id[[1L]]
  h5_path <- file.path(cache_dir, unit_id, row$file_name[[1L]])
  if (!file.exists(h5_path)) stop("missing HCA H5 input for ", unit_id, call. = FALSE)
  observed_bytes <- file.info(h5_path)$size
  observed_sha <- sha_file(h5_path)
  if (!identical(as.numeric(observed_bytes), as.numeric(row$expected_bytes[[1L]])) ||
      !identical(tolower(observed_sha), tolower(row$expected_sha256[[1L]]))) {
    stop("HCA H5 input identity mismatch for ", unit_id, call. = FALSE)
  }

  message("MV8-I loading ", unit_id, flush = TRUE)
  raw <- Seurat::Read10X_h5(h5_path, use.names = FALSE, unique.features = TRUE)
  raw_names <- Seurat::Read10X_h5(h5_path, use.names = TRUE, unique.features = TRUE)
  if (is.list(raw) || is.list(raw_names) || !identical(dim(raw), dim(raw_names))) {
    stop("HCA H5 is not a single Gene Expression matrix for ", unit_id, call. = FALSE)
  }
  counts <- as(raw, "dgCMatrix")
  feature_names <- rownames(raw_names)
  nfeature <- Matrix::colSums(counts > 0)
  total <- Matrix::colSums(counts)
  mito <- grepl("^MT-", feature_names)
  ribo <- grepl("^(RPS|RPL)", feature_names)
  percent_mito <- 100 * Matrix::colSums(counts[mito, , drop = FALSE]) / total
  percent_ribo <- 100 * Matrix::colSums(counts[ribo, , drop = FALSE]) / total
  eligible <- nfeature >= 500 & nfeature <= 9000 & percent_mito <= 20 & percent_ribo > 5 & total > 0
  if (sum(eligible) < 384L) {
    stop("fewer than 384 cells satisfy the frozen QC rule for ", unit_id, call. = FALSE)
  }
  selected <- select_matched_cells(colnames(counts)[eligible], n = 384L, seed = 20260805L)
  selected_idx <- match(selected, colnames(counts))
  raw_selected <- as(counts[, selected_idx, drop = FALSE], "dgCMatrix")
  rownames(raw_selected) <- rownames(raw)
  colnames(raw_selected) <- selected
  rm(raw_names, counts, raw)
  gc(verbose = FALSE)

  feature_ids <- sub("\\.[0-9]+$", "", rownames(raw_selected))
  panel_idx <- match(stable, feature_ids)
  if (anyNA(panel_idx) || anyDuplicated(panel_idx)) {
    stop("common475 stable-ID mapping is incomplete for ", unit_id, call. = FALSE)
  }
  obj <- Seurat::CreateSeuratObject(counts = raw_selected, min.cells = 0, min.features = 0, project = unit_id)
  obj <- Seurat::SCTransform(obj, assay = "RNA", new.assay.name = "SCT", variable.features.n = 3000L,
                              return.only.var.genes = FALSE, verbose = FALSE)
  sct <- Seurat::GetAssayData(obj, assay = "SCT", layer = "data")
  panel_sct_idx <- match(rownames(raw_selected)[panel_idx], rownames(sct))
  if (anyNA(panel_sct_idx) || !all(selected %in% colnames(sct))) {
    stop("SCT output does not retain the common475 panel for ", unit_id, call. = FALSE)
  }
  sct <- as.matrix(sct[panel_sct_idx, selected, drop = FALSE])
  rownames(sct) <- panel$feature_id
  colnames(sct) <- selected
  center <- reference$center
  scale <- reference$scale
  if (!identical(names(center), panel$feature_id) || !identical(names(scale), panel$feature_id)) {
    stop("common475 reference center/scale names do not match for ", unit_id, call. = FALSE)
  }
  standardized <- sweep(sweep(sct, 1L, center, "-"), 1L, scale, "/")
  source_obj <- new_dual_view_source(standardized, sample_id = unit_id, cohort = "mv08i_hca_external",
                                     representation = "sct_global_descriptive_v1", fit_scope_id = reference$pca_model$fit_scope_id,
                                     subsample_seed = 20260805L, standardization_id = reference$pca_model$standardization_id,
                                     contract_profile = "scientific_common475", expected_genes = 475L)
  coords <- construct_frozen_cell_topology_view(source_obj,
    coordinates = t(standardized) %*% reference$pca_model$rotation,
    coordinate_contract_id = "mv08i_immutable_common475_reference_projection_v1",
    coordinate_fit_cache_key = reference$pca_model$cache_key)
  gene_view <- construct_gene_topology_view(source_obj)
  validate_topology_view(coords)
  validate_topology_view(gene_view)
  run_one_view <- function(view) {
    started <- proc.time()[["elapsed"]]
    result <- run_topology_view_ph(view, max_dim = 1L, threshold = -1, field = 2L)
    elapsed <- proc.time()[["elapsed"]] - started
    d <- result$diagram
    rows <- lapply(c("H0", "H1"), function(dim_label) {
      dim_num <- if (dim_label == "H0") 0 else 1
      part <- d[d[, "dimension"] == dim_num, , drop = FALSE]
      finite <- part[is.finite(part[, "death"]), , drop = FALSE]
      positive <- finite[(finite[, "death"] - finite[, "birth"]) > 0, , drop = FALSE]
      data.frame(unit_id = unit_id, view_id = view$view_id, homology_dimension = dim_label,
        point_count = length(view$point_ids), finite_intervals = nrow(finite),
        positive_intervals = nrow(positive), essential_intervals = sum(is.infinite(part[, "death"])),
        max_persistence = if (nrow(positive)) max(positive[, "death"] - positive[, "birth"]) else 0,
        diagram_sha256 = result$provenance$diagram_sha256, elapsed_seconds = elapsed,
        panel_sha256 = panel_digest, selected_cells = 384L,
        selected_cell_sha256 = attr(selected, "selected_cell_sha256"),
        labels_outcomes_opened = FALSE, landscapes_computed = FALSE, fusion_computed = FALSE,
        stringsAsFactors = FALSE)
    })
    list(summary = do.call(rbind, rows), diagram = d)
  }
  cell_result <- run_one_view(coords)
  gene_result <- run_one_view(gene_view)
  saveRDS(list(cell_topology_v1 = cell_result$diagram, gene_topology_v1 = gene_result$diagram,
               unit_id = unit_id, panel_sha256 = panel_digest, selected_cell_sha256 = attr(selected, "selected_cell_sha256")),
          file.path(out_dir, "private", paste0(unit_id, "-topology-diagrams.rds")))
  identity <- data.frame(unit_id = unit_id, file_sha256 = observed_sha, file_bytes = observed_bytes,
    qc_eligible_cells = sum(eligible), selected_cells = 384L, selected_cell_sha256 = attr(selected, "selected_cell_sha256"),
    panel_genes = 475L, panel_sha256 = panel_digest, source_matrix_sha256 = digest::digest(standardized, algo = "sha256"),
    cell_view_payload_sha256 = coords$payload_sha256, gene_view_payload_sha256 = gene_view$payload_sha256,
    outcome_label_state = "closed", labels_outcomes_opened = FALSE, landscapes_computed = FALSE, fusion_computed = FALSE,
    stringsAsFactors = FALSE)
  list(summary = rbind(cell_result$summary, gene_result$summary), identity = identity)
}

results <- lapply(seq_len(nrow(manifest)), function(i) run_one(manifest[i, , drop = FALSE]))
summary <- do.call(rbind, lapply(results, `[[`, "summary"))
identity <- do.call(rbind, lapply(results, `[[`, "identity"))
validation <- data.frame(
  check_id = c("unit_count", "panel_identity", "selected_cells", "separate_views", "separate_h0_h1", "label_firewall", "fusion_closed", "landscapes_deferred"),
  passed = c(nrow(identity) == 8L, all(summary$panel_sha256 == panel_digest), all(summary$selected_cells == 384L),
             identical(sort(unique(summary$view_id)), sort(c("cell_topology_v1", "gene_topology_v1"))),
             nrow(summary) == 32L, all(!summary$labels_outcomes_opened), all(!summary$fusion_computed), all(!summary$landscapes_computed)),
  evidence = c("8 verified HCA units", "ordered common475 stable-ID panel", "384 cells per unit", "cell and gene views", "8 units x 2 views x H0/H1", "labels/outcomes closed", "fusion closed", "landscapes deferred to separate exact pass"),
  stringsAsFactors = FALSE)
atomic_csv(summary, file.path(out_dir, "mv08i-topology-summary.csv"))
atomic_csv(identity, file.path(out_dir, "mv08i-input-identity.csv"))
atomic_csv(validation, file.path(out_dir, "mv08i-topology-validation.csv"))
writeLines(c(
  paste0("# MV8-I common475 HCA topology validation (", format(Sys.Date()), ")"), "",
  "This label-closed run uses the eight verified adult bone-marrow HCA units as an external common475 feasibility/replication gate.", "",
  "- Each unit is independently QC-screened and deterministically reduced to 384 cells.",
  "- Cell and gene topology views are separate; H0 and H1 remain separate.",
  "- The frozen 475-gene reference transform is reused without HCA-derived refitting.",
  "- Landscapes are computed in a separate exact all-active-level pass.",
  "- Labels, outcomes, fusion, exact-500 recovery, and deletion remain closed.", "",
  "This is technical external-validation evidence, not a biological claim."),
  file.path(out_dir, paste0("MV08I_TOPOLOGY_VALIDATION_", format(Sys.Date()), ".md")))
cat("completed MV8-I topology validation for ", nrow(identity), " units\n", sep = "")
