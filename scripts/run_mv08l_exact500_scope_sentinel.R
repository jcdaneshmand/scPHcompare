#!/usr/bin/env Rscript

# Label-closed MV8-L feasibility sentinel. It tests a gene-side all-QC-cell
# source only and deliberately stops before typed topology, persistence, or
# landscapes.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop(
    "usage: run_mv08l_exact500_scope_sentinel.R <filtered.h5> <reference.rds> <exact500.csv> <common475.csv> <output-dir>",
    call. = FALSE
  )
}

filtered_path <- normalizePath(args[[1L]], mustWork = TRUE)
reference_path <- normalizePath(args[[2L]], mustWork = TRUE)
exact_path <- normalizePath(args[[3L]], mustWork = TRUE)
common_path <- normalizePath(args[[4L]], mustWork = TRUE)
output_dir <- normalizePath(args[[5L]], mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
library(Seurat)
library(Matrix)

atomic_csv <- function(x, path) {
  partial <- paste0(path, ".partial")
  utils::write.csv(x, partial, row.names = FALSE, quote = TRUE, na = "")
  if (!file.rename(partial, path)) stop("failed to publish ", basename(path), call. = FALSE)
}
sha_file <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
sha_object <- function(x) digest::digest(x, algo = "sha256", serialize = TRUE)

exact <- utils::read.csv(exact_path, check.names = FALSE, stringsAsFactors = FALSE)
common <- utils::read.csv(common_path, check.names = FALSE, stringsAsFactors = FALSE)
if (nrow(exact) != 500L || !all(c("panel_order", "feature_id", "panel_sha256") %in% names(exact))) {
  stop("expected ordered exact-500 panel", call. = FALSE)
}
if (nrow(common) != 475L || !all(c("panel_order_475", "feature_id", "common475_axis_sha256") %in% names(common))) {
  stop("expected ordered common-475 panel", call. = FALSE)
}
exact <- exact[order(exact$panel_order), , drop = FALSE]
common <- common[order(common$panel_order_475), , drop = FALSE]
exact_sha <- unique(tolower(exact$panel_sha256))
common_sha <- unique(tolower(common$common475_axis_sha256))
if (length(exact_sha) != 1L || length(common_sha) != 1L ||
    !identical(as.integer(exact$panel_order), seq_len(500L)) ||
    !all(common$feature_id %in% exact$feature_id)) {
  stop("panel identities are incompatible", call. = FALSE)
}

reference <- readRDS(reference_path)
if (!identical(reference$panel$feature_id, exact$feature_id) ||
    !identical(tolower(reference$identity$panel_sha256), exact_sha) ||
    !identical(rownames(reference$pca_model$rotation), exact$feature_id) ||
    length(reference$center) != 500L || length(reference$scale) != 500L ||
    any(!is.finite(reference$scale)) || any(reference$scale <= 0)) {
  stop("frozen exact-500 reference is incompatible", call. = FALSE)
}

raw <- Seurat::Read10X_h5(filtered_path, use.names = FALSE, unique.features = TRUE)
raw_names <- Seurat::Read10X_h5(filtered_path, use.names = TRUE, unique.features = TRUE)
if (is.list(raw) || is.list(raw_names) || !identical(dim(raw), dim(raw_names))) {
  stop("input must be a single Gene Expression feature-by-cell matrix", call. = FALSE)
}
counts <- as(raw, "dgCMatrix")
stable_ids <- sub("\\.[0-9]+$", "", sub("^.*-", "", exact$feature_id))
feature_ids <- sub("\\.[0-9]+$", "", rownames(counts))
panel_index <- match(stable_ids, feature_ids)
if (anyNA(panel_index) || anyDuplicated(panel_index)) {
  stop("all ordered exact-500 stable IDs must map exactly once", call. = FALSE)
}
common_index <- match(common$feature_id, exact$feature_id)
if (anyNA(common_index) || anyDuplicated(common_index)) {
  stop("common475 must be an ordered exact500 subset", call. = FALSE)
}

nfeature <- Matrix::colSums(counts > 0)
total <- Matrix::colSums(counts)
names_text <- rownames(raw_names)
mito <- grepl("^MT-", names_text)
ribo <- grepl("^(RPS|RPL)", names_text)
percent_mito <- 100 * Matrix::colSums(counts[mito, , drop = FALSE]) / total
percent_ribo <- 100 * Matrix::colSums(counts[ribo, , drop = FALSE]) / total
eligible <- nfeature >= 500 & nfeature <= 9000 & percent_mito <= 20 & percent_ribo > 5 & total > 0
eligible_counts <- as(counts[, eligible, drop = FALSE], "dgCMatrix")
if (ncol(eligible_counts) < 384L) stop("fewer than 384 cells pass frozen QC", call. = FALSE)
eligible_counts <- eligible_counts[, order(colnames(eligible_counts), method = "radix"), drop = FALSE]
panel_source_ids <- rownames(counts)[panel_index]
raw_panel <- eligible_counts[panel_index, , drop = FALSE]
raw_detected <- Matrix::rowSums(raw_panel > 0)

distance_summary <- function(standardized) {
  gene_sd <- apply(standardized, 1L, stats::sd)
  if (any(!is.finite(gene_sd)) || any(gene_sd <= sqrt(.Machine$double.eps))) {
    return(list(valid = FALSE, zero_variance = sum(!is.finite(gene_sd) | gene_sd <= sqrt(.Machine$double.eps)),
                minimum_sd = suppressWarnings(min(gene_sd, na.rm = TRUE)), summary = NULL))
  }
  correlation <- stats::cor(t(standardized), method = "pearson")
  distance <- sqrt(pmax(0, 2 * (1 - pmin(1, pmax(-1, correlation)))))
  diag(distance) <- 0
  off_diagonal <- distance[row(distance) != col(distance)]
  common_distance <- distance[common_index, common_index, drop = FALSE]
  list(
    valid = isTRUE(all(is.finite(distance))) && isTRUE(all(is.finite(common_distance))) &&
      isTRUE(isSymmetric(distance)) && isTRUE(isSymmetric(common_distance)) &&
      isTRUE(all(diag(distance) == 0)) && isTRUE(all(off_diagonal >= 0)),
    zero_variance = 0L, minimum_sd = min(gene_sd),
    summary = list(
      exact500_distance_sha256 = sha_object(distance),
      common475_distance_sha256 = sha_object(common_distance),
      exact500_offdiagonal_minimum = min(off_diagonal),
      exact500_offdiagonal_median = stats::median(off_diagonal),
      exact500_offdiagonal_maximum = max(off_diagonal),
      common475_offdiagonal_median = stats::median(common_distance[row(common_distance) != col(common_distance)]),
      common475_submatrix_finite = all(is.finite(common_distance))
    )
  )
}

run_configuration <- function(config_id, min_cells) {
  started <- proc.time()[["elapsed"]]
  object <- Seurat::CreateSeuratObject(
    counts = eligible_counts, min.cells = 0, min.features = 0, project = "HCA_BM_002"
  )
  object <- Seurat::SCTransform(
    object, assay = "RNA", new.assay.name = "SCT", variable.features.n = 3000L,
    return.only.var.genes = FALSE, verbose = FALSE, min_cells = min_cells
  )
  sct <- Seurat::GetAssayData(object, assay = "SCT", layer = "data")
  sct_index <- match(panel_source_ids, rownames(sct))
  retained <- !is.na(sct_index)
  result <- list(
    config_id = config_id, min_cells = min_cells,
    elapsed_seconds = proc.time()[["elapsed"]] - started,
    sct_feature_count = nrow(sct), panel_retained = sum(retained),
    missing_panel_count = sum(!retained), exact500_retained = all(retained),
    standardized_finite = FALSE, zero_variance_gene_count = NA_integer_,
    minimum_gene_sd = NA_real_, correlation_chord_valid = FALSE,
    exact500_distance_sha256 = NA_character_, common475_distance_sha256 = NA_character_,
    exact500_offdiagonal_minimum = NA_real_, exact500_offdiagonal_median = NA_real_,
    exact500_offdiagonal_maximum = NA_real_, common475_offdiagonal_median = NA_real_,
    common475_submatrix_finite = FALSE
  )
  if (!all(retained)) return(result)
  sct_panel <- as.matrix(sct[sct_index, colnames(eligible_counts), drop = FALSE])
  rownames(sct_panel) <- exact$feature_id
  standardized <- sweep(sweep(sct_panel, 1L, reference$center, "-"), 1L, reference$scale, "/")
  result$standardized_finite <- all(is.finite(standardized))
  if (!result$standardized_finite) return(result)
  distance_result <- distance_summary(standardized)
  result$zero_variance_gene_count <- distance_result$zero_variance
  result$minimum_gene_sd <- distance_result$minimum_sd
  result$correlation_chord_valid <- distance_result$valid
  if (!is.null(distance_result$summary)) {
    for (name in names(distance_result$summary)) result[[name]] <- distance_result$summary[[name]]
  }
  result
}

default_result <- run_configuration("sct_default_min_cells5_all_qc", 5L)
results <- list(default_result)
if (!isTRUE(default_result$exact500_retained) || !isTRUE(default_result$correlation_chord_valid)) {
  results[[length(results) + 1L]] <- run_configuration("sct_lowcount_min_cells1_all_qc", 1L)
}
summary <- do.call(rbind, lapply(results, function(x) as.data.frame(x, stringsAsFactors = FALSE)))
summary$eligible_cells <- ncol(eligible_counts)
summary$raw_panel_detection_minimum <- min(raw_detected)
summary$raw_panel_detection_median <- stats::median(raw_detected)
summary$raw_panel_detection_maximum <- max(raw_detected)
summary$exact500_panel_sha256 <- exact_sha
summary$common475_panel_sha256 <- common_sha
summary$reference_record_sha256 <- sha_file(reference_path)
summary$labels_outcomes_opened <- FALSE
summary$persistence_computed <- FALSE
summary$landscapes_computed <- FALSE
summary$fusion_computed <- FALSE
atomic_csv(summary, file.path(output_dir, "mv08l-exact500-scope-summary.csv"))

retained_candidate <- summary[summary$exact500_retained, , drop = FALSE]
finite_candidate <- summary[summary$exact500_retained & summary$standardized_finite, , drop = FALSE]
nonzero_candidate <- summary[summary$exact500_retained & summary$standardized_finite &
                              summary$zero_variance_gene_count == 0L, , drop = FALSE]
distance_candidate <- summary[summary$exact500_retained & summary$standardized_finite &
                               summary$zero_variance_gene_count == 0L & summary$correlation_chord_valid, , drop = FALSE]
validation <- data.frame(
  check_id = c("input_panel_mapping", "frozen_reference_binding", "qc_before_panel", "all_qc_cells_used",
               "raw_exact500_detected", "candidate_retains_exact500", "candidate_finite_standardization",
               "candidate_no_effective_zero_variance", "candidate_correlation_chord_valid",
               "candidate_common475_submatrix_finite", "label_firewall", "persistence_deferred"),
  passed = c(TRUE, TRUE, ncol(eligible_counts) == 4614L, identical(colnames(eligible_counts), sort(colnames(eligible_counts), method = "radix")),
             min(raw_detected) > 0L, nrow(retained_candidate) >= 1L, nrow(finite_candidate) >= 1L,
             nrow(nonzero_candidate) >= 1L, nrow(distance_candidate) >= 1L,
             nrow(distance_candidate) >= 1L, TRUE, TRUE),
  evidence = c("all 500 stable IDs map exactly once", "frozen ordered panel, center/scale, and reference record bind",
               paste0("eligible=", ncol(eligible_counts), "; frozen QC rule evaluated before panel diagnostics"),
               "all eligible cell barcodes sorted deterministically; no panel-aware filter",
               paste0("minimum detected cells across panel=", min(raw_detected)),
               if (nrow(retained_candidate)) paste0(retained_candidate$config_id[[1L]], "; retained=500/500") else "no candidate retained 500/500",
               if (nrow(finite_candidate)) paste0(finite_candidate$config_id[[1L]], "; finite=TRUE") else "no finite full-panel candidate",
               if (nrow(nonzero_candidate)) paste0(nonzero_candidate$config_id[[1L]], "; zero_variance=0") else "no nonzero-variance full-panel candidate",
               if (nrow(distance_candidate)) paste0(distance_candidate$config_id[[1L]], "; finite symmetric correlation-chord") else "no valid correlation-chord candidate",
               if (nrow(distance_candidate)) paste0(distance_candidate$config_id[[1L]], "; common475 submatrix finite") else "no valid common475 submatrix",
               "labels/outcomes remain closed", "persistence and landscapes remain deferred"),
  stringsAsFactors = FALSE
)
atomic_csv(validation, file.path(output_dir, "mv08l-exact500-scope-validation.csv"))

identity <- data.frame(
  contract_id = "mv08l_exact500_view_specific_scope_prefreeze_v1", unit_id = "HCA_BM_002",
  filtered_h5_sha256 = sha_file(filtered_path), exact500_panel_sha256 = exact_sha,
  common475_panel_sha256 = common_sha, reference_record_sha256 = sha_file(reference_path),
  eligible_cells = ncol(eligible_counts), labels_outcomes_opened = FALSE,
  persistence_computed = FALSE, landscapes_computed = FALSE, fusion_computed = FALSE,
  stringsAsFactors = FALSE
)
atomic_csv(identity, file.path(output_dir, "mv08l-exact500-scope-identity.csv"))

cat("MV8-L exact500 scope checks=", sum(validation$passed), "/", nrow(validation), "\n", sep = "")
quit(status = if (all(validation$passed)) 0L else 2L)
