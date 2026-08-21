#!/usr/bin/env Rscript

# Label-closed MV8-M sentinel. Compare aggregate geometry from three declared
# gene representations on one frozen HCA unit; stop before PH or biology.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop(
    "usage: run_mv08m_exact500_representation_sentinel.R <filtered.h5> <reference.rds> <exact500.csv> <common475.csv> <output-dir>",
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
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Matrix))

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
    !identical(as.integer(common$panel_order_475), seq_len(475L)) ||
    anyDuplicated(exact$feature_id) || anyDuplicated(common$feature_id) ||
    !all(common$feature_id %in% exact$feature_id)) {
  stop("panel identities are incompatible", call. = FALSE)
}
common_index <- match(common$feature_id, exact$feature_id)

reference <- readRDS(reference_path)
if (!identical(reference$panel$feature_id, exact$feature_id) ||
    !identical(tolower(reference$identity$panel_sha256), exact_sha) ||
    !identical(rownames(reference$pca_model$rotation), exact$feature_id)) {
  stop("frozen exact-500 reference identity is incompatible", call. = FALSE)
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

# Frozen QC is evaluated before panel expression or representation inspection.
nfeature <- Matrix::colSums(counts > 0)
total <- Matrix::colSums(counts)
names_text <- rownames(raw_names)
mito <- grepl("^MT-", names_text)
ribo <- grepl("^(RPS|RPL)", names_text)
percent_mito <- 100 * Matrix::colSums(counts[mito, , drop = FALSE]) / total
percent_ribo <- 100 * Matrix::colSums(counts[ribo, , drop = FALSE]) / total
eligible <- nfeature >= 500 & nfeature <= 9000 & percent_mito <= 20 & percent_ribo > 5 & total > 0
eligible_counts <- as(counts[, eligible, drop = FALSE], "dgCMatrix")
eligible_counts <- eligible_counts[, order(colnames(eligible_counts), method = "radix"), drop = FALSE]
if (ncol(eligible_counts) != 4614L) {
  stop("frozen HCA_BM_002 QC scope changed: expected 4614 cells", call. = FALSE)
}
panel_source_ids <- rownames(counts)[panel_index]
raw_panel <- eligible_counts[panel_index, , drop = FALSE]
if (min(Matrix::rowSums(raw_panel > 0)) <= 0L) {
  stop("at least one exact-panel feature is absent after frozen QC", call. = FALSE)
}

started <- proc.time()[["elapsed"]]
object <- Seurat::CreateSeuratObject(
  counts = eligible_counts, min.cells = 0, min.features = 0, project = "HCA_BM_002"
)
object <- Seurat::SCTransform(
  object, assay = "RNA", new.assay.name = "SCT", variable.features.n = 3000L,
  return.only.var.genes = FALSE, verbose = FALSE, min_cells = 5L
)
sct_seconds <- proc.time()[["elapsed"]] - started

sct_data <- Seurat::GetAssayData(object, assay = "SCT", layer = "data")
sct_index <- match(panel_source_ids, rownames(sct_data))
if (anyNA(sct_index) || anyDuplicated(sct_index)) {
  stop("standard SCTransform did not retain the exact ordered 500 panel", call. = FALSE)
}
sct_data_panel <- as.matrix(sct_data[sct_index, colnames(eligible_counts), drop = FALSE])
rownames(sct_data_panel) <- exact$feature_id

residual_started <- proc.time()[["elapsed"]]
object <- Seurat::GetResidual(
  object = object, features = panel_source_ids, assay = "SCT", umi.assay = "RNA",
  replace.value = FALSE, na.rm = TRUE, verbose = FALSE
)
residual_seconds <- proc.time()[["elapsed"]] - residual_started
residuals <- methods::slot(object[["SCT"]], "scale.data")
residual_index <- match(panel_source_ids, rownames(residuals))
if (anyNA(residual_index) || anyDuplicated(residual_index)) {
  stop("GetResidual did not materialize the exact ordered 500 panel", call. = FALSE)
}
residual_panel <- as.matrix(residuals[residual_index, colnames(eligible_counts), drop = FALSE])
rownames(residual_panel) <- exact$feature_id

normalize_started <- proc.time()[["elapsed"]]
object <- Seurat::NormalizeData(
  object, assay = "RNA", normalization.method = "LogNormalize",
  scale.factor = 10000, verbose = FALSE
)
normalize_seconds <- proc.time()[["elapsed"]] - normalize_started
rna_data <- Seurat::GetAssayData(object, assay = "RNA", layer = "data")
rna_index <- match(panel_source_ids, rownames(rna_data))
if (anyNA(rna_index) || anyDuplicated(rna_index)) {
  stop("LogNormalize did not retain the exact ordered 500 panel", call. = FALSE)
}
lognorm_panel <- as.matrix(rna_data[rna_index, colnames(eligible_counts), drop = FALSE])
rownames(lognorm_panel) <- exact$feature_id

distance_geometry <- function(values, indices) {
  selected <- values[indices, , drop = FALSE]
  gene_sd <- apply(selected, 1L, stats::sd)
  threshold <- sqrt(.Machine$double.eps)
  zero_count <- sum(!is.finite(gene_sd) | gene_sd <= threshold)
  result <- list(
    finite = all(is.finite(selected)), zero_variance_gene_count = zero_count,
    minimum_gene_sd = suppressWarnings(min(gene_sd, na.rm = TRUE)), valid = FALSE,
    distance = NULL, distance_sha256 = NA_character_, offdiagonal_minimum = NA_real_,
    offdiagonal_median = NA_real_, offdiagonal_maximum = NA_real_
  )
  if (!result$finite || zero_count != 0L) return(result)
  correlation <- stats::cor(t(selected), method = "pearson")
  # Subscripted clipping preserves matrix dimensions.  pmax()/pmin() with a
  # scalar first argument can legally drop the dim attribute before diag<-.
  correlation[correlation < -1] <- -1
  correlation[correlation > 1] <- 1
  distance <- sqrt(2 * (1 - correlation))
  distance[distance < 0] <- 0
  diag(distance) <- 0
  off <- distance[upper.tri(distance)]
  result$valid <- all(is.finite(distance)) && isSymmetric(distance) &&
    all(diag(distance) == 0) && all(off >= 0)
  if (!result$valid) return(result)
  result$distance <- distance
  result$distance_sha256 <- sha_object(distance)
  result$offdiagonal_minimum <- min(off)
  result$offdiagonal_median <- stats::median(off)
  result$offdiagonal_maximum <- max(off)
  result
}

representations <- list(
  sct_data_log1p_corrected_umi = sct_data_panel,
  sct_pearson_residual = residual_panel,
  rna_lognormalize_10000 = lognorm_panel
)
representation_semantics <- c(
  sct_data_log1p_corrected_umi = "SCT assay data layer; log1p corrected UMI counts",
  sct_pearson_residual = "SCT assay scale.data; GetResidual from fitted SCT model and RNA counts",
  rna_lognormalize_10000 = "RNA assay data layer; LogNormalize scale.factor=10000"
)
representation_role <- c(
  sct_data_log1p_corrected_umi = "later_production_contract_comparator",
  sct_pearson_residual = "original_corrected_contract_candidate",
  rna_lognormalize_10000 = "diagnostic_only"
)

geometry <- lapply(representations, function(values) {
  list(exact = distance_geometry(values, seq_len(500L)),
       common = distance_geometry(values, common_index))
})

summary_rows <- lapply(names(representations), function(id) {
  values <- representations[[id]]
  exact_result <- geometry[[id]]$exact
  common_result <- geometry[[id]]$common
  data.frame(
    representation_id = id,
    representation_role = unname(representation_role[[id]]),
    semantics = unname(representation_semantics[[id]]),
    feature_count = nrow(values), cell_count = ncol(values), values_finite = all(is.finite(values)),
    exact500_zero_variance_gene_count = exact_result$zero_variance_gene_count,
    exact500_minimum_gene_sd = exact_result$minimum_gene_sd,
    exact500_correlation_chord_valid = exact_result$valid,
    exact500_distance_sha256 = exact_result$distance_sha256,
    exact500_offdiagonal_minimum = exact_result$offdiagonal_minimum,
    exact500_offdiagonal_median = exact_result$offdiagonal_median,
    exact500_offdiagonal_maximum = exact_result$offdiagonal_maximum,
    common475_zero_variance_gene_count = common_result$zero_variance_gene_count,
    common475_minimum_gene_sd = common_result$minimum_gene_sd,
    common475_correlation_chord_valid = common_result$valid,
    common475_distance_sha256 = common_result$distance_sha256,
    common475_offdiagonal_minimum = common_result$offdiagonal_minimum,
    common475_offdiagonal_median = common_result$offdiagonal_median,
    common475_offdiagonal_maximum = common_result$offdiagonal_maximum,
    exact500_viable = nrow(values) == 500L && ncol(values) == 4614L &&
      all(is.finite(values)) && exact_result$zero_variance_gene_count == 0L && exact_result$valid,
    labels_outcomes_opened = FALSE, persistence_computed = FALSE,
    landscapes_computed = FALSE, clustering_computed = FALSE, fusion_computed = FALSE,
    stringsAsFactors = FALSE
  )
})
representation_summary <- do.call(rbind, summary_rows)
atomic_csv(representation_summary, file.path(output_dir, "mv08m-representation-summary.csv"))

top_k_overlap <- function(left, right, k = 10L) {
  overlaps <- vapply(seq_len(nrow(left)), function(i) {
    left_order <- order(left[i, ], seq_len(nrow(left)), method = "radix")
    right_order <- order(right[i, ], seq_len(nrow(right)), method = "radix")
    left_neighbors <- head(left_order[left_order != i], k)
    right_neighbors <- head(right_order[right_order != i], k)
    length(intersect(left_neighbors, right_neighbors)) / k
  }, numeric(1L))
  mean(overlaps)
}

pairs <- utils::combn(names(representations), 2L, simplify = FALSE)
stability_rows <- lapply(pairs, function(pair) {
  left <- geometry[[pair[[1L]]]]$common$distance
  right <- geometry[[pair[[2L]]]]$common$distance
  pair_valid <- !is.null(left) && !is.null(right)
  if (!pair_valid) {
    return(data.frame(
      representation_a = pair[[1L]], representation_b = pair[[2L]], common475_pair_valid = FALSE,
      pair_count = 0L, pearson_distance_correlation = NA_real_, spearman_distance_correlation = NA_real_,
      median_absolute_distance_difference = NA_real_, p95_absolute_distance_difference = NA_real_,
      mean_top10_neighbor_overlap = NA_real_, interpretation = "descriptive_no_threshold",
      stringsAsFactors = FALSE
    ))
  }
  left_values <- left[upper.tri(left)]
  right_values <- right[upper.tri(right)]
  differences <- abs(left_values - right_values)
  data.frame(
    representation_a = pair[[1L]], representation_b = pair[[2L]], common475_pair_valid = TRUE,
    pair_count = length(left_values),
    pearson_distance_correlation = stats::cor(left_values, right_values, method = "pearson"),
    spearman_distance_correlation = stats::cor(left_values, right_values, method = "spearman"),
    median_absolute_distance_difference = stats::median(differences),
    p95_absolute_distance_difference = unname(stats::quantile(differences, 0.95, names = FALSE)),
    mean_top10_neighbor_overlap = top_k_overlap(left, right, 10L),
    interpretation = "descriptive_no_threshold", stringsAsFactors = FALSE
  )
})
stability <- do.call(rbind, stability_rows)
atomic_csv(stability, file.path(output_dir, "mv08m-common475-representation-stability.csv"))

validation <- data.frame(
  check_id = c(
    "input_panel_mapping", "frozen_reference_identity_binding", "qc_before_panel",
    "all_qc_cells_used", "one_standard_sct_fit", "three_representations_complete",
    "pearson_residual_source_bound", "lognormalize_source_bound", "common475_comparisons_complete",
    "stability_is_descriptive", "label_outcome_firewall", "downstream_firewall"
  ),
  passed = c(
    TRUE, TRUE, ncol(eligible_counts) == 4614L,
    identical(colnames(eligible_counts), sort(colnames(eligible_counts), method = "radix")),
    TRUE, nrow(representation_summary) == 3L && all(representation_summary$feature_count == 500L) &&
      all(representation_summary$cell_count == 4614L) && all(representation_summary$values_finite),
    all(panel_source_ids %in% rownames(residuals)),
    all(panel_source_ids %in% rownames(rna_data)),
    nrow(stability) == 3L && all(stability$common475_pair_valid) && all(stability$pair_count == choose(475L, 2L)),
    all(stability$interpretation == "descriptive_no_threshold"), TRUE, TRUE
  ),
  evidence = c(
    "all ordered exact500 stable IDs map exactly once; common475 is an ordered subset",
    "MV07 record binds panel and lineage; its data-layer center/scale is not transferred across representations",
    paste0("eligible=", ncol(eligible_counts), "; frozen QC precedes panel inspection"),
    "all eligible barcodes sorted deterministically; no gene-aware filter",
    paste0("SCTransform min_cells=5; return.only.var.genes=FALSE; seconds=", format(sct_seconds, digits = 12)),
    "SCT data, GetResidual Pearson residual, and RNA LogNormalize each evaluated at 500 x 4614",
    paste0("Seurat::GetResidual -> SCT scale.data; seconds=", format(residual_seconds, digits = 12)),
    paste0("Seurat::NormalizeData LogNormalize scale.factor=10000 -> RNA data; seconds=", format(normalize_seconds, digits = 12)),
    "three common475 pairwise comparisons each contain 112575 unique distances",
    "no MV08-G panel-sensitivity threshold repurposed for representation selection",
    "labels and outcomes remain closed",
    "no PH, landscapes, clustering, fusion, manuscript claim, other unit, or deletion"
  ), stringsAsFactors = FALSE
)
atomic_csv(validation, file.path(output_dir, "mv08m-validation.csv"))

identity <- data.frame(
  contract_id = "mv08m_exact500_gene_representation_prefreeze_v1", unit_id = "HCA_BM_002",
  filtered_h5_sha256 = sha_file(filtered_path), exact500_panel_sha256 = exact_sha,
  common475_panel_sha256 = common_sha, reference_record_sha256 = sha_file(reference_path),
  reference_payload_contract_id = if (!is.null(reference$identity$payload_contract_id))
    as.character(reference$identity$payload_contract_id) else "mv05d0_sct_data_matrix_v1_lineage",
  eligible_cells = ncol(eligible_counts), sct_fit_count = 1L,
  seurat_version = as.character(utils::packageVersion("Seurat")),
  sctransform_version = as.character(utils::packageVersion("sctransform")),
  total_elapsed_seconds = proc.time()[["elapsed"]] - started,
  labels_outcomes_opened = FALSE, persistence_computed = FALSE,
  landscapes_computed = FALSE, clustering_computed = FALSE, fusion_computed = FALSE,
  stringsAsFactors = FALSE
)
atomic_csv(identity, file.path(output_dir, "mv08m-identity.csv"))

cat(
  "MV8-M representation audit checks=", sum(validation$passed), "/", nrow(validation),
  "; viable=", paste(representation_summary$representation_id[representation_summary$exact500_viable], collapse = ","),
  "\n", sep = ""
)
quit(status = if (all(validation$passed)) 0L else 2L)
