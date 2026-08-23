#!/usr/bin/env Rscript

# Construct one current-input selected-384-fit common475 baseline. This child
# opens only the two explicitly bound HCA matrices and never runs PH.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) stop(paste(
  "usage: run_mv08s_same_axis_baseline_entry.R <mv08s-prefreeze>",
  "<unit-id> <filtered-h5> <raw-h5> <reference-rds> <common-panel> <output-rds>"
), call. = FALSE)
prefreeze <- normalizePath(args[[1L]], mustWork = TRUE)
unit_id <- args[[2L]]
filtered_path <- normalizePath(args[[3L]], mustWork = TRUE)
raw_path <- normalizePath(args[[4L]], mustWork = TRUE)
reference_path <- normalizePath(args[[5L]], mustWork = TRUE)
panel_path <- normalizePath(args[[6L]], mustWork = TRUE)
output_path <- normalizePath(args[[7L]], mustWork = FALSE)
if (file.exists(output_path)) stop("refusing to overwrite MV8-S baseline", call. = FALSE)
for (package in c("digest", "Matrix", "Seurat")) {
  if (!requireNamespace(package, quietly = TRUE)) stop(package, " required", call. = FALSE)
}
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
if (requireNamespace("future", quietly = TRUE)) future::plan(future::sequential)
options(warn = 1)
source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/mv07g_sentinel.R")
source("R/mv08g_panel_sensitivity.R")
source("R/mv08s_ph_sentinel.R")

sha_file <- function(path) digest::digest(
  file = path, algo = "sha256", serialize = FALSE
)
binding <- utils::read.csv(
  file.path(prefreeze, "mv08s-external-input-bindings.csv"),
  check.names = FALSE, stringsAsFactors = FALSE
)
binding <- binding[binding$unit_id == unit_id, , drop = FALSE]
if (nrow(binding) != 1L || binding$selected_cells != 384L ||
    binding$selection_seed != 20260805L || binding$panel_genes != 475L ||
    binding$outcome_label_state != "closed" ||
    isTRUE(binding$biological_outcomes_computed) ||
    sha_file(filtered_path) != binding$filtered_h5_sha256 ||
    sha_file(raw_path) != binding$raw_h5_sha256 ||
    sha_file(reference_path) != binding$reference_rds_sha256 ||
    sha_file(panel_path) != binding$panel_file_sha256) {
  stop("MV8-S external baseline input binding drift", call. = FALSE)
}

panel <- utils::read.csv(panel_path, check.names = FALSE, stringsAsFactors = FALSE)
if (nrow(panel) != 475L ||
    !identical(as.integer(panel$panel_order_475), 1:475) ||
    !all(panel$common475_axis_sha256 == binding$panel_sha256)) {
  stop("MV8-S common475 panel drift", call. = FALSE)
}
stable <- sub("\\.[0-9]+$", "", sub("^.*-", "", panel$feature_id))
reference <- readRDS(reference_path)
mv08g_validate_source_record_v1(reference)
if (nrow(reference$panel) != 475L || length(reference$center) != 475L ||
    length(reference$scale) != 475L || reference$pca_model$n_components != 30L ||
    reference$cache_key != binding$reference_cache_key) {
  stop("MV8-S common475 reference drift", call. = FALSE)
}

filtered <- Seurat::Read10X_h5(
  filtered_path, use.names = FALSE, unique.features = TRUE
)
raw <- Seurat::Read10X_h5(raw_path, use.names = FALSE, unique.features = TRUE)
filtered_names <- Seurat::Read10X_h5(
  filtered_path, use.names = TRUE, unique.features = TRUE
)
if (is.list(filtered) || is.list(raw) || is.list(filtered_names) ||
    !identical(rownames(filtered), rownames(raw)) ||
    !identical(dim(filtered_names), dim(filtered))) {
  stop("MV8-S HCA matrix axes are incompatible", call. = FALSE)
}
feature_ids <- sub("\\.[0-9]+$", "", rownames(filtered))
panel_index <- match(stable, feature_ids)
if (anyNA(panel_index) || anyDuplicated(panel_index)) {
  stop("MV8-S panel mapping is incomplete or ambiguous", call. = FALSE)
}
counts <- methods::as(filtered, "dgCMatrix")
nfeature <- Matrix::colSums(counts > 0)
total <- Matrix::colSums(counts)
names_text <- rownames(filtered_names)
mito <- grepl("^MT-", names_text)
ribo <- grepl("^(RPS|RPL)", names_text)
percent_mito <- 100 * Matrix::colSums(counts[mito, , drop = FALSE]) / total
percent_ribo <- 100 * Matrix::colSums(counts[ribo, , drop = FALSE]) / total
eligible <- nfeature >= 500 & nfeature <= 9000 & percent_mito <= 20 &
  percent_ribo > 5
if (sum(eligible) != binding$qc_eligible_cells || sum(eligible) < 384L) {
  stop("MV8-S frozen HCA QC axis drift", call. = FALSE)
}
selected <- select_matched_cells(
  colnames(counts)[eligible], n = 384L, seed = 20260805L
)
if (attr(selected, "selected_cell_sha256") != binding$selected_cell_sha256) {
  stop("MV8-S selected-cell digest drift", call. = FALSE)
}
raw_index <- match(selected, colnames(raw))
if (anyNA(raw_index)) stop("MV8-S selected cells are absent from raw H5", call. = FALSE)
raw_selected <- methods::as(raw[, raw_index, drop = FALSE], "dgCMatrix")
colnames(raw_selected) <- selected
if (any(Matrix::colSums(raw_selected) <= 0)) {
  stop("MV8-S selected raw matrix contains zero-count cells", call. = FALSE)
}

object <- Seurat::CreateSeuratObject(
  counts = raw_selected, min.cells = 0, min.features = 0, project = unit_id
)
object <- Seurat::SCTransform(
  object, assay = "RNA", new.assay.name = "SCT", min_cells = 5L,
  variable.features.n = 3000L, return.only.var.genes = FALSE,
  verbose = FALSE
)
sct <- Seurat::GetAssayData(object, assay = "SCT", layer = "data")
sct_index <- match(rownames(raw)[panel_index], rownames(sct))
if (anyNA(sct_index) || !all(selected %in% colnames(sct))) {
  stop("MV8-S SCT baseline does not retain the frozen axes", call. = FALSE)
}
sct <- as.matrix(sct[sct_index, selected, drop = FALSE])
storage.mode(sct) <- "double"
rownames(sct) <- panel$feature_id
colnames(sct) <- selected
if (any(!is.finite(sct))) stop("MV8-S SCT baseline is nonfinite", call. = FALSE)

standardized <- sweep(
  sweep(sct, 1L, reference$center, "-"), 1L, reference$scale, "/"
)
source_record <- new_dual_view_source(
  standardized,
  sample_id = unit_id,
  cohort = "external8",
  representation = "sct_data_selected384_fit_same_axis",
  fit_scope_id = reference$pca_model$fit_scope_id,
  subsample_seed = 20260805L,
  standardization_id = reference$pca_model$standardization_id,
  contract_profile = "scientific_common475"
)
coordinates <- t(standardized) %*% reference$pca_model$rotation
cell_view <- construct_frozen_cell_topology_view(
  source_record, coordinates,
  coordinate_contract_id = "mv08s_immutable_common475_projection_v1",
  coordinate_fit_cache_key = reference$pca_model$cache_key
)
gene_view <- construct_gene_topology_view(source_record)
identity <- list(
  contract_id = "mv08s_same_axis_external_baseline_identity_v1",
  unit_id = unit_id,
  dataset_scope = "external8",
  filtered_h5_sha256 = binding$filtered_h5_sha256,
  raw_h5_sha256 = binding$raw_h5_sha256,
  filtered_cells = ncol(filtered),
  qc_eligible_cells = sum(eligible),
  selected_cells = 384L,
  selection_seed = 20260805L,
  selected_cell_sha256 = attr(selected, "selected_cell_sha256"),
  fit_scope = "same_frozen_selected384_only",
  representation_id = "sct_data_selected384_fit_same_axis",
  normalization = list(
    method = "SCTransform", min_cells = 5L, variable_features_n = 3000L,
    return_only_var_genes = FALSE, workers = 1L
  ),
  panel_id = "common475",
  panel_sha256 = binding$panel_sha256,
  reference_file_sha256 = binding$reference_rds_sha256,
  reference_cache_key = binding$reference_cache_key,
  standardized_matrix_sha256 = .mv08s_sha_object(standardized),
  cell_view_payload_sha256 = cell_view$payload_sha256,
  gene_view_payload_sha256 = gene_view$payload_sha256,
  outcome_label_state = "closed",
  biological_outcomes_computed = FALSE
)
record <- mv08s_new_baseline_record_v1(
  identity, panel[c("panel_order_475", "feature_id", "gene",
                    "common475_axis_sha256")],
  list(cell_topology_v1 = cell_view, gene_topology_v1 = gene_view)
)
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
partial <- tempfile(pattern = paste0(basename(output_path), "."),
                    tmpdir = dirname(output_path))
on.exit(if (file.exists(partial)) unlink(partial), add = TRUE)
saveRDS(record, partial, compress = FALSE, version = 3)
if (!file.rename(partial, output_path)) {
  stop("failed to atomically publish MV8-S baseline", call. = FALSE)
}
cat("MV8-S baseline unit=", unit_id, "; views=2; selected=384; checks=pass\n",
    sep = "")
