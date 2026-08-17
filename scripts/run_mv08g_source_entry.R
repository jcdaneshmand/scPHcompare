#!/usr/bin/env Rscript

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1")
options(warn = 2)
for (package in c("digest", "Matrix")) {
  if (!requireNamespace(package, quietly = TRUE)) stop(package, " required.")
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop("usage: run_mv08g_source_entry.R PREFREEZE PRIMARY_CACHE ADDED_CACHE OUTPUT SEED")
}
prefreeze <- args[[1L]]
primary_cache <- args[[2L]]
added_cache <- args[[3L]]
output <- args[[4L]]
seed <- as.integer(args[[5L]])
source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/mv03_pilot.R")
source("R/mv05_resource_safe_execution.R")
source("R/mv07g_sentinel.R")
source("R/mv07h_full_topology.R")
source("R/mv08g_panel_sensitivity.R")
sha <- function(path) digest::digest(file = path, algo = "sha256",
                                     serialize = FALSE)

manifest <- read.csv(file.path(prefreeze, "mv08g-cache-manifest.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)
panel <- read.csv(file.path(prefreeze, "mv08g-panel.csv"),
                  stringsAsFactors = FALSE, check.names = FALSE)
mv07h_validate_cache_manifest_v1(manifest)
mv08g_validate_common475_panel_v1(panel)
part <- manifest[manifest$seed == seed, , drop = FALSE]
part <- part[order(part$sample_id, method = "radix"), , drop = FALSE]
if (nrow(part) != 124L || !(seed %in% .mv08g_seeds)) {
  stop("MV8-G source axis differs from the prefreeze.")
}
paths <- ifelse(part$source_tier == "primary90",
  file.path(primary_cache, part$private_cache_file),
  file.path(added_cache, part$private_cache_file))
if (!all(file.exists(paths)) ||
    !identical(unname(vapply(paths, sha, character(1L))),
               unname(part$private_cache_sha256))) {
  stop("MV8-G input cache identity drift.")
}
records <- lapply(paths, readRDS)
invisible(lapply(records, mv05d0_validate_normalization_cache_record_v2))
if (!identical(vapply(records, function(record) record$identity$sample_id,
                       character(1L)), part$sample_id) ||
    !identical(vapply(records, function(record) record$identity$seed,
                       integer(1L)), rep(seed, 124L)) ||
    !identical(vapply(records, `[[`, character(1L), "cache_key"),
               part$normalization_cache_key) ||
    !identical(vapply(records, `[[`, character(1L), "payload_sha256"),
               part$payload_sha256)) {
  stop("MV8-G normalization cache payload identity drift.")
}
matrices <- lapply(records, mv05d0_sct_matrix_from_cache_v1)
names(matrices) <- part$sample_id
selected_cells <- lapply(matrices, colnames)
input_cache_keys <- stats::setNames(
  vapply(records, `[[`, character(1L), "cache_key"), part$sample_id)
selected_cell_sha256 <- stats::setNames(vapply(records, function(record) {
  record$identity$selected_cell_sha256
}, character(1L)), part$sample_id)
rm(records)
invisible(gc())

fit_scope_id <- paste0("mv08g_full124_common475_seed_", seed)
prepared <- prepare_mv03_sources(
  matrices, panel$feature_id,
  cohort = "mv08g_full124_reference_sensitivity",
  representation = "sct_common475_reference_sensitivity_v1",
  fit_scope_id = fit_scope_id, seed = seed, selected_cells = selected_cells,
  contract_profile = "scientific_common475", expected_genes = 475L,
  expected_cells = 384L, expected_pcs = 30L)
rm(matrices)
invisible(gc())
pca_model <- fit_cell_topology_pca(
  prepared$sources, n_components = 30L, pca_seed = seed)
sample_ids <- sort(part$sample_id, method = "radix")
views <- lapply(sample_ids, function(sample_id) {
  source_object <- prepared$sources[[sample_id]]
  list(
    cell_topology_v1 = construct_cell_topology_view(
      source_object, pca_model, n_components = 30L),
    gene_topology_v1 = construct_gene_topology_view(source_object))
})
names(views) <- sample_ids
identity <- list(
  contract_id = "mv08g_source_identity_v1", seed = seed,
  fit_scope_id = fit_scope_id, fit_samples = 124L, fit_cells = 47616L,
  panel_size = 475L,
  common475_axis_sha256 = unique(panel$common475_axis_sha256),
  parent_panel_sha256 = unique(panel$parent_panel_sha256),
  panel_file_sha256 = sha(file.path(prefreeze, "mv08g-panel.csv")),
  cache_manifest_sha256 = sha(file.path(prefreeze, "mv08g-cache-manifest.csv")),
  input_cache_keys = input_cache_keys,
  selected_cell_sha256 = selected_cell_sha256,
  sample_ids = sample_ids,
  standardization_id = prepared$standardization_id,
  pca_model_cache_key = pca_model$cache_key,
  view_cache_keys = lapply(views, function(pair) {
    vapply(pair, `[[`, character(1L), "cache_key")
  }), outcome_label_state = "closed", biological_outcomes_computed = FALSE)
record <- mv08g_new_source_record_v1(
  identity, panel[c("panel_order_475", "panel_order_500", "feature_id", "gene")],
  prepared, pca_model, views)
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
if (file.exists(output)) {
  existing <- readRDS(output)
  mv08g_validate_source_record_v1(existing)
  if (existing$cache_key != record$cache_key) {
    stop("Existing MV8-G source bundle is stale; refusing overwrite.")
  }
  message("Reused MV8-G source bundle: ", basename(output))
} else {
  partial <- tempfile(pattern = basename(output), tmpdir = dirname(output))
  saveRDS(record, partial, compress = FALSE, version = 3)
  if (!file.rename(partial, output)) {
    unlink(partial)
    stop("Failed to atomically publish MV8-G source bundle.")
  }
  message("Built MV8-G source bundle: ", basename(output))
}
