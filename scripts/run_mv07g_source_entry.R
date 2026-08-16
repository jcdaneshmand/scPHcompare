#!/usr/bin/env Rscript

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1")
options(warn = 2)
for (package in c("digest", "Matrix")) {
  if (!requireNamespace(package, quietly = TRUE)) stop(package, " required.")
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop("usage: run_mv07g_source_entry.R PREFREEZE PRIMARY_CACHE ADDED_CACHE OUTPUT SEED")
}
prefreeze <- args[[1]]
primary_cache <- args[[2]]
added_cache <- args[[3]]
output <- args[[4]]
seed <- as.integer(args[[5]])
source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/mv03_pilot.R")
source("R/mv05_resource_safe_execution.R")
source("R/mv07g_sentinel.R")
sha <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}
manifest <- read.csv(file.path(prefreeze, "mv07g-cache-manifest.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)
axis <- read.csv(file.path(prefreeze, "mv07g-sentinel-axis.csv"),
                 stringsAsFactors = FALSE, check.names = FALSE)
panel <- read.csv(file.path(prefreeze, "mv07g-panel.csv"),
                  stringsAsFactors = FALSE, check.names = FALSE)
part <- manifest[manifest$seed == seed, , drop = FALSE]
part <- part[order(part$sample_id, method = "radix"), , drop = FALSE]
sentinel_ids <- sort(axis$sample_id[axis$seed == seed], method = "radix")
if (nrow(part) != 124L || length(sentinel_ids) != 6L ||
    nrow(panel) != 500L ||
    !identical(panel$panel_order, seq_len(500L)) ||
    length(unique(panel$panel_sha256)) != 1L) {
  stop("MV7-G source axes differ from the frozen contract.")
}
paths <- ifelse(
  part$source_tier == "primary90",
  file.path(primary_cache, part$private_cache_file),
  file.path(added_cache, part$private_cache_file)
)
if (!all(file.exists(paths)) ||
    !identical(unname(vapply(paths, sha, character(1L))),
               unname(part$private_cache_sha256))) {
  stop("MV7-G source cache identity drift.")
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
  stop("MV7-G source cache payload identity drift.")
}
matrices <- lapply(records, mv05d0_sct_matrix_from_cache_v1)
names(matrices) <- part$sample_id
selected_cells <- lapply(matrices, colnames)
input_cache_keys <- stats::setNames(
  vapply(records, `[[`, character(1L), "cache_key"), part$sample_id)
selected_cell_sha256 <- stats::setNames(
  vapply(records, function(record) record$identity$selected_cell_sha256,
         character(1L)), part$sample_id)
rm(records)
invisible(gc())
fit_scope_id <- paste0("mv07e_full124_seed_", seed)
prepared <- prepare_mv03_sources(
  matrices, panel$feature_id,
  cohort = "mv07e_full124_global_descriptive",
  representation = "sct_global_descriptive_v1",
  fit_scope_id = fit_scope_id, seed = seed,
  selected_cells = selected_cells
)
rm(matrices)
invisible(gc())
pca_model <- fit_cell_topology_pca(
  prepared$sources, n_components = 30L, pca_seed = seed
)
views <- lapply(sentinel_ids, function(sample_id) {
  source_object <- prepared$sources[[sample_id]]
  list(
    cell_topology_v1 = construct_cell_topology_view(
      source_object, pca_model, n_components = 30L),
    gene_topology_v1 = construct_gene_topology_view(source_object)
  )
})
names(views) <- sentinel_ids
identity <- list(
  contract_id = "mv07g_source_identity_v1",
  seed = seed, fit_scope_id = fit_scope_id,
  fit_samples = 124L, fit_cells = 47616L, panel_size = 500L,
  panel_sha256 = unique(panel$panel_sha256),
  panel_file_sha256 = sha(file.path(prefreeze, "mv07g-panel.csv")),
  cache_manifest_sha256 = sha(file.path(prefreeze, "mv07g-cache-manifest.csv")),
  input_cache_keys = input_cache_keys,
  selected_cell_sha256 = selected_cell_sha256,
  sentinel_ids = sentinel_ids,
  standardization_id = prepared$standardization_id,
  pca_model_cache_key = pca_model$cache_key,
  view_cache_keys = lapply(views, function(pair) {
    vapply(pair, `[[`, character(1L), "cache_key")
  }),
  outcome_label_state = "closed", biological_outcomes_computed = FALSE
)
record <- mv07g_new_source_record_v1(
  identity, panel[c("panel_order", "feature_id", "gene")],
  prepared, pca_model, views
)
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
if (file.exists(output)) {
  existing <- readRDS(output)
  mv07g_validate_source_record_v1(existing)
  if (!identical(existing$cache_key, record$cache_key)) {
    stop("Existing MV7-G source bundle has stale identity; refusing overwrite.")
  }
  message("Reused MV7-G source bundle: ", basename(output))
} else {
  partial <- tempfile(pattern = basename(output), tmpdir = dirname(output))
  saveRDS(record, partial, compress = FALSE, version = 3)
  if (!file.rename(partial, output)) {
    unlink(partial)
    stop("Failed to atomically publish MV7-G source bundle.")
  }
  message("Built MV7-G source bundle: ", basename(output))
}
