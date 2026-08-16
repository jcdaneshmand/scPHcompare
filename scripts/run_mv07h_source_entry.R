#!/usr/bin/env Rscript

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1")
options(warn = 2)
for (package in c("digest", "Matrix")) {
  if (!requireNamespace(package, quietly = TRUE)) stop(package, " required.")
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop("usage: run_mv07h_source_entry.R PREFREEZE PRIMARY_CACHE ADDED_CACHE MV07G_SOURCE OUTPUT SEED")
}
prefreeze <- args[[1L]]
primary_cache <- args[[2L]]
added_cache <- args[[3L]]
parent_path <- args[[4L]]
output <- args[[5L]]
seed <- as.integer(args[[6L]])
source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/mv03_pilot.R")
source("R/mv05_resource_safe_execution.R")
source("R/mv07g_sentinel.R")
source("R/mv07h_full_topology.R")

manifest <- read.csv(file.path(prefreeze, "mv07h-cache-manifest.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)
panel <- read.csv(file.path(prefreeze, "mv07h-panel.csv"),
                  stringsAsFactors = FALSE, check.names = FALSE)
mv07h_validate_cache_manifest_v1(manifest)
part <- manifest[manifest$seed == seed, , drop = FALSE]
part <- part[order(part$sample_id, method = "radix"), , drop = FALSE]
if (nrow(part) != 124L || nrow(panel) != 500L ||
    !identical(as.integer(panel$panel_order), seq_len(500L))) {
  stop("MV7-H source axes differ from the prefreeze.")
}
parent <- readRDS(parent_path)
mv07g_validate_source_record_v1(parent)
if (parent$identity$seed != seed ||
    !identical(parent$identity$panel_sha256, unique(panel$panel_sha256)) ||
    !identical(parent$pca_model$fit_sample_ids, part$sample_id)) {
  stop("MV7-H parent transform differs from the exact seed axis.")
}
paths <- ifelse(
  part$source_tier == "primary90",
  file.path(primary_cache, part$private_cache_file),
  file.path(added_cache, part$private_cache_file)
)
if (!all(file.exists(paths)) ||
    !identical(unname(vapply(paths, .mv07h_sha256, character(1L))),
               unname(part$private_cache_sha256))) {
  stop("MV7-H input cache identity drift.")
}
views <- vector("list", nrow(part))
names(views) <- part$sample_id
for (index in seq_len(nrow(part))) {
  cache <- readRDS(paths[[index]])
  mv05d0_validate_normalization_cache_record_v2(cache)
  if (cache$identity$sample_id != part$sample_id[[index]] ||
      cache$identity$seed != seed ||
      cache$cache_key != part$normalization_cache_key[[index]] ||
      cache$payload_sha256 != part$payload_sha256[[index]] ||
      cache$identity$selected_cell_sha256 !=
        parent$identity$selected_cell_sha256[[part$sample_id[[index]]]]) {
    stop("MV7-H normalization cache payload identity drift.")
  }
  matrix <- mv05d0_sct_matrix_from_cache_v1(cache)
  if (!all(panel$feature_id %in% rownames(matrix)) || ncol(matrix) != 384L) {
    stop("MV7-H cache lacks the exact panel/cell shape.")
  }
  value <- as.matrix(matrix[panel$feature_id, , drop = FALSE])
  value <- sweep(sweep(value, 1L, parent$center, "-"),
                 1L, parent$scale, "/")
  source_object <- new_dual_view_source(
    value, sample_id = part$sample_id[[index]],
    cohort = parent$pca_model$cohort,
    representation = parent$pca_model$representation,
    fit_scope_id = parent$pca_model$fit_scope_id,
    subsample_seed = seed,
    standardization_id = parent$pca_model$standardization_id
  )
  expected_source_key <- parent$pca_model$source_cache_keys[[
    part$sample_id[[index]]
  ]]
  if (!identical(source_object$cache_key, expected_source_key)) {
    stop("MV7-H reconstructed source differs from the accepted PCA fit input.")
  }
  views[[index]] <- list(
    cell_topology_v1 = construct_cell_topology_view(
      source_object, parent$pca_model, n_components = 30L),
    gene_topology_v1 = construct_gene_topology_view(source_object)
  )
  rm(cache, matrix, value, source_object)
  if (index %% 10L == 0L) invisible(gc(FALSE))
}
for (sample_id in names(parent$views)) {
  for (view_id in .mv07h_views) {
    old <- parent$views[[sample_id]][[view_id]]
    new <- views[[sample_id]][[view_id]]
    if (old$cache_key != new$cache_key ||
        old$payload_sha256 != new$payload_sha256) {
      stop("MV7-H sentinel typed view does not reproduce MV7-G exactly.")
    }
  }
}
record <- mv07h_new_source_record_v1(parent, views)
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
if (file.exists(output)) {
  existing <- readRDS(output)
  mv07h_validate_source_record_v1(existing)
  if (existing$cache_key != record$cache_key) {
    stop("Existing MV7-H source bundle is stale; refusing overwrite.")
  }
  message("Reused MV7-H source bundle: ", basename(output))
} else {
  partial <- tempfile(pattern = basename(output), tmpdir = dirname(output))
  saveRDS(record, partial, compress = FALSE, version = 3)
  if (!file.rename(partial, output)) {
    unlink(partial)
    stop("Failed to atomically publish MV7-H source bundle.")
  }
  message("Built MV7-H source bundle: ", basename(output))
}
