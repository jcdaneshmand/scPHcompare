#!/usr/bin/env Rscript

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1")
options(warn = 2)
if (!requireNamespace("TDA", quietly = TRUE)) stop("TDA required.")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("usage: run_mv08g_ph_fallback_entry.R SOURCE SAMPLE_ID OUTPUT")
}
source_path <- args[[1L]]
sample_id <- args[[2L]]
output <- args[[3L]]
source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/mv07g_sentinel.R")
source("R/mv07h_full_topology.R")
source("R/mv08g_panel_sensitivity.R")
source_record <- readRDS(source_path)
mv08g_validate_source_record_v1(source_record)
if (!sample_id %in% names(source_record$views)) {
  stop("MV8-G fallback sample is absent from the source bundle.")
}
view <- source_record$views[[sample_id]][["gene_topology_v1"]]
result <- mv07h_run_topology_view_ph_gudhi_v1(
  view, max_dim = 1L, threshold = -1, field = 2L)
record <- mv08g_new_ph_record_v1(
  source_record, sample_id, "gene_topology_v1", result)
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
if (file.exists(output)) {
  existing <- readRDS(output)
  mv08g_validate_ph_record_v1(existing, view)
  if (existing$cache_key != record$cache_key) {
    stop("Existing MV8-G fallback record is stale; refusing overwrite.")
  }
  message("Reused MV8-G fallback record: ", basename(output))
} else {
  partial <- tempfile(pattern = basename(output), tmpdir = dirname(output))
  saveRDS(record, partial, compress = FALSE, version = 3)
  if (!file.rename(partial, output)) {
    unlink(partial)
    stop("Failed to atomically publish MV8-G fallback record.")
  }
  message("Built MV8-G GUDHI fallback record: ", basename(output))
}
