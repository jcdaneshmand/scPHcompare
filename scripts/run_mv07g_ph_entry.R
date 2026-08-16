#!/usr/bin/env Rscript

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1")
options(warn = 2)
if (!requireNamespace("ripserr", quietly = TRUE)) stop("ripserr required.")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("usage: run_mv07g_ph_entry.R SOURCE SAMPLE_ID VIEW_ID OUTPUT")
}
source_path <- args[[1]]
sample_id <- args[[2]]
view_id <- match.arg(args[[3]], c("cell_topology_v1", "gene_topology_v1"))
output <- args[[4]]
source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/mv07g_sentinel.R")
source_record <- readRDS(source_path)
mv07g_validate_source_record_v1(source_record)
if (!sample_id %in% names(source_record$views)) {
  stop("MV7-G PH sample is absent from the source bundle.")
}
view <- source_record$views[[sample_id]][[view_id]]
result <- run_topology_view_ph(view, max_dim = 1L, threshold = -1, field = 2L)
record <- mv07g_new_ph_record_v1(source_record, sample_id, view_id, result)
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
if (file.exists(output)) {
  existing <- readRDS(output)
  mv07g_validate_ph_record_v1(existing, view)
  if (!identical(existing$cache_key, record$cache_key)) {
    stop("Existing MV7-G PH record has stale identity; refusing overwrite.")
  }
  message("Reused MV7-G PH record: ", basename(output))
} else {
  partial <- tempfile(pattern = basename(output), tmpdir = dirname(output))
  saveRDS(record, partial, compress = FALSE, version = 3)
  if (!file.rename(partial, output)) {
    unlink(partial)
    stop("Failed to atomically publish MV7-G PH record.")
  }
  message("Built MV7-G PH record: ", basename(output))
}
