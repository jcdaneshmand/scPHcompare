#!/usr/bin/env Rscript

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1")
options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("usage: run_mv06d_ph_entry.R SOURCE_RDS ROLE VIEW_ID OUTPUT_RDS",
       call. = FALSE)
}
source_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
role <- match.arg(args[[2L]], c("held_out", "training"))
view_id <- match.arg(args[[3L]], c("cell_topology_v1", "gene_topology_v1"))
output_path <- args[[4L]]

source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/mv06d_matched_profile.R")
record <- readRDS(source_path)
mv06d_validate_source_record_v1(record)
view <- record$payload$views[[role]][[view_id]]
result <- run_topology_view_ph(view, max_dim = 1L, threshold = -1, field = 2L)
ph <- mv06d_new_ph_record_v1(
  record$cache_key, record$identity$sentinel_ids[[role]], role, view, result
)
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
if (file.exists(output_path)) {
  existing <- readRDS(output_path)
  mv06d_validate_ph_record_v1(existing)
  if (!identical(existing$cache_key, ph$cache_key)) {
    stop("Existing MV6-D PH record has a stale identity; refusing overwrite.",
         call. = FALSE)
  }
  message("Reused validated MV6-D PH record: ", basename(output_path))
} else {
  partial <- tempfile(pattern = basename(output_path),
                      tmpdir = dirname(output_path))
  saveRDS(ph, partial, compress = FALSE, version = 3)
  if (!file.rename(partial, output_path)) {
    unlink(partial)
    stop("Failed to atomically publish MV6-D PH record.", call. = FALSE)
  }
  message("Built MV6-D PH record: ", basename(output_path))
}
