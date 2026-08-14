#!/usr/bin/env Rscript

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1")
options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("usage: run_mv06d_landscape_entry.R HELD_OUT_PH TRAINING_PH OUTPUT_RDS",
       call. = FALSE)
}
first_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
second_path <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
output_path <- args[[3L]]

source("R/landscape_contract.R")
source("R/landscape_reference.R")
source("R/landscape_public_api.R")
source("R/mv06d_matched_profile.R")
first <- readRDS(first_path)
second <- readRDS(second_path)
mv06d_validate_ph_record_v1(first)
mv06d_validate_ph_record_v1(second)
result <- persistence_landscape_distance(
  first$topology_result$diagram, second$topology_result$diagram,
  method = "auto", exact_max_intervals = 500L, abs_tol = 1e-8,
  rel_tol = 1e-8, subdivisions = 200L,
  first_id = first$cache_key, second_id = second$cache_key
)
record <- mv06d_new_landscape_record_v1(first, second, result)
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
if (file.exists(output_path)) {
  existing <- readRDS(output_path)
  if (!inherits(existing, "scph_mv06d_landscape_record_v1") ||
      !identical(existing$payload_sha256, digest::digest(
        list(identity = existing$identity, result = existing$result),
        algo = "sha256", serialize = TRUE
      )) || !identical(existing$cache_key, record$cache_key)) {
    stop("Existing MV6-D landscape record is stale; refusing overwrite.",
         call. = FALSE)
  }
  message("Reused validated MV6-D landscape record: ", basename(output_path))
} else {
  partial <- tempfile(pattern = basename(output_path),
                      tmpdir = dirname(output_path))
  saveRDS(record, partial, compress = FALSE, version = 3)
  if (!file.rename(partial, output_path)) {
    unlink(partial)
    stop("Failed to atomically publish MV6-D landscape record.", call. = FALSE)
  }
  message("Built MV6-D landscape record: ", basename(output_path))
}
