#!/usr/bin/env Rscript

Sys.setenv(
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1"
)
options(warn = 2)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(
    "usage: run_mv05d2_ph_entry.R PILOT_CSV ROW_INDEX FOLD_CACHE_DIR ",
    "RESULT_RDS",
    call. = FALSE
  )
}
manifest_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
row_index <- as.integer(args[[2L]])
fold_cache_dir <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
result_path <- args[[4L]]
dir.create(dirname(result_path), recursive = TRUE, showWarnings = FALSE)

source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/mv03_pilot.R")
source("R/mv05_resource_safe_execution.R")
source("R/mv05d2_ph_profiling.R")

file_sha <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}
implementation_files <- c(
  "R/toy_baseline.R", "R/dual_view_topology.R", "R/mv03_pilot.R",
  "R/mv05_resource_safe_execution.R", "R/mv05d2_ph_profiling.R"
)
implementation_sha <- digest::digest(
  stats::setNames(vapply(implementation_files, file_sha, character(1L)),
                  implementation_files),
  algo = "sha256", serialize = TRUE
)
manifest <- utils::read.csv(
  manifest_path, stringsAsFactors = FALSE, check.names = FALSE
)
if (is.na(row_index) || row_index < 1L || row_index > nrow(manifest)) {
  stop("MV5-D2 row index is outside the pilot manifest.", call. = FALSE)
}
job <- manifest[row_index, , drop = FALSE]
if (job$outcome_label_state != "closed" ||
    as.logical(job$biological_outcomes_computed) ||
    job$point_count != 384L || job$coordinate_count != 30L ||
    job$view_id != "cell_topology_v1") {
  stop("Requested MV5-D2 job violates the pilot contract.", call. = FALSE)
}
fold_cache_path <- file.path(fold_cache_dir, job$fold_cache_file)
if (!file.exists(fold_cache_path) ||
    !identical(file_sha(fold_cache_path), job$fold_cache_sha256)) {
  stop("MV5-D1 source cache is missing or differs from the manifest.",
       call. = FALSE)
}
fold_record <- readRDS(fold_cache_path)
mv05d1_validate_cell_fold_record_v1(fold_record)
if (!identical(fold_record$cache_key, job$fold_cache_key)) {
  stop("MV5-D1 source cache identity differs from the pilot manifest.",
       call. = FALSE)
}
view <- fold_record$payload$cell_views[[job$sample_id]]
validate_topology_view(view)
if (is.null(view) || !identical(view$cache_key, job$view_cache_key) ||
    !identical(view$payload_sha256, job$view_payload_sha256)) {
  stop("MV5-D1 typed view differs from the pilot manifest.", call. = FALSE)
}
identity <- mv05d2_ph_identity_v1(
  job_id = job$job_id, fold_id = job$fold_id,
  fit_scope_id = job$fit_scope_id, held_out_study = job$held_out_study,
  seed = job$seed, sample_id = job$sample_id,
  execution_role = job$execution_role,
  missing_feature_count = job$missing_feature_count,
  fold_cache_key = job$fold_cache_key,
  fold_cache_sha256 = job$fold_cache_sha256,
  view_cache_key = job$view_cache_key,
  view_payload_sha256 = job$view_payload_sha256,
  pilot_manifest_sha256 = file_sha(manifest_path),
  implementation_sha256 = implementation_sha,
  runtime = mv05d2_ph_runtime_v1()
)
disposition <- mv05d2_ph_cache_disposition_v1(
  result_path, identity$cache_key
)
if (identical(disposition, "reuse_validated")) {
  message("Reused validated MV5-D2 PH cache: ", job$job_id)
  quit(status = 0L, save = "no")
}
result <- run_topology_view_ph(
  view, max_dim = identity$max_dim, threshold = identity$threshold,
  field = identity$field
)
record <- mv05d2_new_ph_record_v1(identity, view, result)
partial <- tempfile(pattern = paste0(basename(result_path), "."),
                    tmpdir = dirname(result_path))
saveRDS(record, partial, compress = FALSE, version = 3)
if (file.exists(result_path)) {
  unlink(partial)
  stop("Refusing to overwrite an MV5-D2 PH cache.", call. = FALSE)
}
if (!file.rename(partial, result_path)) {
  unlink(partial)
  stop("Failed to atomically publish an MV5-D2 PH cache.", call. = FALSE)
}
message("Built atomic MV5-D2 PH cache: ", job$job_id)

