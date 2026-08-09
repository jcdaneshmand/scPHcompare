#!/usr/bin/env Rscript

Sys.setenv(
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1"
)
options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop(
    "usage: run_mv05d3_ph_group.R FULL_MANIFEST_CSV GROUP_ID ",
    "FOLD_CACHE_DIR RESULT_DIR VIEW_AUDIT_CSV RUN_MODE",
    call. = FALSE
  )
}
manifest_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
group_id <- args[[2L]]
fold_cache_dir <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
result_dir <- args[[4L]]
audit_path <- args[[5L]]
run_mode <- match.arg(args[[6L]], c("build_or_resume", "validate_resume"))
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(audit_path), recursive = TRUE, showWarnings = FALSE)

source("R/provenance_utils.R")
source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/mv03_pilot.R")
source("R/mv05_resource_safe_execution.R")
source("R/mv05d2_ph_profiling.R")
source("R/mv05d3_ph_production.R")
file_sha <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}
implementation_files <- c(
  "R/toy_baseline.R", "R/dual_view_topology.R", "R/mv03_pilot.R",
  "R/mv05_resource_safe_execution.R", "R/mv05d2_ph_profiling.R",
  "R/mv05d3_ph_production.R", "scripts/run_mv05d3_ph_group.R"
)
implementation_sha <- digest::digest(
  stats::setNames(vapply(implementation_files, file_sha, character(1L)),
                  implementation_files),
  algo = "sha256", serialize = TRUE
)
manifest <- utils::read.csv(
  manifest_path, stringsAsFactors = FALSE, check.names = FALSE
)
mv05d3_validate_full_manifest_v1(manifest)
jobs <- manifest[manifest$group_id == group_id, , drop = FALSE]
jobs <- jobs[order(jobs$view_order), , drop = FALSE]
if (nrow(jobs) != 90L || length(unique(jobs$fold_cache_file)) != 1L) {
  stop("Requested MV5-D3 group is absent or incomplete.", call. = FALSE)
}
fold_cache_path <- file.path(fold_cache_dir, jobs$fold_cache_file[[1L]])
if (!file.exists(fold_cache_path) || file_sha(fold_cache_path) !=
    jobs$fold_cache_sha256[[1L]]) {
  stop("MV5-D3 source fold cache is missing or stale.", call. = FALSE)
}
fold_record <- readRDS(fold_cache_path)
mv05d1_validate_cell_fold_record_v1(fold_record)
if (fold_record$cache_key != jobs$fold_cache_key[[1L]]) {
  stop("MV5-D3 source fold-cache identity is stale.", call. = FALSE)
}
existing <- if (file.exists(audit_path)) {
  utils::read.csv(audit_path, stringsAsFactors = FALSE, check.names = FALSE)
} else {
  NULL
}
if (!is.null(existing) &&
    (anyDuplicated(existing$job_id) ||
     any(!existing$job_id %in% jobs$job_id) ||
     any(existing$disposition != "built_atomic"))) {
  stop("Existing MV5-D3 view checkpoint is not safely resumable.",
       call. = FALSE)
}
rows <- if (is.null(existing)) list() else split(existing, seq_len(nrow(existing)))
for (index in seq_len(nrow(jobs))) {
  job <- jobs[index, , drop = FALSE]
  view <- fold_record$payload$cell_views[[job$sample_id]]
  validate_topology_view(view)
  if (is.null(view) || view$cache_key != job$view_cache_key ||
      view$payload_sha256 != job$view_payload_sha256) {
    stop("MV5-D3 source view differs from the full manifest.", call. = FALSE)
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
  result_file <- paste0(
    sprintf("view_%03d", job$view_order), "__",
    gsub("[^A-Za-z0-9_.-]", "_", job$sample_id), ".rds"
  )
  result_path <- file.path(result_dir, result_file)
  prior <- if (is.null(existing)) existing else
    existing[existing$job_id == job$job_id, , drop = FALSE]
  if (file.exists(result_path)) {
    if (is.null(prior) || nrow(prior) != 1L) {
      stop("MV5-D3 result exists without one checkpoint; refusing ambiguity.",
           call. = FALSE)
    }
    record <- readRDS(result_path)
    mv05d3_validate_record_static_v1(record)
    if (record$identity$cache_key != identity$cache_key ||
        file_sha(result_path) != prior$result_file_sha256) {
      stop("Existing MV5-D3 result or checkpoint is stale.", call. = FALSE)
    }
    next
  }
  if (!is.null(prior) && nrow(prior) > 0L) {
    stop("MV5-D3 checkpoint exists without its result; refusing rebuild.",
         call. = FALSE)
  }
  if (run_mode == "validate_resume") {
    stop("Resume validation found a missing production result.", call. = FALSE)
  }
  started <- proc.time()[["elapsed"]]
  result <- run_topology_view_ph(
    view, max_dim = 1L, threshold = -1, field = 2L
  )
  record <- mv05d2_new_ph_record_v1(identity, view, result)
  mv05d3_validate_record_static_v1(record)
  partial <- tempfile(pattern = paste0(result_file, "."), tmpdir = result_dir)
  saveRDS(record, partial, compress = FALSE, version = 3)
  if (file.exists(result_path) || !file.rename(partial, result_path)) {
    unlink(partial)
    stop("Failed to atomically publish an MV5-D3 result.", call. = FALSE)
  }
  operation_seconds <- proc.time()[["elapsed"]] - started
  diagram <- record$topology_result$diagram
  rows[[length(rows) + 1L]] <- data.frame(
    contract_id = "mv05d3_cell_ph_view_metric_v1",
    job_id = job$job_id, group_id = job$group_id,
    group_order = job$group_order, view_order = job$view_order,
    fold_id = job$fold_id, seed = job$seed, sample_id = job$sample_id,
    execution_role = job$execution_role,
    missing_feature_count = job$missing_feature_count,
    mapping_stratum = job$mapping_stratum,
    disposition = "built_atomic", operation_seconds = operation_seconds,
    h0_intervals = sum(diagram[, "dimension"] == 0),
    h1_intervals = sum(diagram[, "dimension"] == 1),
    finite_intervals = record$topology_result$provenance$finite_interval_count,
    infinite_intervals =
      record$topology_result$provenance$infinite_interval_count,
    h0_mst_maximum_absolute_error =
      record$h0_mst_oracle$maximum_absolute_error,
    h0_mst_tolerance = record$h0_mst_oracle$tolerance,
    h0_mst_oracle_passed = record$h0_mst_oracle$passed,
    diagram_sha256 = record$topology_result$provenance$diagram_sha256,
    record_cache_key = record$cache_key,
    implementation_sha256 = identity$implementation_sha256,
    result_file = result_file,
    result_size_bytes = unname(file.info(result_path)$size),
    result_file_sha256 = file_sha(result_path),
    landscape_jobs_executed = 0L, distance_jobs_executed = 0L,
    clustering_jobs_executed = 0L, integration_jobs_executed = 0L,
    gene_view_jobs_executed = 0L, biological_outcomes_computed = FALSE,
    outcome_label_state = "closed", stringsAsFactors = FALSE
  )
  checkpoint <- do.call(rbind, rows)
  checkpoint <- checkpoint[order(checkpoint$view_order), , drop = FALSE]
  write_provenance_csv(checkpoint, audit_path)
}
final <- utils::read.csv(
  audit_path, stringsAsFactors = FALSE, check.names = FALSE
)
mv05d3_validate_view_metrics_v1(final, expected_jobs = 90L)
message("Completed MV5-D3 group: ", group_id, " (", run_mode, ")")
