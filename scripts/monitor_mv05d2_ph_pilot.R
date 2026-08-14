#!/usr/bin/env Rscript

options(warn = 2)
for (package in c("processx", "ps", "digest")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(package, " is required for monitored MV5-D2 execution.",
         call. = FALSE)
  }
}

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% 7:8) {
  stop(
    "usage: monitor_mv05d2_ph_pilot.R PILOT_CSV FOLD_CACHE_DIR ",
    "RESULT_DIR FAILURE_DIR METRICS_CSV TIMEOUT_SECONDS RSS_CAP_BYTES ",
    "[all|repeat]",
    call. = FALSE
  )
}
manifest_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
fold_cache_dir <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
result_dir <- args[[3L]]
failure_dir <- args[[4L]]
metrics_path <- args[[5L]]
timeout_seconds <- as.numeric(args[[6L]])
rss_cap_bytes <- as.numeric(args[[7L]])
execution_mode <- if (length(args) == 8L) args[[8L]] else "all"
if (!execution_mode %in% c("all", "repeat")) {
  stop("MV5-D2 execution mode must be all or repeat.", call. = FALSE)
}
stage_worker_cap_seconds <- 3600
if (!is.finite(timeout_seconds) || timeout_seconds <= 0 ||
    !is.finite(rss_cap_bytes) || rss_cap_bytes <= 0) {
  stop("MV5-D2 timeout and RSS cap must be positive.", call. = FALSE)
}
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(failure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(metrics_path), recursive = TRUE, showWarnings = FALSE)

source("R/provenance_utils.R")
source("R/mv05d2_ph_profiling.R")
jobs <- utils::read.csv(
  manifest_path, stringsAsFactors = FALSE, check.names = FALSE
)
jobs$manifest_row_index <- seq_len(nrow(jobs))
required <- c(
  "job_id", "fold_id", "seed", "sample_id", "execution_role",
  "missing_feature_count", "mapping_stratum", "repeat_required",
  "outcome_label_state", "biological_outcomes_computed"
)
if (!all(required %in% names(jobs)) || nrow(jobs) == 0L ||
    anyDuplicated(jobs$job_id) || any(jobs$outcome_label_state != "closed") ||
    any(as.logical(jobs$biological_outcomes_computed))) {
  stop("MV5-D2 pilot manifest is invalid.", call. = FALSE)
}
if (execution_mode == "repeat") {
  jobs <- jobs[as.logical(jobs$repeat_required), , drop = FALSE]
  if (nrow(jobs) != 5L || length(unique(jobs$seed)) != 5L ||
      length(unique(jobs$execution_role)) != 2L) {
    stop("MV5-D2 repeat subset violates its frozen design.", call. = FALSE)
  }
}

process_tree_rss <- function(pid) {
  root <- tryCatch(ps::ps_handle(pid), error = function(error) NULL)
  if (is.null(root)) return(NA_real_)
  handles <- c(list(root), tryCatch(
    ps::ps_children(root, recursive = TRUE), error = function(error) list()
  ))
  sum(vapply(handles, function(handle) {
    tryCatch(as.numeric(ps::ps_memory_info(handle)[["rss"]]),
             error = function(error) 0)
  }, numeric(1L)))
}

existing <- if (file.exists(metrics_path)) {
  utils::read.csv(metrics_path, stringsAsFactors = FALSE, check.names = FALSE)
} else {
  NULL
}
if (!is.null(existing) &&
    (anyDuplicated(existing$job_id) ||
     any(!existing$job_id %in% jobs$job_id) ||
     any(existing$disposition != "completed"))) {
  stop("Existing MV5-D2 metrics are not safely resumable.", call. = FALSE)
}
rows <- if (is.null(existing)) list() else split(existing, seq_len(nrow(existing)))
completed_ids <- if (is.null(existing)) character() else existing$job_id
stage_started <- Sys.time()
for (index in seq_len(nrow(jobs))) {
  job <- jobs[index, , drop = FALSE]
  if (job$job_id %in% completed_ids) next
  stem <- gsub("[^A-Za-z0-9_.-]", "_", job$job_id)
  result_path <- file.path(result_dir, paste0(stem, ".rds"))
  stdout_path <- file.path(failure_dir, paste0(stem, "__stdout.txt"))
  stderr_path <- file.path(failure_dir, paste0(stem, "__stderr.txt"))
  if (file.exists(result_path)) {
    stop("Result exists without a completed metric; refusing ambiguous resume: ",
         result_path, call. = FALSE)
  }
  started <- Sys.time()
  process <- processx::process$new(
    command = Sys.which("Rscript"),
    args = c(
      "--vanilla", "scripts/run_mv05d2_ph_entry.R", manifest_path,
      as.character(job$manifest_row_index), fold_cache_dir, result_path
    ),
    stdout = "|", stderr = "|", cleanup_tree = TRUE
  )
  peak_rss <- 0
  disposition <- "running"
  while (process$is_alive()) {
    elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    rss <- process_tree_rss(process$get_pid())
    if (is.finite(rss)) peak_rss <- max(peak_rss, rss)
    stage_elapsed <- as.numeric(difftime(
      Sys.time(), stage_started, units = "secs"
    ))
    if (elapsed > timeout_seconds || stage_elapsed > stage_worker_cap_seconds) {
      disposition <- "timeout_guard"
      process$kill_tree()
      break
    }
    if (peak_rss > rss_cap_bytes) {
      disposition <- "rss_guard"
      process$kill_tree()
      break
    }
    Sys.sleep(0.1)
  }
  process$wait(timeout = 5000)
  stdout <- paste(process$read_all_output_lines(), collapse = "\n")
  stderr <- paste(process$read_all_error_lines(), collapse = "\n")
  exit_status <- process$get_exit_status()
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  if (identical(disposition, "running")) {
    disposition <- if (identical(exit_status, 0L) && file.exists(result_path)) {
      "completed"
    } else {
      "failed"
    }
  }
  if (nzchar(stdout)) writeLines(stdout, stdout_path, useBytes = TRUE)
  if (nzchar(stderr)) writeLines(stderr, stderr_path, useBytes = TRUE)
  record <- if (identical(disposition, "completed")) readRDS(result_path) else NULL
  if (!is.null(record)) mv05d2_validate_ph_record_v1(record)
  diagram <- if (is.null(record)) NULL else record$topology_result$diagram
  rows[[length(rows) + 1L]] <- data.frame(
    contract_id = "mv05d2_cell_ph_resource_metric_v1",
    job_id = job$job_id, fold_id = job$fold_id, seed = job$seed,
    sample_id = job$sample_id, execution_role = job$execution_role,
    missing_feature_count = job$missing_feature_count,
    mapping_stratum = job$mapping_stratum,
    repeat_required = job$repeat_required,
    disposition = disposition,
    exit_status = if (is.null(exit_status)) NA_integer_ else exit_status,
    elapsed_seconds = elapsed,
    peak_process_tree_rss_bytes = peak_rss,
    timeout_seconds = timeout_seconds,
    rss_cap_bytes = rss_cap_bytes,
    stage_worker_cap_seconds = stage_worker_cap_seconds,
    maximum_heavy_workers = 1L,
    h0_intervals = if (is.null(diagram)) NA_integer_ else
      sum(diagram[, "dimension"] == 0),
    h1_intervals = if (is.null(diagram)) NA_integer_ else
      sum(diagram[, "dimension"] == 1),
    finite_intervals = if (is.null(record)) NA_integer_ else
      record$topology_result$provenance$finite_interval_count,
    infinite_intervals = if (is.null(record)) NA_integer_ else
      record$topology_result$provenance$infinite_interval_count,
    h0_mst_maximum_absolute_error = if (is.null(record)) NA_real_ else
      record$h0_mst_oracle$maximum_absolute_error,
    h0_mst_tolerance = if (is.null(record)) NA_real_ else
      record$h0_mst_oracle$tolerance,
    h0_mst_oracle_passed = if (is.null(record)) FALSE else
      record$h0_mst_oracle$passed,
    diagram_sha256 = if (is.null(record)) NA_character_ else
      record$topology_result$provenance$diagram_sha256,
    record_cache_key = if (is.null(record)) NA_character_ else
      record$cache_key,
    result_file = if (file.exists(result_path)) basename(result_path) else NA,
    result_size_bytes = if (file.exists(result_path))
      unname(file.info(result_path)$size) else NA_real_,
    result_file_sha256 = if (file.exists(result_path))
      digest::digest(file = result_path, algo = "sha256", serialize = FALSE)
      else NA_character_,
    landscape_jobs_executed = 0L, distance_jobs_executed = 0L,
    clustering_jobs_executed = 0L, integration_jobs_executed = 0L,
    gene_view_jobs_executed = 0L, biological_outcomes_computed = FALSE,
    outcome_label_state = "closed", stringsAsFactors = FALSE
  )
  metrics <- do.call(rbind, rows)
  write_provenance_csv(metrics, metrics_path)
  if (!identical(disposition, "completed")) quit(status = 2L, save = "no")
}
metrics <- do.call(rbind, rows)
mv05d2_validate_resource_metrics_v1(
  metrics, expected_jobs = nrow(jobs), timeout_seconds = timeout_seconds,
  rss_cap_bytes = rss_cap_bytes,
  stage_worker_cap_seconds = stage_worker_cap_seconds
)
message("Completed ", nrow(metrics), " bounded MV5-D2 PH jobs.")
