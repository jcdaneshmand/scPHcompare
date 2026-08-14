#!/usr/bin/env Rscript

options(warn = 2)

for (package in c("processx", "ps", "digest")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(package, " is required for monitored MV-03 PH execution.",
         call. = FALSE)
  }
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) {
  stop(
    "Usage: monitor_mv03_ph_stage.R <job-manifest-csv> <result-dir> ",
    "<failure-dir> <metrics-csv> <timeout-seconds> <rss-cap-bytes> ",
    "<stage-worker-cap-seconds>",
    call. = FALSE
  )
}

manifest_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
result_dir <- args[[2L]]
failure_dir <- args[[3L]]
metrics_path <- args[[4L]]
timeout_seconds <- as.numeric(args[[5L]])
rss_cap_bytes <- as.numeric(args[[6L]])
stage_worker_cap_seconds <- as.numeric(args[[7L]])
if (!is.finite(timeout_seconds) || timeout_seconds <= 0 ||
    !is.finite(rss_cap_bytes) || rss_cap_bytes <= 0 ||
    !is.finite(stage_worker_cap_seconds) || stage_worker_cap_seconds <= 0) {
  stop("Timeout, RSS cap, and stage cap must be positive finite values.",
       call. = FALSE)
}
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(failure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(metrics_path), recursive = TRUE, showWarnings = FALSE)

jobs <- utils::read.csv(
  manifest_path, stringsAsFactors = FALSE, check.names = FALSE
)
required <- c("stage", "cohort", "representation", "sample_id", "seed",
              "view_id", "point_count", "coordinate_count", "view_cache_key",
              "view_rds")
if (!all(required %in% names(jobs)) || nrow(jobs) == 0L) {
  stop("Job manifest has an incompatible schema or no jobs.", call. = FALSE)
}

process_tree_rss <- function(pid) {
  root <- tryCatch(ps::ps_handle(pid), error = function(e) NULL)
  if (is.null(root)) return(NA_real_)
  handles <- c(list(root), tryCatch(
    ps::ps_children(root, recursive = TRUE), error = function(e) list()
  ))
  values <- vapply(handles, function(handle) {
    tryCatch(as.numeric(ps::ps_memory_info(handle)[["rss"]]),
             error = function(e) 0)
  }, numeric(1L))
  sum(values)
}

rows <- vector("list", nrow(jobs))
stage_started <- Sys.time()
for (index in seq_len(nrow(jobs))) {
  job <- jobs[index, , drop = FALSE]
  stem <- paste(job$stage, job$cohort, job$representation, job$sample_id,
                job$seed, job$view_id, sep = "__")
  result_path <- file.path(result_dir, paste0(stem, ".rds"))
  stdout_path <- file.path(failure_dir, paste0(stem, "__stdout.txt"))
  stderr_path <- file.path(failure_dir, paste0(stem, "__stderr.txt"))
  if (file.exists(result_path)) {
    stop("Refusing to overwrite existing result: ", result_path,
         call. = FALSE)
  }
  started <- Sys.time()
  process <- processx::process$new(
    command = Sys.which("Rscript"),
    args = c("--vanilla", "scripts/run_mv03_ph_job.R", job$view_rds,
             result_path),
    stdout = "|", stderr = "|", cleanup_tree = TRUE
  )
  peak_rss <- 0
  disposition <- "running"
  while (process$is_alive()) {
    elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    rss <- process_tree_rss(process$get_pid())
    if (is.finite(rss)) peak_rss <- max(peak_rss, rss)
    if (elapsed > timeout_seconds) {
      disposition <- "timeout"
      process$kill_tree()
      break
    }
    if (peak_rss > rss_cap_bytes) {
      disposition <- "rss_cap_exceeded"
      process$kill_tree()
      break
    }
    Sys.sleep(0.25)
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
  result <- if (identical(disposition, "completed")) readRDS(result_path) else NULL
  rows[[index]] <- data.frame(
    stage = job$stage,
    cohort = job$cohort,
    representation = job$representation,
    sample_id = job$sample_id,
    seed = job$seed,
    view_id = job$view_id,
    point_count = job$point_count,
    coordinate_count = job$coordinate_count,
    view_cache_key = job$view_cache_key,
    disposition = disposition,
    exit_status = if (is.null(exit_status)) NA_integer_ else exit_status,
    elapsed_seconds = elapsed,
    peak_process_tree_rss_bytes = peak_rss,
    timeout_seconds = timeout_seconds,
    rss_cap_bytes = rss_cap_bytes,
    stage_worker_cap_seconds = stage_worker_cap_seconds,
    max_dim = 1L,
    threshold = -1,
    field = 2L,
    h0_intervals = if (is.null(result)) NA_integer_ else
      sum(result$diagram[, "dimension"] == 0),
    h1_intervals = if (is.null(result)) NA_integer_ else
      sum(result$diagram[, "dimension"] == 1),
    finite_intervals = if (is.null(result)) NA_integer_ else
      result$provenance$finite_interval_count,
    infinite_intervals = if (is.null(result)) NA_integer_ else
      result$provenance$infinite_interval_count,
    diagram_sha256 = if (is.null(result)) NA_character_ else
      result$provenance$diagram_sha256,
    result_file = if (file.exists(result_path))
      gsub("\\\\", "/", result_path) else NA_character_,
    result_file_sha256 = if (file.exists(result_path))
      digest::digest(file = result_path, algo = "sha256", serialize = FALSE)
      else NA_character_,
    stdout_file = if (file.exists(stdout_path))
      gsub("\\\\", "/", stdout_path) else NA_character_,
    stderr_file = if (file.exists(stderr_path))
      gsub("\\\\", "/", stderr_path) else NA_character_,
    stringsAsFactors = FALSE
  )
  utils::write.csv(do.call(rbind, rows[seq_len(index)]), metrics_path,
                   row.names = FALSE)
  if (!identical(disposition, "completed")) break
}
metrics <- do.call(rbind, Filter(Negate(is.null), rows))
stage_elapsed <- as.numeric(difftime(Sys.time(), stage_started, units = "secs"))
if (!nrow(metrics) || any(metrics$disposition != "completed") ||
    sum(metrics$elapsed_seconds) > stage_worker_cap_seconds ||
    stage_elapsed > stage_worker_cap_seconds) {
  quit(status = 2L, save = "no")
}
message("Completed ", nrow(metrics), " Stage ",
        paste(unique(metrics$stage), collapse = ";"),
        " jobs within frozen caps.")
