#!/usr/bin/env Rscript

options(warn = 2)
for (package in c("processx", "ps", "digest")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(package, " is required for monitored MV5-H execution.", call. = FALSE)
  }
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 11L) {
  stop(
    "usage: monitor_mv05h_ph_groups.R MANIFEST SOURCE_ROOT RESULT_ROOT ",
    "VIEW_AUDIT_ROOT LOG_ROOT GROUP_METRICS MAX_GROUPS MAX_WORKERS ",
    "GROUP_TIMEOUT_SECONDS RSS_CAP_BYTES STAGE_CAP_SECONDS", call. = FALSE
  )
}
manifest_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
source_root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
result_root <- args[[3L]]
audit_root <- args[[4L]]
log_root <- args[[5L]]
metrics_path <- args[[6L]]
max_groups <- as.integer(args[[7L]])
max_workers <- as.integer(args[[8L]])
timeout_seconds <- as.numeric(args[[9L]])
rss_cap_bytes <- as.numeric(args[[10L]])
stage_cap_seconds <- as.numeric(args[[11L]])
if (is.na(max_groups) || max_groups < 1L || max_groups > 75L ||
    is.na(max_workers) || max_workers < 1L || max_workers > 2L ||
    !is.finite(timeout_seconds) || timeout_seconds <= 0 ||
    !is.finite(rss_cap_bytes) || rss_cap_bytes <= 0 ||
    !is.finite(stage_cap_seconds) || stage_cap_seconds <= 0) {
  stop("MV5-H monitor guards are invalid.", call. = FALSE)
}
for (path in c(result_root, audit_root, log_root, dirname(metrics_path))) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}
source("R/provenance_utils.R")
source("R/mv05d2_ph_profiling.R")
source("R/mv05h_integrated_ph_production.R")
manifest <- utils::read.csv(
  manifest_path, stringsAsFactors = FALSE, check.names = FALSE
)
mv05h_validate_manifest_v1(manifest)
groups <- unique(manifest[c("group_id", "group_order", "fold_id", "seed")])
groups <- groups[order(groups$group_order), , drop = FALSE]
groups <- groups[groups$group_order <= max_groups, , drop = FALSE]
existing <- if (file.exists(metrics_path)) {
  utils::read.csv(metrics_path, stringsAsFactors = FALSE, check.names = FALSE)
} else {
  NULL
}
if (!is.null(existing) &&
    (anyDuplicated(existing$group_id) ||
     any(!existing$group_id %in% groups$group_id) ||
     any(existing$disposition != "completed"))) {
  stop("Existing MV5-H group metrics are not safely resumable.", call. = FALSE)
}
rows <- if (is.null(existing)) list() else split(existing, seq_len(nrow(existing)))
completed <- if (is.null(existing)) character() else existing$group_id
pending <- groups[!groups$group_id %in% completed, , drop = FALSE]
running <- list()
failed <- FALSE
stage_started <- Sys.time()

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

launch_group <- function(group) {
  stem <- gsub("[^A-Za-z0-9_.-]", "_", group$group_id)
  result_dir <- file.path(result_root, stem)
  audit_path <- file.path(audit_root, paste0(stem, "__views.csv"))
  process <- processx::process$new(
    command = Sys.which("Rscript"),
    args = c(
      "--vanilla", "scripts/run_mv05h_ph_group.R", manifest_path,
      group$group_id, source_root, result_dir, audit_path, "build_or_resume"
    ),
    stdout = "|", stderr = "|", cleanup_tree = TRUE
  )
  list(
    process = process, group = group, started = Sys.time(), peak_rss = 0,
    result_dir = result_dir, audit_path = audit_path,
    stdout_path = file.path(log_root, paste0(stem, "__stdout.txt")),
    stderr_path = file.path(log_root, paste0(stem, "__stderr.txt")),
    disposition = "running"
  )
}

while (nrow(pending) > 0L || length(running) > 0L) {
  while (!failed && length(running) < max_workers && nrow(pending) > 0L) {
    group <- pending[1L, , drop = FALSE]
    pending <- pending[-1L, , drop = FALSE]
    running[[group$group_id]] <- launch_group(group)
  }
  Sys.sleep(0.25)
  for (id in names(running)) {
    state <- running[[id]]
    process <- state$process
    elapsed <- as.numeric(difftime(Sys.time(), state$started, units = "secs"))
    if (process$is_alive()) {
      rss <- process_tree_rss(process$get_pid())
      if (is.finite(rss)) state$peak_rss <- max(state$peak_rss, rss)
      stage_elapsed <- as.numeric(difftime(
        Sys.time(), stage_started, units = "secs"
      ))
      if (elapsed > timeout_seconds || stage_elapsed > stage_cap_seconds) {
        state$disposition <- "timeout_guard"
        process$kill_tree()
      } else if (state$peak_rss > rss_cap_bytes) {
        state$disposition <- "rss_guard"
        process$kill_tree()
      }
      running[[id]] <- state
      next
    }
    process$wait(timeout = 5000)
    stdout <- paste(process$read_all_output_lines(), collapse = "\n")
    stderr <- paste(process$read_all_error_lines(), collapse = "\n")
    if (nzchar(stdout)) writeLines(stdout, state$stdout_path, useBytes = TRUE)
    if (nzchar(stderr)) writeLines(stderr, state$stderr_path, useBytes = TRUE)
    exit_status <- process$get_exit_status()
    disposition <- state$disposition
    if (disposition == "running") {
      disposition <- if (identical(exit_status, 0L) &&
                         file.exists(state$audit_path)) "completed" else "failed"
    }
    view_metrics <- if (disposition == "completed") {
      utils::read.csv(
        state$audit_path, stringsAsFactors = FALSE, check.names = FALSE
      )
    } else {
      NULL
    }
    if (!is.null(view_metrics)) {
      mv05h_validate_view_metrics_v1(view_metrics, expected_jobs = 90L)
    }
    rows[[length(rows) + 1L]] <- data.frame(
      contract_id = "mv05h_group_resource_metric_v1",
      group_id = state$group$group_id,
      group_order = state$group$group_order,
      fold_id = state$group$fold_id, seed = state$group$seed,
      disposition = disposition,
      exit_status = if (is.null(exit_status)) NA_integer_ else exit_status,
      completed_views = if (is.null(view_metrics)) 0L else nrow(view_metrics),
      elapsed_seconds = elapsed,
      peak_process_tree_rss_bytes = state$peak_rss,
      private_result_bytes = if (is.null(view_metrics)) 0 else
        sum(view_metrics$result_size_bytes),
      group_timeout_seconds = timeout_seconds,
      rss_cap_bytes = rss_cap_bytes, stage_cap_seconds = stage_cap_seconds,
      maximum_heavy_workers = max_workers,
      landscape_jobs_executed = 0L, distance_jobs_executed = 0L,
      retrieval_jobs_executed = 0L, clustering_jobs_executed = 0L,
      gene_view_jobs_executed = 0L, fusion_jobs_executed = 0L,
      new_data_jobs_executed = 0L, biological_outcomes_computed = FALSE,
      outcome_label_state = "closed", stringsAsFactors = FALSE
    )
    metrics <- do.call(rbind, rows)
    metrics <- metrics[order(metrics$group_order), , drop = FALSE]
    write_provenance_csv(metrics, metrics_path)
    if (disposition != "completed") failed <- TRUE
    running[[id]] <- NULL
  }
  if (failed && length(running) == 0L) break
}
metrics <- do.call(rbind, rows)
if (failed || nrow(metrics) != nrow(groups)) quit(status = 2L, save = "no")
mv05h_validate_group_metrics_v1(
  metrics, expected_groups = nrow(groups),
  elapsed_cap_seconds = timeout_seconds, rss_cap_bytes = rss_cap_bytes,
  stage_cap_seconds = stage_cap_seconds
)
message("Completed ", nrow(metrics), " MV5-H integrated PH groups.")
