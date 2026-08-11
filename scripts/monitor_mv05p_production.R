#!/usr/bin/env Rscript

options(warn = 2)
for (package in c("processx", "ps", "digest")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(package, " is required for monitored MV5-P production.",
         call. = FALSE)
  }
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("usage: monitor_mv05p_production.R OUTPUT_ROOT GROUP_METRICS",
       call. = FALSE)
}
source("R/provenance_utils.R")
source("R/mv05p_distance_production.R")
output_root <- args[[1L]]
metrics_path <- args[[2L]]
group_queue_path <- "docs/audits/mv05o-production-group-queue-2026-08-10.csv"
source_freeze_path <- "docs/audits/mv05o-source-freeze-2026-08-10.csv"
combined_projection_path <-
  "docs/audits/mv05n-combined-resource-projection-2026-08-10.csv"
landscape_projection_path <-
  "docs/audits/mv05n-full-matrix-resource-projection-2026-08-10.csv"
python <- "tmp/mv05p/python-runtime/bin/python"
if (!file.exists(python)) {
  stop("MV5-P Python executable is missing.", call. = FALSE)
}
runtime_freeze <- normalizePath(
  "tmp/mv05p/python-runtime-freeze.txt", mustWork = TRUE)
persim_files <- normalizePath(
  "tmp/mv05p/persim-installed-files.csv", mustWork = TRUE)
file_sha <- function(path) digest::digest(file = path, algo = "sha256",
                                          serialize = FALSE)
helper_sha <-
  "f97ea2cc64e8eff0b92070ff563a21d75524ae5b3aff2a909ed42eec80422146"
worker_sha <-
  "013f9194677e6b02fcdd96aca9f07d1db6a65d956215155d73fd2042c260cfdf"
read_public <- function(path) utils::read.csv(
  normalizePath(path, mustWork = TRUE), stringsAsFactors = FALSE,
  check.names = FALSE
)
write_atomic <- function(value, path) {
  temporary <- paste0(path, ".partial.", Sys.getpid())
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  write_provenance_csv(value, temporary)
  if (!file.rename(temporary, path)) {
    stop("MV5-P metrics atomic rename failed.", call. = FALSE)
  }
}
directory_bytes <- function(path) {
  files <- list.files(path, recursive = TRUE, full.names = TRUE,
                      all.files = TRUE, no.. = TRUE)
  files <- files[file.exists(files) & !dir.exists(files)]
  if (!length(files)) return(0)
  sum(file.info(files)$size, na.rm = TRUE)
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

if (file_sha("R/mv05p_distance_production.R") != helper_sha ||
    file_sha("scripts/run_mv05p_production_group.R") != worker_sha) {
  stop("MV5-P orchestration source drifted.", call. = FALSE)
}
if (file_sha(runtime_freeze) !=
      "470355c653dffd102f6a0abfd7dcb1978f2b696e8df870a6171af9855a95f1f6" ||
    file_sha(persim_files) !=
      "8f7ee101c08c4f50b044d102553c3a3a1375d4c8f60a1acdbd52cf2e00050988") {
  stop("MV5-P Python runtime provenance drifted.", call. = FALSE)
}
runtime_check <- processx::run(
  python, c("-c", "import persim; assert persim.__version__ == '0.3.8'"),
  error_on_status = FALSE
)
if (!identical(runtime_check$status, 0L)) {
  stop("MV5-P frozen Persim 0.3.8 runtime is unavailable.", call. = FALSE)
}

source_freeze <- read_public(source_freeze_path)
if (nrow(source_freeze) != 18L ||
    any(!file.exists(source_freeze$artifact_path)) ||
    any(vapply(source_freeze$artifact_path, file_sha, character(1L)) !=
          source_freeze$sha256) ||
    length(unique(source_freeze$source_freeze_sha256)) != 1L ||
    unique(source_freeze$source_freeze_sha256) !=
      "541e7d3aa8acce5d512bbb4819c034735eef47387e91a63abccfa259f53d6de1") {
  stop("MV5-P committed source freeze validation failed.", call. = FALSE)
}
groups <- read_public(group_queue_path)
mv05p_validate_group_queue_v1(groups)
groups <- groups[order(groups$execution_order), , drop = FALSE]
projection <- mv05p_group_projection_v1(
  groups, read_public(combined_projection_path),
  read_public(landscape_projection_path)
)
if (abs(sum(projection$projected_worker_seconds) / 3600 -
        16.1170472529584) > 1e-10 ||
    abs(sum(projection$projected_private_bytes) - 1277893355) > 1) {
  stop("MV5-P launch projection differs from the prefreeze.", call. = FALSE)
}
for (path in c(output_root, file.path(output_root, "logs"),
               dirname(metrics_path))) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}
partial <- list.files(output_root, pattern = "[.]partial[.]", recursive = TRUE,
                      full.names = TRUE, all.files = TRUE)
if (length(partial)) {
  stop("MV5-P found a partial artifact; automatic recovery is prohibited.",
       call. = FALSE)
}
existing <- if (file.exists(metrics_path)) read_public(metrics_path) else NULL
if (!is.null(existing) &&
    (anyDuplicated(existing$production_group_id) ||
     any(!existing$production_group_id %in% groups$production_group_id) ||
     any(existing$disposition != "completed") ||
     any(existing$outcome_label_state != "closed") ||
     any(as.logical(existing$biological_outcomes_computed)) ||
     any(as.integer(existing$clustering_jobs_executed) != 0L))) {
  stop("MV5-P existing metrics are not safely resumable.", call. = FALSE)
}
rows <- if (is.null(existing)) list() else split(existing, seq_len(nrow(existing)))
completed <- if (is.null(existing)) character() else existing$production_group_id
pending <- groups[!groups$production_group_id %in% completed, , drop = FALSE]
running <- list()
failed <- FALSE
last_storage_check <- Sys.time() - 60

budget <- function() {
  completed_seconds <- if (is.null(existing) && !length(rows)) 0 else {
    current <- if (length(rows)) do.call(rbind, rows) else existing
    sum(current$elapsed_seconds[current$disposition == "completed"])
  }
  running_elapsed <- if (!length(running)) numeric() else vapply(
    running, function(state) as.numeric(difftime(
      Sys.time(), state$started, units = "secs")), numeric(1L))
  running_ids <- names(running)
  incomplete_ids <- c(running_ids, pending$production_group_id)
  projected <- projection[
    projection$production_group_id %in% incomplete_ids, , drop = FALSE]
  remaining_seconds <- sum(projected$projected_worker_seconds)
  if (length(running_elapsed)) {
    expected_running <- projection$projected_worker_seconds[
      match(running_ids, projection$production_group_id)]
    remaining_seconds <- remaining_seconds - sum(pmin(
      running_elapsed, expected_running
    ))
  }
  completed_bytes <- if (!length(rows)) 0 else sum(vapply(rows, function(row) {
    if (row$disposition == "completed") as.numeric(row$private_artifact_bytes) else 0
  }, numeric(1L)))
  mv05p_launch_budget_v1(
    completed_seconds + sum(running_elapsed), remaining_seconds,
    completed_bytes, sum(projected$projected_private_bytes)
  )
}

launch_group <- function(group) {
  if (file_sha("scripts/run_mv05p_production_group.R") != worker_sha) {
    stop("MV5-P group worker changed during production.", call. = FALSE)
  }
  id <- group$production_group_id
  stem <- gsub("[^A-Za-z0-9_.-]", "_", id)
  process <- processx::process$new(
    command = Sys.which("Rscript"),
    args = c("--vanilla", "scripts/run_mv05p_production_group.R",
             id, output_root, python),
    stdout = "|", stderr = "|", cleanup_tree = TRUE
  )
  list(
    process = process, group = group, started = Sys.time(), peak_rss = 0,
    stdout_path = file.path(output_root, "logs", paste0(stem, "__stdout.txt")),
    stderr_path = file.path(output_root, "logs", paste0(stem, "__stderr.txt")),
    disposition = "running"
  )
}

while (nrow(pending) > 0L || length(running) > 0L) {
  while (!failed && length(running) < 2L && nrow(pending) > 0L) {
    current_budget <- budget()
    if (!current_budget$launch_authorized) {
      failed <- TRUE
      message("MV5-P launch budget guard closed the queue.")
      break
    }
    group <- pending[1L, , drop = FALSE]
    pending <- pending[-1L, , drop = FALSE]
    running[[group$production_group_id]] <- launch_group(group)
  }
  Sys.sleep(0.5)
  now <- Sys.time()
  check_storage <- as.numeric(difftime(now, last_storage_check,
                                      units = "secs")) >= 10
  if (check_storage) {
    last_storage_check <- now
    if (directory_bytes(output_root) > 10737418240) {
      for (id in names(running)) {
        running[[id]]$disposition <- "storage_guard"
        running[[id]]$process$kill_tree()
      }
      failed <- TRUE
    }
  }
  for (id in names(running)) {
    state <- running[[id]]
    process <- state$process
    elapsed <- as.numeric(difftime(Sys.time(), state$started, units = "secs"))
    if (process$is_alive()) {
      rss <- process_tree_rss(process$get_pid())
      if (is.finite(rss)) state$peak_rss <- max(state$peak_rss, rss)
      if (elapsed > 900) {
        state$disposition <- "timeout_guard"
        process$kill_tree()
      } else if (state$peak_rss > 4294967296) {
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
    group_status_path <- file.path(
      output_root, "groups", gsub("[^A-Za-z0-9_.-]", "_", id),
      "group-status.csv"
    )
    disposition <- state$disposition
    if (disposition == "running") {
      disposition <- if (identical(exit_status, 0L) &&
                         file.exists(group_status_path)) "completed" else "failed"
    }
    group_status <- if (disposition == "completed") {
      read_public(group_status_path)
    } else NULL
    if (!is.null(group_status) &&
        (nrow(group_status) != 1L ||
         group_status$production_group_id != id ||
         group_status$outcome_label_state != "closed" ||
         as.logical(group_status$biological_outcomes_computed) ||
         as.integer(group_status$clustering_jobs_executed) != 0L)) {
      disposition <- "failed_validation"
    }
    expected <- projection[projection$production_group_id == id, ]
    rows[[length(rows) + 1L]] <- data.frame(
      contract_id = "mv05p_group_resource_metric_v1",
      production_group_id = id,
      execution_order = state$group$execution_order,
      representation = state$group$representation,
      disposition = disposition,
      exit_status = if (is.null(exit_status)) NA_integer_ else exit_status,
      landscape_chunks = if (is.null(group_status)) 0L else
        group_status$landscape_chunks,
      landscape_rows = if (is.null(group_status)) 0L else
        group_status$landscape_rows,
      baseline_units = if (is.null(group_status)) 0L else
        group_status$baseline_units,
      baseline_rows = if (is.null(group_status)) 0L else
        group_status$baseline_rows,
      elapsed_seconds = elapsed,
      projected_worker_seconds = expected$projected_worker_seconds,
      peak_process_tree_rss_bytes = state$peak_rss,
      peak_rss_provenance = "live_process_tree",
      private_artifact_count = if (is.null(group_status)) 0L else
        group_status$private_artifact_count,
      private_artifact_bytes = if (is.null(group_status)) 0 else
        group_status$private_artifact_bytes,
      projected_private_bytes = expected$projected_private_bytes,
      source_freeze_sha256 = unique(source_freeze$source_freeze_sha256),
      runtime_freeze_sha256 = file_sha(runtime_freeze),
      outcome_label_state = "closed", biological_outcomes_computed = FALSE,
      clustering_jobs_executed = 0L, stringsAsFactors = FALSE
    )
    metrics <- do.call(rbind, rows)
    metrics <- metrics[order(metrics$execution_order), , drop = FALSE]
    write_atomic(metrics, metrics_path)
    if (disposition != "completed") failed <- TRUE
    running[[id]] <- NULL
  }
  if (failed && length(running) == 0L) break
}
metrics <- if (length(rows)) do.call(rbind, rows) else data.frame()
if (failed || nrow(metrics) != 150L || any(metrics$disposition != "completed")) {
  quit(status = 2L, save = "no")
}
final_budget <- budget()
if (sum(metrics$landscape_chunks) != 4340L ||
    sum(metrics$landscape_rows) != 1050700L ||
    sum(metrics$baseline_units) != 225L ||
    sum(metrics$baseline_rows) != 788025L ||
    sum(metrics$elapsed_seconds) / 3600 > 21.6 ||
    max(metrics$peak_process_tree_rss_bytes) > 4294967296 ||
    directory_bytes(output_root) > 10737418240 ||
    any(metrics$outcome_label_state != "closed") ||
    any(as.logical(metrics$biological_outcomes_computed)) ||
    any(as.integer(metrics$clustering_jobs_executed) != 0L)) {
  stop("MV5-P final production accounting failed.", call. = FALSE)
}
message("Completed all 150 frozen MV5-P production groups.")
