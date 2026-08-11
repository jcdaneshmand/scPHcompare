#!/usr/bin/env Rscript

options(warn = 2)
for (package in c("processx", "ps", "digest")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(package, " is required for monitored MV5-X execution.", call. = FALSE)
  }
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 13L) {
  stop(paste(
    "usage: monitor_mv05x_configuration.R EXECUTION_QUEUE PRIVATE_INVENTORY",
    "D1_ROOT G_ROOT RESULT_ROOT LOG_ROOT RESOURCE_OUTPUT EXPECTED_HEAD",
    "IMPLEMENTATION_SHA256 PYTHON_EXECUTABLE PYTHON_SCRIPT_SHA256",
    "RUN_MODE MAX_UNITS"
  ), call. = FALSE)
}

source("R/provenance_utils.R")
source("R/mv05v_robustness_prefreeze.R")
source("R/mv05x_configuration_execution.R")
queue_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
inventory_path <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
d1_root <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
g_root <- normalizePath(args[[4L]], winslash = "/", mustWork = TRUE)
result_root <- args[[5L]]
log_root <- args[[6L]]
resource_output <- args[[7L]]
expected_head <- args[[8L]]
implementation_sha <- args[[9L]]
python_executable <- normalizePath(args[[10L]], winslash = "/", mustWork = TRUE)
python_script_sha <- args[[11L]]
run_mode <- match.arg(args[[12L]], c("build_or_resume", "validate_resume"))
max_units <- as.integer(args[[13L]])
if (is.na(max_units) || max_units != 150L) {
  stop("MV5-X configuration monitor requires exactly 150 groups.", call. = FALSE)
}
head <- trimws(system2("git", c("rev-parse", "HEAD"), stdout = TRUE))
if (!identical(head, expected_head)) {
  stop("MV5-X monitor is not at its bound prospective HEAD.", call. = FALSE)
}

queue <- utils::read.csv(queue_path, stringsAsFactors = FALSE,
                         check.names = FALSE)
mv05x_validate_configuration_queue_v1(queue)
queue <- queue[as.logical(queue$execution_authorized), , drop = FALSE]
queue <- queue[order(queue$configuration_execution_order), , drop = FALSE]
if (nrow(queue) != max_units ||
    any(queue$implementation_sha256 != implementation_sha)) {
  stop("MV5-X bound configuration identity is stale.", call. = FALSE)
}
dir.create(result_root, recursive = TRUE, showWarnings = FALSE)
dir.create(log_root, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(resource_output), recursive = TRUE, showWarnings = FALSE)

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
tree_bytes <- function(path) {
  files <- list.files(path, recursive = TRUE, full.names = TRUE,
                      all.files = TRUE, no.. = TRUE)
  if (!length(files)) return(0)
  info <- file.info(files)
  sum(info$size[!info$isdir], na.rm = TRUE)
}
safe_name <- function(value) gsub("[^A-Za-z0-9_.-]", "_", value)
write_checkpoint <- function(value) {
  temporary <- tempfile(pattern = ".mv05x_resources_",
                        tmpdir = dirname(resource_output))
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  write_provenance_csv(value, temporary)
  if (!file.rename(temporary, resource_output)) {
    stop("MV5-X resource checkpoint publication failed.", call. = FALSE)
  }
}

required_metrics <- c(
  "robustness_group_id", "configuration_execution_order", "disposition",
  "elapsed_seconds", "labels_opened", "outcomes_computed"
)
metrics <- NULL
if (file.exists(resource_output)) {
  metrics <- utils::read.csv(resource_output, stringsAsFactors = FALSE,
                             check.names = FALSE)
  expected_prefix <- queue$robustness_group_id[seq_len(nrow(metrics))]
  if (!all(required_metrics %in% names(metrics)) ||
      nrow(metrics) > nrow(queue) ||
      !identical(metrics$robustness_group_id, expected_prefix) ||
      any(!metrics$disposition %in% c("completed", "reused_validated")) ||
      any(as.logical(metrics$labels_opened)) ||
      any(as.logical(metrics$outcomes_computed))) {
    stop("MV5-X existing resource checkpoint is invalid.", call. = FALSE)
  }
}
completed_count <- if (is.null(metrics)) 0L else nrow(metrics)
prior_elapsed <- if (is.null(metrics)) 0 else sum(metrics$elapsed_seconds)
stage_started <- Sys.time()
failed <- FALSE

if (completed_count < nrow(queue)) for (index in seq.int(completed_count + 1L,
                                                         nrow(queue))) {
  unit <- queue[index, , drop = FALSE]
  stem <- safe_name(unit$robustness_group_id)
  final_dir <- file.path(result_root, stem)
  preexisting <- dir.exists(final_dir)
  started <- Sys.time()
  process <- processx::process$new(
    command = Sys.which("Rscript"),
    args = c(
      "--vanilla", "scripts/run_mv05w_full_group.R", queue_path,
      inventory_path, d1_root, g_root, result_root,
      unit$robustness_group_id, expected_head, implementation_sha,
      python_executable, python_script_sha, run_mode
    ),
    stdout = "|", stderr = "|", cleanup_tree = TRUE
  )
  peak_rss <- 0
  guard <- "running"
  while (process$is_alive()) {
    Sys.sleep(0.25)
    elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    stage_elapsed <- prior_elapsed + as.numeric(difftime(
      Sys.time(), stage_started, units = "secs"
    ))
    rss <- process_tree_rss(process$get_pid())
    if (is.finite(rss)) peak_rss <- max(peak_rss, rss)
    storage <- tree_bytes(result_root)
    if (elapsed > 600 || stage_elapsed > 30 * 60 * 60) {
      guard <- "timeout_guard"
      process$kill_tree()
    } else if (peak_rss > 4 * 1024^3) {
      guard <- "rss_guard"
      process$kill_tree()
    } else if (storage > 16 * 1024^3) {
      guard <- "storage_guard"
      process$kill_tree()
    }
  }
  process$wait(timeout = 5000)
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  stdout <- paste(process$read_all_output_lines(), collapse = "\n")
  stderr <- paste(process$read_all_error_lines(), collapse = "\n")
  if (nzchar(stdout)) writeLines(
    stdout, file.path(log_root, paste0(stem, "__stdout.txt")), useBytes = TRUE
  )
  if (nzchar(stderr)) writeLines(
    stderr, file.path(log_root, paste0(stem, "__stderr.txt")), useBytes = TRUE
  )
  exit_status <- process$get_exit_status()
  status_path <- file.path(final_dir, "status.csv")
  status <- if (file.exists(status_path)) utils::read.csv(
    status_path, stringsAsFactors = FALSE, check.names = FALSE
  ) else NULL
  disposition <- guard
  if (guard == "running") {
    disposition <- if (identical(exit_status, 0L) && !is.null(status) &&
                       nrow(status) == 1L && status$status == "completed") {
      if (preexisting) "reused_validated" else "completed"
    } else "failed"
  }
  row <- data.frame(
    contract_id = "mv05x_pc20_resource_metric_v1",
    robustness_group_id = unit$robustness_group_id,
    configuration_execution_order = unit$configuration_execution_order,
    execution_order = unit$execution_order, fold_id = unit$fold_id,
    seed = unit$seed, representation = unit$representation,
    configuration_id = unit$configuration_id, disposition = disposition,
    exit_status = if (is.null(exit_status)) NA_integer_ else exit_status,
    completed_views = if (is.null(status)) 0L else status$completed_views,
    landscape_summary_rows = if (is.null(status)) 0L else
      status$landscape_summary_rows,
    landscape_pair_rows = if (is.null(status)) 0L else
      status$landscape_pair_rows,
    energy_pair_rows = if (is.null(status)) 0L else status$energy_pair_rows,
    method_rows = if (is.null(status)) 0L else status$method_rows,
    elapsed_seconds = elapsed,
    peak_process_tree_rss_bytes = peak_rss,
    cumulative_private_bytes = tree_bytes(result_root),
    unit_elapsed_cap_seconds = 600,
    rss_cap_bytes = 4 * 1024^3,
    program_elapsed_cap_seconds = 30 * 60 * 60,
    storage_cap_bytes = 16 * 1024^3, maximum_workers = 1L,
    labels_opened = FALSE, outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
  metrics <- if (is.null(metrics)) row else rbind(metrics, row)
  write_checkpoint(metrics)
  if (disposition %in% c("failed", "timeout_guard", "rss_guard",
                          "storage_guard")) {
    failed <- TRUE
    break
  }
}

if (failed || is.null(metrics) || nrow(metrics) != nrow(queue)) {
  quit(status = 2L, save = "no")
}
if (any(!metrics$disposition %in% c("completed", "reused_validated")) ||
    any(metrics$completed_views != 90L) ||
    any(metrics$landscape_summary_rows != 180L) ||
    any(metrics$elapsed_seconds > 600) ||
    sum(metrics$elapsed_seconds) > 30 * 60 * 60 ||
    any(metrics$peak_process_tree_rss_bytes > 4 * 1024^3) ||
    max(metrics$cumulative_private_bytes) > 16 * 1024^3 ||
    any(metrics$landscape_pair_rows != queue$landscape_request_rows) ||
    any(metrics$energy_pair_rows != queue$energy_request_rows) ||
    any(metrics$method_rows != queue$assembled_method_rows) ||
    any(as.logical(metrics$labels_opened)) ||
    any(as.logical(metrics$outcomes_computed))) {
  stop("MV5-X PC20 execution completed outside its contract.", call. = FALSE)
}
message("Completed 150 MV5-X PC20 groups (", run_mode, ").")
