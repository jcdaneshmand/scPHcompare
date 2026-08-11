#!/usr/bin/env Rscript

options(warn = 2)
for (package in c("processx", "ps", "digest")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(package, " is required for monitored MV5-W smoke.",
         call. = FALSE)
  }
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 13L) {
  stop(paste(
    "usage: monitor_mv05w_smoke.R EXECUTION_QUEUE PRIVATE_INVENTORY",
    "D1_ROOT G_ROOT RESULT_ROOT LOG_ROOT RESOURCE_OUTPUT EXPECTED_HEAD",
    "IMPLEMENTATION_SHA256 PYTHON_EXECUTABLE PYTHON_SCRIPT_SHA256",
    "RUN_MODE MAX_UNITS"
  ), call. = FALSE)
}

source("R/provenance_utils.R")
source("R/mv05t_robustness_gate.R")
source("R/mv05u_robustness_admission.R")
source("R/mv05v_robustness_prefreeze.R")
source("R/mv05w_launch_readiness.R")
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
if (is.na(max_units) || max_units != 1L) {
  stop("MV5-W smoke monitor permits exactly one group.", call. = FALSE)
}

queue <- utils::read.csv(
  queue_path, stringsAsFactors = FALSE, check.names = FALSE
)
mv05w_validate_smoke_queue_v1(queue)
queue <- queue[queue$execution_order <= max_units, , drop = FALSE]
dir.create(result_root, recursive = TRUE, showWarnings = FALSE)
dir.create(log_root, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(resource_output), recursive = TRUE, showWarnings = FALSE)
if (file.exists(resource_output)) {
  stop("MV5-W monitor refuses to overwrite a resource output.",
       call. = FALSE)
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
tree_bytes <- function(path) {
  files <- list.files(path, recursive = TRUE, full.names = TRUE,
                      all.files = TRUE, no.. = TRUE)
  if (!length(files)) return(0)
  info <- file.info(files)
  sum(info$size[!info$isdir], na.rm = TRUE)
}
safe_name <- function(value) gsub("[^A-Za-z0-9_.-]", "_", value)

rows <- vector("list", nrow(queue))
stage_started <- Sys.time()
failed <- FALSE
for (index in seq_len(nrow(queue))) {
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
    stage_elapsed <- as.numeric(difftime(
      Sys.time(), stage_started, units = "secs"
    ))
    rss <- process_tree_rss(process$get_pid())
    if (is.finite(rss)) peak_rss <- max(peak_rss, rss)
    storage <- tree_bytes(result_root)
    if (elapsed > 600 || stage_elapsed > 1200) {
      guard <- "timeout_guard"
      process$kill_tree()
    } else if (peak_rss > 4 * 1024^3) {
      guard <- "rss_guard"
      process$kill_tree()
    } else if (storage > 4 * 1024^3) {
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
  rows[[index]] <- data.frame(
    contract_id = "mv05w_smoke_resource_metric_v1",
    robustness_group_id = unit$robustness_group_id,
    execution_order = unit$execution_order, fold_id = unit$fold_id,
    seed = unit$seed,
    representation = unit$representation,
    configuration_id = unit$configuration_id,
    disposition = disposition,
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
    stage_elapsed_cap_seconds = 1200,
    storage_cap_bytes = 4 * 1024^3,
    maximum_workers = 1L,
    labels_opened = FALSE, outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
  if (disposition %in% c("failed", "timeout_guard", "rss_guard",
                          "storage_guard")) {
    failed <- TRUE
    break
  }
}
metrics <- do.call(rbind, rows[!vapply(rows, is.null, logical(1L))])
write_provenance_csv(metrics, resource_output)
if (failed || nrow(metrics) != nrow(queue)) {
  quit(status = 2L, save = "no")
}
if (any(!metrics$disposition %in% c("completed", "reused_validated")) ||
    any(metrics$completed_views != 90L) ||
    any(metrics$landscape_summary_rows != 180L) ||
    any(metrics$elapsed_seconds > 600) ||
    sum(metrics$elapsed_seconds) > 1200 ||
    any(metrics$peak_process_tree_rss_bytes > 4 * 1024^3) ||
    max(metrics$cumulative_private_bytes) > 4 * 1024^3 ||
    any(metrics$landscape_pair_rows != queue$landscape_request_rows) ||
    any(metrics$energy_pair_rows != queue$energy_request_rows) ||
    any(metrics$method_rows != queue$assembled_method_rows)) {
  stop("MV5-W smoke completed outside its resource contract.",
       call. = FALSE)
}
message("Completed one MV5-W full-pair smoke group (", run_mode, ").")
