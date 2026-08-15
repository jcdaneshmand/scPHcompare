#!/usr/bin/env Rscript

for (package in c("processx", "ps")) if (!requireNamespace(
  package, quietly = TRUE
)) stop(package, " is required for MV6-G completion.", call. = FALSE)
options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 12L) stop(
  "usage: run_mv06g_completion_monitor.R COMPLETION_POLICY QUEUE PARENT ",
  "SOURCE_GROUPS REBIND_POLICY REBIND_SOURCES GROUP_ROOT RUST_LIBRARY ",
  "PRIVATE_ROOT METRIC_DIR GROUP_ID COMPLETION_SOURCES", call. = FALSE
)
source("R/mv06f_production.R")
source("R/mv06g_production.R")
source("R/mv06g_completion.R")
paths <- vapply(args[c(1:8, 12)], normalizePath, character(1L), winslash = "/",
                mustWork = TRUE)
read_csv <- function(path) utils::read.csv(path, stringsAsFactors = FALSE,
                                            check.names = FALSE)
policy <- read_csv(paths[[1L]]); queue <- read_csv(paths[[2L]])
parent <- read_csv(paths[[3L]]); source_groups <- read_csv(paths[[4L]])
rebind <- read_csv(paths[[5L]]); rebind_sources <- read_csv(paths[[6L]])
completion_sources <- read_csv(paths[[9L]])
row <- queue[queue$group_id == args[[11L]], , drop = FALSE]
source_group <- source_groups[source_groups$group_id == args[[11L]],
                              , drop = FALSE]
if (nrow(row) != 1L || nrow(source_group) != 1L ||
    !identical(as.integer(row$biological_pairs),
               as.integer(row$query_biological_pairs)) ||
    !as.logical(policy$remaining_production_authorized) ||
    policy$scientific_implementation_root_sha256 !=
      rebind$production_implementation_root_sha256 ||
    policy$execution_implementation_root_sha256 !=
      mv06g_completion_root_v1(completion_sources) ||
    !identical(rebind_sources$path, mv06g_production_source_paths_v1()) ||
    !identical(unname(vapply(rebind_sources$path, .mv06f_sha256, character(1L))),
               unname(rebind_sources$sha256)) ||
    .mv06f_sha256(paths[[8L]]) != policy$rust_library_sha256) {
  stop("MV6-G completion monitor identities are stale.", call. = FALSE)
}
private_root <- args[[9L]]; metric_dir <- args[[10L]]
group_output <- file.path(private_root, "groups")
dir.create(group_output, recursive = TRUE, showWarnings = FALSE)
dir.create(metric_dir, recursive = TRUE, showWarnings = FALSE)
safe <- mv06g_safe_group_name_v1(row$group_id)
final_dir <- file.path(group_output, safe)
metric_path <- file.path(metric_dir, paste0(safe, "__resource.csv"))
prior_paths <- list.files(metric_dir, pattern = "__resource\\.csv$",
                          full.names = TRUE)
prior <- if (length(prior_paths)) {
  do.call(rbind, lapply(prior_paths, read_csv))
} else data.frame()
prior_worker <- if (nrow(prior)) sum(prior$charged_worker_seconds) else 0
if (dir.exists(final_dir)) {
  if (!file.exists(metric_path)) stop(
    "MV6-G completion group lacks its resource metric.", call. = FALSE
  )
  mv06g_validate_production_group_v1(final_dir, row, parent, rebind,
                                     source_group)
  mv06g_validate_completion_metric_v1(read_csv(metric_path), row, policy)
  message("Reused validated MV6-G completion group: ", row$group_id)
  quit(status = 0L, save = "no")
}
if (file.exists(metric_path)) stop(
  "MV6-G completion metric exists without its group.", call. = FALSE
)
partials <- list.files(group_output, pattern = paste0("^", safe,
  "\\.partial\\."), full.names = TRUE)
if (length(partials)) stop(
  "MV6-G completion partial state requires quarantine.", call. = FALSE
)
log_dir <- file.path(private_root, "runner-logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
stdout_log <- file.path(log_dir, paste0(safe, ".stdout.log"))
stderr_log <- file.path(log_dir, paste0(safe, ".stderr.log"))
if (file.exists(stdout_log) || file.exists(stderr_log)) stop(
  "MV6-G completion runner logs require quarantine.", call. = FALSE
)
tree_rss <- function(pid) {
  root <- tryCatch(ps::ps_handle(pid), error = function(error) NULL)
  if (is.null(root)) return(0)
  handles <- c(list(root), tryCatch(ps::ps_children(root, recursive = TRUE),
                                   error = function(error) list()))
  sum(vapply(handles, function(handle) tryCatch(
    as.numeric(ps::ps_memory_info(handle)[["rss"]]),
    error = function(error) 0
  ), numeric(1L)))
}
private_bytes <- function() {
  files <- list.files(private_root, recursive = TRUE, full.names = TRUE)
  files <- files[file.info(files)$isdir %in% FALSE]
  if (!length(files)) 0 else sum(file.info(files)$size, na.rm = TRUE)
}
runner_args <- c(
  "--vanilla", "scripts/run_mv06g_group.R", paths[[2L]], paths[[3L]],
  paths[[4L]], paths[[5L]], paths[[6L]], paths[[7L]], paths[[8L]],
  row$group_id, group_output
)
started <- Sys.time(); peak <- 0; failure <- NA_character_
process <- processx::process$new(Sys.which("Rscript"), runner_args,
  stdout = stdout_log, stderr = stderr_log, cleanup_tree = TRUE
)
while (process$is_alive()) {
  Sys.sleep(0.25)
  peak <- max(peak, tree_rss(process$get_pid()))
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  guard <- mv06g_completion_guard_v1(
    prior_worker, elapsed, peak, private_bytes(), policy
  )
  if (!guard$pass) {
    failure <- guard$disposition
    process$kill_tree()
  }
}
process$wait(timeout = 5000); status <- process$get_exit_status()
elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
complete <- identical(status, 0L) && is.na(failure) && dir.exists(final_dir)
if (complete) mv06g_validate_production_group_v1(
  final_dir, row, parent, rebind, source_group
)
metric <- data.frame(
  contract_id = "mv06g_serial_completion_metric_v1", group_id = row$group_id,
  execution_order = as.integer(row$execution_order),
  disposition = if (complete) "completed" else if (!is.na(failure)) failure
    else "failed", exit_status = status, elapsed_seconds = elapsed,
  charged_worker_seconds = elapsed, peak_process_tree_rss_bytes = peak,
  cumulative_worker_seconds = prior_worker + elapsed,
  cumulative_private_bytes = private_bytes(), scientific_group_complete = complete,
  scientific_implementation_root_sha256 =
    policy$scientific_implementation_root_sha256,
  execution_implementation_root_sha256 =
    policy$execution_implementation_root_sha256,
  rust_library_sha256 = policy$rust_library_sha256,
  runner_stdout_sha256 = .mv06f_sha256(stdout_log),
  runner_stderr_sha256 = .mv06f_sha256(stderr_log),
  runner_stderr_tail = paste(tail(readLines(stderr_log, warn = FALSE), 20L),
                             collapse = " | "),
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  fusion_evaluations = 0L, outcome_jobs = 0L, stringsAsFactors = FALSE
)
utils::write.csv(metric, metric_path, row.names = FALSE, na = "")
if (!complete) stop("MV6-G completion stopped at ", row$group_id,
                    call. = FALSE)
message("Completed monitored MV6-G group: ", row$group_id)
