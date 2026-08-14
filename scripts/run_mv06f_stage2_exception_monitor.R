#!/usr/bin/env Rscript

options(warn = 2)
for (package in c("processx", "ps", "digest")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(package, " is required for MV6-F resource-exception diagnosis.",
         call. = FALSE)
  }
}
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1", RAYON_NUM_THREADS = "1")
args <- getOption("mv06f.args", commandArgs(trailingOnly = TRUE))
if (length(args) != 13L) {
  stop("usage: run_mv06f_stage2_exception_monitor.R AUTHORIZATION QUEUE ",
       "CONTRACT SOURCES CANDIDATE FOLDS RESOURCES PANEL CACHE_DIR ",
       "RUST_LIBRARY PRIVATE_ROOT PUBLIC_METRIC GROUP_ID", call. = FALSE)
}
paths <- vapply(args[1:10], normalizePath, character(1L), winslash = "/",
                mustWork = TRUE)
private_root <- args[[11L]]; metric_path <- args[[12L]]; group_id <- args[[13L]]
source("R/mv06f_production.R")
queue <- utils::read.csv(paths[[2L]], stringsAsFactors = FALSE,
                         check.names = FALSE)
contract <- utils::read.csv(paths[[3L]], stringsAsFactors = FALSE,
                            check.names = FALSE)
sources <- utils::read.csv(paths[[4L]], stringsAsFactors = FALSE,
                           check.names = FALSE)
authorization <- utils::read.csv(paths[[1L]], stringsAsFactors = FALSE,
                                 check.names = FALSE)
row <- queue[queue$group_id == group_id & queue$stage == "stage_2",
             , drop = FALSE]
source_root <- digest::digest(stats::setNames(sources$sha256, sources$path),
                              algo = "sha256", serialize = TRUE)
if (nrow(row) != 1L || nrow(contract) != 1L || nrow(authorization) != 1L ||
    authorization$contract_id != "mv06f_stage2_resource_exception_v1" ||
    authorization$group_id != group_id ||
    authorization$diagnostic_cap_bytes != 12 * 1024^3 ||
    authorization$automatic_retry != FALSE ||
    authorization$monitor_sha256 != .mv06f_sha256(
      "scripts/run_mv06f_stage2_exception_monitor.R"
    ) || contract$queue_root_sha256 != mv06f_queue_root_v1(queue) ||
    contract$implementation_root_sha256 != source_root ||
    .mv06f_sha256(paths[[10L]]) != contract$rust_library_sha256 ||
    authorization$implementation_root_sha256 != source_root ||
    authorization$rust_library_sha256 != contract$rust_library_sha256) {
  stop("MV6-F resource-exception authorization is stale.", call. = FALSE)
}
group_root <- file.path(private_root, "groups")
log_root <- file.path(private_root, "logs")
dir.create(group_root, recursive = TRUE, showWarnings = FALSE)
dir.create(log_root, recursive = TRUE, showWarnings = FALSE)
safe_id <- gsub("[^A-Za-z0-9_.-]", "_", group_id)
final_dir <- file.path(group_root, safe_id)
if (dir.exists(final_dir) || file.exists(metric_path) ||
    length(list.files(group_root, pattern = paste0(
      "^", safe_id, "\\.partial\\."
    )))) {
  stop("MV6-F diagnostic requires quarantined prior state and no result.",
       call. = FALSE)
}
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
tree_bytes <- function() {
  files <- list.files(private_root, recursive = TRUE, full.names = TRUE,
                      all.files = TRUE, no.. = TRUE)
  if (!length(files)) return(0)
  info <- file.info(files); sum(info$size[!info$isdir], na.rm = TRUE)
}
stdout_path <- file.path(log_root, paste0(safe_id, "__diagnostic_stdout.txt"))
stderr_path <- file.path(log_root, paste0(safe_id, "__diagnostic_stderr.txt"))
runner_args <- c(
  "--vanilla", "scripts/run_mv06f_group.R", paths[[2L]], paths[[3L]],
  paths[[4L]], paths[[5L]], paths[[6L]], paths[[7L]], paths[[8L]], paths[[9L]],
  paths[[10L]], group_id, group_root
)
started <- Sys.time()
process <- processx::process$new(
  Sys.which("Rscript"), runner_args, stdout = stdout_path,
  stderr = stderr_path, cleanup_tree = TRUE
)
peak <- 0; failure <- NA_character_
while (process$is_alive()) {
  Sys.sleep(0.25)
  peak <- max(peak, tree_rss(process$get_pid()))
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  if (elapsed > 1800) failure <- "diagnostic_elapsed_cap_exceeded"
  if (peak > 12 * 1024^3) failure <- "diagnostic_rss_cap_exceeded"
  if (tree_bytes() > 10 * 1024^3) failure <- "private_storage_cap_exceeded"
  if (!is.na(failure)) process$kill_tree()
}
process$wait(timeout = 5000)
elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
exit_status <- process$get_exit_status()
complete <- is.na(failure) && identical(exit_status, 0L) && dir.exists(final_dir)
if (complete) mv06f_validate_group_directory_v1(
  final_dir, row, contract$queue_root_sha256,
  contract$implementation_root_sha256, contract$rust_library_sha256
)
disposition <- if (!is.na(failure)) failure else if (complete) {
  "diagnostic_completed"
} else "diagnostic_failed"
metric <- data.frame(
  contract_id = "mv06f_stage2_resource_exception_metric_v1",
  group_id = group_id, fold_id = row$fold_id, seed = row$seed,
  prior_peak_process_tree_rss_bytes = authorization$prior_peak_rss_bytes,
  prior_cap_bytes = authorization$prior_cap_bytes,
  elapsed_seconds = elapsed, peak_process_tree_rss_bytes = peak,
  diagnostic_cap_bytes = 12 * 1024^3, disposition = disposition,
  exit_status = if (is.null(exit_status)) NA_integer_ else exit_status,
  group_directory_complete = complete, private_root_bytes = tree_bytes(),
  queue_root_sha256 = contract$queue_root_sha256,
  implementation_root_sha256 = contract$implementation_root_sha256,
  rust_library_sha256 = contract$rust_library_sha256,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  fusion_jobs = 0L, clustering_jobs = 0L, outcome_jobs = 0L,
  stringsAsFactors = FALSE
)
dir.create(dirname(metric_path), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(metric, metric_path, row.names = FALSE, na = "")
if (!complete) stop("MV6-F resource diagnosis stopped: ", disposition,
                    call. = FALSE)
message("Completed MV6-F 12-GiB resource diagnosis for ", group_id, ".")
