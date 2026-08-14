#!/usr/bin/env Rscript

for (package in c("processx", "ps", "digest")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(package, " is required for monitored MV6-F stage 1.",
         call. = FALSE)
  }
}
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1", RAYON_NUM_THREADS = "1")
options(warn = 2)
args <- getOption("mv06f.args", commandArgs(trailingOnly = TRUE))
if (length(args) != 11L) {
  stop("usage: run_mv06f_stage1_monitor.R QUEUE CONTRACT SOURCES CANDIDATE ",
       "FOLDS RESOURCES PANEL CACHE_DIR RUST_LIBRARY PRIVATE_ROOT ",
       "PUBLIC_METRIC_CSV", call. = FALSE)
}
paths <- vapply(args[1:9], normalizePath, character(1L), winslash = "/",
                mustWork = TRUE)
private_root <- args[[10L]]
metric_path <- args[[11L]]
source("R/mv06f_production.R")
queue <- utils::read.csv(paths[[1L]], stringsAsFactors = FALSE,
                         check.names = FALSE)
contract <- utils::read.csv(paths[[2L]], stringsAsFactors = FALSE,
                            check.names = FALSE)
sources <- utils::read.csv(paths[[3L]], stringsAsFactors = FALSE,
                           check.names = FALSE)
mv06f_validate_queue_v1(queue)
stage <- queue[queue$stage == "stage_1_maximum", , drop = FALSE]
if (nrow(stage) != 1L || nrow(contract) != 1L ||
    contract$queue_root_sha256 != mv06f_queue_root_v1(queue) ||
    !all(file.exists(sources$path)) ||
    !identical(tolower(unname(vapply(
      sources$path, .mv06f_sha256, character(1L)
    ))), tolower(unname(sources$sha256))) ||
    .mv06f_sha256(paths[[9L]]) != contract$rust_library_sha256 ||
    contract$implementation_root_sha256 != digest::digest(
      stats::setNames(sources$sha256, sources$path),
      algo = "sha256", serialize = TRUE
    )) {
  stop("MV6-F stage-1 monitor preflight identity validation failed.",
       call. = FALSE)
}
dir.create(private_root, recursive = TRUE, showWarnings = FALSE)
group_root <- file.path(private_root, "groups")
log_root <- file.path(private_root, "logs")
dir.create(group_root, recursive = TRUE, showWarnings = FALSE)
dir.create(log_root, recursive = TRUE, showWarnings = FALSE)
safe_id <- gsub("[^A-Za-z0-9_.-]", "_", stage$group_id)
final_dir <- file.path(group_root, safe_id)
partials <- list.files(group_root, pattern = paste0(
  "^", safe_id, "\\.partial\\."
), full.names = TRUE)
if (length(partials)) {
  stop("MV6-F stage 1 found partial state; quarantine is required.",
       call. = FALSE)
}
if (file.exists(metric_path)) {
  stop("MV6-F stage-1 metric already exists; refusing overwrite.",
       call. = FALSE)
}

process_tree_rss <- function(pid) {
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
  values <- list.files(private_root, recursive = TRUE, full.names = TRUE,
                       all.files = TRUE, no.. = TRUE)
  values <- values[file.info(values)$isdir %in% FALSE]
  if (!length(values)) 0 else sum(file.info(values)$size, na.rm = TRUE)
}
stdout_path <- file.path(log_root, paste0(safe_id, "__stdout.txt"))
stderr_path <- file.path(log_root, paste0(safe_id, "__stderr.txt"))
runner_args <- c(
  "--vanilla", "scripts/run_mv06f_group.R", paths[[1L]], paths[[2L]],
  paths[[3L]], paths[[4L]], paths[[5L]], paths[[6L]], paths[[7L]], paths[[8L]],
  paths[[9L]], stage$group_id, group_root
)
started <- Sys.time()
process <- processx::process$new(
  command = Sys.which("Rscript"), args = runner_args,
  stdout = stdout_path, stderr = stderr_path, cleanup_tree = TRUE
)
peak <- 0
cap_failure <- ""
while (process$is_alive()) {
  Sys.sleep(0.25)
  peak <- max(peak, process_tree_rss(process$get_pid()))
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  if (elapsed > 1200) {
    cap_failure <- "stage1_elapsed_cap_exceeded"
    process$kill_tree()
  } else if (peak > 6 * 1024^3) {
    cap_failure <- "stage1_rss_cap_exceeded"
    process$kill_tree()
  } else if (private_bytes() > 10 * 1024^3) {
    cap_failure <- "private_storage_cap_exceeded"
    process$kill_tree()
  }
}
process$wait(timeout = 5000)
elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
exit_status <- process$get_exit_status()
completed <- identical(exit_status, 0L) && dir.exists(final_dir) &&
  !nzchar(cap_failure)
if (completed) {
  mv06f_validate_group_directory_v1(
    final_dir, stage, contract$queue_root_sha256,
    contract$implementation_root_sha256, contract$rust_library_sha256
  )
}
disposition <- if (nzchar(cap_failure)) cap_failure else if (completed) {
  "completed"
} else "failed"
metric <- data.frame(
  contract_id = "mv06f_stage1_resource_metric_v1",
  group_id = stage$group_id, fold_id = stage$fold_id, seed = stage$seed,
  biological_pairs = stage$biological_pairs,
  landscape_component_rows = stage$landscape_component_rows,
  disposition = disposition,
  exit_status = if (is.null(exit_status)) NA_integer_ else exit_status,
  elapsed_seconds = elapsed, peak_process_tree_rss_bytes = peak,
  private_root_bytes = private_bytes(), elapsed_cap_seconds = 1200,
  rss_cap_bytes = 6 * 1024^3, private_storage_cap_bytes = 10 * 1024^3,
  group_directory_complete = completed,
  queue_root_sha256 = contract$queue_root_sha256,
  implementation_root_sha256 = contract$implementation_root_sha256,
  rust_library_sha256 = contract$rust_library_sha256,
  monitor_sha256 = .mv06f_sha256("scripts/run_mv06f_stage1_monitor.R"),
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  fusion_jobs = 0L, clustering_jobs = 0L, outcome_jobs = 0L,
  stringsAsFactors = FALSE
)
dir.create(dirname(metric_path), recursive = TRUE, showWarnings = FALSE)
partial_metric <- tempfile(pattern = basename(metric_path),
                           tmpdir = dirname(metric_path))
utils::write.csv(metric, partial_metric, row.names = FALSE, na = "")
if (!file.rename(partial_metric, metric_path)) {
  stop("MV6-F failed to atomically publish its stage-1 metric.",
       call. = FALSE)
}
if (!completed) {
  stop("MV6-F stage 1 stopped with disposition: ", disposition,
       call. = FALSE)
}
message("Completed monitored MV6-F stage 1 in ", format(elapsed, digits = 7),
        " seconds with peak process-tree RSS ", peak, " bytes.")
