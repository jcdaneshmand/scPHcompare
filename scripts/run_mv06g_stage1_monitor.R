#!/usr/bin/env Rscript

for (package in c("processx", "ps", "digest")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(package, " is required for monitored MV6-G stage one.",
         call. = FALSE)
  }
}
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1", RAYON_NUM_THREADS = "1")
options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 10L) {
  stop("usage: run_mv06g_stage1_monitor.R QUEUE PARENT_CONTRACT ",
       "SOURCE_GROUPS LAUNCH SOURCES GROUP_ROOT RUST_LIBRARY PRIVATE_ROOT ",
       "PUBLIC_METRIC_CSV RUN_ID", call. = FALSE)
}
source("R/mv06f_production.R")
source("R/mv06g_fusion_prefreeze.R")
source("R/mv06g_stage1.R")
paths <- vapply(args[1:7], normalizePath, character(1L), winslash = "/",
                mustWork = TRUE)
queue <- utils::read.csv(paths[[1L]], stringsAsFactors = FALSE,
                         check.names = FALSE)
parent <- utils::read.csv(paths[[2L]], stringsAsFactors = FALSE,
                          check.names = FALSE)
source_groups <- utils::read.csv(paths[[3L]], stringsAsFactors = FALSE,
                                 check.names = FALSE)
launch <- utils::read.csv(paths[[4L]], stringsAsFactors = FALSE,
                          check.names = FALSE)
sources <- utils::read.csv(paths[[5L]], stringsAsFactors = FALSE,
                           check.names = FALSE)
stage <- queue[queue$stage == "stage_1_maximum", , drop = FALSE]
source_group <- source_groups[source_groups$group_id == stage$group_id,
                              , drop = FALSE]
if (nrow(stage) != 1L || nrow(parent) != 1L || nrow(launch) != 1L ||
    nrow(source_group) != 1L ||
    launch$parent_contract_sha256 != .mv06f_sha256(paths[[2L]]) ||
    launch$rust_library_sha256 != .mv06f_sha256(paths[[7L]]) ||
    !identical(sources$path, mv06g_stage1_source_paths_v1()) ||
    !identical(unname(vapply(sources$path, .mv06f_sha256, character(1L))),
               unname(sources$sha256)) ||
    launch$stage1_implementation_root_sha256 != digest::digest(
      stats::setNames(sources$sha256, sources$path), algo = "sha256",
      serialize = TRUE
    )) {
  stop("MV6-G stage-one monitor preflight failed.", call. = FALSE)
}
private_root <- args[[8L]]
metric_path <- args[[9L]]
run_id <- args[[10L]]
if (!run_id %in% c("primary", "repeat")) {
  stop("MV6-G stage-one run ID must be primary or repeat.", call. = FALSE)
}
dir.create(private_root, recursive = TRUE, showWarnings = FALSE)
group_root <- file.path(private_root, "groups")
log_root <- file.path(private_root, "logs")
dir.create(group_root, recursive = TRUE, showWarnings = FALSE)
dir.create(log_root, recursive = TRUE, showWarnings = FALSE)
safe <- gsub("[^A-Za-z0-9_.-]", "_", stage$group_id)
final_dir <- file.path(group_root, safe)
partials <- list.files(group_root, pattern = paste0("^", safe, "\\.partial\\."),
                       full.names = TRUE)
if (length(partials) || dir.exists(final_dir) || file.exists(metric_path)) {
  stop("MV6-G monitored run requires clean private/output state.", call. = FALSE)
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
stdout_path <- file.path(log_root, paste0(safe, "__", run_id, "__stdout.txt"))
stderr_path <- file.path(log_root, paste0(safe, "__", run_id, "__stderr.txt"))
runner_args <- c("--vanilla", "scripts/run_mv06g_stage1_group.R", paths,
                 group_root)
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
  if (elapsed > launch$elapsed_cap_seconds) {
    cap_failure <- "stage1_elapsed_cap_exceeded"; process$kill_tree()
  } else if (peak > launch$rss_cap_bytes) {
    cap_failure <- "stage1_rss_cap_exceeded"; process$kill_tree()
  } else if (private_bytes() > launch$private_storage_cap_bytes) {
    cap_failure <- "private_storage_cap_exceeded"; process$kill_tree()
  }
}
process$wait(timeout = 5000)
elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
exit_status <- process$get_exit_status()
completed <- identical(exit_status, 0L) && dir.exists(final_dir) &&
  !nzchar(cap_failure)
if (completed) {
  mv06g_validate_stage1_group_v1(final_dir, stage, parent, launch,
                                 source_group)
}
disposition <- if (nzchar(cap_failure)) cap_failure else if (completed) {
  "completed"
} else "failed"
metric <- data.frame(
  contract_id = "mv06g_stage1_resource_metric_v1", run_id = run_id,
  group_id = stage$group_id, fold_id = stage$fold_id, seed = stage$seed,
  training_biological_pairs = 2080L, training_component_rows = 8320L,
  query_ranking_rows = 14625L, disposition = disposition,
  exit_status = if (is.null(exit_status)) NA_integer_ else exit_status,
  elapsed_seconds = elapsed, peak_process_tree_rss_bytes = peak,
  private_root_bytes = private_bytes(),
  elapsed_cap_seconds = launch$elapsed_cap_seconds,
  rss_cap_bytes = launch$rss_cap_bytes,
  private_storage_cap_bytes = launch$private_storage_cap_bytes,
  group_directory_complete = completed,
  parent_contract_sha256 = launch$parent_contract_sha256,
  stage1_implementation_root_sha256 =
    launch$stage1_implementation_root_sha256,
  rust_library_sha256 = launch$rust_library_sha256,
  monitor_sha256 = .mv06f_sha256("scripts/run_mv06g_stage1_monitor.R"),
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  fusion_evaluations = 0L, outcome_jobs = 0L, stringsAsFactors = FALSE
)
dir.create(dirname(metric_path), recursive = TRUE, showWarnings = FALSE)
partial <- tempfile(pattern = basename(metric_path), tmpdir = dirname(metric_path))
utils::write.csv(metric, partial, row.names = FALSE, na = "")
if (!file.rename(partial, metric_path)) {
  stop("MV6-G monitor failed atomic metric publication.", call. = FALSE)
}
if (!completed) {
  stop("MV6-G stage one stopped with disposition: ", disposition,
       call. = FALSE)
}
message("Completed monitored MV6-G stage one in ", format(elapsed, digits = 7),
        " seconds with peak process-tree RSS ", peak, " bytes.")
