#!/usr/bin/env Rscript

for (package in c("processx", "ps")) if (!requireNamespace(
  package, quietly = TRUE
)) stop(package, " is required for MV6-G production rebind.", call. = FALSE)
options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 10L) {
  stop("usage: run_mv06g_production_rebind_monitor.R QUEUE PARENT ",
       "SOURCE_GROUPS POLICY SOURCES GROUP_ROOT RUST_LIBRARY PRIVATE_ROOT ",
       "METRIC GROUP_ID", call. = FALSE)
}
source("R/mv06f_production.R")
source("R/mv06g_production.R")
paths <- vapply(args[1:7], normalizePath, character(1L), winslash = "/",
                mustWork = TRUE)
queue <- utils::read.csv(paths[[1L]], stringsAsFactors = FALSE,
                         check.names = FALSE)
parent <- utils::read.csv(paths[[2L]], stringsAsFactors = FALSE,
                          check.names = FALSE)
source_groups <- utils::read.csv(paths[[3L]], stringsAsFactors = FALSE,
                                 check.names = FALSE)
policy <- utils::read.csv(paths[[4L]], stringsAsFactors = FALSE,
                          check.names = FALSE)
row <- queue[queue$group_id == args[[10L]], , drop = FALSE]
source_group <- source_groups[source_groups$group_id == args[[10L]],
                              , drop = FALSE]
if (nrow(row) != 1L || row$group_id != policy$stage1_equivalence_group_id ||
    !as.logical(policy$rebind_equivalence_authorized)) {
  stop("MV6-G rebind monitor is limited to the stage-one group.", call. = FALSE)
}
private_root <- args[[8L]]; metric_path <- args[[9L]]
dir.create(private_root, recursive = TRUE, showWarnings = FALSE)
group_root <- file.path(private_root, "groups"); dir.create(group_root,
  recursive = TRUE, showWarnings = FALSE)
safe <- gsub("[^A-Za-z0-9_.-]", "_", row$group_id)
if (dir.exists(file.path(group_root, safe)) || file.exists(metric_path)) {
  stop("MV6-G production rebind requires clean state.", call. = FALSE)
}
rss <- function(pid) {
  root <- tryCatch(ps::ps_handle(pid), error = function(error) NULL)
  if (is.null(root)) return(0)
  handles <- c(list(root), tryCatch(ps::ps_children(root, recursive = TRUE),
                                   error = function(error) list()))
  sum(vapply(handles, function(handle) tryCatch(
    as.numeric(ps::ps_memory_info(handle)[["rss"]]),
    error = function(error) 0
  ), numeric(1L)))
}
bytes <- function() {
  files <- list.files(private_root, recursive = TRUE, full.names = TRUE)
  files <- files[file.info(files)$isdir %in% FALSE]
  if (!length(files)) 0 else sum(file.info(files)$size, na.rm = TRUE)
}
runner_args <- c("--vanilla", "scripts/run_mv06g_group.R", paths,
                 row$group_id, group_root)
started <- Sys.time(); peak <- 0; failure <- ""
process <- processx::process$new(Sys.which("Rscript"), runner_args,
  stdout = "|", stderr = "|", cleanup_tree = TRUE)
while (process$is_alive()) {
  Sys.sleep(0.25); peak <- max(peak, rss(process$get_pid()))
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  if (elapsed > policy$elapsed_cap_seconds_per_group) {
    failure <- "elapsed_cap"; process$kill_tree()
  } else if (peak > policy$rss_cap_bytes_per_group) {
    failure <- "rss_cap"; process$kill_tree()
  } else if (bytes() > policy$private_storage_cap_bytes) {
    failure <- "storage_cap"; process$kill_tree()
  }
}
process$wait(timeout = 5000); status <- process$get_exit_status()
elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
complete <- identical(status, 0L) && !nzchar(failure) &&
  dir.exists(file.path(group_root, safe))
if (complete) mv06g_validate_production_group_v1(
  file.path(group_root, safe), row, parent, policy, source_group
)
metric <- data.frame(
  contract_id = "mv06g_production_rebind_metric_v1", group_id = row$group_id,
  disposition = if (complete) "completed" else if (nzchar(failure)) failure else
    "failed", exit_status = status, elapsed_seconds = elapsed,
  peak_process_tree_rss_bytes = peak, private_root_bytes = bytes(),
  complete = complete,
  production_implementation_root_sha256 =
    policy$production_implementation_root_sha256,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  fusion_evaluations = 0L, outcome_jobs = 0L, stringsAsFactors = FALSE
)
dir.create(dirname(metric_path), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(metric, metric_path, row.names = FALSE, na = "")
if (!complete) stop("MV6-G production rebind failed.", call. = FALSE)
message("Completed monitored MV6-G production rebind.")
