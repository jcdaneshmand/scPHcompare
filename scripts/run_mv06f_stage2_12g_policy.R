#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 12L) {
  stop("usage: run_mv06f_stage2_12g_policy.R POLICY QUEUE CONTRACT SOURCES ",
       "CANDIDATE FOLDS RESOURCES PANEL CACHE_DIR RUST_LIBRARY PRIVATE_ROOT ",
       "PUBLIC_METRIC_DIR", call. = FALSE)
}
source("R/mv06f_production.R")
policy <- utils::read.csv(args[[1L]], stringsAsFactors = FALSE,
                          check.names = FALSE)
queue <- utils::read.csv(args[[2L]], stringsAsFactors = FALSE,
                         check.names = FALSE)
contract <- utils::read.csv(args[[3L]], stringsAsFactors = FALSE,
                            check.names = FALSE)
sources <- utils::read.csv(args[[4L]], stringsAsFactors = FALSE,
                           check.names = FALSE)
stage2 <- queue[queue$stage == "stage_2", , drop = FALSE]
stage2 <- stage2[order(stage2$execution_order), , drop = FALSE]
source_root <- digest::digest(stats::setNames(sources$sha256, sources$path),
                              algo = "sha256", serialize = TRUE)
if (nrow(policy) != 74L || !identical(policy$group_id, stage2$group_id) ||
    any(policy$diagnostic_cap_bytes != 12 * 1024^3) ||
    any(as.logical(policy$automatic_retry)) ||
    any(policy$monitor_sha256 != .mv06f_sha256(
      "scripts/run_mv06f_stage2_exception_monitor.R"
    )) || any(policy$driver_sha256 != .mv06f_sha256(
      "scripts/run_mv06f_stage2_12g_policy.R"
    )) || contract$implementation_root_sha256 != source_root ||
    contract$queue_root_sha256 != mv06f_queue_root_v1(queue) ||
    .mv06f_sha256(args[[10L]]) != contract$rust_library_sha256) {
  stop("MV6-F 12-GiB execution policy is stale.", call. = FALSE)
}
private_root <- args[[11L]]; metric_dir <- args[[12L]]
group_root <- file.path(private_root, "groups")
dir.create(group_root, recursive = TRUE, showWarnings = FALSE)
dir.create(metric_dir, recursive = TRUE, showWarnings = FALSE)
auth_dir <- file.path(private_root, "policy-authorizations")
dir.create(auth_dir, recursive = TRUE, showWarnings = FALSE)
safe_name <- function(value) gsub("[^A-Za-z0-9_.-]", "_", value)
for (index in seq_len(nrow(stage2))) {
  unit <- stage2[index, , drop = FALSE]
  safe <- safe_name(unit$group_id)
  final_dir <- file.path(group_root, safe)
  metric_path <- file.path(metric_dir, paste0(safe, "__resource.csv"))
  if (dir.exists(final_dir)) {
    mv06f_validate_group_directory_v1(
      final_dir, unit, contract$queue_root_sha256,
      contract$implementation_root_sha256, contract$rust_library_sha256
    )
    if (unit$execution_order == 2L) {
      message("Reused admitted 12-GiB diagnostic group: ", unit$group_id)
      next
    }
    if (!file.exists(metric_path)) {
      stop("MV6-F policy found a group without its resource metric.",
           call. = FALSE)
    }
    metric <- utils::read.csv(metric_path, stringsAsFactors = FALSE)
    if (metric$disposition != "diagnostic_completed") {
      stop("MV6-F prior policy metric is not successful.", call. = FALSE)
    }
    message("Reused validated 12-GiB policy group: ", unit$group_id)
    next
  }
  auth_path <- file.path(auth_dir, paste0(safe, "__authorization.csv"))
  utils::write.csv(policy[index, , drop = FALSE], auth_path,
                   row.names = FALSE, na = "")
  child_args <- c(
    "--vanilla", "scripts/run_mv06f_stage2_exception_monitor.R",
    auth_path, args[[2L]], args[[3L]], args[[4L]], args[[5L]], args[[6L]],
    args[[7L]], args[[8L]], args[[9L]], args[[10L]], private_root,
    metric_path, unit$group_id
  )
  status <- system2(
    Sys.which("Rscript"),
    c("--vanilla", "scripts/run_mv06f_stage2_exception_monitor.R",
      shQuote(child_args[-c(1L, 2L)]))
  )
  if (status != 0L || !dir.exists(final_dir) || !file.exists(metric_path)) {
    stop("MV6-F 12-GiB policy stopped at ", unit$group_id, call. = FALSE)
  }
  metric <- utils::read.csv(metric_path, stringsAsFactors = FALSE)
  if (metric$disposition != "diagnostic_completed") {
    stop("MV6-F 12-GiB group did not complete.", call. = FALSE)
  }
  message("Completed MV6-F 12-GiB policy group ", index, "/74: ",
          unit$group_id)
}
message("Completed all MV6-F stage-two groups under the 12-GiB policy.")
