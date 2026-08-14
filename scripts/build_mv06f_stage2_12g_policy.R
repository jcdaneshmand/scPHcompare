#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 9L) {
  stop("usage: build_mv06f_stage2_12g_policy.R QUEUE CONTRACT SOURCES ",
       "FIRST_FAILURE DIAGNOSTIC_RESULT EXCEPTION_MONITOR DRIVER ",
       "RUST_LIBRARY OUTPUT", call. = FALSE)
}
source("R/mv06f_production.R")
queue <- utils::read.csv(args[[1L]], stringsAsFactors = FALSE,
                         check.names = FALSE)
contract <- utils::read.csv(args[[2L]], stringsAsFactors = FALSE,
                            check.names = FALSE)
sources <- utils::read.csv(args[[3L]], stringsAsFactors = FALSE,
                           check.names = FALSE)
failure <- utils::read.csv(args[[4L]], stringsAsFactors = FALSE,
                           check.names = FALSE)
diagnostic <- utils::read.csv(args[[5L]], stringsAsFactors = FALSE,
                              check.names = FALSE)
stage2 <- queue[queue$stage == "stage_2", , drop = FALSE]
stage2 <- stage2[order(stage2$execution_order), , drop = FALSE]
source_root <- digest::digest(stats::setNames(sources$sha256, sources$path),
                              algo = "sha256", serialize = TRUE)
if (nrow(stage2) != 74L || nrow(failure) != 1L || nrow(diagnostic) != 1L ||
    failure$group_id != diagnostic$group_id ||
    failure$disposition != "group_rss_cap_exceeded" ||
    diagnostic$disposition != "diagnostic_completed" ||
    !as.logical(diagnostic$group_directory_complete) ||
    diagnostic$peak_process_tree_rss_bytes <= 8 * 1024^3 ||
    diagnostic$peak_process_tree_rss_bytes >= 12 * 1024^3 ||
    contract$queue_root_sha256 != mv06f_queue_root_v1(queue) ||
    contract$implementation_root_sha256 != source_root ||
    .mv06f_sha256(args[[8L]]) != contract$rust_library_sha256) {
  stop("MV6-F 12-GiB policy inputs are invalid.", call. = FALSE)
}
policy <- data.frame(
  contract_id = "mv06f_stage2_resource_exception_v1",
  group_id = stage2$group_id, fold_id = stage2$fold_id, seed = stage2$seed,
  execution_order = stage2$execution_order,
  biological_pairs = stage2$biological_pairs,
  landscape_component_rows = stage2$landscape_component_rows,
  parent_failure_sha256 = .mv06f_sha256(args[[4L]]),
  diagnostic_result_sha256 = .mv06f_sha256(args[[5L]]),
  prior_peak_rss_bytes = failure$peak_process_tree_rss_bytes,
  prior_cap_bytes = 8 * 1024^3, diagnostic_cap_bytes = 12 * 1024^3,
  elapsed_cap_seconds = 1800, private_storage_cap_bytes = 10 * 1024^3,
  monitor_sha256 = .mv06f_sha256(args[[6L]]),
  driver_sha256 = .mv06f_sha256(args[[7L]]),
  scientific_runner_sha256 = sources$sha256[
    sources$path == "scripts/run_mv06f_group.R"
  ], queue_root_sha256 = contract$queue_root_sha256,
  implementation_root_sha256 = contract$implementation_root_sha256,
  rust_library_sha256 = contract$rust_library_sha256,
  maximum_workers = 1L, automatic_retry = FALSE,
  authorization_scope = "serial_stage2_12g_after_observed_diagnostic",
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  fusion_jobs = 0L, clustering_jobs = 0L, outcome_jobs = 0L,
  stringsAsFactors = FALSE
)
dir.create(dirname(args[[9L]]), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(policy, args[[9L]], row.names = FALSE, na = "")
message("Prefroze serial 12-GiB policy for 74 MV6-F stage-two groups.")
