#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) {
  stop("usage: build_mv06f_stage2_exception_prefreeze.R QUEUE CONTRACT ",
       "SOURCES FAILED_METRIC MONITOR RUST_LIBRARY OUTPUT", call. = FALSE)
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
row <- queue[queue$group_id == failure$group_id, , drop = FALSE]
source_root <- digest::digest(stats::setNames(sources$sha256, sources$path),
                              algo = "sha256", serialize = TRUE)
if (nrow(failure) != 1L || nrow(row) != 1L || row$stage != "stage_2" ||
    failure$disposition != "group_rss_cap_exceeded" ||
    failure$peak_process_tree_rss_bytes <= 8 * 1024^3 ||
    failure$peak_process_tree_rss_bytes >= 12 * 1024^3 ||
    contract$queue_root_sha256 != mv06f_queue_root_v1(queue) ||
    contract$implementation_root_sha256 != source_root ||
    .mv06f_sha256(args[[6L]]) != contract$rust_library_sha256) {
  stop("MV6-F stage-two exception prefreeze inputs are invalid.",
       call. = FALSE)
}
authorization <- data.frame(
  contract_id = "mv06f_stage2_resource_exception_v1",
  group_id = failure$group_id, fold_id = failure$fold_id,
  seed = failure$seed,
  parent_failure_sha256 = .mv06f_sha256(args[[4L]]),
  prior_peak_rss_bytes = failure$peak_process_tree_rss_bytes,
  prior_cap_bytes = 8 * 1024^3, diagnostic_cap_bytes = 12 * 1024^3,
  elapsed_cap_seconds = 1800, private_storage_cap_bytes = 10 * 1024^3,
  monitor_sha256 = .mv06f_sha256(args[[5L]]),
  scientific_runner_sha256 = sources$sha256[
    sources$path == "scripts/run_mv06f_group.R"
  ], queue_root_sha256 = contract$queue_root_sha256,
  implementation_root_sha256 = contract$implementation_root_sha256,
  rust_library_sha256 = contract$rust_library_sha256,
  maximum_workers = 1L, automatic_retry = FALSE,
  authorization_scope = "one_group_resource_diagnosis_only",
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  fusion_jobs = 0L, clustering_jobs = 0L, outcome_jobs = 0L,
  stringsAsFactors = FALSE
)
dir.create(dirname(args[[7L]]), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(authorization, args[[7L]], row.names = FALSE, na = "")
message("Prefroze one MV6-F 12-GiB resource diagnosis.")
