#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8L) {
  stop("usage: build_mv06g_production_prefreeze.R QUEUE PARENT SOURCE_GROUPS ",
       "WORKLOAD STAGE1_VALIDATION STAGE1_LAUNCH RUST_LIBRARY OUTPUT_DIR",
       call. = FALSE)
}
source("R/mv06f_production.R")
source("R/mv06g_fusion_prefreeze.R")
source("R/mv06g_stage1.R")
source("R/mv06g_production.R")
queue <- utils::read.csv(args[[1L]], stringsAsFactors = FALSE,
                         check.names = FALSE)
parent <- utils::read.csv(args[[2L]], stringsAsFactors = FALSE,
                          check.names = FALSE)
source_groups <- utils::read.csv(args[[3L]], stringsAsFactors = FALSE,
                                 check.names = FALSE)
workload <- utils::read.csv(args[[4L]], stringsAsFactors = FALSE,
                            check.names = FALSE)
stage1_validation <- utils::read.csv(args[[5L]], stringsAsFactors = FALSE,
                                     check.names = FALSE)
stage1_launch <- utils::read.csv(args[[6L]], stringsAsFactors = FALSE,
                                 check.names = FALSE)
rust <- normalizePath(args[[7L]], winslash = "/", mustWork = TRUE)
output <- args[[8L]]; dir.create(output, recursive = TRUE, showWarnings = FALSE)
mv06f_validate_queue_v1(queue)
remaining <- workload[workload$execution_order != 1L, , drop = FALSE]
remaining <- remaining[order(remaining$execution_order), , drop = FALSE]
rownames(remaining) <- NULL
paths <- mv06g_production_source_paths_v1()
if (!all(file.exists(paths))) stop("MV6-G production sources are incomplete.",
                                   call. = FALSE)
sources <- data.frame(
  contract_id = "mv06g_production_source_inventory_v1", path = paths,
  sha256 = unname(vapply(paths, .mv06f_sha256, character(1L))),
  stringsAsFactors = FALSE
)
root <- digest::digest(stats::setNames(sources$sha256, sources$path),
                       algo = "sha256", serialize = TRUE)
if (nrow(remaining) != 74L || sum(remaining$training_biological_pairs) !=
      260595L || sum(remaining$training_component_rows) != 1042380L ||
    sum(remaining$query_biological_pairs) != 33725L ||
    sum(remaining$query_ranking_rows) != 303525L ||
    nrow(stage1_validation) != 10L || any(!stage1_validation$passed) ||
    stage1_launch$stage1_implementation_root_sha256 !=
      "6a76a11d1b211fcf89fddcd67b7591161619950023420c0bfccbbdc65b76ce82" ||
    parent$rust_library_sha256 != .mv06f_sha256(rust)) {
  stop("MV6-G complete-production prefreeze inputs are stale.", call. = FALSE)
}
policy <- data.frame(
  contract_id = "mv06g_production_rebind_prefreeze_v1",
  parent_stage1_commit = "416d01d",
  parent_contract_sha256 = .mv06f_sha256(args[[2L]]),
  queue_root_sha256 = parent$queue_root_sha256,
  production_implementation_root_sha256 = root,
  rust_library_sha256 = parent$rust_library_sha256,
  remaining_groups = 74L, training_biological_pairs = 260595L,
  training_component_rows = 1042380L,
  query_biological_pairs = 33725L, query_ranking_rows = 303525L,
  stage1_equivalence_group_id = workload$group_id[workload$execution_order == 1L],
  required_equivalence_artifacts = 3L,
  elapsed_cap_seconds_per_group = 1800L,
  rss_cap_bytes_per_group = 12 * 1024^3,
  private_storage_cap_bytes = 5 * 1024^3,
  total_worker_cap_seconds = 12 * 3600,
  maximum_workers = 1L, automatic_retry = FALSE,
  rebind_equivalence_authorized = TRUE,
  remaining_production_authorized = FALSE,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  fusion_evaluations = 0L, outcome_jobs = 0L,
  disposition = "prefreeze_pass_rebind_equivalence_only",
  stringsAsFactors = FALSE
)
write_csv <- function(value, name) utils::write.csv(
  value, file.path(output, name), row.names = FALSE, na = ""
)
write_csv(policy, "mv06g-production-policy.csv")
write_csv(sources, "mv06g-production-sources.csv")
write_csv(remaining, "mv06g-production-queue.csv")
message("Prefroze MV6-G production rebind equivalence only.")
