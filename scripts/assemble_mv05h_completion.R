#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8L) {
  stop(
    "usage: assemble_mv05h_completion.R VIEW_AUDIT_ROOT GROUP_METRICS ",
    "INDEPENDENT_SUMMARY REPEAT_SUMMARY RESUME_VALIDATION VIEW_OUT ",
    "COMPLETION_OUT FAILURE_OUT", call. = FALSE
  )
}
audit_root <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
read_public <- function(path) utils::read.csv(
  normalizePath(path, winslash = "/", mustWork = TRUE),
  stringsAsFactors = FALSE, check.names = FALSE
)
source("R/provenance_utils.R")
source("R/mv05h_integrated_ph_production.R")
group_metrics <- read_public(args[[2L]])
independent <- read_public(args[[3L]])
repeat_summary <- read_public(args[[4L]])
resume <- read_public(args[[5L]])
audit_files <- sort(
  list.files(audit_root, pattern = "__views[.]csv$", full.names = TRUE),
  method = "radix"
)
if (length(audit_files) != 75L) {
  stop("MV5-H completion requires 75 view-audit files.", call. = FALSE)
}
views <- do.call(rbind, lapply(audit_files, function(path) {
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}))
views <- views[order(views$group_order, views$view_order), , drop = FALSE]
rownames(views) <- NULL
mv05h_validate_view_metrics_v1(views, expected_jobs = 6750L)
if (nrow(group_metrics) != 75L ||
    any(group_metrics$disposition != "completed") ||
    !isTRUE(independent$all_checks_pass[[1L]]) ||
    !isTRUE(repeat_summary$all_pass[[1L]]) ||
    !isTRUE(resume$complete_snapshot_identical[[1L]])) {
  stop("MV5-H completion evidence is incomplete.", call. = FALSE)
}
failures <- data.frame(
  contract_id = character(), group_id = character(), job_id = character(),
  failure_code = character(), failure_detail = character(),
  outcome_label_state = character(), biological_outcomes_computed = logical(),
  stringsAsFactors = FALSE
)
completion <- data.frame(
  contract_id = "mv05h_integrated_cell_ph_completion_v1",
  groups = nrow(group_metrics), views = nrow(views),
  held_out_views = sum(views$execution_role == "held_out"),
  training_views = sum(views$execution_role == "training"),
  completed_views = sum(views$disposition == "built_atomic"),
  failed_views = nrow(failures), h0_intervals = sum(views$h0_intervals),
  h1_intervals = sum(views$h1_intervals),
  stored_mst_checks_passed = sum(views$h0_mst_oracle_passed),
  independently_validated_views = independent$views_validated[[1L]],
  independently_recomputed_mst_checks = independent$fresh_mst_checks[[1L]],
  exact_repeat_views = repeat_summary$views_compared[[1L]],
  exact_repeat_views_passed = repeat_summary$byte_matches[[1L]],
  resume_groups = resume$groups[[1L]],
  resume_groups_rebuilt = resume$resumed_groups_rebuilt[[1L]],
  group_worker_seconds = sum(group_metrics$elapsed_seconds),
  view_operation_seconds = sum(views$operation_seconds),
  median_group_seconds = stats::median(group_metrics$elapsed_seconds),
  maximum_group_seconds = max(group_metrics$elapsed_seconds),
  peak_process_tree_rss_bytes =
    max(group_metrics$peak_process_tree_rss_bytes),
  private_result_bytes = sum(views$result_size_bytes),
  all_records_independently_validated = TRUE,
  full_integrated_cell_ph_complete = TRUE,
  landscape_jobs_executed = 0L, distance_jobs_executed = 0L,
  retrieval_jobs_executed = 0L, clustering_jobs_executed = 0L,
  gene_view_jobs_executed = 0L, fusion_jobs_executed = 0L,
  new_data_jobs_executed = 0L, biological_outcomes_computed = FALSE,
  outcome_label_state = "closed", stringsAsFactors = FALSE
)
write_provenance_csv(views, args[[6L]])
write_provenance_csv(completion, args[[7L]])
write_provenance_csv(failures, args[[8L]])
message("Assembled complete public MV5-H view and completion evidence.")
