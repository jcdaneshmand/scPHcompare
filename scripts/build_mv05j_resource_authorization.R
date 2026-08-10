#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 9L) {
  stop(
    "usage: build_mv05j_resource_authorization.R METRICS_CSV ",
    "VALIDATION_SUMMARY_CSV REPRESENTATIVE_ADMISSION_CSV ",
    "MAXIMUM_ADMISSION_CSV RESUME_CSV PUBLIC_REPEAT_CSV RESOURCE_OUTPUT_CSV ",
    "AUTHORIZATION_OUTPUT_CSV OBSERVED_WALL_SECONDS", call. = FALSE
  )
}
source("R/provenance_utils.R")
read_one <- function(path) {
  value <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(value) != 1L) stop("Expected one-row evidence: ", path,
                              call. = FALSE)
  value
}
metrics <- utils::read.csv(args[[1L]], stringsAsFactors = FALSE,
                           check.names = FALSE)
validation <- read_one(args[[2L]])
representative <- read_one(args[[3L]])
maximum <- read_one(args[[4L]])
resume <- read_one(args[[5L]])
repeat_evidence <- utils::read.csv(
  args[[6L]], stringsAsFactors = FALSE, check.names = FALSE
)
wall_seconds <- as.numeric(args[[9L]])

zero_columns <- c(
  "retrieval_evaluation_jobs_executed", "clustering_jobs_executed",
  "integration_jobs_executed", "gene_topology_jobs_executed",
  "fusion_jobs_executed", "new_data_jobs_executed",
  "held_out_scale_fit_jobs_executed"
)
preflight <- nrow(metrics) == 75L &&
  all(metrics$disposition == "completed") &&
  all(metrics$exit_status == 0L) &&
  sum(metrics$biological_pairs) == 35350L &&
  sum(metrics$retrieval_rows) == 176750L &&
  sum(metrics$completed_methods) == 375L &&
  all(metrics$failed_methods == 0L) &&
  all(metrics$elapsed_seconds <= 600) &&
  max(metrics$peak_process_tree_rss_bytes) <= 4294967296 &&
  sum(metrics$elapsed_seconds) <= 7200 &&
  sum(metrics$private_result_bytes) <= 2147483648 &&
  all(unlist(metrics[zero_columns], use.names = FALSE) == 0) &&
  all(metrics$outcome_label_state == "closed") &&
  !any(as.logical(metrics$biological_outcomes_computed)) &&
  validation$groups == 75L && validation$biological_pairs == 35350L &&
  validation$retrieval_rows == 176750L &&
  validation$failed_method_groups == 0L &&
  validation$maximum_numeric_difference <= 1e-12 &&
  representative$repeat_byte_identical && maximum$repeat_byte_identical &&
  max(representative$energy_max_absolute_difference,
      representative$pseudobulk_max_absolute_difference,
      maximum$energy_max_absolute_difference,
      maximum$pseudobulk_max_absolute_difference) <= 1e-12 &&
  resume$exact_snapshot_unchanged && resume$unchanged_files == 150L &&
  nrow(repeat_evidence) == 6L && all(repeat_evidence$byte_identical) &&
  is.finite(wall_seconds) && wall_seconds > 0

resource <- metrics
resource$contract_id <- "mv05j_public_group_resource_metric_v1"
write_provenance_csv(resource, args[[7L]])
authorization <- data.frame(
  contract_id = "mv05j_integrated_retrieval_evaluation_authorization_v1",
  decision = if (preflight) {
    "approve_separate_prediction_locked_integrated_retrieval_evaluation"
  } else "hold",
  authorized = preflight,
  groups = nrow(metrics), seeds = length(unique(metrics$seed)),
  biological_pairs = sum(metrics$biological_pairs),
  retrieval_rows = sum(metrics$retrieval_rows),
  completed_method_groups = sum(metrics$completed_methods),
  failed_method_groups = sum(metrics$failed_methods),
  aggregate_worker_seconds = sum(metrics$elapsed_seconds),
  aggregate_scientific_operation_seconds = sum(metrics$operation_seconds),
  observed_queue_wall_seconds = wall_seconds,
  median_group_seconds = stats::median(metrics$elapsed_seconds),
  maximum_group_seconds = max(metrics$elapsed_seconds),
  peak_process_tree_rss_bytes = max(metrics$peak_process_tree_rss_bytes),
  private_result_bytes = sum(metrics$private_result_bytes),
  maximum_numeric_difference = validation$maximum_numeric_difference,
  immutable_resume_files = resume$unchanged_files,
  byte_identical_public_artifacts = sum(repeat_evidence$byte_identical),
  retrieval_evaluation_jobs_executed = 0L,
  clustering_jobs_executed = 0L, integration_jobs_executed = 0L,
  gene_topology_jobs_executed = 0L, fusion_jobs_executed = 0L,
  new_data_jobs_executed = 0L, held_out_scale_fit_jobs_executed = 0L,
  biological_outcomes_computed = FALSE, outcome_label_state = "closed",
  authorized_scope = paste(
    "later separately specified prediction-locked integrated retrieval",
    "endpoint evaluation after public-hash verification"
  ),
  prohibited_scope = paste(
    "method tuning or selection; SCT comparison; clustering; gene topology;",
    "fusion; new data; manuscript claims"
  ),
  stringsAsFactors = FALSE
)
write_provenance_csv(authorization, args[[8L]])
if (!preflight) stop("MV5-J authorization preflight failed.", call. = FALSE)
message("MV5-J authorizes only a later prediction-locked evaluation sprint.")
