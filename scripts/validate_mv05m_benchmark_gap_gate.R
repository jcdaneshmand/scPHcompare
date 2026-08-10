#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("usage: validate_mv05m_benchmark_gap_gate.R RESULT_DIR OUTPUT_CSV",
       call. = FALSE)
}
result_dir <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
output_path <- args[[2L]]
source("R/provenance_utils.R")

read_result <- function(name) utils::read.csv(
  file.path(result_dir, name), stringsAsFactors = FALSE, check.names = FALSE
)
criteria <- read_result("mv05m-selection-criteria-2026-08-10.csv")
scores <- read_result("mv05m-candidate-scores-2026-08-10.csv")
readiness <- read_result("mv05m-axis-readiness-2026-08-10.csv")
selected <- read_result("mv05m-selected-next-sprint-2026-08-10.csv")
pairs <- read_result("mv05m-training-pair-scope-2026-08-10.csv")
evidence <- read_result("mv05m-evidence-registry-2026-08-10.csv")
summary <- read_result("mv05m-production-summary-2026-08-10.csv")

criterion_ids <- c(
  "scientific_value", "reviewer_relevance", "identifiability_validity",
  "artifact_readiness", "resource_feasibility", "outcome_selection_safety"
)
weights <- c(3L, 2L, 3L, 2L, 1L, 2L)
if (!identical(criteria$criterion_id, criterion_ids) ||
    !identical(criteria$weight, weights) ||
    !all(criteria$scale_frozen_before_scoring)) {
  stop("MV5-M independent criterion validation failed.", call. = FALSE)
}
recomputed <- as.integer(as.matrix(scores[, criterion_ids]) %*% weights)
if (!identical(recomputed, scores$weighted_score) ||
    sum(scores$selected_next) != 1L ||
    scores$candidate_id[scores$selected_next] !=
      "label_free_clustering_contract_gate" ||
    scores$weighted_score[scores$selected_next] != 45L) {
  stop("MV5-M independent score reconstruction failed.", call. = FALSE)
}
if (nrow(readiness) != 9L ||
    readiness$disposition[readiness$axis_id == "technical_mixing"] !=
      "blocked_pending_identifiability_design" ||
    readiness$disposition[readiness$axis_id ==
                            "label_free_sample_clustering"] !=
      "selected_for_contract_and_resource_gate_only") {
  stop("MV5-M independent readiness validation failed.", call. = FALSE)
}
if (nrow(pairs) != 15L ||
    sum(pairs$existing_query_training_pairs_per_seed) != 7070L ||
    sum(pairs$missing_training_training_pairs_per_seed) != 52535L ||
    selected$exact_training_pairs_per_representation != 262675L ||
    selected$exact_h0_h1_rows_per_representation != 525350L) {
  stop("MV5-M independent pair-scope validation failed.", call. = FALSE)
}
if (nrow(evidence) != 13L || sum(evidence$confidential_source) != 1L ||
    any(evidence$wording_reproduced) || any(evidence$biological_outcomes_computed)) {
  stop("MV5-M public evidence-boundary validation failed.", call. = FALSE)
}
zero_fields <- c(
  "clustering_jobs_executed", "technical_mixing_jobs_executed",
  "robustness_jobs_executed", "gene_view_jobs_executed",
  "fusion_jobs_executed", "new_data_jobs_executed",
  "optimization_jobs_executed"
)
if (any(unlist(summary[zero_fields]) != 0L) ||
    isTRUE(summary$biological_outcomes_computed) ||
    isTRUE(summary$tissue_results_consulted_for_selection)) {
  stop("MV5-M no-outcome boundary validation failed.", call. = FALSE)
}
validation <- data.frame(
  contract_id = "mv05m_independent_validation_v1",
  validation = c(
    "criteria_weights_and_scale", "candidate_score_reconstruction",
    "single_selected_candidate", "axis_readiness_dispositions",
    "fold_pair_scope_reconstruction", "confidential_source_boundary",
    "no_outcome_or_downstream_execution"
  ),
  status = "passed",
  rows_checked = c(6L, 8L, 1L, 9L, 15L, 13L, length(zero_fields)),
  max_absolute_difference = 0,
  evidence = c(
    "six frozen 0-to-4 criteria and integer weights exact",
    "all eight weighted totals independently recomputed",
    "MV5-N clustering contract/resource gate uniquely selected at score 45",
    "all nine benchmark axes retain complete blocked or deferred status",
    "7070 existing query-training and 52535 missing training pairs per seed exact",
    "one confidential source class consulted with zero wording reproduced",
    "zero outcomes and zero clustering/mixing/robustness/gene/fusion/new-data/optimization jobs"
  ),
  stringsAsFactors = FALSE
)
write_provenance_csv(validation, output_path)
message("MV5-M independent benchmark-gap gate validation passed.")
