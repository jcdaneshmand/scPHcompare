#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("usage: build_mv05m_benchmark_gap_gate.R OUTPUT_DIR", call. = FALSE)
}
output_dir <- args[[1L]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
source("R/provenance_utils.R")
source("R/mv05m_benchmark_gap_gate.R")

criteria <- .mv05m_criteria
criteria$contract_id <- "mv05m_selection_criteria_v1"
criteria$scale_frozen_before_scoring <- TRUE
criteria <- criteria[, c(
  "contract_id", "criterion_id", "weight", "minimum", "maximum",
  "direction", "scale_frozen_before_scoring"
)]

scores <- mv05m_score_candidates_v1()
scores$contract_id <- "mv05m_candidate_score_v1"
scores$biological_outcomes_computed <- FALSE
scores$tissue_results_consulted_for_selection <- FALSE
scores <- scores[, c(
  "contract_id", "candidate_id", .mv05m_criteria$criterion_id,
  "weighted_score", "validity_gate", "selection_eligible",
  "selection_rank", "selected_next", "preliminary_disposition",
  "biological_outcomes_computed", "tissue_results_consulted_for_selection"
)]
readiness <- mv05m_axis_readiness_v1()
readiness$contract_id <- "mv05m_axis_readiness_v1"
readiness <- readiness[, c(
  "contract_id", "axis_id", "scientific_estimand", "primary_unit",
  "current_evidence", "dominant_validity_risk", "missing_computation",
  "disposition"
)]
selected <- mv05m_selected_sprint_v1(scores)

endpoints <- utils::read.csv(
  "docs/audits/mv05e-query-endpoints-2026-08-08.csv",
  stringsAsFactors = FALSE, check.names = FALSE
)
one_axis <- endpoints[
  endpoints$method_id == "cell_landscape_h0_v1" &
    endpoints$seed == 20260805L, , drop = FALSE
]
folds <- split(one_axis, one_axis$fold_id)
pair_scope <- do.call(rbind, lapply(sort(names(folds), method = "radix"),
  function(fold_id) {
    part <- folds[[fold_id]]
    query_samples <- length(unique(part$query_sample_id))
    training_samples <- unique(part$training_samples)
    if (length(training_samples) != 1L ||
        query_samples + training_samples != 90L) {
      stop("MV5-M fold pair-scope reconstruction failed.", call. = FALSE)
    }
    data.frame(
      contract_id = "mv05m_training_pair_scope_v1",
      fold_id = fold_id,
      held_out_samples = query_samples,
      training_samples = training_samples,
      existing_query_training_pairs_per_seed =
        query_samples * training_samples,
      missing_training_training_pairs_per_seed =
        training_samples * (training_samples - 1L) / 2L,
      query_query_pairs_not_required_for_inductive_assignment =
        query_samples * (query_samples - 1L) / 2L,
      stringsAsFactors = FALSE
    )
  }))
if (sum(pair_scope$existing_query_training_pairs_per_seed) != 7070L ||
    sum(pair_scope$missing_training_training_pairs_per_seed) != 52535L) {
  stop("MV5-M aggregate pair scope does not match the frozen fold geometry.",
       call. = FALSE)
}

evidence <- data.frame(
  contract_id = "mv05m_evidence_registry_v1",
  evidence_id = c(
    "dissertation_aims", "preprint_aims", "project_evidence",
    "reviewer_workstreams", "statistical_benchmark_contract",
    "mv05e_sct_retrieval", "mv05k_integrated_retrieval",
    "mv05l_representation_did", "sct_distance_pair_scope",
    "integrated_distance_pair_scope", "clustering_contract_helpers",
    "gene_view_feasibility", "optimization_gate"
  ),
  source = c(
    "local_untracked_dissertation_pdf_via_existing_evidence_ledger",
    "local_untracked_preprint_pdf_via_existing_evidence_ledger",
    "docs/PROJECT_EVIDENCE.md",
    "confidential_reports_via_git_ignored_response_matrix",
    "docs/specifications/MV05_STATISTICAL_BENCHMARK_PLAN_V1.md",
    "docs/audits/MV05E_PREDICTION_LOCKED_RETRIEVAL_EVALUATION_2026-08-08.md",
    "docs/audits/MV05K_PREDICTION_LOCKED_INTEGRATED_RETRIEVAL_EVALUATION_2026-08-10.md",
    "docs/audits/MV05L_LOCKED_REPRESENTATION_COMPARISON_2026-08-10.md",
    "R/mv05d4_landscape_production.R",
    "R/mv05i_integrated_landscape_production.R",
    "R/mv05_benchmark_contract.R",
    "docs/audits/MV05C_EXISTING_DATA_FEASIBILITY_2026-08-06.md",
    "docs/audits/LANDSCAPE_ORACLE_AND_DIAGRAM_ELIGIBILITY_2026-08-05.md"
  ),
  evidence_use = c(
    "sample_level_topological_clustering_and_future_direction aims",
    "integration impact and sample clustering aims",
    "scientific risks page mappings and public reviewer themes",
    "generalized relevance ranking only no wording published",
    "frozen endpoints label firewall clustering hierarchy and multiplicity",
    "accepted biological conservation evidence",
    "accepted integrated biological conservation evidence",
    "accepted paired representation effect evidence",
    "query_to_training artifacts explicitly exclude full clustering matrices",
    "query_to_training artifacts explicitly exclude full clustering matrices",
    "tested label_free stable_k helper available",
    "documents partial SCT and unavailable integrated gene folds",
    "documents current negative Rust decision and measured feasibility"
  ),
  confidential_source = c(FALSE, FALSE, FALSE, TRUE, rep(FALSE, 9L)),
  wording_reproduced = FALSE,
  biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)

summary <- data.frame(
  contract_id = "mv05m_production_summary_v1",
  source_commit = "a3a98895d47b99f8a12691b10a2320b82d195433",
  axes_audited = nrow(readiness), candidates_scored = nrow(scores),
  criteria = nrow(criteria), selected_candidates = sum(scores$selected_next),
  selected_next_sprint = selected$next_sprint_id,
  selected_candidate_id = selected$selected_candidate_id,
  exact_training_pairs_per_representation =
    selected$exact_training_pairs_per_representation,
  exact_h0_h1_rows_per_representation =
    selected$exact_h0_h1_rows_per_representation,
  biological_outcomes_computed = FALSE,
  tissue_results_consulted_for_selection = FALSE,
  clustering_jobs_executed = 0L, technical_mixing_jobs_executed = 0L,
  robustness_jobs_executed = 0L, gene_view_jobs_executed = 0L,
  fusion_jobs_executed = 0L, new_data_jobs_executed = 0L,
  optimization_jobs_executed = 0L,
  decision_state = "mv5n_contract_and_resource_gate_selected_not_executed",
  stringsAsFactors = FALSE
)

outputs <- list(
  "mv05m-selection-criteria-2026-08-10.csv" = criteria,
  "mv05m-candidate-scores-2026-08-10.csv" = scores,
  "mv05m-axis-readiness-2026-08-10.csv" = readiness,
  "mv05m-selected-next-sprint-2026-08-10.csv" = selected,
  "mv05m-training-pair-scope-2026-08-10.csv" = pair_scope,
  "mv05m-evidence-registry-2026-08-10.csv" = evidence,
  "mv05m-production-summary-2026-08-10.csv" = summary
)
for (name in names(outputs)) {
  write_provenance_csv(outputs[[name]], file.path(output_dir, name))
}
message("MV5-M benchmark-gap gate audit artifacts written.")
