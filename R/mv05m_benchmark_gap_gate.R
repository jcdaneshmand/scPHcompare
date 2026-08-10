# Deterministic, outcome-free benchmark-gap scoring for MV5-M.

.mv05m_criteria <- data.frame(
  criterion_id = c(
    "scientific_value", "reviewer_relevance", "identifiability_validity",
    "artifact_readiness", "resource_feasibility", "outcome_selection_safety"
  ),
  weight = c(3L, 2L, 3L, 2L, 1L, 2L),
  minimum = 0L,
  maximum = 4L,
  direction = "higher_is_better",
  stringsAsFactors = FALSE
)

.mv05m_candidates <- data.frame(
  candidate_id = c(
    "label_free_clustering_contract_gate",
    "retrieval_robustness_sensitivity",
    "external_validation_new_data",
    "technical_mixing_evaluation",
    "integration_method_expansion",
    "optimization_or_rust",
    "gene_topology",
    "cell_gene_fusion"
  ),
  scientific_value = c(4L, 4L, 4L, 4L, 4L, 1L, 3L, 2L),
  reviewer_relevance = c(4L, 3L, 4L, 4L, 4L, 1L, 2L, 2L),
  identifiability_validity = c(3L, 4L, 4L, 1L, 3L, 4L, 2L, 1L),
  artifact_readiness = c(3L, 2L, 0L, 2L, 1L, 3L, 1L, 0L),
  resource_feasibility = c(2L, 3L, 0L, 3L, 0L, 4L, 1L, 1L),
  outcome_selection_safety = c(4L, 3L, 4L, 4L, 3L, 4L, 4L, 2L),
  validity_gate = c(TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE),
  preliminary_disposition = c(
    "eligible_contract_and_resource_gate",
    "eligible_after_clustering_gate",
    "deferred_existing_data_first",
    "blocked_current_sample_level_identifiability",
    "deferred_until_single_integration_baseline_complete",
    "deferred_no_measured_bottleneck_need",
    "blocked_incomplete_view_eligibility",
    "blocked_component_views_not_both_eligible"
  ),
  stringsAsFactors = FALSE
)

mv05m_score_candidates_v1 <- function(
    candidates = .mv05m_candidates, criteria = .mv05m_criteria) {
  score_columns <- criteria$criterion_id
  if (!is.data.frame(candidates) || !is.data.frame(criteria) ||
      !all(score_columns %in% names(candidates)) ||
      anyDuplicated(candidates$candidate_id) ||
      anyDuplicated(criteria$criterion_id)) {
    stop("MV5-M candidate or criterion registry is malformed.", call. = FALSE)
  }
  scores <- as.matrix(candidates[, score_columns, drop = FALSE])
  storage.mode(scores) <- "integer"
  lower <- matrix(criteria$minimum, nrow(scores), ncol(scores), byrow = TRUE)
  upper <- matrix(criteria$maximum, nrow(scores), ncol(scores), byrow = TRUE)
  if (anyNA(scores) || any(scores < lower) || any(scores > upper)) {
    stop("MV5-M criterion scores are outside the frozen 0-to-4 scale.",
         call. = FALSE)
  }
  candidates$weighted_score <- as.integer(scores %*% criteria$weight)
  candidates$selection_eligible <- as.logical(candidates$validity_gate)
  eligible <- candidates[candidates$selection_eligible, , drop = FALSE]
  if (!nrow(eligible)) {
    stop("No MV5-M candidate passes the validity gate.", call. = FALSE)
  }
  eligible <- eligible[order(
    -eligible$weighted_score,
    -eligible$reviewer_relevance,
    -eligible$artifact_readiness,
    eligible$candidate_id,
    method = "radix"
  ), , drop = FALSE]
  selected <- eligible$candidate_id[[1L]]
  candidates$selected_next <- candidates$candidate_id == selected
  candidates$selection_rank <- NA_integer_
  candidates$selection_rank[match(eligible$candidate_id,
                                  candidates$candidate_id)] <-
    seq_len(nrow(eligible))
  candidates <- candidates[order(
    ifelse(candidates$selected_next, 0L, 1L),
    ifelse(is.na(candidates$selection_rank), .Machine$integer.max,
           candidates$selection_rank),
    candidates$candidate_id,
    method = "radix"
  ), , drop = FALSE]
  rownames(candidates) <- NULL
  candidates
}

mv05m_axis_readiness_v1 <- function() {
  data.frame(
    axis_id = c(
      "biological_conservation_retrieval", "technical_mixing",
      "label_free_sample_clustering", "robustness_sensitivity",
      "integration_method_expansion", "gene_topology", "cell_gene_fusion",
      "external_validation_new_data", "optimization_rust"
    ),
    scientific_estimand = c(
      "cross-study tissue MRR and fixed 1-NN; topology minus matched energy; SCT versus integrated DID",
      "residual study or approach separation conditional on biology, reported separately from conservation",
      "training-only stable sample partitions and held-out assignment agreement with biological and technical labels opened post prediction",
      "change in frozen retrieval effects across prespecified cell-count, feature, dimension, and distance perturbations",
      "incremental conservation and technical-mixing behavior across additional inductive integration methods",
      "gene-coexpression topology versus matched gene-correlation baseline",
      "incremental value of transparent cell/gene distance fusion over both eligible components",
      "replication of frozen primary effects in a separately sourced heterogeneous collection",
      "measured runtime or memory reduction at numerical and scientific equivalence"
    ),
    primary_unit = c(
      "held-out biological sample blocked by study",
      "study-pair summaries or cell neighborhoods conditional on biological state",
      "held-out biological sample assigned from training-only sample clusters",
      "held-out biological sample blocked by study",
      "held-out biological sample blocked by study",
      "biological sample with common named genes",
      "biological sample",
      "held-out biological sample blocked by external study",
      "validated computational operation"
    ),
    current_evidence = c(
      "complete_mv05e_mv05k_mv05l",
      "not_executed_and_current_query_training_pairs_exclude_same_study",
      "pilot_only_and_current_production_distance_artifacts_are_query_to_training",
      "five_cell_subsample_seeds_only_no_frozen_sensitivity_grid",
      "one_inductive_seurat_reference_mapping_only",
      "sct_partial_and_integrated_structurally_unavailable_in_mv05c",
      "not_eligible_until_both_views_pass",
      "no_external_collection_accepted",
      "landscape_engine_profiled_and_current_rust_gate_negative"
    ),
    dominant_validity_risk = c(
      "none_beyond_stated_existing_data_scope",
      "eligible_studies_are_nested_within_tissue_and_loso_query_pairs_have_no_same_study_comparator",
      "training_matrices_missing_and_out_of_sample_assignment_must_be_frozen_before_labels",
      "post_outcome_sensitivity_cannot_rescue_or_replace_null_primary_results",
      "large_method_search_and_native_spaces_risk_unstructured_multiplicity",
      "fold_specific_feature_eligibility_and_query_active_gene_sets",
      "fusion_weight_selection_would_be_outcome_sensitive",
      "dataset_shift_and acquisition_or harmonization burden",
      "optimization_without_a_current_bottleneck_would_displace_scientific_work"
    ),
    missing_computation = c(
      "none_for_frozen_retrieval_axis",
      "identifiable conditional technical endpoint and compatible pair or neighborhood artifacts",
      "262675 training_training pairs per representation across folds and seeds plus frozen clustering prediction machinery",
      "prospective sensitivity grid and bounded cache_reuse plan",
      "method shortlist native induction contracts and resource gates",
      "eligible fold_complete gene views and matched baselines",
      "component eligibility normalization and outcome_closed weight rule",
      "data discovery harmonization and preregistered external protocol",
      "fresh profile demonstrating threshold crossing"
    ),
    disposition = c(
      "complete_confirmatory_with_null_negative_result",
      "blocked_pending_identifiability_design",
      "selected_for_contract_and_resource_gate_only",
      "next_after_clustering_gate_if_still_warranted",
      "deferred",
      "blocked",
      "blocked",
      "deferred_existing_data_first",
      "deferred"
    ),
    stringsAsFactors = FALSE
  )
}

mv05m_selected_sprint_v1 <- function(scores = mv05m_score_candidates_v1()) {
  selected <- scores[scores$selected_next, , drop = FALSE]
  if (nrow(selected) != 1L ||
      selected$candidate_id != "label_free_clustering_contract_gate") {
    stop("MV5-M deterministic selection does not resolve to MV5-N.",
         call. = FALSE)
  }
  data.frame(
    contract_id = "mv05m_selected_next_sprint_v1",
    selected_candidate_id = selected$candidate_id,
    selected_weighted_score = selected$weighted_score,
    next_sprint_id = "MV5-N",
    next_sprint_name = "label_closed_clustering_contract_and_complete_matrix_resource_gate",
    permitted_work = paste(
      "freeze training-only PAM and held-out-medoid assignment;",
      "freeze k=2:10 five-seed one-SE stability;",
      "inventory and profile missing training-training matrices;",
      "validate deterministic synthetic and bounded real label-closed admission"
    ),
    prohibited_work = paste(
      "no full matrix production; no biological or technical label opening;",
      "no clustering outcome; no spectral promotion; no method or tissue selection"
    ),
    exact_training_pairs_per_representation = 262675L,
    exact_h0_h1_rows_per_representation = 525350L,
    existing_query_training_pairs_per_representation = 35350L,
    landscape_only_sct_lower_bound_hours = 8.655,
    landscape_only_integrated_lower_bound_hours = 4.280,
    lower_bound_scope =
      "linear_projection_from_accepted_component_row_throughput_excludes_baselines_validation_io_and_repeat",
    biological_outcomes_computed = FALSE,
    tissue_results_consulted_for_selection = FALSE,
    stringsAsFactors = FALSE
  )
}
