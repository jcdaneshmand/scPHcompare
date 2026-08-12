# MV5-AM selection-resistant four-panel synthesis helpers.

.mv05am_true <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  tolower(trimws(as.character(x))) == "true"
}

#' Return the immutable four-configuration robustness sequence.
mv05am_panel_registry_v1 <- function() {
  data.frame(
    contract_id = "mv05am_panel_registry_v1",
    panel_order = 1:4,
    panel_id = c("pc20_vs_pc30", "cosine_chord_vs_euclidean",
                 "nested192_vs_384_cells", "nested256_vs_384_cells"),
    accepted_sprint = c("MV5-Z", "MV5-AD", "MV5-AH", "MV5-AL"),
    accepted_completion = c("d4c332a", "0b32d76", "fbda3eb", "d5c4b7e"),
    changed_factor = c("coordinate_count", "point_metric",
                       "cell_realization_depth", "cell_realization_depth"),
    candidate_state = c("20_coordinates", "cosine_chord", "192_cells",
                        "256_cells"),
    reference_state = c("30_coordinates", "euclidean", "384_cells",
                        "384_cells"),
    configuration_sequence_exhausted_after_panel = c(FALSE, FALSE, FALSE, TRUE),
    stringsAsFactors = FALSE)
}

#' Validate and bind one complete, unsliced accepted result panel.
mv05am_bind_complete_panel_v1 <- function(panel_id, production, macro,
                                           intervals, primary) {
  registry <- mv05am_panel_registry_v1()
  if (length(panel_id) != 1L || !panel_id %in% registry$panel_id ||
      !is.data.frame(production) || nrow(production) != 1L ||
      !is.data.frame(macro) || nrow(macro) != 24L ||
      !is.data.frame(intervals) || nrow(intervals) != 24L ||
      !is.data.frame(primary) || nrow(primary) != 4L) {
    stop("MV5-AM requires one complete 24-estimand, four-test panel.",
         call. = FALSE)
  }
  expected <- c(prediction_groups = 150L, ranking_rows = 282800L,
    outcome_groups = 150L, query_method_rows = 3600L,
    long_query_endpoint_rows = 7200L, macro_estimands = 24L,
    intervals = 24L, primary_tests = 4L, representations = 2L,
    methods = 4L, endpoints = 2L, samples = 90L, studies = 15L,
    tissues = 5L, seeds = 5L)
  if (!all(names(expected) %in% names(production)) ||
      any(vapply(names(expected), function(name) {
        as.integer(production[[name]][[1L]]) != expected[[name]]
      }, logical(1L))) ||
      anyDuplicated(macro$estimand_id) || anyDuplicated(intervals$estimand_id) ||
      !setequal(macro$estimand_id, intervals$estimand_id) ||
      !all(primary$estimand_id %in% intervals$estimand_id) ||
      any(!is.finite(macro$estimate)) ||
      any(!is.finite(intervals$estimate)) ||
      any(!is.finite(intervals$ci_lower)) ||
      any(!is.finite(intervals$ci_upper)) ||
      any(!is.finite(primary$raw_p_value)) ||
      any(!is.finite(primary$holm_p_value))) {
    stop("MV5-AM panel counts, identities, or complete reporting drifted.",
         call. = FALSE)
  }
  prohibited <- grep(
    "clustering_executed|other_configuration.*executed|nested.*executed|method_selection_executed|equivalence_claim_authorized|refit_or_rerank_executed",
    names(production), value = TRUE)
  if (length(prohibited) &&
      any(.mv05am_true(unlist(production[1L, prohibited], use.names = FALSE)))) {
    stop("MV5-AM source panel reports a prohibited operation.", call. = FALSE)
  }
  row <- registry[registry$panel_id == panel_id, , drop = FALSE]
  data.frame(
    contract_id = "mv05am_complete_panel_binding_v1",
    panel_order = row$panel_order,
    panel_id = panel_id,
    accepted_sprint = row$accepted_sprint,
    complete_panel_bound = TRUE,
    macro_estimands = nrow(macro), intervals = nrow(intervals),
    primary_tests = nrow(primary), result_rows_read = 52L,
    result_rows_used_for_continuation_decision = 0L,
    subgroup_or_result_selection_used = FALSE,
    stringsAsFactors = FALSE)
}

#' Add immutable panel identity to a complete result table.
mv05am_stack_panel_v1 <- function(panel_id, table, table_type) {
  registry <- mv05am_panel_registry_v1()
  expected <- c(macro = 24L, interval = 24L, primary = 4L)
  if (!panel_id %in% registry$panel_id || !table_type %in% names(expected) ||
      !is.data.frame(table) || nrow(table) != expected[[table_type]] ||
      !"estimand_id" %in% names(table)) {
    stop("MV5-AM cannot stack an incomplete panel table.", call. = FALSE)
  }
  row <- registry[registry$panel_id == panel_id, , drop = FALSE]
  cbind(data.frame(
    synthesis_contract_id = paste0("mv05am_complete_", table_type, "_panel_v1"),
    panel_order = row$panel_order, panel_id = panel_id,
    accepted_sprint = row$accepted_sprint,
    changed_factor = row$changed_factor,
    candidate_state = row$candidate_state,
    reference_state = row$reference_state,
    stringsAsFactors = FALSE), table, stringsAsFactors = FALSE)
}

#' Structural comparability contract; contains no estimates or p-values.
mv05am_comparability_v1 <- function() {
  registry <- mv05am_panel_registry_v1()
  data.frame(
    contract_id = "mv05am_cross_panel_comparability_v1",
    panel_order = registry$panel_order,
    panel_id = registry$panel_id,
    same_fold_axis = TRUE, same_seed_axis = TRUE,
    same_representation_axis = TRUE, same_method_family_axis = TRUE,
    same_endpoint_axis = TRUE, same_sample_study_tissue_axis = TRUE,
    same_complete_24_estimand_shape = TRUE,
    one_named_factor_contrast = TRUE,
    estimates_directly_poolable_as_one_effect = FALSE,
    reason_not_pooled = "distinct_interventions_with_separate_prespecified_estimands",
    stringsAsFactors = FALSE)
}

#' Preserve the corrected landscape definition through synthesis.
mv05am_landscape_contract_v1 <- function() {
  data.frame(
    contract_id = "mv05am_landscape_definition_v1",
    contract_item = c("finite_intervals", "essential_h0", "level_policy",
      "integration", "dimension_policy", "composite_policy", "grid_policy",
      "normalization_policy"),
    required_state = c("all_finite_intervals", "exclude_infinite_interval",
      "all_consecutive_active_levels", "exact_linear_critical_pair_squared_l2",
      "h0_h1_separate_primary_components",
      "unweighted_euclidean_h0_h1_composite_descriptive_only",
      "no_fixed_uniform_grid", "no_posthoc_scale_or_level_weighting"),
    preserved_in_all_four_panels = TRUE,
    changed_by_synthesis = FALSE,
    stringsAsFactors = FALSE)
}

#' Fixed evidence gaps ordered without consulting result values.
mv05am_evidence_gaps_v1 <- function() {
  data.frame(
    contract_id = "mv05am_evidence_gap_registry_v1",
    prerequisite_order = 1:5,
    gap_id = c("public_landscape_contract", "complete_clustering_pairs",
      "gene_view_execution", "cell_gene_fusion", "new_data_external_validation"),
    observed_state = c(
      "correct_exact_engine_exists_but_general_reference_not_scientific_default",
      "directed_retrieval_pairs_do_not_identify_complete_within_training_matrices",
      "planned_not_executed", "planned_not_executed",
      "deferred_until_existing_data_evidence_gap_requires_it"),
    next_action_eligible = c(TRUE, FALSE, FALSE, FALSE, FALSE),
    eligibility_reason = c(
      "must_make_scientific_definition_and_public_API_unambiguous_before_expansion",
      "requires_separate_prospective_distance_scope_after_public_contract",
      "requires_cell_view_contract_stability_first",
      "requires_both_views_validated_first",
      "existing_data_sequence_not_yet_exhausted"),
    result_value_used_for_priority = FALSE,
    stringsAsFactors = FALSE)
}

#' Mandatory structural continuation criteria.
mv05am_continuation_criteria_v1 <- function() {
  data.frame(
    contract_id = "mv05am_continuation_criteria_v1",
    criterion_order = 1:8,
    criterion_id = c("four_panels_complete", "configuration_sequence_exhausted",
      "selection_firewall", "landscape_definition_preserved",
      "no_effect_pooling", "clustering_nonidentifiable",
      "public_contract_gap_precedes_expansion", "stop_before_new_execution"),
    mandatory = TRUE, prospective = TRUE,
    stringsAsFactors = FALSE)
}

#' Make the sole structural next-action decision.
mv05am_decide_v1 <- function(registry, bindings, comparability, gaps,
                             criterion_pass) {
  forbidden <- c("estimate", "ci_lower", "ci_upper", "raw_p_value",
                 "holm_p_value", "tissue", "winner", "rank")
  inputs <- list(registry, bindings, comparability, gaps, criterion_pass)
  if (any(vapply(inputs, function(x) any(forbidden %in% names(x)), logical(1L)))) {
    stop("MV5-AM decision interface received a prohibited result input.",
         call. = FALSE)
  }
  criteria <- mv05am_continuation_criteria_v1()
  if (!is.data.frame(registry) || nrow(registry) != 4L ||
      !identical(registry$panel_order, 1:4) ||
      !isTRUE(registry$configuration_sequence_exhausted_after_panel[[4L]]) ||
      !is.data.frame(bindings) || nrow(bindings) != 4L ||
      !all(.mv05am_true(bindings$complete_panel_bound)) ||
      any(bindings$result_rows_used_for_continuation_decision != 0L) ||
      !is.data.frame(comparability) || nrow(comparability) != 4L ||
      any(.mv05am_true(comparability$estimates_directly_poolable_as_one_effect)) ||
      !is.data.frame(gaps) || nrow(gaps) != 5L ||
      !identical(gaps$gap_id[[1L]], "public_landscape_contract") ||
      !identical(which(.mv05am_true(gaps$next_action_eligible)), 1L) ||
      !is.data.frame(criterion_pass) ||
      !identical(criterion_pass$criterion_id, criteria$criterion_id) ||
      !all(.mv05am_true(criterion_pass$passed))) {
    stop("MV5-AM continuation decision cannot be authorized.", call. = FALSE)
  }
  data.frame(
    contract_id = "mv05am_continuation_decision_v1",
    decision = "authorize_public_landscape_contract_reconciliation_prefreeze_only",
    authorized_next_sprint = "MV5-AN",
    decision_basis = "structural_prerequisite_order_only",
    numerical_result_rows_consumed_by_decision = 0L,
    favorable_result_selection_used = FALSE,
    fifth_robustness_configuration_authorized = FALSE,
    new_calculation_authorized = FALSE, clustering_authorized = FALSE,
    gene_or_fusion_authorized = FALSE, new_data_authorized = FALSE,
    rust_or_optimization_authorized = FALSE,
    public_default_change_authorized = FALSE,
    prefreeze_only = TRUE, stringsAsFactors = FALSE)
}
