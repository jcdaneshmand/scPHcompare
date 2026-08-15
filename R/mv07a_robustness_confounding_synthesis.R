# MV7-A selection-resistant robustness and confounding evidence-map helpers.

.mv07a_true <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  tolower(trimws(as.character(x))) == "true"
}

#' Return the dissertation-aligned landscape contract carried into MV-07.
mv07a_landscape_contract_v1 <- function() {
  data.frame(
    contract_id = "mv07a_landscape_contract_v1",
    item = c("finite_intervals", "essential_h0", "level_policy",
      "integration", "dimension_policy", "grid_policy", "level_cap_policy",
      "composite_policy"),
    required_state = c("all_finite_positive_persistence_intervals",
      "exclude_infinite_interval", "all_consecutive_active_levels",
      "exact_or_error_controlled_squared_l2_on_dimension_support",
      "h0_h1_separate", "no_universal_fixed_grid", "no_universal_level_cap",
      "unweighted_h0_h1_euclidean_composite_descriptive_only"),
    preserved = TRUE, changed_by_mv07a = FALSE,
    stringsAsFactors = FALSE)
}

#' Return the complete, result-independent robustness coverage registry.
mv07a_robustness_coverage_v1 <- function() {
  data.frame(
    contract_id = "mv07a_robustness_coverage_v1",
    axis_order = 1:14,
    axis_id = c("cell_realization_depth", "gene_panel_size", "cell_pca_dimension",
      "cell_point_metric", "gene_point_metric", "sampling_seed",
      "homology_dimension", "landscape_grid_or_level_cap", "integration_state",
      "clustering_algorithm", "outlier_influence", "cell_type_composition",
      "filtration_threshold", "external_dataset"),
    cell_view_state = c("complete_192_256_384", "not_applicable",
      "complete_20_30", "complete_euclidean_cosine_chord", "not_applicable",
      "complete_five_seeds", "complete_h0_h1_separate", "not_applicable",
      "complete_sct_integrated", "complete_pam_average_linkage_secondary",
      "not_run", "not_run", "not_applicable_complete_vr_primary", "not_run"),
    gene_view_state = c("fixed_384_not_perturbed", "fixed_500_not_perturbed",
      "fixed_30_cell_pcs_not_perturbed", "not_applicable",
      "fixed_correlation_chord_not_perturbed", "complete_five_seeds",
      "complete_h0_h1_separate", "not_applicable", "sct_only",
      "not_run", "not_run", "not_run", "not_applicable_complete_vr_primary",
      "not_run"),
    coverage = c("partial", "gap", "partial", "partial", "gap", "complete",
      "complete", "complete", "partial", "partial", "gap", "unavailable",
      "complete", "deferred"),
    new_ph_required = c(FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE,
      TRUE, FALSE, FALSE, FALSE, FALSE, TRUE),
    interpretation = c(
      "cell-view depth sensitivity complete; gene view remains fixed at 384 cells",
      "500-gene global core is stable across seeds but panel-size sensitivity is unknown",
      "cell PC20 versus PC30 complete; matched gene view uses the fixed 30-PC cell transform",
      "cell Euclidean versus cosine-chord complete; gene metric is a separate geometry",
      "correlation-chord is the only accepted gene geometry",
      "all accepted cell and gene results average the same five fixed seeds",
      "H0 and H1 are retained and reported separately in every corrected analysis",
      "exact or error-controlled integration replaces arbitrary grids and caps",
      "SCT versus integrated is evaluated for cell topology; gene topology is SCT only",
      "PAM and average linkage are secondary clustering sensitivities for cell representations",
      "study/sample influence has not been quantified across the final separate views",
      "no harmonized cell-type annotations exist in the accepted metadata",
      "the primary contract uses complete Vietoris-Rips rather than a tuned sparse threshold",
      "existing-data sequence must identify a named unresolved estimand before new data"),
    stringsAsFactors = FALSE)
}

#' Return the result-independent confounding coverage registry.
mv07a_confounding_coverage_v1 <- function() {
  data.frame(
    contract_id = "mv07a_confounding_coverage_v1",
    axis_order = 1:10,
    axis_id = c("tissue", "study", "sequencing_approach", "retained_cell_count",
      "library_size", "preprocessing_representation", "study_influence",
      "tissue_heterogeneity", "cell_type_composition", "gene_panel_transductivity"),
    observed_design = c("five_tissue_equal_macro", "fifteen_study_loso_and_blocking",
      "approach_recorded_but_partly_nested", "metadata_and_matched_depths_available",
      "not_available_in_accepted_metadata", "sct_and_inductive_integrated_cell_views",
      "not_quantified_for_final_cell_gene_results", "reported_descriptively",
      "annotations_unavailable", "panel_availability_seen_across_existing_samples"),
    coverage = c("complete", "complete", "partial", "partial", "unavailable",
      "complete_cell_only", "gap", "complete_descriptive", "unavailable",
      "complete_boundary"),
    admitted_next_diagnostic = c(FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE,
      TRUE, FALSE, FALSE),
    causal_adjustment_authorized = FALSE,
    note = c(
      "equal-tissue macro prevents large tissues from dominating but does not remove tissue-study nesting",
      "LOSO and study-block resampling are complete; leave-one-study influence remains useful",
      "technology effects are not separately identifiable where approach is nested in study or tissue",
      "association and stratified sensitivity can be assessed without new PH",
      "do not invent a library-size proxy from retained cells",
      "integration comparison is cell-view evidence and cannot be generalized to the gene view",
      "prospectively compute leave-one-study and leave-one-tissue influence",
      "retain all five tissues; do not select favorable tissues",
      "record as unavailable and limit biological mechanism claims",
      "the fixed panel is technically transductive and not external validation"),
    stringsAsFactors = FALSE)
}

#' Return the fixed evidence-gap order without numerical outcomes.
mv07a_gap_registry_v1 <- function() {
  data.frame(
    contract_id = "mv07a_gap_registry_v1", prerequisite_order = 1:6,
    gap_id = c("existing_outcome_influence", "retained_cell_count_association",
      "approach_stratification", "gene_panel_size_sensitivity",
      "gene_metric_sensitivity", "external_validation"),
    next_sprint_eligible = c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE),
    requires_new_ph = c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE),
    decision_basis = c(
      "reuse_locked_sample_level_outcomes_and_block_structure",
      "reuse_authoritative_metadata_without_recomputing_topology",
      "reuse_authoritative_metadata_with_identifiability_limits",
      "defer_until_no_new_ph_diagnostics_and_claim_need_are known",
      "defer_until_no_new_ph_diagnostics_and_claim_need_are known",
      "trigger_only_for_named_estimand_unresolved_by_existing_data"),
    result_value_used_for_priority = FALSE,
    stringsAsFactors = FALSE)
}

#' Return claim boundaries that follow from study design, not favorable results.
mv07a_claim_boundaries_v1 <- function() {
  data.frame(
    contract_id = "mv07a_claim_boundaries_v1", claim_order = 1:9,
    claim_family = c("landscape_definition", "cell_topology", "gene_topology",
      "fusion", "integration", "tissue_generalization", "technology_effect",
      "biological_mechanism", "external_generalization"),
    current_status = c("supported_method_contract", "existing_data_secondary",
      "existing_data_secondary", "negative", "existing_data_negative_or_null",
      "blocked_descriptive", "not_identifiable_causally", "prohibited",
      "prohibited"),
    allowed_language = c(
      "all-active-level exact-or-error-controlled landscapes with separate H0 and H1",
      "cell-view topology was evaluated under fixed blocked existing-data conditions",
      "gene-view topology was evaluated conditionally on a fixed transductive 500-gene panel",
      "equal-weight fusion did not reliably outperform both component views",
      "integration did not improve the prespecified topology increment in this corpus",
      "report complete tissue heterogeneity without tissue selection",
      "describe associations only and state nesting limitations",
      "do not infer pathways cell states or mechanisms from sample topology alone",
      "do not claim performance on unseen datasets"),
    stringsAsFactors = FALSE)
}

#' Decide the next sprint using only structural coverage and gap order.
mv07a_decide_v1 <- function(robustness, confounding, gaps, landscape) {
  forbidden <- c("estimate", "ci_lower", "ci_upper", "raw_p_value",
    "holm_p_value", "mrr", "balanced_accuracy", "winner", "rank")
  inputs <- list(robustness, confounding, gaps, landscape)
  if (any(vapply(inputs, function(x) any(tolower(names(x)) %in% forbidden),
      logical(1L)))) {
    stop("MV7-A decision interface received a prohibited numerical result.", call. = FALSE)
  }
  if (!is.data.frame(robustness) || nrow(robustness) != 14L ||
      !is.data.frame(confounding) || nrow(confounding) != 10L ||
      !is.data.frame(gaps) || nrow(gaps) != 6L ||
      !identical(which(.mv07a_true(gaps$next_sprint_eligible)), 1:3) ||
      any(gaps$requires_new_ph[1:3]) ||
      !is.data.frame(landscape) || nrow(landscape) != 8L ||
      !all(.mv07a_true(landscape$preserved))) {
    stop("MV7-A structural continuation criteria failed.", call. = FALSE)
  }
  data.frame(
    contract_id = "mv07a_continuation_decision_v1",
    decision = "authorize_no_new_ph_confounding_prefreeze_only",
    authorized_next_sprint = "MV7-B",
    selection_basis = "prerequisite_order_and_existing_artifact_reuse_only",
    numerical_result_rows_consumed = 0L,
    new_ph_authorized = FALSE, new_data_authorized = FALSE,
    method_or_weight_selection_authorized = FALSE,
    advanced_fusion_authorized = FALSE, default_change_authorized = FALSE,
    manuscript_claim_promotion_authorized = FALSE,
    author_team_gate_reached = FALSE, prefreeze_only = TRUE,
    stringsAsFactors = FALSE)
}
