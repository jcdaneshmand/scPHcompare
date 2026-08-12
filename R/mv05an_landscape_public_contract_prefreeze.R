# MV5-AN public landscape contract reconciliation helpers.

mv05an_classify_landscape_function_v1 <- function(file, function_name) {
  if (length(file) != 1L || length(function_name) != 1L) {
    stop("MV5-AN classification requires scalar file and function names.")
  }
  class <- if (file == "landscape_reference.R") {
    "corrected_exact_or_error_controlled_engine"
  } else if (file == "landscape_convergence.R") {
    "internal_diagnostic"
  } else if (file == "cross_iteration_functions.R") {
    "legacy_curve_approximation"
  } else if (file == "PH_PostProcessing_andAnalysis.R") {
    "legacy_grid_artifact_workflow"
  } else if (file == "mv05am_four_panel_synthesis_gate.R") {
    "audit_contract_only"
  } else if (file == "landscape_contract.R") {
    if (function_name %in% c("validate_landscape_grid",
        "validate_landscape_levels", "canonicalize_landscape_matrix",
        "finite_landscape_diagram", "normalize_landscape_grids",
        "normalize_landscape_level_sets", "trapezoidal_l2")) {
      "compatible_validation_or_numeric_primitive"
    } else if (function_name == "landscape_specification_registry") {
      "versioned_semantic_registry"
    } else {
      "versioned_fixed_grid_contract"
    }
  } else {
    "unclassified"
  }
  data.frame(file = file, function_name = function_name,
    pathway_class = class, ambiguous = class == "unclassified",
    public_export = FALSE,
    stringsAsFactors = FALSE)
}

mv05an_target_contract_v1 <- function() {
  data.frame(
    contract_id = "mv05an_target_public_landscape_contract_v1",
    item_order = 1:12,
    contract_item = c("diagram_input", "finite_interval_policy",
      "essential_h0_policy", "level_policy", "dimension_policy",
      "distance_integrand", "integration_policy", "adaptive_fallback",
      "composite_policy", "grid_and_cap_policy", "provenance",
      "failure_policy"),
    required_state = c("typed_dimension_birth_death_with_validated_finite_values",
      "all_positive_persistence_finite_intervals", "exclude_infinite_h0",
      "all_consecutive_active_levels_zero_pad_missing_depth",
      "compute_and_return_h0_h1_separately",
      "sum_squared_level_differences_over_filtration",
      "exact_linear_critical_pair_segments_when_within_explicit_guard",
      "partitioned_error_controlled_quadrature_with_refinement_and_hard_failure",
      "unweighted_euclidean_h0_h1_norm_secondary_descriptive",
      "no_universal_uniform_grid_or_level_cap",
      "version_method_tolerances_interval_counts_hashes_runtime_and_error",
      "no_silent_fallback_partial_result_or_nonfinite_distance"),
    immutable_scientific_semantics = TRUE,
    stringsAsFactors = FALSE)
}

mv05an_public_api_decision_v1 <- function() {
  data.frame(
    contract_id = "mv05an_public_api_decision_v1",
    api_id = c("persistence_landscape_distance",
               "persistence_landscape_distance_matrix"),
    export_in_later_sprint = TRUE,
    versioned_result_class = c("scph_landscape_distance_v1",
                               "scph_landscape_distance_matrix_v1"),
    scientific_contract = "full_l2_error_controlled_v1",
    legacy_function_redirected = FALSE,
    legacy_artifact_overwritten = FALSE,
    workflow_default_changed = FALSE,
    reason = c("safe_pairwise_public_entrypoint_without_legacy_shape_collision",
      "complete_symmetric_matrix_with_per_dimension_and_provenance_sidecars"),
    stringsAsFactors = FALSE)
}

mv05an_migration_plan_v1 <- function() {
  data.frame(
    contract_id = "mv05an_migration_plan_v1",
    migration_order = 1:6,
    stage_id = c("add_versioned_pair_api", "add_versioned_matrix_api",
      "freeze_legacy_artifact_reader", "add_explicit_legacy_mode",
      "evaluate_workflow_default_migration", "retire_ambiguous_artifacts"),
    authorized_in_mv05ao = c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE),
    behavior_change = c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE),
    prerequisite = c("MV5_AN_complete", "pair_api_validated",
      "legacy_schema_manifest", "versioned_APIs_available",
      "realistic_compatibility_and_resource_evidence",
      "documented_major_or_deprecation_release"),
    stringsAsFactors = FALSE)
}

mv05an_implementation_scope_v1 <- function() {
  data.frame(
    contract_id = "mv05an_mv05ao_scope_v1",
    scope_item = c("new_pair_api", "new_matrix_api", "result_classes",
      "provenance_schema", "analytic_oracles", "exact_adaptive_agreement",
      "legacy_schema_detection", "explicit_legacy_mode", "documentation",
      "resource_smoke", "independent_validation", "repeat_and_resume"),
    required = TRUE,
    default_change_authorized = FALSE,
    legacy_artifact_rewrite_authorized = FALSE,
    stringsAsFactors = FALSE)
}

mv05an_decide_v1 <- function(function_inventory, entrypoints, artifacts,
                             target, migration, validation_pass) {
  if (!is.data.frame(function_inventory) || nrow(function_inventory) < 40L ||
      any(function_inventory$ambiguous) || !is.data.frame(entrypoints) ||
      nrow(entrypoints) != 6L || !is.data.frame(artifacts) ||
      nrow(artifacts) != 6L || !is.data.frame(target) || nrow(target) != 12L ||
      !all(target$immutable_scientific_semantics) ||
      !is.data.frame(migration) || nrow(migration) != 6L ||
      !identical(which(migration$authorized_in_mv05ao), 1:4) ||
      !is.data.frame(validation_pass) || !all(validation_pass$passed)) {
    stop("MV5-AN cannot authorize implementation from incomplete inventory.",
         call. = FALSE)
  }
  data.frame(
    contract_id = "mv05an_continuation_decision_v1",
    decision = "authorize_mv05ao_versioned_public_landscape_api_implementation",
    authorized_next_sprint = "MV5-AO",
    new_versioned_api_authorized = TRUE,
    current_workflow_default_change_authorized = FALSE,
    existing_function_redirection_authorized = FALSE,
    legacy_artifact_overwrite_authorized = FALSE,
    new_scientific_calculation_on_project_data_authorized = FALSE,
    clustering_gene_fusion_new_data_rust_claims_authorized = FALSE,
    stringsAsFactors = FALSE)
}
