# MV5-AR opt-in corrected-landscape integration prefreeze helpers.

mv05ar_prefreeze_tables_v1 <- function() {
  controls <- data.frame(
    contract_id = "mv05ar_control_contract_v1",
    field = c("contract_id", "enabled", "max_wall_seconds", "max_pairs",
              "max_rss_bytes", "workers", "existing", "downstream_use",
              "method", "exact_max_intervals", "abs_tol", "rel_tol"),
    value = c("scph_corrected_landscape_workflow_control_v1", "TRUE",
              "required_positive", "required_positive", "default_2147483648",
              "1", "resume_or_fail", "artifacts_only", "auto", "500",
              "1e-8", "1e-8"),
    caller_override = c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE,
                        FALSE, FALSE, FALSE, FALSE),
    stringsAsFactors = FALSE
  )
  boundaries <- data.frame(
    contract_id = "mv05ar_entrypoint_boundary_v1",
    entrypoint = c("run_unified_pipeline", "run_postprocessing_pipeline",
      "process_iteration_calculate_matrices", "run_modular_analysis",
      "run_cross_iteration", "apply_custom_iteration_overrides"),
    later_change = c("add_null_default_pass_through", "add_null_default_orchestrator",
      "none_legacy_unchanged", "none_no_corrected_consumption",
      "none_no_corrected_consumption", "none_legacy_keys_unchanged"),
    default_behavior_changed = FALSE,
    corrected_downstream_consumption = FALSE,
    stringsAsFactors = FALSE
  )
  artifacts <- data.frame(
    contract_id = "mv05ar_artifact_contract_v1",
    order = 1:7,
    artifact = c("resource-plan-v1.csv", "input-manifest-v1.csv",
      "pairs/<pair-id>--<pair-cache-key>.rds", "pair-index-v1.csv",
      "distance-matrix-v1.rds", "provenance-v1.csv", "completion-v1.csv"),
    write_policy = c("atomic_create", "atomic_create", "atomic_sharded_create",
      "atomic_create_after_pairs", "atomic_create_after_pairs", "atomic_create",
      "atomic_create_last"),
    overwrite_allowed = FALSE,
    legacy_filename_reused = FALSE,
    stringsAsFactors = FALSE
  )
  resources <- data.frame(
    contract_id = "mv05ar_resource_admission_v1",
    policy = c("exact_exact_pair_seconds", "adaptive_pair_seconds",
      "iteration_overhead_seconds", "minimum_rss_bytes", "default_rss_bytes",
      "observed_h0_min", "observed_h0_max", "observed_h1_min",
      "observed_h1_max", "workers"),
    value = c(30, 240, 30, 1610612736, 2147483648, 383, 499, 79, 2802, 1),
    outside_envelope_action = c(rep("not_applicable", 5),
                                rep("profiling_required", 4), "hard_fail"),
    stringsAsFactors = FALSE
  )
  stages <- data.frame(
    contract_id = "mv05ar_migration_stage_v1",
    stage = 1:5,
    stage_id = c("prefreeze_additive_artifacts", "implement_additive_opt_in",
      "evaluate_realistic_workflow_smoke", "prefreeze_downstream_consumption",
      "consider_default_migration"),
    authorized_now = c(TRUE, FALSE, FALSE, FALSE, FALSE),
    changes_workflow_behavior = c(FALSE, TRUE, FALSE, TRUE, TRUE),
    prerequisite = c("MV5_AP_R1", "MV5_AR_complete", "implementation_validated",
      "smoke_accepted", "major_release_and_author_decision"),
    stringsAsFactors = FALSE
  )
  list(controls = controls, boundaries = boundaries, artifacts = artifacts,
       resources = resources, stages = stages)
}

mv05ar_decide_v1 <- function(pathways_reconciled, schemas_unambiguous,
                             legacy_coexists, resource_policy_bounded,
                             default_changes, behavior_changes,
                             artifact_writes, project_calculations) {
  flags <- c(pathways_reconciled, schemas_unambiguous, legacy_coexists,
             resource_policy_bounded)
  counts <- c(default_changes, behavior_changes, artifact_writes,
              project_calculations)
  if (!is.logical(flags) || anyNA(flags) || !is.numeric(counts) ||
      anyNA(counts) || any(counts < 0)) {
    stop("Invalid MV5-AR decision inputs.", call. = FALSE)
  }
  passed <- all(flags) && all(counts == 0)
  data.frame(
    contract_id = "mv05ar_continuation_decision_v1",
    decision = if (passed) {
      "authorize_additive_opt_in_artifact_implementation_only"
    } else "stop_and_reconcile_integration_prefreeze",
    additive_implementation_authorized = passed,
    corrected_downstream_consumption_authorized = FALSE,
    workflow_default_change_authorized = FALSE,
    legacy_artifact_rewrite_authorized = FALSE,
    optimization_authorized = FALSE,
    next_sprint = if (passed) "MV5-AS" else "none",
    stringsAsFactors = FALSE
  )
}
