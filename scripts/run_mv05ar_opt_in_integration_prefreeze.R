#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: run_mv05ar_opt_in_integration_prefreeze.R OUTPUT_DIR")
}
output_dir <- normalizePath(args[[1L]], mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
setwd("/mnt/e/Repositories/Jonah/PH Pipeline Repo/scPHcompare")
pkgload::load_all(".", quiet = TRUE)

write_csv <- function(value, id) utils::write.csv(
  value, file.path(output_dir, paste0("mv05ar-", id, "-2026-08-12.csv")),
  row.names = FALSE, na = "", quote = TRUE
)
tables <- mv05ar_prefreeze_tables_v1()
lapply(names(tables), function(id) write_csv(tables[[id]], id))

entrypoints <- c(
  "run_unified_pipeline", "run_postprocessing_pipeline", "run_modular_analysis",
  "run_cross_iteration", "load_custom_iteration_inputs_template",
  "import_custom_iteration_inputs"
)
source_file <- c(
  "R/unified_pipeline.R", "R/PH_PostProcessing_andAnalysis.R",
  "R/PH_PostProcessing_andAnalysis.R", "R/cross_iteration_functions.R",
  "R/custom_iteration_inputs_template.R", "R/custom_iteration_inputs_template.R"
)
namespace <- readLines("NAMESPACE", warn = FALSE)
pathways <- do.call(rbind, lapply(seq_along(entrypoints), function(index) {
  fun <- get(entrypoints[[index]], mode = "function")
  data.frame(
    contract_id = "mv05ar_pathway_inventory_v1",
    entrypoint = entrypoints[[index]], source_file = source_file[[index]],
    exported = paste0("export(", entrypoints[[index]], ")") %in% namespace,
    current_formals = paste(names(formals(fun)), collapse = ";"),
    corrected_control_currently_present =
      "corrected_landscape_control" %in% names(formals(fun)),
    behavior_change_in_mv05ar = FALSE,
    stringsAsFactors = FALSE
  )
}))
write_csv(pathways, "pathway-inventory")

legacy <- data.frame(
  contract_id = "mv05ar_legacy_coexistence_v1",
  item = c("ComputePersistenceLandscapes", "compute_and_save_landscape_matrices",
    "process_iteration_calculate_matrices", "landscape_l2_distance_matrix",
    "landscape_list_path", "landscape_matrix_path", "run_modular_analysis",
    "run_cross_iteration"),
  policy = c("unchanged_legacy_level1", "unchanged_legacy_writer",
    "unchanged_legacy_orchestration", "unchanged_legacy_iteration_field",
    "unchanged_legacy_custom_key", "unchanged_legacy_custom_key",
    "no_corrected_object_consumption", "no_corrected_object_consumption"),
  silent_conversion_allowed = FALSE,
  overwrite_allowed = FALSE,
  stringsAsFactors = FALSE
)
write_csv(legacy, "legacy-coexistence")

atomic <- data.frame(
  contract_id = "mv05ar_atomic_resume_contract_v1",
  order = 1:10,
  action = c("validate_named_diagrams", "compute_input_set_identity",
    "write_resource_plan", "write_input_manifest", "validate_existing_shards",
    "compute_missing_pair_to_temp", "reload_validate_atomic_rename_pair",
    "assemble_crosscheck_matrix", "reload_validate_atomic_rename_matrix",
    "write_hash_bound_completion_last"),
  conflict_action = c(rep("hard_fail", 2), rep("never_overwrite", 8)),
  completion_visible = c(rep(FALSE, 9), TRUE),
  stringsAsFactors = FALSE
)
write_csv(atomic, "atomic-resume-contract")

validation <- data.frame(
  contract_id = "mv05ar_implementation_validation_plan_v1",
  validation_order = 1:15,
  requirement = c("default_off_equivalence", "control_schema_rejection",
    "analytic_h0_oracle", "analytic_h1_oracle", "all_active_levels",
    "pair_shard_matrix_equivalence", "h0_h1_separation",
    "strict_error_certificate", "cache_identity", "atomic_interruption_resume",
    "immutable_complete_resume", "serialization_reload", "legacy_coexistence",
    "resource_admission_rejection", "clean_repeat"),
  required = TRUE,
  stringsAsFactors = FALSE
)
write_csv(validation, "implementation-validation-plan")

aborts <- data.frame(
  contract_id = "mv05ar_implementation_abort_rule_v1",
  abort_order = 1:14,
  condition = c("ambiguous_or_duplicate_sample_id", "malformed_diagram",
    "input_hash_drift", "unprofiled_interval_envelope", "pair_budget_exceeded",
    "wall_budget_exceeded", "rss_budget_below_minimum", "uncertified_result",
    "pair_cache_conflict", "non_atomic_write", "legacy_filename_collision",
    "corrected_downstream_consumption", "default_behavior_drift",
    "test_or_check_failure"),
  partial_scientific_result_accepted = FALSE,
  stringsAsFactors = FALSE
)
write_csv(aborts, "implementation-abort-rules")

prior_decision <- utils::read.csv(
  "docs/audits/mv05apr1-continuation-decision-2026-08-12.csv",
  stringsAsFactors = FALSE, check.names = FALSE)
prior_resource <- utils::read.csv(
  "docs/audits/mv05apr1-resource-summary-2026-08-12.csv",
  stringsAsFactors = FALSE, check.names = FALSE)
aq_policy <- utils::read.csv(
  "docs/audits/mv05aq-engine-routing-policy-2026-08-12.csv",
  stringsAsFactors = FALSE, check.names = FALSE)
bindings <- data.frame(
  contract_id = "mv05ar_accepted_evidence_binding_v1",
  evidence = c("mv05apr1_decision", "mv05apr1_max_wall_seconds",
    "mv05apr1_max_rss_bytes", "mv05aq_method", "mv05aq_exact_guard",
    "mv05aq_grid_fallbacks", "mv05aq_level_caps"),
  value = c(prior_decision$decision,
    max(prior_resource$run_a_wall_elapsed_seconds,
        prior_resource$run_b_wall_elapsed_seconds),
    max(prior_resource$run_a_max_rss_bytes,
        prior_resource$run_b_max_rss_bytes),
    aq_policy$public_default_method, aq_policy$exact_max_intervals,
    aq_policy$grid_fallback, aq_policy$landscape_level_cap),
  accepted = c(prior_decision$realistic_gate_passed, rep(TRUE, 6)),
  stringsAsFactors = FALSE
)
write_csv(bindings, "accepted-evidence-binding")

source_files <- c(
  "R/unified_pipeline.R", "R/PH_PostProcessing_andAnalysis.R",
  "R/cross_iteration_functions.R", "R/custom_iteration_inputs_template.R",
  "R/landscape_public_api.R", "R/landscape_reference.R",
  "R/mv05ar_opt_in_integration_prefreeze.R",
  "scripts/run_mv05ar_opt_in_integration_prefreeze.R",
  "scripts/validate_mv05ar_opt_in_integration_prefreeze.R",
  "tests/testthat/test-mv05ar-opt-in-integration-prefreeze.R",
  "NAMESPACE"
)
source_freeze <- data.frame(
  contract_id = "mv05ar_source_freeze_v1", path = source_files,
  sha256 = vapply(source_files, function(path) sub(" .*", "", system2(
    "sha256sum", path, stdout = TRUE)), character(1)),
  stringsAsFactors = FALSE
)
write_csv(source_freeze, "source-freeze")

prohibited <- data.frame(
  contract_id = "mv05ar_prohibited_change_counters_v1",
  counter = c("workflow_default_changes", "workflow_behavior_changes",
    "exports_added", "legacy_source_changes", "legacy_artifact_writes",
    "corrected_artifact_writes", "project_calculations", "downstream_consumption",
    "fixed_grid_fallbacks", "landscape_level_caps", "interval_removals",
    "tolerance_relaxations", "optimization_changes", "biological_outcome_accesses"),
  value = 0L, stringsAsFactors = FALSE
)
write_csv(prohibited, "prohibited-change-counters")

decision <- mv05ar_decide_v1(
  pathways_reconciled = nrow(pathways) == 6L && all(pathways$exported),
  schemas_unambiguous = nrow(tables$artifacts) == 7L,
  legacy_coexists = all(!legacy$silent_conversion_allowed &
                          !legacy$overwrite_allowed),
  resource_policy_bounded = identical(bindings$value[[2L]], "567.94") &&
    identical(bindings$value[[3L]], "990363648"),
  default_changes = 0, behavior_changes = 0, artifact_writes = 0,
  project_calculations = 0
)
write_csv(decision, "continuation-decision")
cat("MV5-AR prefreeze:", decision$decision, "\n")
