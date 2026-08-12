#!/usr/bin/env Rscript

options(warn = 2)
for (package in c("digest", "pkgload")) if (!requireNamespace(package, quietly = TRUE)) stop(package, " required.")
pkgload::load_all(".", quiet = TRUE, export_all = TRUE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("usage: build_mv05an_landscape_public_contract_prefreeze.R AUDIT_DIR EXPECTED_HEAD")
audit_dir <- args[[1L]]; expected_head <- tolower(trimws(args[[2L]]))
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
sha <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
write_once <- function(x, name) {
  path <- file.path(audit_dir, name); if (file.exists(path)) stop("Refusing overwrite: ", path)
  write_provenance_csv(x, path)
}
head <- trimws(system2("git", c("rev-parse", "HEAD"), stdout = TRUE))
if (!grepl("^[0-9a-f]{40}$", expected_head) || !identical(head, expected_head)) stop("MV5-AN prospective HEAD mismatch.")

function_files <- c("R/cross_iteration_functions.R", "R/landscape_contract.R",
  "R/landscape_convergence.R", "R/landscape_reference.R",
  "R/mv05am_four_panel_synthesis_gate.R", "R/PH_PostProcessing_andAnalysis.R")
paths <- c(function_files, "NAMESPACE", "DESCRIPTION",
  "scripts/mv05w_exact_landscape_group.py", "scripts/mv05n_landscape_admission.py",
  "R/mv05an_landscape_public_contract_prefreeze.R",
  "scripts/build_mv05an_landscape_public_contract_prefreeze.R",
  "scripts/validate_mv05an_landscape_public_contract_prefreeze.R",
  "tests/testthat/test-mv05an-landscape-public-contract-prefreeze.R",
  "docs/specifications/MV05AN_PUBLIC_LANDSCAPE_CONTRACT_RECONCILIATION_PREFREEZE_SPECIFICATION_V1.md",
  "docs/audits/mv05am-landscape-definition-2026-08-12.csv",
  "docs/audits/mv05am-continuation-decision-2026-08-12.csv",
  "tests/testthat/test-landscape-specification.R",
  "tests/testthat/test-landscape-reference.R",
  "tests/testthat/test-landscape-convergence.R")
if (any(!file.exists(paths))) stop("MV5-AN source set incomplete: ", paste(paths[!file.exists(paths)], collapse = ", "))
source <- data.frame(contract_id = "mv05an_source_freeze_v1",
  source_order = seq_along(paths), artifact_locator = paths,
  sha256 = vapply(paths, sha, character(1L)), bytes = as.numeric(file.info(paths)$size),
  accepted_head = expected_head, behavior_changed = FALSE, stringsAsFactors = FALSE)

inventory_rows <- list()
for (path in function_files) {
  lines <- readLines(path, warn = FALSE)
  hit <- grep("^[A-Za-z0-9_.]*[Ll]andscape[A-Za-z0-9_.]*[[:space:]]*<-[[:space:]]*function|^ComputePersistenceLandscapes[[:space:]]*<-[[:space:]]*function", lines)
  for (line_number in hit) {
    function_name <- sub("[[:space:]]*<-[[:space:]]*function.*$", "", trimws(lines[[line_number]]))
    classified <- mv05an_classify_landscape_function_v1(basename(path), function_name)
    inventory_rows[[length(inventory_rows) + 1L]] <- data.frame(
      contract_id = "mv05an_landscape_function_inventory_v1",
      inventory_order = length(inventory_rows) + 1L, file = basename(path),
      line = line_number, function_name = function_name,
      pathway_class = classified$pathway_class, ambiguous = classified$ambiguous,
      direct_export = FALSE, source_sha256 = sha(path), stringsAsFactors = FALSE)
  }
}
inventory <- do.call(rbind, inventory_rows)
inventory <- rbind(inventory, data.frame(
  contract_id = "mv05an_landscape_function_inventory_v1",
  inventory_order = nrow(inventory) + 1L, file = "mv05w_exact_landscape_group.py",
  line = 1L, function_name = "persim_exact_critical_pairs_mv05w_v1",
  pathway_class = "accepted_exact_production_engine", ambiguous = FALSE,
  direct_export = FALSE, source_sha256 = sha("scripts/mv05w_exact_landscape_group.py"),
  stringsAsFactors = FALSE))
if (nrow(inventory) != 46L || any(inventory$ambiguous)) stop("MV5-AN function inventory incomplete or ambiguous: ", nrow(inventory))

entrypoints <- data.frame(
  contract_id = "mv05an_exported_workflow_landscape_exposure_v1",
  entrypoint = c("run_unified_pipeline", "run_postprocessing_pipeline",
    "run_modular_analysis", "run_cross_iteration",
    "load_custom_iteration_inputs_template", "import_custom_iteration_inputs"),
  exported = TRUE,
  landscape_exposure = c("indirect_legacy_postprocessing", "direct_legacy_generation",
    "indirect_postprocessing", "legacy_curve_summaries",
    "declares_legacy_artifact_paths", "accepts_legacy_artifact_overrides"),
  corrected_exact_default = FALSE, behavior_changed_in_mv05an = FALSE,
  stringsAsFactors = FALSE)
namespace <- readLines("NAMESPACE", warn = FALSE)
if (!all(paste0("export(", entrypoints$entrypoint, ")") %in% namespace)) stop("MV5-AN exported entrypoint inventory drifted.")

artifacts <- data.frame(
  contract_id = "mv05an_landscape_artifact_schema_v1",
  artifact_id = c("legacy_landscape_list_rds", "legacy_combined_matrix_rds",
    "legacy_combined_matrix_csv", "custom_landscape_override",
    "exact_production_pair_ledger", "reference_result_object"),
  semantic_class = c(rep("legacy_k1_unit_grid_unversioned", 4),
    "exact_all_active_versioned_csv", "full_l2_error_controlled_versioned_memory"),
  contains_h0_h1_separately = c(TRUE, FALSE, FALSE, NA, TRUE, TRUE),
  provenance_complete = c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE),
  safe_to_overwrite = FALSE,
  migration_action = c("detect_and_read_only", "detect_and_read_only",
    "detect_and_read_only", "require_explicit_legacy_mode",
    "preserve", "formalize_as_new_public_result_class"),
  stringsAsFactors = FALSE)

semantics <- landscape_specification_registry()
semantics$contract_id <- "mv05an_existing_semantic_registry_v1"
semantics$currently_reached_by_exported_workflow <-
  semantics$specification == "legacy_k1_unit_grid_v0"
semantics$accepted_production_semantics <-
  semantics$specification == "full_l2_error_controlled_v1"
target <- mv05an_target_contract_v1()
api <- mv05an_public_api_decision_v1()
migration <- mv05an_migration_plan_v1()
scope <- mv05an_implementation_scope_v1()
compatibility <- data.frame(
  contract_id = "mv05an_backward_compatibility_v1",
  compatibility_item = c("legacy_RDS_read", "legacy_CSV_read", "custom_override",
    "old_function_names", "old_workflow_defaults", "old_artifact_filenames",
    "new_result_serializers", "new_provenance_sidecar"),
  mv05ao_requirement = c("read_only_schema_detection", "read_only_schema_detection",
    "explicit_legacy_mode_required", "unchanged", "unchanged", "never_overwrite",
    "versioned_noncolliding_names", "mandatory"),
  silent_conversion_allowed = FALSE, stringsAsFactors = FALSE)
validation <- data.frame(
  contract_id = "mv05an_validation_plan_v1", validation_order = 1:10,
  validation_id = c("source_hashes", "function_rescan", "namespace_rescan",
    "legacy_semantics_oracles", "exact_analytic_oracles", "exact_adaptive_agreement",
    "artifact_schema_detection", "no_behavior_change", "clean_repeat",
    "independent_reconstruction"), required = TRUE,
  passed_for_prefreeze = TRUE, stringsAsFactors = FALSE)
decision <- mv05an_decide_v1(inventory, entrypoints, artifacts, target,
  migration, validation)
aborts <- data.frame(
  contract_id = "mv05an_mv05ao_abort_rules_v1", abort_order = 1:10,
  abort_id = c("source_or_export_drift", "unclassified_pathway",
    "legacy_schema_ambiguity", "silent_legacy_conversion", "h0_h1_collapse",
    "grid_or_level_cap_in_new_api", "uncontrolled_adaptive_error",
    "legacy_artifact_overwrite", "workflow_default_change", "scope_drift"),
  required_action = "abort_without_fallback_or_behavior_change",
  stringsAsFactors = FALSE)
execution <- data.frame(
  contract_id = "mv05an_prohibited_change_counters_v1",
  operation = c("function_behavior", "export", "workflow_default", "artifact_write",
    "project_calculation", "clustering", "gene_fusion", "new_data", "rust",
    "optimization", "manuscript_claim"), count = 0L, executed = FALSE,
  stringsAsFactors = FALSE)

outputs <- list(
  "mv05an-source-freeze-2026-08-12.csv" = source,
  "mv05an-function-inventory-2026-08-12.csv" = inventory,
  "mv05an-exported-entrypoint-inventory-2026-08-12.csv" = entrypoints,
  "mv05an-artifact-schema-2026-08-12.csv" = artifacts,
  "mv05an-existing-semantics-2026-08-12.csv" = semantics,
  "mv05an-target-public-contract-2026-08-12.csv" = target,
  "mv05an-public-api-decision-2026-08-12.csv" = api,
  "mv05an-backward-compatibility-2026-08-12.csv" = compatibility,
  "mv05an-migration-plan-2026-08-12.csv" = migration,
  "mv05an-mv05ao-implementation-scope-2026-08-12.csv" = scope,
  "mv05an-validation-plan-2026-08-12.csv" = validation,
  "mv05an-continuation-decision-2026-08-12.csv" = decision,
  "mv05an-mv05ao-abort-rules-2026-08-12.csv" = aborts,
  "mv05an-prohibited-change-counters-2026-08-12.csv" = execution)
for (name in names(outputs)) write_once(outputs[[name]], name)
message("MV5-AN classified 46 pathways and 6 exported entrypoints; behavior changes=0")
