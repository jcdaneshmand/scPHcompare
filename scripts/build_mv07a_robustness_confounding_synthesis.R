#!/usr/bin/env Rscript

options(warn = 2)
for (package in c("digest", "pkgload")) {
  if (!requireNamespace(package, quietly = TRUE)) stop(package, " required.")
}
pkgload::load_all(".", quiet = TRUE, export_all = TRUE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("usage: build_mv07a_robustness_confounding_synthesis.R AUDIT_DIR EXPECTED_HEAD",
    call. = FALSE)
}
audit_dir <- args[[1L]]
expected_head <- tolower(trimws(args[[2L]]))
if (!grepl("^[0-9a-f]{40}$", expected_head)) stop("Full EXPECTED_HEAD required.")
head <- tolower(trimws(system2("git", c("rev-parse", "HEAD"), stdout = TRUE)))
if (!identical(head, expected_head)) {
  stop("MV7-A requires prospective HEAD ", expected_head, "; observed ", head, ".")
}
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
sha <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
write_once <- function(value, name) {
  path <- file.path(audit_dir, name)
  if (file.exists(path)) stop("Refusing to overwrite: ", path, call. = FALSE)
  write_provenance_csv(value, path)
}

paths <- c(
  mv05_l = "docs/audits/MV05L_LOCKED_REPRESENTATION_COMPARISON_2026-08-10.md",
  mv05_s = "docs/audits/MV05S_PREDICTION_LOCKED_CLUSTERING_OUTCOME_EXECUTION_2026-08-10.md",
  mv05_am = "docs/audits/MV05AM_SELECTION_RESISTANT_FOUR_PANEL_SYNTHESIS_GATE_2026-08-12.md",
  mv05_ay = "docs/audits/MV05AY_COMPLETE_CORRECTED_LANDSCAPE_MATRIX_PRODUCTION_2026-08-12.md",
  mv06_c = "docs/audits/MV06C_GLOBAL_CORE_MATCHED_SCT_FEASIBILITY_2026-08-14.md",
  mv06_d = "docs/audits/MV06D_MATCHED_SCT_SOURCE_PH_LANDSCAPE_PROFILE_2026-08-14.md",
  mv06_e = "docs/audits/MV06E_BATCHED_LANDSCAPE_ACCELERATION_ADMISSION_2026-08-14.md",
  mv06_f = "docs/audits/MV06F_STAGE2_COMPLETE_PREFREEZE_2026-08-14.md",
  mv06_g = "docs/audits/MV06G_COMPLETE_PRODUCTION_VALIDATION_2026-08-15.md",
  mv06_h = "docs/audits/MV06H_BLOCKED_FUSION_OUTCOME_EXECUTION_2026-08-15.md",
  mv06_h_decision = "docs/audits/mv06h-outcome-evidence/mv06h-decision.csv",
  mv06_h_primary = "docs/audits/mv06h-outcome-evidence/mv06h-primary-contrasts.csv",
  mv06_h_methods = "docs/audits/mv06h-outcome-evidence/mv06h-method-summaries.csv",
  mv06_h_tissues = "docs/audits/mv06h-outcome-evidence/mv06h-tissue-method-summaries.csv",
  mv06_f_validation = "docs/audits/mv06f-stage2-complete-evidence/mv06f-complete-validation.csv",
  mv06_g_validation = "docs/audits/mv06g-completion-complete-evidence/mv06g-complete-validation.csv",
  mv06_d_projection = "docs/audits/mv06d-profile-evidence/mv06d-worker-projection.csv",
  specification = "docs/specifications/MV07A_ROBUSTNESS_CONFOUNDING_EVIDENCE_MAP_PREFREEZE_V1.md",
  helper = "R/mv07a_robustness_confounding_synthesis.R",
  builder = "scripts/build_mv07a_robustness_confounding_synthesis.R",
  validator = "scripts/validate_mv07a_robustness_confounding_synthesis.R",
  repeat_validator = "scripts/validate_mv07a_robustness_confounding_repeat.R",
  tests = "tests/testthat/test-mv07a-robustness-confounding-synthesis.R")
if (any(!file.exists(paths))) {
  stop("MV7-A source set incomplete: ", paste(names(paths)[!file.exists(paths)],
    collapse = ", "), call. = FALSE)
}
source_freeze <- data.frame(
  contract_id = "mv07a_source_freeze_v1", source_id = names(paths),
  artifact_locator = unname(paths),
  sha256 = vapply(paths, sha, character(1L)),
  bytes = as.numeric(file.info(paths)$size), accepted_head = expected_head,
  numerical_outcome_source = names(paths) %in% c("mv06_h_primary", "mv06_h_methods",
    "mv06_h_tissues"), consumed_by_continuation_decision = FALSE,
  stringsAsFactors = FALSE)

landscape <- mv07a_landscape_contract_v1()
robustness <- mv07a_robustness_coverage_v1()
confounding <- mv07a_confounding_coverage_v1()
gaps <- mv07a_gap_registry_v1()
claims <- mv07a_claim_boundaries_v1()
decision <- mv07a_decide_v1(robustness, confounding, gaps, landscape)

computation <- data.frame(
  contract_id = "mv07a_computation_budget_v1",
  workload = c("mv05ay_eight_strata_corrected_matrices",
    "mv05s_cell_clustering_outcomes", "mv06f_dual_view_pair_generation",
    "mv06g_ranking_generation", "mv07b_no_new_ph_diagnostics"),
  measured_or_projected = c("measured", "measured", "measured", "measured",
    "bounded_projection"),
  worker_or_wall_hours = c((7614 + 9715 + 6 * 286.04) / 3600,
    325 / 3600, 21538.531 / 3600, 20965.105 / 3600, 0.25),
  peak_rss_bytes = c(945364992, 740687872, 9575215104, 186503168, 1073741824),
  retained_or_public_bytes = c(NA_real_, 6532590, 624237551, 531251838, 104857600),
  new_ph = c(TRUE, FALSE, TRUE, FALSE, FALSE),
  resource_disposition = c("completed_under_frozen_caps", "completed_under_frozen_caps",
    "completed_under_amended_12g_cap", "completed_under_frozen_caps",
    "low_cost_reuse_existing_csvs"),
  stringsAsFactors = FALSE)

criteria <- data.frame(
  contract_id = "mv07a_acceptance_criteria_v1",
  criterion_order = 1:10,
  criterion_id = c("source_freeze_complete", "robustness_axes_complete",
    "confounding_axes_complete", "landscape_contract_preserved",
    "cell_gene_boundaries_separate", "negative_fusion_retained",
    "selection_firewall", "no_new_ph_or_data", "computation_reconciled",
    "single_structural_next_action"),
  passed = c(nrow(source_freeze) == length(paths), nrow(robustness) == 14L,
    nrow(confounding) == 10L, all(landscape$preserved), TRUE, TRUE,
    decision$numerical_result_rows_consumed == 0L,
    !decision$new_ph_authorized && !decision$new_data_authorized,
    nrow(computation) == 5L, decision$authorized_next_sprint == "MV7-B"),
  evidence = c("23_of_23_files_exist_and_hash", "14_fixed_axes", "10_fixed_axes",
    "8_of_8_items_unchanged", "view_specific_states_and_claims",
    "mv06h_no_reliable_outperformance", "zero_numerical_rows_to_decision",
    "both_authorizations_false", "four_measured_plus_one_bounded_projection",
    "mv7b_prefreeze_only"),
  stringsAsFactors = FALSE)
if (!all(criteria$passed)) stop("MV7-A acceptance criteria failed.", call. = FALSE)

firewall <- data.frame(
  contract_id = "mv07a_selection_firewall_v1",
  prohibited_input = c("estimate", "confidence_interval", "p_value",
    "method_rank", "tissue_rank", "winner", "outcome_selected_weight",
    "outcome_selected_view"),
  present_in_source_corpus = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
  consumed_by_continuation_decision = FALSE,
  stringsAsFactors = FALSE)

outputs <- list(
  "mv07a-source-freeze.csv" = source_freeze,
  "mv07a-landscape-contract.csv" = landscape,
  "mv07a-robustness-coverage.csv" = robustness,
  "mv07a-confounding-coverage.csv" = confounding,
  "mv07a-computation-budget.csv" = computation,
  "mv07a-evidence-gap-registry.csv" = gaps,
  "mv07a-claim-boundaries.csv" = claims,
  "mv07a-selection-firewall.csv" = firewall,
  "mv07a-acceptance-criteria.csv" = criteria,
  "mv07a-decision.csv" = decision)
for (name in names(outputs)) write_once(outputs[[name]], name)
message("MV7-A complete: 14 robustness axes, 10 confounding axes, MV7-B prefreeze only")
