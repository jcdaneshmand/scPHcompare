#!/usr/bin/env Rscript

options(warn = 2)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required.")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("usage: validate_mv05am_four_panel_synthesis_gate.R AUDIT_DIR OUTPUT_CSV")
}
audit_dir <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
output <- args[[2L]]
readc <- function(name) read.csv(file.path(audit_dir, name),
  stringsAsFactors = FALSE, check.names = FALSE)
sha <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
truth <- function(x) tolower(as.character(x)) == "true"
checks <- list(); add <- function(id, passed, evidence) {
  checks[[length(checks) + 1L]] <<- data.frame(
    contract_id = "mv05am_independent_validation_v1",
    validation_id = id, passed = isTRUE(passed), evidence = evidence,
    production_synthesis_or_decision_helper_called = FALSE,
    new_calculation_executed = FALSE, labels_opened = FALSE,
    outcomes_recomputed = FALSE, stringsAsFactors = FALSE)
}

source <- readc("mv05am-source-freeze-2026-08-12.csv")
source_ok <- nrow(source) == 38L && all(file.exists(source$artifact_locator)) &&
  all(vapply(source$artifact_locator, sha, character(1L)) == source$sha256) &&
  !any(truth(source$consumed_by_continuation_helper)) &&
  sum(truth(source$numerical_result_source)) == 12L
add("source_hashes_and_decision_boundary", source_ok,
    paste0(nrow(source), "_sources_12_numerical_zero_decision_consumption"))

registry <- readc("mv05am-panel-registry-2026-08-12.csv")
expected_panels <- c("pc20_vs_pc30", "cosine_chord_vs_euclidean",
  "nested192_vs_384_cells", "nested256_vs_384_cells")
registry_ok <- nrow(registry) == 4L &&
  identical(as.integer(registry$panel_order), 1:4) &&
  identical(registry$panel_id, expected_panels) &&
  !any(truth(registry$configuration_sequence_exhausted_after_panel[1:3])) &&
  truth(registry$configuration_sequence_exhausted_after_panel[[4L]])
add("canonical_four_panel_registry", registry_ok,
    "four_positions_complete_sequence_exhausted_after_four")

binding <- readc("mv05am-complete-panel-binding-2026-08-12.csv")
binding_ok <- nrow(binding) == 4L &&
  identical(binding$panel_id, expected_panels) &&
  all(truth(binding$complete_panel_bound)) &&
  all(binding$macro_estimands == 24L) && all(binding$intervals == 24L) &&
  all(binding$primary_tests == 4L) && all(binding$result_rows_read == 52L) &&
  all(binding$result_rows_used_for_continuation_decision == 0L) &&
  !any(truth(binding$subgroup_or_result_selection_used))
add("complete_unsliced_bindings", binding_ok,
    "4_panels_24_macro_24_intervals_4_tests_each")

macro <- readc("mv05am-complete-macro-synthesis-2026-08-12.csv")
interval <- readc("mv05am-complete-interval-synthesis-2026-08-12.csv")
primary <- readc("mv05am-complete-primary-synthesis-2026-08-12.csv")
shape_ok <- nrow(macro) == 96L && nrow(interval) == 96L && nrow(primary) == 16L &&
  all(table(macro$panel_id) == 24L) && all(table(interval$panel_id) == 24L) &&
  all(table(primary$panel_id) == 4L) &&
  !anyDuplicated(paste(macro$panel_id, macro$estimand_id)) &&
  !anyDuplicated(paste(interval$panel_id, interval$estimand_id)) &&
  !anyDuplicated(paste(primary$panel_id, primary$estimand_id))
add("complete_synthesis_shape", shape_ok, "96_macro_96_intervals_16_primary")

prefix <- c(pc20_vs_pc30 = "mv05z", cosine_chord_vs_euclidean = "mv05ad",
  nested192_vs_384_cells = "mv05ah", nested256_vs_384_cells = "mv05al")
exact <- TRUE; reconstructed <- c(macro = 0L, interval = 0L, primary = 0L)
for (panel_id in expected_panels) {
  p <- prefix[[panel_id]]
  original_macro <- read.csv(paste0("docs/audits/", p,
    "-macro-estimands-2026-08-11.csv"), stringsAsFactors = FALSE, check.names = FALSE)
  original_interval <- read.csv(paste0("docs/audits/", p,
    "-estimand-intervals-2026-08-11.csv"), stringsAsFactors = FALSE, check.names = FALSE)
  original_primary <- read.csv(paste0("docs/audits/", p,
    "-primary-contrasts-2026-08-11.csv"), stringsAsFactors = FALSE, check.names = FALSE)
  observed_macro <- macro[macro$panel_id == panel_id, ]
  observed_interval <- interval[interval$panel_id == panel_id, ]
  observed_primary <- primary[primary$panel_id == panel_id, ]
  match_exact <- function(original, observed, columns) {
    observed <- observed[match(original$estimand_id, observed$estimand_id), ]
    !anyNA(observed$estimand_id) && all(vapply(columns, function(column) {
      if (is.numeric(original[[column]]) && is.numeric(observed[[column]])) {
        identical(as.numeric(original[[column]]), as.numeric(observed[[column]]))
      } else {
        identical(original[[column]], observed[[column]])
      }
    }, logical(1L)))
  }
  exact <- exact &&
    match_exact(original_macro, observed_macro,
      c("estimand_id", "estimand_order", "estimand_type", "estimand_role",
        "endpoint_id", "representation", "family_id", "estimate")) &&
    match_exact(original_interval, observed_interval,
      c("estimand_id", "estimate", "ci_lower", "ci_upper",
        "bootstrap_replicates", "bootstrap_seed", "inference_status")) &&
    match_exact(original_primary, observed_primary,
      c("estimand_id", "estimate", "ci_lower", "ci_upper", "raw_p_value",
        "holm_p_value", "randomization_replicates", "randomization_seed"))
  reconstructed <- reconstructed + c(nrow(original_macro), nrow(original_interval),
                                      nrow(original_primary))
}
add("independent_numerical_reconstruction", exact &&
      identical(unname(reconstructed), c(96L, 96L, 16L)),
    "all_estimates_intervals_pvalues_and_identities_exact_from_accepted_sources")

comparability <- readc("mv05am-cross-panel-comparability-2026-08-12.csv")
comparability_ok <- nrow(comparability) == 4L &&
  all(truth(comparability$same_fold_axis)) &&
  all(truth(comparability$same_seed_axis)) &&
  all(truth(comparability$same_representation_axis)) &&
  all(truth(comparability$same_method_family_axis)) &&
  all(truth(comparability$same_endpoint_axis)) &&
  all(truth(comparability$same_sample_study_tissue_axis)) &&
  all(truth(comparability$same_complete_24_estimand_shape)) &&
  !any(truth(comparability$estimates_directly_poolable_as_one_effect))
add("cross_panel_comparability_without_pooling", comparability_ok,
    "common_axes_distinct_interventions_no_single_effect_pooling")

landscape <- readc("mv05am-landscape-definition-2026-08-12.csv")
states <- c("all_finite_intervals", "exclude_infinite_interval",
  "all_consecutive_active_levels", "exact_linear_critical_pair_squared_l2",
  "h0_h1_separate_primary_components",
  "unweighted_euclidean_h0_h1_composite_descriptive_only",
  "no_fixed_uniform_grid", "no_posthoc_scale_or_level_weighting")
landscape_ok <- nrow(landscape) == 8L && setequal(landscape$required_state, states) &&
  all(truth(landscape$preserved_in_all_four_panels)) &&
  !any(truth(landscape$changed_by_synthesis))
add("landscape_definition_preserved", landscape_ok,
    "all_levels_exact_l2_h0_h1_separate_no_grid_or_cap")

clustering <- readc("mv05am-clustering-identifiability-2026-08-12.csv")
clustering_ok <- nrow(clustering) == 4L &&
  all(truth(clustering$directed_retrieval_pairs_complete)) &&
  !any(truth(clustering$complete_within_training_pair_matrix)) &&
  !any(truth(clustering$clustering_identifiable_from_accepted_rows)) &&
  !any(truth(clustering$clustering_authorized))
add("clustering_nonidentifiability", clustering_ok,
    "4_of_4_directed_only_panels_no_complete_training_matrices")

gaps <- readc("mv05am-evidence-gap-registry-2026-08-12.csv")
gaps_ok <- nrow(gaps) == 5L && identical(gaps$prerequisite_order, 1:5) &&
  gaps$gap_id[[1L]] == "public_landscape_contract" &&
  identical(which(truth(gaps$next_action_eligible)), 1L) &&
  !any(truth(gaps$result_value_used_for_priority))
add("structural_evidence_gap_order", gaps_ok,
    "public_landscape_contract_first_zero_result_value_priority")

criteria <- readc("mv05am-continuation-criteria-2026-08-12.csv")
firewall <- readc("mv05am-selection-firewall-2026-08-12.csv")
firewall_ok <- nrow(criteria) == 8L && all(truth(criteria$mandatory)) &&
  all(truth(criteria$prospective)) && all(truth(criteria$passed)) &&
  nrow(firewall) == 10L && all(truth(firewall$present_in_public_synthesis)) &&
  !any(truth(firewall$consumed_by_continuation_helper)) &&
  !any(truth(firewall$privileged_or_filtered))
add("selection_firewall", firewall_ok,
    "10_result_or_subgroup_inputs_published_but_zero_decision_consumption")

decision <- readc("mv05am-continuation-decision-2026-08-12.csv")
decision_ok <- nrow(decision) == 1L && decision$authorized_next_sprint == "MV5-AN" &&
  decision$numerical_result_rows_consumed_by_decision == 0L &&
  !truth(decision$favorable_result_selection_used) &&
  !truth(decision$fifth_robustness_configuration_authorized) &&
  !truth(decision$new_calculation_authorized) &&
  !truth(decision$clustering_authorized) &&
  !truth(decision$public_default_change_authorized) && truth(decision$prefreeze_only)
add("rule_bound_continuation_decision", decision_ok,
    "MV5_AN_public_contract_prefreeze_only_no_default_change")

execution <- readc("mv05am-prohibited-execution-counters-2026-08-12.csv")
execution_ok <- nrow(execution) == 14L && !any(truth(execution$executed)) &&
  !any(truth(execution$authorized)) && all(execution$count == 0L)
add("zero_prohibited_execution", execution_ok,
    "14_operations_executed_0_authorized_0")

result <- do.call(rbind, checks)
if (nrow(result) != 12L || !all(result$passed)) {
  stop("MV5-AM validation failed: ",
       paste(result$validation_id[!result$passed], collapse = ", "), call. = FALSE)
}
if (file.exists(output)) stop("Refusing to overwrite validation output.")
write.csv(result, output, row.names = FALSE, na = "")
message("MV5-AM independent validation passed: 12 categories; 96/96/16 rows")
