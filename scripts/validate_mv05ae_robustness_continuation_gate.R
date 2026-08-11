#!/usr/bin/env Rscript

options(warn = 2)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest is required.")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("usage: validate_mv05ae_robustness_continuation_gate.R AUDIT_DIR OUTPUT_CSV",
       call. = FALSE)
}
audit_dir <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
output <- args[[2L]]
readc <- function(name) read.csv(file.path(audit_dir, name),
  stringsAsFactors = FALSE, check.names = FALSE)
sha <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
truth <- function(x) tolower(as.character(x)) == "true"
checks <- list(); add <- function(id, passed, evidence) {
  checks[[length(checks) + 1L]] <<- data.frame(
    contract_id = "mv05ae_independent_validation_v1",
    validation_id = id, passed = isTRUE(passed), evidence = evidence,
    decision_helper_called = FALSE, new_calculation_executed = FALSE,
    labels_opened = FALSE, outcomes_computed = FALSE,
    stringsAsFactors = FALSE)
}

source <- readc("mv05ae-source-freeze-2026-08-11.csv")
source_ok <- nrow(source) == 30L && all(file.exists(source$artifact_locator)) &&
  all(vapply(source$artifact_locator, sha, character(1L)) == source$sha256)
add("source_hashes", source_ok, "30_of_30_exact")

order <- readc("mv05ae-configuration-order-2026-08-11.csv")
expected_ids <- c("cells384_pc20_euclidean_v1",
  "cells384_pc30_cosine_chord_v1", "nested_cells_192_pc30_euclidean_v1",
  "nested_cells_256_pc30_euclidean_v1")
order_ok <- nrow(order) == 4L && identical(order$configuration_id, expected_ids) &&
  identical(as.integer(order$configuration_order), 1:4) &&
  identical(order$position_state,
            c("complete", "complete", "next_eligible", "later_closed"))
add("canonical_order", order_ok, "pc20_cosine_complete_nested192_next_nested256_closed")

evidence <- readc("mv05ae-complete-evidence-binding-2026-08-11.csv")
evidence_ok <- nrow(evidence) == 2L &&
  setequal(evidence$analysis_id, c("pc20_vs_pc30", "cosine_chord_vs_euclidean")) &&
  all(evidence$macro_estimands == 24L) && all(evidence$intervals == 24L) &&
  all(evidence$primary_tests == 4L) && all(truth(evidence$complete_evidence_bound)) &&
  !any(truth(evidence$cell_count_sensitivity_answered)) &&
  !any(truth(evidence$subgroup_or_method_selection_used))
add("complete_unsliced_evidence", evidence_ok,
    "two_complete_24_estimand_4_test_panels")

criteria <- readc("mv05ae-continuation-criteria-2026-08-11.csv")
firewall <- readc("mv05ae-selection-firewall-2026-08-11.csv")
criteria_ok <- nrow(criteria) == 9L && identical(criteria$criterion_order, 1:9) &&
  all(truth(criteria$mandatory)) && all(truth(criteria$frozen_before_new_calculation)) &&
  all(truth(criteria$passed)) && nrow(firewall) == 9L &&
  !any(truth(firewall$consumed_by_decision_helper)) &&
  all(truth(firewall$complete_prior_panels_bound))
add("selection_firewall", criteria_ok,
    "9_required_criteria_9_prohibited_selection_inputs")

decision <- readc("mv05ae-continuation-decision-2026-08-11.csv")
decision_ok <- nrow(decision) == 1L &&
  decision$authorized_configuration_id == "nested_cells_192_pc30_euclidean_v1" &&
  decision$authorized_groups == 150L &&
  !truth(decision$favorable_subgroup_or_estimate_selection_used) &&
  !truth(decision$nested_192_calculation_executed) &&
  !truth(decision$labels_opened) && !truth(decision$rankings_computed) &&
  !truth(decision$outcomes_computed) && !truth(decision$clustering_authorized) &&
  !truth(decision$nested_256_authorized)
add("single_configuration_decision", decision_ok,
    "nested192_only_calculation_0_nested256_0")

queue <- readc("mv05ae-execution-queue-2026-08-11.csv")
queue_ok <- nrow(queue) == 150L && !anyDuplicated(queue$robustness_group_id) &&
  length(unique(queue$fold_id)) == 15L && length(unique(queue$seed)) == 5L &&
  setequal(queue$representation, c("sct_whole", "inductive_integrated")) &&
  all(queue$configuration_id == "nested_cells_192_pc30_euclidean_v1") &&
  all(queue$cells == 192L) && all(queue$coordinates == 30L) &&
  all(queue$point_metric == "euclidean") &&
  all(truth(queue$later_label_closed_calculation_authorized)) &&
  !any(truth(queue$labels_opened)) && !any(truth(queue$rankings_computed)) &&
  !any(truth(queue$outcomes_computed)) && !any(truth(queue$execution_completed))
add("queue_axes", queue_ok, "150_groups_15_folds_5_seeds_2_representations")

scope <- readc("mv05ae-execution-scope-2026-08-11.csv")
scope_expected <- c(groups = 150L, views = 13500L, biological_pairs = 70700L,
  landscape_request_rows = 141400L, landscape_subchunks = 720L,
  energy_request_rows = 70700L, assembled_method_rows = 282800L)
scope_ok <- nrow(scope) == 1L && all(vapply(names(scope_expected), function(name) {
  as.integer(scope[[name]]) == scope_expected[[name]]
}, logical(1L))) && scope$cells == 192L && scope$coordinates == 30L &&
  scope$point_metric == "euclidean" && !truth(scope$execution_completed)
add("exact_scope", scope_ok,
    "13500_views_70700_pairs_141400_landscapes_70700_energy_282800_methods")

resource <- readc("mv05ae-resource-envelope-2026-08-11.csv")
resource_ok <- nrow(resource) == 1L && resource$admission_units == 6L &&
  resource$later_max_workers == 1L && resource$later_group_cap_seconds == 600 &&
  resource$later_group_rss_cap_bytes == 4294967296 &&
  resource$later_configuration_cap_worker_hours == 6 &&
  resource$later_storage_cap_bytes == 4294967296 && truth(resource$resource_gate_pass)
add("resource_envelope", resource_ok, "6_real_admissions_two_full_precedents_caps_pass")

validation <- readc("mv05ae-validation-plan-2026-08-11.csv")
abort <- readc("mv05ae-abort-rules-2026-08-11.csv")
guards_ok <- nrow(validation) == 12L && all(truth(validation$required)) &&
  nrow(abort) == 10L && all(abort$required_action ==
    "abort_without_substitution_or_automatic_repair")
add("validation_abort_contract", guards_ok, "12_validations_10_abort_rules")

no_execution <- all(vapply(list(queue, scope, decision), function(x) {
  for (name in intersect(c("labels_opened", "rankings_computed",
                            "outcomes_computed", "execution_completed",
                            "nested_192_calculation_executed"), names(x))) {
    if (any(truth(x[[name]]))) return(FALSE)
  }
  TRUE
}, logical(1L)))
add("zero_new_calculation_or_outcome", no_execution,
    "calculation_0_rankings_0_labels_0_outcomes_0")

result <- do.call(rbind, checks)
if (!all(result$passed)) stop("MV5-AE validation failed: ",
  paste(result$validation_id[!result$passed], collapse = ", "), call. = FALSE)
if (file.exists(output)) stop("Refusing to overwrite validation output.", call. = FALSE)
write.csv(result, output, row.names = FALSE, na = "")
message("MV5-AE independent validation passed: ", nrow(result), " categories")
