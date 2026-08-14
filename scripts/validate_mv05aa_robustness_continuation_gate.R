#!/usr/bin/env Rscript

options(warn = 2)
if (!requireNamespace("digest", quietly = TRUE)) {
  stop("digest is required for MV5-AA validation.", call. = FALSE)
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("usage: validate_mv05aa_robustness_continuation_gate.R AUDIT_DIR",
       call. = FALSE)
}
audit_dir <- args[[1L]]
read_audit <- function(name) utils::read.csv(
  file.path(audit_dir, name), stringsAsFactors = FALSE, check.names = FALSE)
truth <- function(x) !is.na(x) & tolower(trimws(as.character(x))) == "true"
file_sha <- function(path) digest::digest(
  file = path, algo = "sha256", serialize = FALSE)
checks <- list()
record <- function(id, passed, evidence) {
  checks[[length(checks) + 1L]] <<- data.frame(
    contract_id = "mv05aa_independent_validation_v1",
    validation_id = id, passed = isTRUE(passed), evidence = evidence,
    production_scientific_helpers_called = FALSE,
    stringsAsFactors = FALSE)
  if (!isTRUE(passed)) stop("Independent MV5-AA validation failed: ", id,
                            call. = FALSE)
}

sources <- read_audit("mv05aa-source-freeze-2026-08-11.csv")
record("source_freeze", nrow(sources) == 19L &&
  all(file.exists(sources$path)) &&
  identical(unname(vapply(sources$path, file_sha, character(1L))),
            sources$sha256) &&
  identical(unname(file.info(sources$path)$size), as.numeric(sources$bytes)),
  "19 paths independently rehashed and resized")

order <- read_audit("mv05aa-configuration-order-2026-08-11.csv")
expected_ids <- c("cells384_pc20_euclidean_v1",
  "cells384_pc30_cosine_chord_v1",
  "nested_cells_192_pc30_euclidean_v1",
  "nested_cells_256_pc30_euclidean_v1")
record("canonical_order", identical(order$configuration_id, expected_ids) &&
  identical(as.integer(order$configuration_order), 1:4) &&
  identical(order$position_state,
    c("complete", "next_eligible", "later_closed", "later_closed")),
  "PC20 complete; cosine is the unique next position")

evidence <- read_audit("mv05aa-complete-pc20-evidence-binding-2026-08-11.csv")
record("complete_pc20_binding", nrow(evidence) == 1L &&
  evidence$prediction_groups == 150L && evidence$outcome_groups == 150L &&
  evidence$macro_estimands == 24L && evidence$primary_tests == 4L &&
  truth(evidence$complete_evidence_bound) &&
  !truth(evidence$cosine_geometry_answered_by_pc20) &&
  !truth(evidence$subgroup_or_method_selection_used),
  "all 24 estimands and four primary tests are bound unsliced")

criteria <- read_audit("mv05aa-continuation-criteria-2026-08-11.csv")
record("criterion_registry", nrow(criteria) == 8L &&
  identical(as.integer(criteria$criterion_order), 1:8) &&
  all(truth(criteria$mandatory)) && all(truth(criteria$passed)) &&
  all(truth(criteria$frozen_before_new_calculation)),
  "all eight prospectively committed criteria pass")

decision <- read_audit("mv05aa-continuation-decision-2026-08-11.csv")
record("single_decision", nrow(decision) == 1L &&
  decision$decision == "authorize_later_label_closed_cosine_calculation_only" &&
  decision$authorized_configuration_id == expected_ids[[2L]] &&
  decision$authorized_groups == 150L &&
  !truth(decision$favorable_subgroup_selection_used) &&
  !truth(decision$cosine_calculation_executed) &&
  !truth(decision$labels_opened) && !truth(decision$outcomes_computed) &&
  !truth(decision$clustering_authorized) &&
  !truth(decision$nested_configurations_authorized),
  "only later label-closed cosine calculation is authorized")
prohibited_decision_columns <- c("representation", "homology_dimension",
  "tissue", "endpoint", "seed", "estimate", "p_value", "adjusted_p_value")
record("selection_firewall",
  !any(prohibited_decision_columns %in% names(decision)),
  "decision schema has no subgroup, estimate, interval, or p-value input")

queue <- read_audit("mv05aa-cosine-execution-queue-2026-08-11.csv")
record("cosine_queue_axis", nrow(queue) == 150L &&
  !anyDuplicated(queue$robustness_group_id) &&
  all(queue$configuration_id == expected_ids[[2L]]) &&
  length(unique(queue$fold_id)) == 15L &&
  length(unique(queue$seed)) == 5L &&
  setequal(queue$representation, c("sct_whole", "inductive_integrated")) &&
  identical(as.integer(queue$mv05aa_execution_order), 1:150),
  "15 folds x five seeds x two representations")
record("cosine_isolation", all(as.integer(queue$cells) == 384L) &&
  all(as.integer(queue$coordinates) == 30L) &&
  all(queue$point_metric ==
      "euclidean_chord_after_row_unit_normalization") &&
  all(truth(queue$later_label_closed_calculation_authorized)),
  "only point geometry changes from the accepted reference")
record("closed_state", !any(truth(queue$labels_opened)) &&
  !any(truth(queue$outcomes_computed)) &&
  !any(truth(queue$rankings_computed)) &&
  !any(truth(queue$execution_completed)),
  "all 150 groups remain unexecuted and label/outcome closed")

scope <- read_audit("mv05aa-cosine-execution-scope-2026-08-11.csv")
expected_scope <- c(groups = 150L, folds = 15L, seeds = 5L,
  representations = 2L, views = 13500L, biological_pairs = 70700L,
  landscape_request_rows = 141400L, landscape_subchunks = 720L,
  energy_request_rows = 70700L, assembled_method_rows = 282800L)
record("scope_reconstruction", all(vapply(names(expected_scope), function(name) {
  as.integer(scope[[name]]) == expected_scope[[name]]
}, logical(1L))) && scope$views == sum(as.integer(queue$view_count)) &&
  scope$biological_pairs == sum(as.integer(queue$biological_pairs)) &&
  scope$landscape_request_rows == sum(as.integer(queue$landscape_request_rows)),
  "all scope totals independently reconstructed from 150 queue rows")

resource <- read_audit("mv05aa-cosine-resource-envelope-2026-08-11.csv")
x_resources <- utils::read.csv(
  "docs/audits/mv05x-pc20-production-resources-2026-08-11.csv",
  stringsAsFactors = FALSE, check.names = FALSE)
record("resource_precedent", nrow(x_resources) == 150L &&
  isTRUE(all.equal(resource$precedent_worker_seconds,
                   sum(as.numeric(x_resources$elapsed_seconds)))) &&
  isTRUE(all.equal(resource$precedent_max_group_seconds,
                   max(as.numeric(x_resources$elapsed_seconds)))) &&
  isTRUE(all.equal(resource$precedent_peak_rss_bytes,
                   max(as.numeric(x_resources$peak_process_tree_rss_bytes)))) &&
  resource$cosine_max_workers == 1L &&
  resource$cosine_configuration_cap_worker_hours == 8 &&
  resource$cosine_storage_cap_bytes == 4294967296,
  "measured PC20 precedent reconstructed; prospective cosine caps retained")

validation <- read_audit("mv05aa-cosine-validation-plan-2026-08-11.csv")
aborts <- read_audit("mv05aa-cosine-abort-rules-2026-08-11.csv")
record("later_execution_gates", nrow(validation) == 12L &&
  all(truth(validation$required)) && nrow(aborts) == 10L &&
  all(aborts$required_action ==
      "abort_without_substitution_or_automatic_repair"),
  "12 validation classes and 10 hard-abort classes frozen")

result <- do.call(rbind, checks)
path <- file.path(
  audit_dir, "mv05aa-independent-validation-2026-08-11.csv")
if (file.exists(path)) stop("Refusing to overwrite: ", path, call. = FALSE)
utils::write.csv(result, path, row.names = FALSE, na = "")
message("MV5-AA independent validation passed ", nrow(result), " categories.")
