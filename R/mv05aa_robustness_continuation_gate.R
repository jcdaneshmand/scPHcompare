# MV5-AA selection-resistant robustness continuation helpers.

.mv05aa_true <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  tolower(trimws(as.character(x))) == "true"
}

#' Recover and validate the frozen MV5-V configuration order.
mv05aa_configuration_order_v1 <- function(queue) {
  required <- c("configuration_id", "configuration_order", "execution_order")
  if (!is.data.frame(queue) || !all(required %in% names(queue)) ||
      nrow(queue) != 600L) {
    stop("MV5-AA requires the complete 600-row MV5-V queue.", call. = FALSE)
  }
  mapping <- unique(queue[c("configuration_id", "configuration_order")])
  mapping$configuration_order <- as.integer(mapping$configuration_order)
  mapping <- mapping[order(mapping$configuration_order), , drop = FALSE]
  expected <- c(
    "cells384_pc20_euclidean_v1",
    "cells384_pc30_cosine_chord_v1",
    "nested_cells_192_pc30_euclidean_v1",
    "nested_cells_256_pc30_euclidean_v1")
  counts <- table(queue$configuration_id)
  if (nrow(mapping) != 4L || anyDuplicated(mapping$configuration_order) ||
      !identical(mapping$configuration_order, seq_len(4L)) ||
      !identical(mapping$configuration_id, expected) ||
      !identical(unname(as.integer(counts[expected])), rep(150L, 4L)) ||
      !identical(as.integer(queue$execution_order), seq_len(600L))) {
    stop("The frozen MV5-V configuration order is ambiguous or drifted.",
         call. = FALSE)
  }
  mapping$contract_id <- "mv05aa_configuration_order_v1"
  mapping$position_state <- c("complete", "next_eligible", "later_closed",
                              "later_closed")
  mapping[c("contract_id", "configuration_order", "configuration_id",
            "position_state")]
}

#' Validate the complete, unsliced MV5-Z evidence bundle.
mv05aa_bind_pc20_evidence_v1 <- function(production, primary, macro, intervals) {
  required_counts <- c(prediction_groups = 150L, ranking_rows = 282800L,
    outcome_groups = 150L, query_method_rows = 3600L,
    long_query_endpoint_rows = 7200L, macro_estimands = 24L,
    intervals = 24L, primary_tests = 4L)
  if (!is.data.frame(production) || nrow(production) != 1L ||
      !all(names(required_counts) %in% names(production)) ||
      any(vapply(names(required_counts), function(name) {
        as.integer(production[[name]][[1L]]) != required_counts[[name]]
      }, logical(1L))) ||
      nrow(primary) != 4L || nrow(macro) != 24L || nrow(intervals) != 24L ||
      !all(c("clustering_executed", "other_configurations_executed",
             "method_selection_executed", "equivalence_claim_authorized") %in%
           names(production)) ||
      any(.mv05aa_true(unlist(production[1L, c(
        "clustering_executed", "other_configurations_executed",
        "method_selection_executed", "equivalence_claim_authorized")])))) {
    stop("The complete MV5-Z evidence bundle is missing or drifted.",
         call. = FALSE)
  }
  data.frame(
    contract_id = "mv05aa_pc20_evidence_binding_v1",
    evidence_scope = "all_24_estimands_and_all_4_primary_tests_unsliced",
    prediction_groups = 150L, outcome_groups = 150L,
    macro_estimands = 24L, primary_tests = 4L,
    uniform_robustness_established = FALSE,
    equivalence_established = FALSE,
    cosine_geometry_answered_by_pc20 = FALSE,
    subgroup_or_method_selection_used = FALSE,
    complete_evidence_bound = TRUE,
    stringsAsFactors = FALSE)
}

#' Frozen whole-analysis continuation criteria.
mv05aa_continuation_criteria_v1 <- function() {
  data.frame(
    contract_id = "mv05aa_continuation_criteria_v1",
    criterion_order = seq_len(8L),
    criterion_id = c(
      "canonical_next_configuration", "complete_pc20_evidence",
      "distinct_named_alternative", "one_factor_identifiability",
      "artifact_readiness", "bounded_resource_envelope",
      "selection_firewall", "sequential_stop_boundary"),
    required_state = c(
      "cosine_is_position_2_after_complete_pc20",
      "all_24_estimands_and_4_primary_tests_bound_unsliced",
      "radial_scale_geometry_not_answered_by_pc_count",
      "384_cells_30pc_only_point_metric_changes",
      "all_150_coordinate_sources_already_frozen_and_zero_norm_checked",
      "one_worker_8_hours_4_gib_storage_and_per_group_caps",
      "no_representation_dimension_tissue_endpoint_seed_or_pvalue_selection",
      "stop_before_labels_outcomes_clustering_or_nested_cell_configs"),
    mandatory = TRUE,
    frozen_before_new_calculation = TRUE,
    stringsAsFactors = FALSE)
}

#' Make the single allowed MV5-AA continuation decision.
mv05aa_decide_v1 <- function(order, evidence, criterion_pass) {
  criteria <- mv05aa_continuation_criteria_v1()
  if (!is.data.frame(order) || nrow(order) != 4L ||
      order$configuration_id[[2L]] != "cells384_pc30_cosine_chord_v1" ||
      order$position_state[[1L]] != "complete" ||
      order$position_state[[2L]] != "next_eligible" ||
      !is.data.frame(evidence) || nrow(evidence) != 1L ||
      !isTRUE(evidence$complete_evidence_bound[[1L]]) ||
      !isFALSE(evidence$cosine_geometry_answered_by_pc20[[1L]]) ||
      !is.data.frame(criterion_pass) ||
      !identical(criterion_pass$criterion_id, criteria$criterion_id) ||
      !all(.mv05aa_true(criterion_pass$passed))) {
    stop("MV5-AA continuation cannot be authorized.", call. = FALSE)
  }
  data.frame(
    contract_id = "mv05aa_continuation_decision_v1",
    decision = "authorize_later_label_closed_cosine_calculation_only",
    authorized_configuration_id = order$configuration_id[[2L]],
    authorized_groups = 150L,
    pc20_result_used_as = "complete_reporting_and_need_for_prespecified_sensitivity_only",
    favorable_subgroup_selection_used = FALSE,
    cosine_calculation_executed = FALSE,
    labels_opened = FALSE, outcomes_computed = FALSE,
    clustering_authorized = FALSE, nested_configurations_authorized = FALSE,
    stringsAsFactors = FALSE)
}

#' Construct the exact later cosine calculation queue.
mv05aa_cosine_queue_v1 <- function(queue, decision) {
  id <- "cells384_pc30_cosine_chord_v1"
  if (!is.data.frame(decision) || nrow(decision) != 1L ||
      decision$authorized_configuration_id[[1L]] != id ||
      .mv05aa_true(decision$cosine_calculation_executed[[1L]])) {
    stop("MV5-AA decision does not authorize queue construction.",
         call. = FALSE)
  }
  result <- queue[queue$configuration_id == id, , drop = FALSE]
  result <- result[order(result$execution_order), , drop = FALSE]
  if (nrow(result) != 150L || anyDuplicated(result$robustness_group_id) ||
      length(unique(result$fold_id)) != 15L ||
      length(unique(result$seed)) != 5L ||
      !setequal(result$representation,
                c("sct_whole", "inductive_integrated")) ||
      any(as.integer(result$cells) != 384L) ||
      any(as.integer(result$coordinates) != 30L) ||
      any(result$point_metric !=
          "euclidean_chord_after_row_unit_normalization") ||
      any(.mv05aa_true(result$outcomes_computed))) {
    stop("The frozen cosine queue is incomplete or contaminated.",
         call. = FALSE)
  }
  result$mv05aa_execution_order <- seq_len(nrow(result))
  result$later_label_closed_calculation_authorized <- TRUE
  result$labels_opened <- FALSE
  result$outcomes_computed <- FALSE
  result$rankings_computed <- FALSE
  result$execution_completed <- FALSE
  result
}
