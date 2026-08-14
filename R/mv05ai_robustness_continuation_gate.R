# MV5-AI selection-resistant nested-cell continuation helpers.

.mv05ai_true <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  tolower(trimws(as.character(x))) == "true"
}

#' Recover the immutable MV5-V configuration order after three completions.
mv05ai_configuration_order_v1 <- function(queue) {
  required <- c("configuration_id", "configuration_order", "execution_order")
  if (!is.data.frame(queue) || !all(required %in% names(queue)) ||
      nrow(queue) != 600L) {
    stop("MV5-AI requires the complete 600-row MV5-V queue.", call. = FALSE)
  }
  expected <- c(
    "cells384_pc20_euclidean_v1",
    "cells384_pc30_cosine_chord_v1",
    "nested_cells_192_pc30_euclidean_v1",
    "nested_cells_256_pc30_euclidean_v1")
  mapping <- unique(queue[c("configuration_id", "configuration_order")])
  mapping$configuration_order <- as.integer(mapping$configuration_order)
  mapping <- mapping[order(mapping$configuration_order), , drop = FALSE]
  counts <- table(queue$configuration_id)
  if (nrow(mapping) != 4L || anyDuplicated(mapping$configuration_order) ||
      !identical(mapping$configuration_order, 1:4) ||
      !identical(mapping$configuration_id, expected) ||
      !identical(unname(as.integer(counts[expected])), rep(150L, 4L)) ||
      !identical(as.integer(queue$execution_order), seq_len(600L))) {
    stop("The frozen MV5-V configuration order is ambiguous or drifted.",
         call. = FALSE)
  }
  mapping$contract_id <- "mv05ai_configuration_order_v1"
  mapping$position_state <- c("complete", "complete", "complete",
                              "next_eligible")
  mapping[c("contract_id", "configuration_order", "configuration_id",
            "position_state")]
}

#' Bind one complete unsliced 24-estimand robustness result.
mv05ai_bind_complete_evidence_v1 <- function(
    analysis_id, production, primary, macro, intervals) {
  required_counts <- c(
    prediction_groups = 150L, ranking_rows = 282800L,
    outcome_groups = 150L, query_method_rows = 3600L,
    long_query_endpoint_rows = 7200L, macro_estimands = 24L,
    intervals = 24L, primary_tests = 4L)
  if (length(analysis_id) != 1L ||
      !analysis_id %in% c("pc20_vs_pc30", "cosine_chord_vs_euclidean",
                           "nested192_vs_384_cells") ||
      !is.data.frame(production) || nrow(production) != 1L ||
      !all(names(required_counts) %in% names(production)) ||
      any(vapply(names(required_counts), function(name) {
        as.integer(production[[name]][[1L]]) != required_counts[[name]]
      }, logical(1L))) ||
      nrow(primary) != 4L || nrow(macro) != 24L || nrow(intervals) != 24L ||
      anyDuplicated(intervals$estimand_id) ||
      !setequal(primary$estimand_id,
                intervals$estimand_id[intervals$estimand_role %in%
                  c("confirmatory_pc20_sensitivity",
                    "confirmatory_cosine_sensitivity",
                    "confirmatory_nested192_sensitivity")])) {
    stop("A complete robustness evidence bundle is missing or drifted.",
         call. = FALSE)
  }
  prohibited <- intersect(c(
    "clustering_executed", "other_configurations_executed",
    "nested_configurations_executed", "method_selection_executed",
    "equivalence_claim_authorized"), names(production))
  if (length(prohibited) && any(.mv05ai_true(unlist(production[1L, prohibited])))) {
    stop("A robustness evidence bundle reports a prohibited operation.",
         call. = FALSE)
  }
  data.frame(
    contract_id = "mv05ai_complete_evidence_binding_v1",
    analysis_id = analysis_id,
    evidence_scope = "all_24_estimands_intervals_and_4_primary_tests_unsliced",
    prediction_groups = 150L, outcome_groups = 150L,
    macro_estimands = 24L, intervals = 24L, primary_tests = 4L,
    nested256_sensitivity_answered = FALSE,
    subgroup_or_method_selection_used = FALSE,
    complete_evidence_bound = TRUE,
    stringsAsFactors = FALSE)
}

#' Frozen whole-analysis post-nested-192 continuation criteria.
mv05ai_continuation_criteria_v1 <- function() {
  data.frame(
    contract_id = "mv05ai_continuation_criteria_v1",
    criterion_order = seq_len(9L),
    criterion_id = c(
      "canonical_next_configuration", "complete_prior_evidence",
      "distinct_named_alternative", "one_factor_identifiability",
      "nested_source_readiness", "bounded_resource_envelope",
      "selection_firewall", "complete_later_reporting",
      "sequential_stop_boundary"),
    required_state = c(
      "nested_256_is_position_4_after_pc20_cosine_and_nested192_complete",
      "three_complete_24_estimand_4_test_panels_bound_unsliced",
      "nested256_is_prespecified_resolution_between_192_and_384_not_selected_by_results",
      "30pc_euclidean_fixed_only_nested_cell_count_changes_384_to_256",
      "all_13500_views_use_same_sha256_order_with_192_strictly_nested_in_256",
      "observed_nested192_full_execution_supports_conservative_one_worker_caps",
      "no_representation_dimension_tissue_endpoint_seed_estimate_interval_or_pvalue_selection",
      "later_nested_256_outcome_requires_same_24_estimand_complete_panel",
      "stop_before_calculation_ranking_labels_outcomes_or_clustering"),
    mandatory = TRUE, frozen_before_new_calculation = TRUE,
    stringsAsFactors = FALSE)
}

#' Make the single allowed MV5-AI continuation decision.
mv05ai_decide_v1 <- function(order, evidence, criterion_pass) {
  criteria <- mv05ai_continuation_criteria_v1()
  if (!is.data.frame(order) || nrow(order) != 4L ||
      order$configuration_id[[4L]] !=
        "nested_cells_256_pc30_euclidean_v1" ||
      !identical(order$position_state,
                 c("complete", "complete", "complete", "next_eligible")) ||
      !is.data.frame(evidence) || nrow(evidence) != 3L ||
      !setequal(evidence$analysis_id,
                c("pc20_vs_pc30", "cosine_chord_vs_euclidean",
                  "nested192_vs_384_cells")) ||
      !all(.mv05ai_true(evidence$complete_evidence_bound)) ||
      any(.mv05ai_true(evidence$nested256_sensitivity_answered)) ||
      !is.data.frame(criterion_pass) ||
      !identical(criterion_pass$criterion_id, criteria$criterion_id) ||
      !all(.mv05ai_true(criterion_pass$passed))) {
    stop("MV5-AI continuation cannot be authorized.", call. = FALSE)
  }
  data.frame(
    contract_id = "mv05ai_continuation_decision_v1",
    decision = "authorize_later_label_closed_nested_256_calculation_only",
    authorized_configuration_id = order$configuration_id[[4L]],
    authorized_groups = 150L,
    prior_results_used_as =
      "complete_reporting_and_need_for_prespecified_distinct_sensitivity_only",
    favorable_subgroup_or_estimate_selection_used = FALSE,
    nested_256_calculation_executed = FALSE,
    labels_opened = FALSE, rankings_computed = FALSE,
    outcomes_computed = FALSE, clustering_authorized = FALSE,
    nested_256_calculation_authorized = TRUE, stringsAsFactors = FALSE)
}

#' Construct the exact later nested-256 calculation queue.
mv05ai_nested_256_queue_v1 <- function(queue, decision) {
  id <- "nested_cells_256_pc30_euclidean_v1"
  if (!is.data.frame(decision) || nrow(decision) != 1L ||
      decision$authorized_configuration_id[[1L]] != id ||
      .mv05ai_true(decision$nested_256_calculation_executed[[1L]])) {
    stop("MV5-AI decision does not authorize queue construction.", call. = FALSE)
  }
  result <- queue[queue$configuration_id == id, , drop = FALSE]
  result <- result[order(result$execution_order), , drop = FALSE]
  if (nrow(result) != 150L || anyDuplicated(result$robustness_group_id) ||
      length(unique(result$fold_id)) != 15L ||
      length(unique(result$seed)) != 5L ||
      !setequal(result$representation,
                c("sct_whole", "inductive_integrated")) ||
      any(as.integer(result$cells) != 256L) ||
      any(as.integer(result$coordinates) != 30L) ||
      any(result$point_metric != "euclidean") ||
      any(.mv05ai_true(result$outcomes_computed))) {
    stop("The frozen nested-256 queue is incomplete or contaminated.",
         call. = FALSE)
  }
  result$mv05ai_execution_order <- seq_len(nrow(result))
  result$later_label_closed_calculation_authorized <- TRUE
  result$labels_opened <- FALSE; result$rankings_computed <- FALSE
  result$outcomes_computed <- FALSE; result$execution_completed <- FALSE
  result
}
