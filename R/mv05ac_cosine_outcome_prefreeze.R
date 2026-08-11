# MV5-AC prospective cosine-chord robustness-outcome contracts.

.mv05ac_digest <- function(value) {
  digest::digest(value, algo = "sha256", serialize = TRUE)
}

.mv05ac_is_true <- function(value) {
  tolower(as.character(value)) == "true"
}

.mv05ac_hash_ok <- function(value, width = 64L) {
  length(value) > 0L && all(grepl(sprintf("^[0-9a-f]{%d}$", width), value))
}

.mv05ac_assert_preoutcome <- function(value, label = "input") {
  if (!is.data.frame(value)) stop(label, " must be a data frame.", call. = FALSE)
  prohibited <- c(
    "reciprocal_rank", "one_nn_correct", "first_same_tissue_rank",
    "estimate", "effect", "confidence_interval", "p_value",
    "adjusted_p_value", "winner", "method_rank")
  found <- intersect(tolower(names(value)), prohibited)
  if (length(found)) {
    stop(label, " contains result columns before authorization: ",
         paste(found, collapse = ", "), ".", call. = FALSE)
  }
  for (column in intersect(c("outcomes_computed", "evaluation_executed",
                              "ranking_executed", "method_selection_executed"),
                           names(value))) {
    if (any(.mv05ac_is_true(value[[column]]))) {
      stop(label, " reports a prohibited completed operation in ", column,
           ".", call. = FALSE)
    }
  }
  invisible(TRUE)
}

#' Frozen Euclidean-to-cosine-chord method map.
#'
#' @return A data frame defining the two representations and four paired
#'   distance families.
mv05ac_method_map_v1 <- function() {
  representation <- rep(c("sct_whole", "inductive_integrated"), each = 4L)
  family <- rep(c("h0", "h1", "raw_composite", "energy"), 2L)
  baseline <- c(
    "cell_landscape_h0_v1", "cell_landscape_h1_v1",
    "cell_landscape_h0_h1_raw_euclidean_v1",
    "cell_distribution_energy_shared_pca_v1",
    "integrated_cell_landscape_h0_v1",
    "integrated_cell_landscape_h1_v1",
    "integrated_cell_landscape_h0_h1_raw_euclidean_v1",
    "integrated_cell_distribution_energy_v1")
  cosine <- rep(c(
    "cell_landscape_h0_v1", "cell_landscape_h1_v1",
    "cell_landscape_h0_h1_raw_euclidean_v1",
    "cell_distribution_energy_shared_pca_v1"), 2L)
  data.frame(
    contract_id = "mv05ac_method_map_v1",
    representation = representation,
    family_id = family,
    baseline_method_id = baseline,
    cosine_method_id = cosine,
    method_role = rep(c("confirmatory_topology", "confirmatory_topology",
                        "descriptive_uncalibrated_composite",
                        "matched_same_coordinate_baseline"), 2L),
    baseline_coordinate_count = 30L,
    cosine_coordinate_count = 30L,
    cells = 384L,
    baseline_point_metric = "euclidean",
    cosine_point_metric = "euclidean_chord_after_row_unit_normalization",
    pair_scope = "heldout_query_to_training_reference",
    outcomes_computed = FALSE,
    evaluation_executed = FALSE,
    stringsAsFactors = FALSE)
}

#' Frozen cosine-chord robustness endpoint registry.
#'
#' @return A two-row endpoint registry.
mv05ac_endpoint_registry_v1 <- function() {
  data.frame(
    contract_id = "mv05ac_endpoint_registry_v1",
    endpoint_id = c("cross_study_tissue_mrr_v1",
                    "cross_study_tissue_1nn_balanced_accuracy_v1"),
    query_value = c("reciprocal_rank", "fixed_rank_one_correctness"),
    endpoint_role = c("primary", "supportive"),
    aggregation_policy =
      "mean_five_seeds_within_sample_then_samples_within_tissue_then_equal_five_tissue_macro",
    independent_unit = "heldout_biological_sample",
    blocking_unit = "heldout_study_within_tissue",
    direction_policy = c(
      "positive_signed_change_means_higher_cosine_mrr_zero_means_stability",
      "positive_signed_change_means_higher_cosine_balanced_accuracy_zero_means_stability"),
    outcomes_computed = FALSE,
    evaluation_executed = FALSE,
    stringsAsFactors = FALSE)
}

#' Frozen cosine-chord robustness estimands.
#'
#' @return The complete 24-row estimand registry.
mv05ac_estimand_registry_v1 <- function(
    method_map = mv05ac_method_map_v1(),
    endpoints = mv05ac_endpoint_registry_v1()) {
  direct <- merge(
    method_map[c("representation", "family_id", "method_role")],
    endpoints[c("endpoint_id", "endpoint_role")], all = TRUE)
  direct$estimand_type <- "direct_cosine_chord_minus_euclidean"
  direct$estimand_role <- ifelse(
    direct$endpoint_role == "primary", "secondary_sensitivity",
    "supportive_sensitivity")
  direct$formula <- "macro_Y_cosine_chord_method_minus_macro_Y_euclidean_method"
  direct$matched_baseline_family <- NA_character_

  did_map <- method_map[method_map$family_id %in% c("h0", "h1"),
                        c("representation", "family_id", "method_role")]
  did <- merge(did_map,
               endpoints[c("endpoint_id", "endpoint_role")], all = TRUE)
  did$estimand_type <- "topology_increment_cosine_chord_minus_euclidean_difference_in_differences"
  did$estimand_role <- ifelse(
    did$endpoint_role == "primary", "confirmatory_cosine_sensitivity",
    "supportive_cosine_sensitivity")
  did$formula <- paste0(
    "(macro_Y_cosine_chord_topology_minus_macro_Y_cosine_chord_energy)-",
    "(macro_Y_euclidean_topology_minus_macro_Y_euclidean_energy)")
  did$matched_baseline_family <- "energy"

  result <- rbind(direct, did)
  result <- result[order(result$estimand_type, result$representation,
                         result$family_id, result$endpoint_id,
                         method = "radix"), , drop = FALSE]
  result$contract_id <- "mv05ac_estimand_registry_v1"
  result$estimand_order <- seq_len(nrow(result))
  result$estimand_id <- vapply(seq_len(nrow(result)), function(index) {
    paste0("mv05ac_estimand_v1:", .mv05ac_digest(list(
      type = result$estimand_type[[index]],
      representation = result$representation[[index]],
      family = result$family_id[[index]],
      endpoint = result$endpoint_id[[index]])))
  }, character(1L))
  result$point_direction <-
    "positive_means_cosine_chord_increased_the_named_quantity_negative_means_decreased_zero_means_unchanged"
  result$interval_policy <-
    "2000_paired_tissue_stratified_heldout_study_block_bootstrap_type7_95_percentile"
  result$p_value_policy <- ifelse(
    result$estimand_role == "confirmatory_cosine_sensitivity",
    "9999_two_sided_paired_study_block_sign_flips_then_holm_across_four_mrr_did_tests",
    "none_complete_interval_reporting")
  result$multiplicity_family <- ifelse(
    result$estimand_role == "confirmatory_cosine_sensitivity",
    "F1_cosine_topology_increment_mrr_four_tests", "none")
  result$equivalence_or_noninferiority_claim_authorized <- FALSE
  result$outcomes_computed <- FALSE
  result$evaluation_executed <- FALSE
  result <- result[c(
    "contract_id", "estimand_id", "estimand_order", "estimand_type",
    "estimand_role", "representation", "family_id", "method_role",
    "matched_baseline_family", "endpoint_id", "endpoint_role", "formula",
    "point_direction", "interval_policy", "p_value_policy",
    "multiplicity_family", "equivalence_or_noninferiority_claim_authorized",
    "outcomes_computed", "evaluation_executed")]
  mv05ac_validate_estimand_registry_v1(result)
  result
}

mv05ac_validate_estimand_registry_v1 <- function(value) {
  if (!is.data.frame(value) || nrow(value) != 24L ||
      anyDuplicated(value$estimand_id) ||
      sum(value$estimand_type == "direct_cosine_chord_minus_euclidean") != 16L ||
      sum(value$estimand_type ==
            "topology_increment_cosine_chord_minus_euclidean_difference_in_differences") != 8L ||
      sum(value$multiplicity_family ==
            "F1_cosine_topology_increment_mrr_four_tests") != 4L ||
      any(value$equivalence_or_noninferiority_claim_authorized) ||
      any(.mv05ac_is_true(value$outcomes_computed)) ||
      any(.mv05ac_is_true(value$evaluation_executed))) {
    stop("MV5-AC estimand registry violates its complete 24-estimand contract.",
         call. = FALSE)
  }
  .mv05ac_assert_preoutcome(value, "MV5-AC estimand registry")
  invisible(value)
}

#' Build the frozen 150-group cosine-chord outcome-evaluation queue.
#'
#' @param group_scope Structural group rows with exact pairing evidence.
#' @param source_freeze_sha256 SHA-256 identity of the source freeze.
#' @return A validated 150-row queue with execution disabled.
mv05ac_build_evaluation_queue_v1 <- function(group_scope,
                                             source_freeze_sha256) {
  required <- c(
    "robustness_group_id", "fold_id", "seed", "representation",
    "configuration_id", "query_samples", "training_samples",
    "biological_pairs", "method_rows", "baseline_group_id",
    "baseline_group_sha256", "cosine_group_manifest_sha256",
    "private_result_locator")
  if (!is.data.frame(group_scope) || nrow(group_scope) != 150L ||
      !all(required %in% names(group_scope)) ||
      anyDuplicated(group_scope$robustness_group_id) ||
      length(source_freeze_sha256) != 1L ||
      !.mv05ac_hash_ok(source_freeze_sha256)) {
    stop("MV5-AC structural group scope is invalid.", call. = FALSE)
  }
  queue <- group_scope[required]
  queue <- queue[order(queue$representation, queue$fold_id, queue$seed,
                       method = "radix"), , drop = FALSE]
  queue$contract_id <- "mv05ac_cosine_evaluation_queue_v1"
  queue$evaluation_order <- seq_len(nrow(queue))
  queue$evaluation_unit_id <- vapply(seq_len(nrow(queue)), function(index) {
    paste0("mv05ac_eval_v1:", .mv05ac_digest(list(
      robustness_group_id = queue$robustness_group_id[[index]],
      baseline_group_id = queue$baseline_group_id[[index]],
      source_freeze_sha256 = source_freeze_sha256)))
  }, character(1L))
  queue$source_freeze_sha256 <- source_freeze_sha256
  queue$method_families <- 4L
  queue$endpoints <- 2L
  queue$expected_query_method_rows <- queue$query_samples * 4L
  queue$expected_query_endpoint_rows <- queue$query_samples * 4L * 2L
  queue$label_state <- "external_hash_bound_not_joined_to_cosine_predictions"
  queue$ranking_executed <- FALSE
  queue$outcomes_computed <- FALSE
  queue$evaluation_executed <- FALSE
  queue$method_selection_executed <- FALSE
  queue$execution_authorized <- FALSE
  queue <- queue[c(
    "contract_id", "evaluation_unit_id", "evaluation_order",
    "robustness_group_id", "fold_id", "seed", "representation",
    "configuration_id", "query_samples", "training_samples",
    "biological_pairs", "method_rows", "method_families", "endpoints",
    "expected_query_method_rows", "expected_query_endpoint_rows",
    "baseline_group_id", "baseline_group_sha256",
    "cosine_group_manifest_sha256", "private_result_locator",
    "source_freeze_sha256", "label_state", "ranking_executed",
    "outcomes_computed", "evaluation_executed",
    "method_selection_executed", "execution_authorized")]
  mv05ac_validate_evaluation_queue_v1(queue)
  queue
}

mv05ac_validate_evaluation_queue_v1 <- function(queue) {
  if (!is.data.frame(queue) || nrow(queue) != 150L ||
      anyDuplicated(queue$evaluation_unit_id) ||
      length(unique(queue$fold_id)) != 15L ||
      !identical(sort(unique(as.integer(queue$seed))), 20260805:20260809) ||
      !setequal(queue$representation, c("sct_whole", "inductive_integrated")) ||
      any(table(queue$representation) != 75L) ||
      any(table(paste(queue$fold_id, queue$seed, sep = "\r")) != 2L) ||
      sum(queue$biological_pairs) != 70700L ||
      sum(queue$method_rows) != 282800L ||
      sum(queue$expected_query_method_rows) != 3600L ||
      sum(queue$expected_query_endpoint_rows) != 7200L ||
      any(.mv05ac_is_true(queue$ranking_executed)) ||
      any(.mv05ac_is_true(queue$outcomes_computed)) ||
      any(.mv05ac_is_true(queue$evaluation_executed)) ||
      any(.mv05ac_is_true(queue$method_selection_executed)) ||
      any(.mv05ac_is_true(queue$execution_authorized))) {
    stop("MV5-AC queue violates its prospective 150-group contract.",
         call. = FALSE)
  }
  .mv05ac_assert_preoutcome(queue, "MV5-AC evaluation queue")
  invisible(queue)
}
