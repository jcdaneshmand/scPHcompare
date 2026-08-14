# MV5-R prospective clustering-outcome evaluation contracts.

.mv05r_digest <- function(value) {
  digest::digest(value, algo = "sha256", serialize = TRUE)
}

.mv05r_is_true <- function(value) {
  tolower(as.character(value)) == "true"
}

.mv05r_assert_preoutcome <- function(value, label = "input") {
  if (!is.data.frame(value)) stop(label, " must be a data frame.", call. = FALSE)
  prohibited <- c("ari", "nmi", "balanced_accuracy", "correct", "prediction",
                  "estimate", "p_value", "adjusted_p_value", "confidence_interval")
  found <- intersect(tolower(names(value)), prohibited)
  if (length(found)) {
    stop(label, " contains outcome-result columns before authorization: ",
         paste(found, collapse = ", "), ".", call. = FALSE)
  }
  if ("outcomes_computed" %in% names(value) &&
      any(.mv05r_is_true(value$outcomes_computed))) {
    stop(label, " reports computed outcomes before authorization.", call. = FALSE)
  }
  if ("evaluation_executed" %in% names(value) &&
      any(.mv05r_is_true(value$evaluation_executed))) {
    stop(label, " reports executed evaluation before authorization.", call. = FALSE)
  }
  invisible(TRUE)
}

mv05r_algorithm_registry_v1 <- function() {
  data.frame(
    contract_id = "mv05r_algorithm_registry_v1",
    algorithm_id = c("pam_stability_k_v1", "hclust_average_v1"),
    role = c("primary", "sensitivity"),
    selected_k_source = "mv05q_pam_five_seed_one_se_v1",
    refit_authorized = FALSE, oracle_k_authorized = FALSE,
    outcome_driven_tuning_authorized = FALSE,
    outcomes_computed = FALSE, evaluation_executed = FALSE,
    stringsAsFactors = FALSE)
}

mv05r_endpoint_registry_v1 <- function() {
  data.frame(
    contract_id = "mv05r_endpoint_registry_v1",
    endpoint_id = c(
      "training_tissue_ari_v1", "training_tissue_nmi_v1",
      "training_study_ari_v1", "training_study_nmi_v1",
      "training_approach_ari_v1", "training_approach_nmi_v1",
      "heldout_tissue_plurality_balanced_accuracy_v1",
      "heldout_approach_plurality_balanced_accuracy_v1"),
    evaluation_scope = c(rep("overlapping_training_partition_alignment", 6L),
                         rep("heldout_label_prediction_from_frozen_training_cluster", 2L)),
    label_axis = c("tissue", "tissue", "study", "study", "approach",
                   "approach", "tissue", "approach"),
    metric_id = c("adjusted_rand_index", "normalized_mutual_information",
                  "adjusted_rand_index", "normalized_mutual_information",
                  "adjusted_rand_index", "normalized_mutual_information",
                  "balanced_accuracy", "balanced_accuracy"),
    role = c("secondary_biological_alignment", "supportive_biological_alignment",
             "secondary_technical_alignment", "supportive_technical_alignment",
             "secondary_technical_alignment", "supportive_technical_alignment",
             "secondary_heldout_biological_generalization",
             "secondary_heldout_technical_generalization"),
    seed_policy = "five_technical_realizations_report_all_then_average_within_sample_or_fold",
    fold_policy = c(rep("descriptive_foldwise_no_independence_claim", 6L),
                    rep("each_sample_scored_only_in_its_heldout_study_fold", 2L)),
    aggregation_policy = c(
      rep("seed_rows_plus_seed_mean_per_fold_then_complete_fold_distribution", 6L),
      "seed_correctness_mean_per_sample_then_equal_tissue_macro_average",
      "seed_correctness_mean_per_sample_then_equal_approach_macro_average"),
    uncertainty_policy = c(
      rep("descriptive_seed_jackknife_only_no_p_value", 6L),
      "2000_tissue_stratified_study_block_bootstrap_95_percentile",
      "2000_global_study_block_bootstrap_for_approach_95_percentile"),
    multiplicity_family = "none_secondary_clustering_complete_reporting",
    p_value_authorized = FALSE,
    outcome_driven_selection_authorized = FALSE,
    outcomes_computed = FALSE, evaluation_executed = FALSE,
    stringsAsFactors = FALSE)
}

mv05r_plurality_map_v1 <- function(training_partition, labels, label_axis) {
  required_partition <- c("sample_id", "cluster")
  if (!is.data.frame(training_partition) ||
      !all(required_partition %in% names(training_partition)) ||
      !is.data.frame(labels) || !all(c("sample_id", label_axis) %in% names(labels)) ||
      anyDuplicated(training_partition$sample_id) || anyDuplicated(labels$sample_id) ||
      !setequal(training_partition$sample_id, labels$sample_id)) {
    stop("MV5-R plurality-map inputs are malformed or misaligned.", call. = FALSE)
  }
  joined <- merge(training_partition[required_partition],
                  labels[c("sample_id", label_axis)], by = "sample_id", sort = FALSE)
  names(joined)[names(joined) == label_axis] <- "label_value"
  counts <- as.data.frame(table(cluster = joined$cluster,
                                label_value = joined$label_value),
                          stringsAsFactors = FALSE)
  counts <- counts[counts$Freq > 0L, , drop = FALSE]
  rows <- lapply(split(counts, counts$cluster), function(part) {
    part <- part[order(-part$Freq, part$label_value, method = "radix"), , drop = FALSE]
    data.frame(cluster = part$cluster[[1L]], predicted_label = part$label_value[[1L]],
               plurality_count = part$Freq[[1L]],
               plurality_tie_size = sum(part$Freq == max(part$Freq)),
               tie_rule = "maximum_count_then_lexicographic_label_v1",
               stringsAsFactors = FALSE)
  })
  result <- do.call(rbind, rows)
  result <- result[order(as.character(result$cluster), method = "radix"), ]
  rownames(result) <- NULL
  result
}

mv05r_build_evaluation_queue_v1 <- function(analysis_queue,
                                             endpoint_registry =
                                               mv05r_endpoint_registry_v1(),
                                             algorithm_registry =
                                               mv05r_algorithm_registry_v1(),
                                             source_freeze_sha256) {
  required <- c("analysis_group_id", "fold_id", "representation", "distance_id",
                "training_samples")
  if (!is.data.frame(analysis_queue) || nrow(analysis_queue) != 150L ||
      !all(required %in% names(analysis_queue)) ||
      length(source_freeze_sha256) != 1L ||
      !grepl("^[0-9a-f]{64}$", source_freeze_sha256)) {
    stop("MV5-R analysis queue or source freeze is invalid.", call. = FALSE)
  }
  base <- merge(analysis_queue[required],
                algorithm_registry[c("algorithm_id", "role")], all = TRUE)
  names(base)[names(base) == "role"] <- "algorithm_role"
  queue <- merge(base, endpoint_registry[c(
    "endpoint_id", "evaluation_scope", "label_axis", "metric_id", "role",
    "aggregation_policy", "uncertainty_policy", "multiplicity_family")], all = TRUE)
  names(queue)[names(queue) == "role"] <- "endpoint_role"
  queue <- queue[order(queue$fold_id, queue$representation, queue$distance_id,
                       queue$algorithm_id, queue$endpoint_id, method = "radix"), ]
  queue$contract_id <- "mv05r_evaluation_queue_v1"
  queue$execution_order <- seq_len(nrow(queue))
  queue$evaluation_unit_id <- vapply(seq_len(nrow(queue)), function(index) {
    paste0("mv05r_eval_v1:", .mv05r_digest(list(
      analysis_group_id = queue$analysis_group_id[[index]],
      algorithm_id = queue$algorithm_id[[index]],
      endpoint_id = queue$endpoint_id[[index]],
      source_freeze_sha256 = source_freeze_sha256)))
  }, character(1L))
  queue$source_freeze_sha256 <- source_freeze_sha256
  queue$label_source_state <- "frozen_external_not_opened_for_outcomes"
  queue$outcomes_computed <- FALSE
  queue$evaluation_executed <- FALSE
  queue$method_selection_executed <- FALSE
  queue <- queue[c("contract_id", "evaluation_unit_id", "execution_order",
                   "analysis_group_id", "fold_id", "representation", "distance_id",
                   "training_samples", "algorithm_id", "algorithm_role",
                   "endpoint_id", "evaluation_scope", "label_axis", "metric_id",
                   "endpoint_role", "aggregation_policy", "uncertainty_policy",
                   "multiplicity_family", "source_freeze_sha256",
                   "label_source_state", "outcomes_computed",
                   "evaluation_executed", "method_selection_executed")]
  mv05r_validate_evaluation_queue_v1(queue)
  queue
}

mv05r_validate_evaluation_queue_v1 <- function(queue) {
  if (!is.data.frame(queue) || nrow(queue) != 2400L ||
      anyDuplicated(queue$evaluation_unit_id) ||
      length(unique(queue$analysis_group_id)) != 150L ||
      !all(table(queue$analysis_group_id) == 16L) ||
      !setequal(queue$algorithm_id,
                c("pam_stability_k_v1", "hclust_average_v1")) ||
      length(unique(queue$endpoint_id)) != 8L ||
      any(.mv05r_is_true(queue$outcomes_computed)) ||
      any(.mv05r_is_true(queue$evaluation_executed)) ||
      any(.mv05r_is_true(queue$method_selection_executed))) {
    stop("MV5-R evaluation queue violates its prospective 2,400-unit contract.",
         call. = FALSE)
  }
  .mv05r_assert_preoutcome(queue, "MV5-R evaluation queue")
  invisible(queue)
}
