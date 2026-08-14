# MV5-T selection-resistant robustness gate helpers.

.mv05t_digest <- function(value) {
  digest::digest(value, algo = "sha256", serialize = TRUE)
}

mv05t_criteria_v1 <- function() {
  data.frame(
    contract_id = "mv05t_criteria_v1",
    criterion_id = c("scientific_value", "identifiability_validity",
                     "outcome_selection_safety", "artifact_readiness",
                     "resource_feasibility", "reviewer_relevance"),
    weight = c(3L, 3L, 3L, 2L, 1L, 2L), minimum = 0L, maximum = 4L,
    scale_frozen_before_scoring = TRUE, stringsAsFactors = FALSE)
}

mv05t_candidate_registry_v1 <- function() {
  criteria <- mv05t_criteria_v1()
  scores <- data.frame(
    candidate_id = c(
      "nested_cell_count_192_256", "pc20_truncation",
      "cosine_chord_geometry", "existing_component_distance_panel",
      "feature_panel_250_1000", "pc50_refit", "additional_cell_seeds",
      "clustering_k_algorithm_expansion", "small_study_eligibility_rerun",
      "integration_method_expansion"),
    scientific_value = c(4, 4, 4, 2, 4, 3, 1, 2, 3, 4),
    identifiability_validity = c(4, 4, 4, 4, 4, 4, 2, 3, 2, 3),
    outcome_selection_safety = c(4, 4, 4, 4, 3, 4, 4, 1, 2, 3),
    artifact_readiness = c(4, 4, 4, 4, 1, 0, 4, 4, 1, 1),
    resource_feasibility = c(3, 3, 3, 4, 1, 1, 3, 4, 1, 0),
    reviewer_relevance = c(4, 3, 3, 3, 3, 3, 2, 2, 3, 4),
    disposition = c(
      "admit_bounded_resource_pilot", "admit_bounded_resource_pilot",
      "admit_bounded_resource_pilot", "already_complete_no_new_execution",
      "defer_requires_training_only_panel_refit", "blocked_no_50pc_cache",
      "reject_redundant_five_seeds_complete",
      "reject_or_defer_postoutcome_search_or_already_complete",
      "defer_requires_complete_fold_rerun",
      "defer_requires_native_inductive_method_gate"),
    named_alternative = c(
      "topology_or_matched_energy_is_driven_by_fixed_cell_count",
      "geometry_or_topology_depends_on_30pc_dimension",
      "euclidean_radial_scale_drives_geometry",
      "H0_H1_or_distance_family_changes_the_result",
      "training_selected_feature_panel_drives_geometry",
      "higher_coordinate_dimension_changes_geometry",
      "five_cell_realizations_understate_sampling_variability",
      "clustering_result_depends_on_algorithm_or_k_search",
      "single_sample_studies_drive_heldout_behavior",
      "Seurat_specific_induction_drives_the_result"),
    mv05s_values_used_for_candidate_choice = FALSE,
    new_outcomes_computed = FALSE, stringsAsFactors = FALSE)
  score_columns <- criteria$criterion_id
  scores$weighted_score <- rowSums(sweep(
    as.matrix(scores[score_columns]), 2L, criteria$weight, `*`))
  scores$selection_eligible <- scores$disposition ==
    "admit_bounded_resource_pilot"
  scores$selected_for_admission <- scores$selection_eligible
  scores$contract_id <- "mv05t_candidate_registry_v1"
  scores[c("contract_id", "candidate_id", score_columns, "weighted_score",
           "selection_eligible", "selected_for_admission", "disposition",
           "named_alternative", "mv05s_values_used_for_candidate_choice",
           "new_outcomes_computed")]
}

mv05t_configuration_registry_v1 <- function() {
  data.frame(
    contract_id = "mv05t_configuration_registry_v1",
    configuration_id = c(
      "nested_cells_192_pc30_euclidean_v1",
      "nested_cells_256_pc30_euclidean_v1",
      "cells384_pc20_euclidean_v1",
      "cells384_pc30_cosine_chord_v1"),
    candidate_id = c("nested_cell_count_192_256", "nested_cell_count_192_256",
                     "pc20_truncation", "cosine_chord_geometry"),
    cells = c(192L, 256L, 384L, 384L),
    coordinates = c(30L, 30L, 20L, 30L),
    point_metric = c("euclidean", "euclidean", "euclidean",
                     "euclidean_chord_after_row_unit_normalization"),
    comparison_design = "one_factor_at_a_time_against_384cell_30pc_euclidean",
    baseline_reused = TRUE,
    role = "postoutcome_secondary_sensitivity",
    outcomes_authorized = FALSE,
    stringsAsFactors = FALSE)
}

mv05t_nested_point_ids_v1 <- function(sample_id, seed, point_ids, cells) {
  if (length(sample_id) != 1L || length(seed) != 1L || anyNA(point_ids) ||
      anyDuplicated(point_ids) || !cells %in% c(192L, 256L, 384L) ||
      length(point_ids) != 384L) {
    stop("MV5-T nested point-selection input is invalid.", call. = FALSE)
  }
  hashes <- vapply(point_ids, function(point_id) {
    .mv05t_digest(list(contract_id = "mv05t_nested_point_order_v1",
                       sample_id = sample_id, seed = as.integer(seed),
                       point_id = point_id))
  }, character(1L))
  point_ids[order(hashes, point_ids, method = "radix")][seq_len(cells)]
}

mv05t_validate_coordinate_pair_v1 <- function(sct_view, integrated_coordinates) {
  if (!is.list(sct_view) || !is.matrix(sct_view$payload) ||
      !is.matrix(integrated_coordinates) ||
      !identical(dim(sct_view$payload), c(384L, 30L)) ||
      !identical(dim(integrated_coordinates), c(384L, 30L)) ||
      !identical(sct_view$point_ids, rownames(sct_view$payload)) ||
      !identical(sct_view$coordinate_ids, colnames(sct_view$payload)) ||
      !identical(sct_view$point_ids, rownames(integrated_coordinates)) ||
      anyNA(sct_view$payload) || anyNA(integrated_coordinates) ||
      any(!is.finite(sct_view$payload)) ||
      any(!is.finite(integrated_coordinates))) {
    stop("MV5-T paired coordinate identity or shape failed.", call. = FALSE)
  }
  invisible(TRUE)
}

mv05t_build_admission_queue_v1 <- function(source_freeze_sha256) {
  if (length(source_freeze_sha256) != 1L ||
      !grepl("^[0-9a-f]{64}$", source_freeze_sha256)) {
    stop("MV5-T source-freeze identity is invalid.", call. = FALSE)
  }
  folds <- data.frame(
    fold_id = c("large_loso_v1:SRA779509", "large_loso_v1:SRA703206",
                "large_loso_v1:SRA713577"),
    fold_role = c("minimum_training_65", "median_training_86",
                  "maximum_training_89"),
    training_samples = c(65L, 86L, 89L), stringsAsFactors = FALSE)
  configurations <- mv05t_configuration_registry_v1()
  queue <- merge(folds, data.frame(
    representation = c("sct_whole", "inductive_integrated"),
    stringsAsFactors = FALSE), all = TRUE)
  queue <- merge(queue, configurations[c(
    "configuration_id", "candidate_id", "cells", "coordinates",
    "point_metric")], all = TRUE)
  queue$seed <- 20260805L
  queue <- queue[order(queue$fold_id, queue$representation,
                       queue$configuration_id, method = "radix"), ]
  queue$contract_id <- "mv05t_admission_queue_v1"
  queue$execution_order <- seq_len(nrow(queue))
  queue$admission_unit_id <- vapply(seq_len(nrow(queue)), function(index) {
    paste0("mv05t_admission_v1:", .mv05t_digest(list(
      fold_id = queue$fold_id[[index]], seed = queue$seed[[index]],
      representation = queue$representation[[index]],
      configuration_id = queue$configuration_id[[index]],
      source_freeze_sha256 = source_freeze_sha256)))
  }, character(1L))
  queue$source_freeze_sha256 <- source_freeze_sha256
  queue$labels_opened <- FALSE
  queue$outcomes_computed <- FALSE
  queue$admission_executed <- FALSE
  queue <- queue[c(
    "contract_id", "admission_unit_id", "execution_order", "fold_id",
    "fold_role", "training_samples", "seed", "representation",
    "configuration_id", "candidate_id", "cells", "coordinates",
    "point_metric", "source_freeze_sha256", "labels_opened",
    "outcomes_computed", "admission_executed")]
  mv05t_validate_admission_queue_v1(queue)
  queue
}

mv05t_validate_admission_queue_v1 <- function(queue) {
  if (!is.data.frame(queue) || nrow(queue) != 24L ||
      anyDuplicated(queue$admission_unit_id) ||
      length(unique(queue$fold_id)) != 3L ||
      length(unique(queue$representation)) != 2L ||
      length(unique(queue$configuration_id)) != 4L ||
      !all(table(queue$fold_id) == 8L) ||
      any(queue$labels_opened) || any(queue$outcomes_computed) ||
      any(queue$admission_executed)) {
    stop("MV5-T admission queue violates the 24-unit no-outcome gate.",
         call. = FALSE)
  }
  invisible(queue)
}
