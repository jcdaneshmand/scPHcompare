# MV5-Q label-closed clustering-artifact production contracts.

.mv05q_required_seeds <- 20260805:20260809

.mv05q_digest <- function(value) {
  digest::digest(value, algo = "sha256", serialize = TRUE)
}

.mv05q_hash_ok <- function(value) {
  length(value) > 0L && all(grepl("^[0-9a-f]{64}$", value))
}

.mv05q_is_true <- function(value) {
  tolower(as.character(value)) == "true"
}

.mv05q_assert_label_closed <- function(value, label = "input") {
  if (!is.data.frame(value)) {
    stop(label, " must be a data frame.", call. = FALSE)
  }
  .mv05n_assert_no_outcome_columns(value, label)
  if ("outcome_label_state" %in% names(value) &&
      any(value$outcome_label_state != "closed")) {
    stop(label, " opens the outcome-label firewall.", call. = FALSE)
  }
  if ("biological_outcomes_computed" %in% names(value) &&
      any(.mv05q_is_true(value$biological_outcomes_computed))) {
    stop(label, " reports biological outcome computation.", call. = FALSE)
  }
  invisible(TRUE)
}

mv05q_method_alias_registry_v1 <- function() {
  rows <- data.frame(
    representation = rep(c("sct_whole", "inductive_integrated"), each = 5L),
    distance_id = rep(c(
      "cell_landscape_h0_v1", "cell_landscape_h1_v1",
      "cell_landscape_h0_h1_raw_euclidean_v1",
      "cell_distribution_energy_v1",
      "pseudobulk_training_standardized_panel_v1"
    ), 2L),
    training_representation = c(rep("sct_whole", 5L),
                                rep("inductive_integrated", 4L), "sct_whole"),
    training_method_id = rep(c(
      "persistence_landscape_l2_exact_v1",
      "persistence_landscape_l2_exact_v1",
      "derived_from_exact_h0_h1_v1",
      "cell_distribution_energy_v1",
      "pseudobulk_training_standardized_panel_v1"
    ), 2L),
    training_component = rep(c("H0", "H1", "H0_H1_raw", "distance",
                               "distance"), 2L),
    query_bundle_id = c(rep("mv05d5_sct_query_bundle_v1", 5L),
                        rep("mv05j_integrated_query_bundle_v1", 5L)),
    query_representation = c(rep("sct_fold", 5L),
                             rep("inductive_integrated", 4L),
                             "training_standardized_pseudobulk"),
    query_method_id = c(
      "cell_landscape_h0_v1", "cell_landscape_h1_v1",
      "cell_landscape_h0_h1_raw_euclidean_v1",
      "cell_distribution_energy_shared_pca_v1",
      "pseudobulk_shared_panel_euclidean_v1",
      "integrated_cell_landscape_h0_v1",
      "integrated_cell_landscape_h1_v1",
      "integrated_cell_landscape_h0_h1_raw_euclidean_v1",
      "integrated_cell_distribution_energy_v1",
      "pseudobulk_training_standardized_panel_v1"
    ),
    training_source_policy = c(rep("representation_matched", 4L),
                               "shared_sct_training_pseudobulk",
                               rep("representation_matched", 4L),
                               "shared_sct_training_pseudobulk"),
    stringsAsFactors = FALSE
  )
  rows$contract_id <- "mv05q_method_alias_registry_v1"
  rows$outcome_label_state <- "closed"
  rows$biological_outcomes_computed <- FALSE
  rows[c("contract_id", setdiff(names(rows), "contract_id"))]
}

mv05q_validate_alias_registry_v1 <- function(registry) {
  required <- names(mv05q_method_alias_registry_v1())
  expected <- mv05q_method_alias_registry_v1()
  if (!is.data.frame(registry) || !identical(names(registry), required) ||
      nrow(registry) != 10L ||
      !isTRUE(all.equal(registry, expected, check.attributes = FALSE))) {
    stop("MV5-Q method aliases differ from the frozen ten-analysis registry.",
         call. = FALSE)
  }
  .mv05q_assert_label_closed(registry, "alias registry")
  invisible(registry)
}

mv05q_reconstruct_distance_matrix_v1 <- function(pair_rows, sample_ids) {
  required <- c("first_sample_id", "second_sample_id", "distance")
  if (!is.data.frame(pair_rows) || !all(required %in% names(pair_rows))) {
    stop("MV5-Q pair rows are malformed.", call. = FALSE)
  }
  .mv05q_assert_label_closed(pair_rows, "training pair rows")
  ids <- sort(unique(as.character(sample_ids)), method = "radix")
  pair_rows$first_sample_id <- as.character(pair_rows$first_sample_id)
  pair_rows$second_sample_id <- as.character(pair_rows$second_sample_id)
  pair_rows$distance <- as.numeric(pair_rows$distance)
  pair_rows <- pair_rows[order(pair_rows$first_sample_id,
                               pair_rows$second_sample_id,
                               method = "radix"), , drop = FALSE]
  observed_key <- paste(pair_rows$first_sample_id, pair_rows$second_sample_id,
                        sep = "\r")
  expected_pairs <- utils::combn(ids, 2L)
  expected_key <- paste(expected_pairs[1L, ], expected_pairs[2L, ], sep = "\r")
  if (length(ids) < 3L || anyDuplicated(ids) ||
      any(pair_rows$first_sample_id >= pair_rows$second_sample_id) ||
      anyDuplicated(observed_key) || !setequal(observed_key, expected_key) ||
      any(!is.finite(pair_rows$distance)) || any(pair_rows$distance < 0)) {
    stop("MV5-Q pair rows do not form one complete finite training matrix.",
         call. = FALSE)
  }
  matrix <- matrix(0, nrow = length(ids), ncol = length(ids),
                   dimnames = list(ids, ids))
  index <- cbind(match(pair_rows$first_sample_id, ids),
                 match(pair_rows$second_sample_id, ids))
  matrix[index] <- pair_rows$distance
  matrix[cbind(index[, 2L], index[, 1L])] <- pair_rows$distance
  .mv05n_validate_distance_matrix(matrix)
}

mv05q_combine_h0_h1_v1 <- function(h0, h1) {
  h0 <- .mv05n_validate_distance_matrix(h0)
  h1 <- .mv05n_validate_distance_matrix(h1)
  if (!identical(dimnames(h0), dimnames(h1))) {
    stop("MV5-Q H0 and H1 matrices must share the exact ordered axis.",
         call. = FALSE)
  }
  result <- sqrt(h0 ^ 2 + h1 ^ 2)
  dimnames(result) <- dimnames(h0)
  .mv05n_validate_distance_matrix(result)
}

mv05q_validate_query_slice_v1 <- function(query_rows, training_ids) {
  required <- c("fold_id", "seed", "representation", "method_id",
                "query_sample_id", "training_sample_id", "distance")
  if (!is.data.frame(query_rows) || !all(required %in% names(query_rows)) ||
      !nrow(query_rows) || length(unique(query_rows$fold_id)) != 1L ||
      length(unique(query_rows$seed)) != 1L ||
      length(unique(query_rows$representation)) != 1L ||
      length(unique(query_rows$method_id)) != 1L) {
    stop("MV5-Q query slice does not have one frozen identity.", call. = FALSE)
  }
  .mv05q_assert_label_closed(query_rows, "query rows")
  validated <- .mv05n_validate_query_distances(
    query_rows[c("query_sample_id", "training_sample_id", "distance")],
    training_ids
  )
  if (anyDuplicated(validated$query_sample_id) == 0L &&
      nrow(validated) == length(unique(validated$query_sample_id))) {
    stop("MV5-Q query slice unexpectedly lacks training-axis expansion.",
         call. = FALSE)
  }
  validated
}

mv05q_fit_analysis_group_v1 <- function(distance_matrices, query_distances,
                                         fold_id, representation, distance_id,
                                         analysis_group_id = NULL) {
  if (!is.list(distance_matrices) || !is.list(query_distances) ||
      !identical(names(distance_matrices), as.character(.mv05q_required_seeds)) ||
      !identical(names(query_distances), as.character(.mv05q_required_seeds))) {
    stop("MV5-Q requires ordered matrices and query slices for five seeds.",
         call. = FALSE)
  }
  matrices <- lapply(distance_matrices, .mv05n_validate_distance_matrix)
  axes <- lapply(matrices, rownames)
  if (!all(vapply(axes, identical, logical(1L), axes[[1L]]))) {
    stop("MV5-Q five-seed training axes drifted.", call. = FALSE)
  }
  candidates <- mv05n_fit_five_seed_partitions_v1(matrices, method = "pam")
  selection <- mv05n_select_k_v1(candidates)
  selected_k <- as.integer(selection$selected_k)
  group_identity <- list(contract_id = "mv05q_analysis_group_identity_v1",
                         fold_id = fold_id, representation = representation,
                         distance_id = distance_id,
                         ordered_training_ids = axes[[1L]],
                         seeds = .mv05q_required_seeds)
  derived_group_id <- paste0("mv05q_analysis_v1:", .mv05q_digest(group_identity))
  if (is.null(analysis_group_id)) analysis_group_id <- derived_group_id
  if (length(analysis_group_id) != 1L || is.na(analysis_group_id) ||
      !grepl("^mv05q_analysis_v1:[0-9a-f]{64}$", analysis_group_id)) {
    stop("MV5-Q analysis group ID is malformed.", call. = FALSE)
  }
  stability <- selection$summary
  stability$contract_id <- "mv05q_stability_summary_v1"
  stability$analysis_group_id <- analysis_group_id
  stability$fold_id <- fold_id
  stability$representation <- representation
  stability$distance_id <- distance_id
  stability$selected_k <- selected_k
  stability$one_se_threshold <- selection$threshold
  stability$outcome_label_state <- "closed"
  stability$biological_outcomes_computed <- FALSE
  stability <- stability[c("contract_id", "analysis_group_id", "fold_id",
                           "representation", "distance_id", "k",
                           "mean_stability", "monte_carlo_se", "pair_count",
                           "selected_k", "one_se_threshold",
                           "outcome_label_state", "biological_outcomes_computed")]
  partition_rows <- list()
  heldout_rows <- list()
  cursor <- 0L
  for (seed in .mv05q_required_seeds) {
    key <- as.character(seed)
    query <- mv05q_validate_query_slice_v1(query_distances[[key]], axes[[1L]])
    pam <- mv05n_pam_partition_v1(matrices[[key]], selected_k)
    average <- mv05n_average_partition_v1(matrices[[key]], selected_k)
    for (algorithm in c("pam_stability_k_v1", "hclust_average_v1")) {
      fit <- if (algorithm == "pam_stability_k_v1") pam else average
      cursor <- cursor + 1L
      partition_rows[[cursor]] <- data.frame(
        contract_id = "mv05q_selected_training_partition_v1",
        analysis_group_id = analysis_group_id, fold_id = fold_id,
        seed = seed, representation = representation, distance_id = distance_id,
        algorithm_id = algorithm, selected_k = selected_k,
        sample_id = fit$sample_id, cluster = fit$cluster,
        is_medoid = if (algorithm == "pam_stability_k_v1") fit$is_medoid else FALSE,
        outcome_label_state = "closed", biological_outcomes_computed = FALSE,
        stringsAsFactors = FALSE
      )
      assigned <- if (algorithm == "pam_stability_k_v1") {
        mv05n_assign_pam_heldout_v1(query, pam)
      } else {
        mv05n_assign_average_heldout_v1(query, average)
      }
      heldout_rows[[cursor]] <- data.frame(
        contract_id = "mv05q_heldout_assignment_v1",
        analysis_group_id = analysis_group_id, fold_id = fold_id,
        seed = seed, representation = representation, distance_id = distance_id,
        algorithm_id = algorithm, selected_k = selected_k,
        query_sample_id = assigned$query_sample_id, cluster = assigned$cluster,
        assignment_distance = assigned$assignment_distance,
        assignment_reference = assigned$assignment_reference,
        assignment_rule = assigned$assignment_rule,
        outcome_label_state = "closed", biological_outcomes_computed = FALSE,
        stringsAsFactors = FALSE
      )
    }
  }
  candidates$contract_id <- "mv05q_candidate_pam_partition_v1"
  candidates$analysis_group_id <- analysis_group_id
  candidates$fold_id <- fold_id
  candidates$representation <- representation
  candidates$distance_id <- distance_id
  candidates$outcome_label_state <- "closed"
  candidates$biological_outcomes_computed <- FALSE
  candidates <- candidates[c("contract_id", "analysis_group_id", "fold_id",
                             "seed", "representation", "distance_id", "k",
                             "sample_id", "cluster", "outcome_label_state",
                             "biological_outcomes_computed")]
  list(
    analysis_group_id = analysis_group_id,
    selected_k = selected_k,
    candidate_partitions = candidates,
    stability = stability,
    selected_partitions = do.call(rbind, partition_rows),
    heldout_assignments = do.call(rbind, heldout_rows)
  )
}

mv05q_validate_group_result_v1 <- function(result, training_samples,
                                            heldout_samples) {
  required <- c("candidate_partitions", "stability", "selected_partitions",
                "heldout_assignments")
  if (!is.list(result) || !all(required %in% names(result))) {
    stop("MV5-Q group result is incomplete.", call. = FALSE)
  }
  lapply(result[required], .mv05q_assert_label_closed)
  candidates <- result$candidate_partitions
  stability <- result$stability
  partitions <- result$selected_partitions
  heldout <- result$heldout_assignments
  expected_candidate_rows <- 5L * 9L * as.integer(training_samples)
  expected_partition_rows <- 5L * 2L * as.integer(training_samples)
  expected_heldout_rows <- 5L * 2L * as.integer(heldout_samples)
  if (nrow(candidates) != expected_candidate_rows || nrow(stability) != 9L ||
      nrow(partitions) != expected_partition_rows ||
      nrow(heldout) != expected_heldout_rows ||
      !identical(sort(unique(candidates$seed)), .mv05q_required_seeds) ||
      !identical(sort(unique(candidates$k)), 2:10) ||
      !setequal(partitions$algorithm_id,
                c("pam_stability_k_v1", "hclust_average_v1")) ||
      !setequal(heldout$algorithm_id,
                c("pam_stability_k_v1", "hclust_average_v1")) ||
      any(partitions$selected_k != result$selected_k) ||
      any(heldout$selected_k != result$selected_k)) {
    stop("MV5-Q group result violates complete-axis expectations.",
         call. = FALSE)
  }
  invisible(result)
}
