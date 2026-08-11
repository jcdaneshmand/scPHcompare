# MV5-U bounded label-closed robustness admission helpers.

.mv05u_digest <- function(value) {
  digest::digest(value, algo = "sha256", serialize = TRUE)
}

.mv05u_file_sha <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

.mv05u_one_string <- function(value, label) {
  value <- as.character(value)
  if (length(value) != 1L || is.na(value) || !nzchar(value)) {
    stop(label, " must be one non-empty string.", call. = FALSE)
  }
  value
}

.mv05u_source_from_view <- function(view) {
  list(
    contract = list(
      profile = view$contract_profile,
      scientific_eligible = view$scientific_eligible
    ),
    cache_key = view$source_cache_key,
    sample_id = view$sample_id,
    cohort = view$cohort,
    representation = view$representation,
    fit_scope_id = view$fit_scope_id,
    subsample_seed = view$subsample_seed
  )
}

#' Construct one frozen MV5-U robustness view
#'
#' @param view Accepted 384-cell by 30-coordinate cell-topology view.
#' @param configuration_id One MV5-T admitted configuration identifier.
#'
#' @return A typed topology view with exactly one admitted factor changed.
#' @keywords internal
mv05u_transform_view_v1 <- function(view, configuration_id) {
  validate_topology_view(view)
  configuration_id <- .mv05u_one_string(configuration_id, "configuration_id")
  configurations <- mv05t_configuration_registry_v1()
  configuration <- configurations[
    configurations$configuration_id == configuration_id, , drop = FALSE
  ]
  if (nrow(configuration) != 1L ||
      !identical(dim(view$payload), c(384L, 30L)) ||
      !identical(view$point_ids, rownames(view$payload)) ||
      !identical(view$coordinate_ids, colnames(view$payload))) {
    stop("MV5-U source view or configuration is invalid.", call. = FALSE)
  }

  selected_ids <- if (as.integer(configuration$cells) == 384L) {
    view$point_ids
  } else {
    mv05t_nested_point_ids_v1(
      view$sample_id, view$subsample_seed, view$point_ids,
      as.integer(configuration$cells)
    )
  }
  coordinate_ids <- view$coordinate_ids[
    seq_len(as.integer(configuration$coordinates))
  ]
  payload <- view$payload[selected_ids, coordinate_ids, drop = FALSE]
  normalization <- "none"
  if (configuration$candidate_id == "cosine_chord_geometry") {
    norms <- sqrt(rowSums(payload^2))
    if (any(!is.finite(norms)) || any(norms <= 0)) {
      stop("MV5-U cosine-chord view contains a zero-norm cell.",
           call. = FALSE)
    }
    payload <- payload / norms
    normalization <- "row_l2_unit_norm_v1"
  }
  if (any(!is.finite(payload))) {
    stop("MV5-U transformed coordinates are not finite.", call. = FALSE)
  }

  transformations <- view$transformations
  transformations$mv05u_robustness <- list(
    contract_id = "mv05u_one_factor_coordinate_transform_v1",
    configuration_id = configuration_id,
    candidate_id = configuration$candidate_id,
    source_view_cache_key = view$cache_key,
    source_payload_sha256 = view$payload_sha256,
    point_selection = if (nrow(payload) == 384L) {
      "accepted_384_cell_order_preserved_v1"
    } else {
      "sha256_sample_seed_cell_nested_v1"
    },
    selected_cells = nrow(payload),
    leading_coordinates = ncol(payload),
    normalization = normalization,
    comparison_design = configuration$comparison_design,
    outcomes_authorized = FALSE
  )
  result <- .new_topology_view(
    view_id = "cell_topology_v1",
    source = .mv05u_source_from_view(view),
    point_metric = as.character(configuration$point_metric),
    payload = payload,
    point_ids = rownames(payload),
    coordinate_ids = colnames(payload),
    transformations = transformations,
    payload_sha256 = .scientific_digest(payload),
    diagnostics = list(
      duplicated_point_rows = sum(duplicated(payload)),
      mv05u_configuration_id = configuration_id
    )
  )
  validate_topology_view(result)
  result
}

.mv05u_order_pairs <- function(pairs, fold_id, seed, scope) {
  hashes <- vapply(seq_len(nrow(pairs)), function(index) {
    .mv05u_digest(list(
      contract_id = "mv05u_pair_coverage_order_v1",
      fold_id = fold_id,
      seed = as.integer(seed),
      pair_scope = scope,
      first_sample_id = pairs$first_sample_id[[index]],
      second_sample_id = pairs$second_sample_id[[index]]
    ))
  }, character(1L))
  pairs[order(hashes, pairs$first_sample_id, pairs$second_sample_id,
              method = "radix"), , drop = FALSE]
}

#' Freeze the 32-pair coverage set for one admission fold
#'
#' @param fold_record Accepted MV5-D1 fold record.
#' @param pairs_per_scope Number of training-training and query-training pairs.
#'
#' @return Canonically ordered 32-row label-closed pair table.
#' @keywords internal
mv05u_pair_coverage_v1 <- function(fold_record, pairs_per_scope = 16L) {
  mv05d1_validate_cell_fold_record_v1(fold_record)
  mv05u_pair_coverage_from_ids_v1(
    training_ids = fold_record$identity$training_ids,
    query_ids = fold_record$identity$query_ids,
    fold_id = fold_record$identity$fold_id,
    seed = fold_record$identity$seed,
    pairs_per_scope = pairs_per_scope
  )
}

#' Build MV5-U pair coverage from frozen fold identities
#'
#' @param training_ids,query_ids Disjoint sample identifiers.
#' @param fold_id Fold identity.
#' @param seed Frozen technical seed.
#' @param pairs_per_scope Must be 16.
#'
#' @return Canonically ordered 32-row pair table.
#' @keywords internal
mv05u_pair_coverage_from_ids_v1 <- function(
    training_ids, query_ids, fold_id, seed, pairs_per_scope = 16L) {
  pairs_per_scope <- as.integer(pairs_per_scope)
  if (!identical(pairs_per_scope, 16L)) {
    stop("MV5-U requires exactly 16 pairs per retrieval scope.",
         call. = FALSE)
  }
  training <- sort(unique(as.character(training_ids)), method = "radix")
  query <- sort(unique(as.character(query_ids)), method = "radix")
  fold_id <- .mv05u_one_string(fold_id, "fold_id")
  seed <- as.integer(seed)
  if (length(training) < 7L || length(query) < 1L ||
      length(intersect(training, query)) || length(seed) != 1L || is.na(seed)) {
    stop("MV5-U training/query identities are invalid.", call. = FALSE)
  }
  training_grid <- utils::combn(training, 2L)
  training_pairs <- data.frame(
    first_sample_id = training_grid[1L, ],
    second_sample_id = training_grid[2L, ],
    pair_scope = "training_training_unordered",
    stringsAsFactors = FALSE
  )
  query_pairs <- expand.grid(
    first_sample_id = query,
    second_sample_id = training,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  query_pairs$pair_scope <- "heldout_training_directed"
  training_pairs <- head(.mv05u_order_pairs(
    training_pairs, fold_id, seed, "training_training_unordered"
  ), pairs_per_scope)
  query_pairs <- head(.mv05u_order_pairs(
    query_pairs, fold_id, seed, "heldout_training_directed"
  ), pairs_per_scope)
  result <- rbind(training_pairs, query_pairs)
  result$pair_ordinal <- seq_len(nrow(result))
  result$pair_request_id <- vapply(seq_len(nrow(result)), function(index) {
    paste0("mv05u_pair_v1:", .mv05u_digest(list(
      fold_id = fold_id, seed = as.integer(seed),
      pair_scope = result$pair_scope[[index]],
      first_sample_id = result$first_sample_id[[index]],
      second_sample_id = result$second_sample_id[[index]]
    )))
  }, character(1L))
  result$outcome_label_state <- "closed"
  result$biological_outcomes_computed <- FALSE
  result <- result[c(
    "pair_request_id", "pair_ordinal", "pair_scope", "first_sample_id",
    "second_sample_id", "outcome_label_state",
    "biological_outcomes_computed"
  )]
  if (nrow(result) != 32L || anyDuplicated(result$pair_request_id) ||
      any(table(result$pair_scope) != 16L)) {
    stop("MV5-U pair coverage is incomplete.", call. = FALSE)
  }
  result
}

#' Validate a corrected PH result against a variable-size MST oracle
#'
#' @param result Result from `run_topology_view_ph()`.
#' @param view Matching MV5-U topology view.
#' @param tolerance Optional maximum absolute H0 death error.
#'
#' @return H0/H1 oracle metrics.
#' @keywords internal
mv05u_validate_ph_result_v1 <- function(result, view, tolerance = NULL) {
  validate_topology_view(view)
  diagram <- result$diagram
  point_count <- nrow(view$payload)
  if (!inherits(result, "scph_topology_result_v1") || !is.matrix(diagram) ||
      !identical(result$provenance$view_cache_key, view$cache_key) ||
      result$provenance$point_count != point_count ||
      result$provenance$max_dim != 1L ||
      result$provenance$threshold != -1 || result$provenance$field != 2L ||
      result$provenance$invalid_interval_count != 0L ||
      result$provenance$zero_persistence_count != 0L ||
      result$provenance$essential_h0_count != 1L) {
    stop("MV5-U PH result violates the corrected typed contract.",
         call. = FALSE)
  }
  h0 <- diagram[diagram[, "dimension"] == 0, , drop = FALSE]
  h1 <- diagram[diagram[, "dimension"] == 1, , drop = FALSE]
  finite_h0 <- sort(h0[is.finite(h0[, "death"]), "death"], method = "radix")
  oracle <- mv05d2_h0_mst_deaths_v1(view$payload)
  if (length(finite_h0) != point_count - 1L || nrow(h0) != point_count ||
      any(!is.finite(h1[, "death"])) || any(h1[, "death"] <= h1[, "birth"])) {
    stop("MV5-U diagram has invalid H0/H1 interval structure.",
         call. = FALSE)
  }
  if (is.null(tolerance)) tolerance <- max(1e-7, max(oracle) * 1e-7)
  maximum_error <- max(abs(finite_h0 - oracle))
  if (!is.finite(maximum_error) || maximum_error > tolerance) {
    stop("MV5-U finite H0 deaths disagree with the MST oracle.",
         call. = FALSE)
  }
  list(
    contract_id = "mv05u_variable_size_h0_mst_oracle_v1",
    point_count = point_count,
    finite_h0_intervals = length(finite_h0),
    finite_h1_intervals = nrow(h1),
    maximum_absolute_error = maximum_error,
    tolerance = tolerance,
    passed = TRUE
  )
}

#' Validate the bound MV5-U execution queue
#'
#' @param queue Data frame produced by the prospective MV5-U binder.
#'
#' @return `queue`, invisibly.
#' @keywords internal
mv05u_validate_execution_queue_v1 <- function(queue) {
  required <- c(
    "contract_id", "admission_unit_id", "execution_order", "fold_id",
    "fold_role", "training_samples", "seed", "representation",
    "configuration_id", "candidate_id", "cells", "coordinates",
    "point_metric", "mv05t_source_freeze_sha256", "source_freeze_sha256",
    "mv05t_queue_sha256",
    "implementation_sha256", "prospective_head",
    "python_executable_sha256", "python_version", "persim_version",
    "numpy_version", "scipy_version", "pair_coverage_per_scope",
    "view_count", "landscape_dimensions", "labels_opened",
    "outcomes_computed", "admission_executed"
  )
  configurations <- mv05t_configuration_registry_v1()
  if (!is.data.frame(queue) || !all(required %in% names(queue)) ||
      nrow(queue) != 24L || anyDuplicated(queue$admission_unit_id) ||
      !identical(as.integer(queue$execution_order), 1:24) ||
      length(unique(queue$fold_id)) != 3L ||
      length(unique(queue$representation)) != 2L ||
      length(unique(queue$configuration_id)) != 4L ||
      any(table(queue$fold_id) != 8L) ||
      !setequal(queue$configuration_id, configurations$configuration_id) ||
      any(queue$seed != 20260805L) ||
      any(queue$pair_coverage_per_scope != 16L) ||
      any(queue$view_count != 90L) ||
      any(queue$landscape_dimensions != "H0|H1") ||
      length(unique(queue$mv05t_source_freeze_sha256)) != 1L ||
      length(unique(queue$source_freeze_sha256)) != 1L ||
      length(unique(queue$mv05t_queue_sha256)) != 1L ||
      length(unique(queue$implementation_sha256)) != 1L ||
      length(unique(queue$prospective_head)) != 1L ||
      length(unique(queue$python_executable_sha256)) != 1L ||
      length(unique(queue$python_version)) != 1L ||
      length(unique(queue$persim_version)) != 1L ||
      length(unique(queue$numpy_version)) != 1L ||
      length(unique(queue$scipy_version)) != 1L ||
      any(!grepl("^[0-9a-f]{64}$", queue$mv05t_source_freeze_sha256)) ||
      any(!grepl("^[0-9a-f]{64}$", queue$source_freeze_sha256)) ||
      any(!grepl("^[0-9a-f]{64}$", queue$mv05t_queue_sha256)) ||
      any(!grepl("^[0-9a-f]{64}$", queue$implementation_sha256)) ||
      any(!grepl("^[0-9a-f]{64}$", queue$python_executable_sha256)) ||
      any(!grepl("^[0-9a-f]{40}$", queue$prospective_head)) ||
      any(queue$persim_version != "0.3.8") ||
      any(queue$labels_opened) || any(queue$outcomes_computed) ||
      any(queue$admission_executed)) {
    stop("MV5-U execution queue violates the 24-unit admission contract.",
         call. = FALSE)
  }
  invisible(queue)
}

#' Derive the resource-only MV5-U continuation gate
#'
#' @param metrics Complete 24-row admission resource table.
#' @param new_private_bytes Total newly published private bytes.
#'
#' @return One-row no-outcome decision table.
#' @keywords internal
mv05u_resource_decision_v1 <- function(metrics, new_private_bytes) {
  required <- c(
    "admission_unit_id", "disposition", "elapsed_seconds",
    "peak_process_tree_rss_bytes", "completed_views",
    "landscape_pair_rows", "energy_pair_rows", "labels_opened",
    "outcomes_computed"
  )
  new_private_bytes <- as.numeric(new_private_bytes)
  valid <- is.data.frame(metrics) && all(required %in% names(metrics)) &&
    nrow(metrics) == 24L && !anyDuplicated(metrics$admission_unit_id) &&
    all(metrics$disposition == "completed") &&
    all(metrics$completed_views == 90L) &&
    all(metrics$landscape_pair_rows == 64L) &&
    all(metrics$energy_pair_rows == 32L) &&
    all(is.finite(metrics$elapsed_seconds)) &&
    all(metrics$elapsed_seconds <= 600) &&
    sum(metrics$elapsed_seconds) <= 7200 &&
    all(metrics$peak_process_tree_rss_bytes <= 4 * 1024^3) &&
    is.finite(new_private_bytes) && new_private_bytes <= 2 * 1024^3 &&
    !any(metrics$labels_opened) && !any(metrics$outcomes_computed)
  data.frame(
    contract_id = "mv05u_resource_decision_v1",
    admission_complete = valid,
    full_robustness_authorized = FALSE,
    next_action = if (valid) {
      "prospectively_freeze_streamed_full_robustness_execution_gate"
    } else {
      "stop_and_review_admission_failure"
    },
    total_elapsed_seconds = if (is.data.frame(metrics)) {
      sum(metrics$elapsed_seconds)
    } else NA_real_,
    maximum_unit_elapsed_seconds = if (is.data.frame(metrics)) {
      max(metrics$elapsed_seconds)
    } else NA_real_,
    maximum_process_tree_rss_bytes = if (is.data.frame(metrics)) {
      max(metrics$peak_process_tree_rss_bytes)
    } else NA_real_,
    new_private_bytes = new_private_bytes,
    labels_opened = FALSE,
    outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
}
