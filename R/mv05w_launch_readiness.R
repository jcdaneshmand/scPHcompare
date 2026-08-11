# MV5-W full-robustness launch-readiness helpers.

.mv05w_digest <- function(value) {
  digest::digest(value, algo = "sha256", serialize = TRUE)
}

#' Validate an MV5-W queue with exactly one prospectively authorized smoke.
#'
#' @param queue Rebound MV5-V group queue with runtime identities.
#' @return `invisible(queue)` on success.
#' @keywords internal
mv05w_validate_smoke_queue_v1 <- function(queue) {
  runtime <- c(
    "python_executable_sha256", "python_version", "persim_version",
    "numpy_version", "scipy_version"
  )
  if (!is.data.frame(queue) || !all(runtime %in% names(queue)) ||
      nrow(queue) != 600L ||
      sum(as.logical(queue$execution_authorized)) != 1L ||
      which(as.logical(queue$execution_authorized)) != 1L ||
      any(as.logical(queue$execution_completed)) ||
      length(unique(queue$python_executable_sha256)) != 1L ||
      !grepl("^[0-9a-f]{64}$", unique(queue$python_executable_sha256))) {
    stop("MV5-W requires exactly its first bound group as the smoke queue.",
         call. = FALSE)
  }
  frozen <- queue
  frozen$execution_authorized <- FALSE
  mv05v_validate_group_queue_v1(frozen)
  invisible(queue)
}

#' Build the complete heldout-to-training pair coverage for one group.
#'
#' @param fold_record Accepted label-closed fold record.
#' @param robustness_group_id Bound MV5-W group identity.
#' @param expected_axis_sha256 Expected accepted biological pair-axis hash.
#' @param expected_view_count Expected complete sample count; 90 in production.
#' @return Complete directed biological-pair requests.
#' @keywords internal
mv05w_full_pair_coverage_v1 <- function(
    fold_record, robustness_group_id, expected_axis_sha256,
    expected_view_count = 90L) {
  training <- sort(unique(as.character(fold_record$identity$training_ids)),
                   method = "radix")
  query <- sort(setdiff(
    names(fold_record$payload$cell_views), training
  ), method = "radix")
  if (!length(training) || !length(query) ||
      length(intersect(training, query)) ||
      length(c(training, query)) != as.integer(expected_view_count)) {
    stop("MV5-W fold does not expose its complete sample query axis.",
         call. = FALSE)
  }
  axis <- expand.grid(
    first_sample_id = query, second_sample_id = training,
    stringsAsFactors = FALSE
  )
  axis_key <- sort(paste(
    axis$first_sample_id, axis$second_sample_id, sep = "\r"
  ), method = "radix")
  if (!identical(.mv05w_digest(axis_key), expected_axis_sha256)) {
    stop("MV5-W complete biological pair axis is stale.", call. = FALSE)
  }
  axis$pair_request_id <- paste0(
    "mv05w_pair_v1:", vapply(seq_len(nrow(axis)), function(index) {
      .mv05w_digest(list(
        contract_id = "mv05w_full_pair_v1",
        robustness_group_id = robustness_group_id,
        first_sample_id = axis$first_sample_id[[index]],
        second_sample_id = axis$second_sample_id[[index]]
      ))
    }, character(1L))
  )
  axis <- axis[order(axis$pair_request_id, method = "radix"), , drop = FALSE]
  axis$pair_ordinal <- seq_len(nrow(axis))
  axis$pair_scope <- "held_out_query_to_training_reference"
  axis$outcome_label_state <- "closed"
  axis$biological_outcomes_computed <- FALSE
  axis <- axis[c(
    "pair_request_id", "pair_ordinal", "pair_scope", "first_sample_id",
    "second_sample_id", "outcome_label_state",
    "biological_outcomes_computed"
  )]
  rownames(axis) <- NULL
  axis
}

#' Assemble the four label-closed robustness method rows.
#'
#' @param landscape_pairs Complete exact H0/H1 landscape rows.
#' @param energy_pairs Complete matched energy rows.
#' @param robustness_group_id Bound group identity.
#' @return Four rows per biological pair, without labels or outcomes.
#' @keywords internal
mv05w_assemble_method_rows_v1 <- function(
    landscape_pairs, energy_pairs, robustness_group_id) {
  required_landscape <- c(
    "pair_request_id", "first_sample_id", "second_sample_id",
    "homology_dimension", "distance", "exact", "all_active_levels",
    "level_cap_applied"
  )
  required_energy <- c(
    "pair_request_id", "first_sample_id", "second_sample_id", "distance"
  )
  if (!is.data.frame(landscape_pairs) ||
      !all(required_landscape %in% names(landscape_pairs)) ||
      !is.data.frame(energy_pairs) ||
      !all(required_energy %in% names(energy_pairs)) ||
      anyDuplicated(paste(landscape_pairs$pair_request_id,
                          landscape_pairs$homology_dimension)) ||
      anyDuplicated(energy_pairs$pair_request_id) ||
      !setequal(landscape_pairs$homology_dimension, c("H0", "H1"))) {
    stop("MV5-W method assembly inputs are incomplete.", call. = FALSE)
  }
  strict <- function(values) tolower(as.character(values)) == "true"
  if (any(!strict(landscape_pairs$exact)) ||
      any(!strict(landscape_pairs$all_active_levels)) ||
      any(strict(landscape_pairs$level_cap_applied))) {
    stop("MV5-W method assembly requires exact all-active landscapes.",
         call. = FALSE)
  }
  h0 <- landscape_pairs[landscape_pairs$homology_dimension == "H0", ]
  h1 <- landscape_pairs[landscape_pairs$homology_dimension == "H1", ]
  h1 <- h1[match(h0$pair_request_id, h1$pair_request_id), ]
  energy <- energy_pairs[match(h0$pair_request_id, energy_pairs$pair_request_id), ]
  if (anyNA(h1$pair_request_id) || anyNA(energy$pair_request_id) ||
      any(h0$first_sample_id != h1$first_sample_id) ||
      any(h0$second_sample_id != h1$second_sample_id) ||
      any(h0$first_sample_id != energy$first_sample_id) ||
      any(h0$second_sample_id != energy$second_sample_id)) {
    stop("MV5-W method assembly pair axes disagree.", call. = FALSE)
  }
  base <- h0[c("pair_request_id", "first_sample_id", "second_sample_id")]
  make <- function(method_id, distance) data.frame(
    contract_id = "mv05w_label_closed_method_row_v1",
    robustness_group_id = robustness_group_id,
    pair_request_id = base$pair_request_id,
    query_sample_id = base$first_sample_id,
    training_sample_id = base$second_sample_id,
    method_id = method_id, distance = as.numeric(distance),
    outcome_label_state = "closed", outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
  result <- do.call(rbind, list(
    make("cell_landscape_h0_v1", h0$distance),
    make("cell_landscape_h1_v1", h1$distance),
    make("cell_landscape_h0_h1_raw_euclidean_v1",
         sqrt(h0$distance^2 + h1$distance^2)),
    make("cell_distribution_energy_shared_pca_v1", energy$distance)
  ))
  result <- result[order(
    result$method_id, result$pair_request_id, method = "radix"
  ), , drop = FALSE]
  rownames(result) <- NULL
  if (nrow(result) != 4L * nrow(energy_pairs) ||
      any(!is.finite(result$distance) | result$distance < 0)) {
    stop("MV5-W method assembly cardinality is invalid.", call. = FALSE)
  }
  result
}
