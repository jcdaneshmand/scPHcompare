# MV5-V streamed full-robustness production prefreeze contracts.

.mv05v_digest <- function(value) {
  digest::digest(value, algo = "sha256", serialize = TRUE)
}

.mv05v_hash_ok <- function(values) {
  is.character(values) && length(values) > 0L &&
    all(!is.na(values) & grepl("^[0-9a-f]{64}$", values))
}

.mv05v_assert_closed <- function(values, label) {
  if (any(tolower(as.character(values)) != "false")) {
    stop("MV5-V requires label-closed `", label, "`.", call. = FALSE)
  }
}

#' Validate the four selection-resistant robustness configurations.
#'
#' @param configurations The committed MV5-T configuration registry.
#' @return `configurations`, ordered canonically.
#' @keywords internal
mv05v_validate_configurations_v1 <- function(configurations) {
  required <- c(
    "configuration_id", "candidate_id", "cells", "coordinates",
    "point_metric", "comparison_design", "baseline_reused", "role",
    "outcomes_authorized"
  )
  expected <- c(
    "cells384_pc20_euclidean_v1", "cells384_pc30_cosine_chord_v1",
    "nested_cells_192_pc30_euclidean_v1",
    "nested_cells_256_pc30_euclidean_v1"
  )
  if (!is.data.frame(configurations) ||
      !all(required %in% names(configurations)) ||
      nrow(configurations) != 4L ||
      anyDuplicated(configurations$configuration_id) ||
      !setequal(configurations$configuration_id, expected) ||
      any(as.logical(configurations$outcomes_authorized)) ||
      any(!as.logical(configurations$baseline_reused)) ||
      any(configurations$role != "postoutcome_secondary_sensitivity")) {
    stop("MV5-V requires exactly the four committed MV5-T configurations.",
         call. = FALSE)
  }
  configurations <- configurations[
    order(configurations$configuration_id, method = "radix"), , drop = FALSE
  ]
  configurations$configuration_order <- seq_len(nrow(configurations))
  rownames(configurations) <- NULL
  configurations
}

#' Summarize the accepted held-out-query pair axes without copying labels.
#'
#' @param sct_pairs Accepted MV5-D4 pair manifest.
#' @param integrated_pairs Accepted MV5-I pair manifest.
#' @param max_rows_per_subchunk Maximum exact-landscape requests per subchunk.
#' @return One label-closed row per fold/seed/representation group.
#' @keywords internal
mv05v_build_pair_scope_v1 <- function(sct_pairs, integrated_pairs,
                                       max_rows_per_subchunk = 250L) {
  max_rows_per_subchunk <- as.integer(max_rows_per_subchunk)
  required <- c(
    "group_id", "group_order", "fold_id", "seed", "homology_dimension",
    "query_sample_id", "training_sample_id", "representation", "pair_scope",
    "outcome_label_state", "biological_outcomes_computed"
  )
  inputs <- list(sct_whole = sct_pairs,
                 inductive_integrated = integrated_pairs)
  expected_rows <- 70700L
  rows <- list()
  cursor <- 0L
  for (representation in names(inputs)) {
    pairs <- inputs[[representation]]
    if (!is.data.frame(pairs) || !all(required %in% names(pairs)) ||
        nrow(pairs) != expected_rows ||
        any(pairs$representation != representation) ||
        any(pairs$pair_scope != "held_out_query_to_training_reference") ||
        any(pairs$outcome_label_state != "closed") ||
        any(as.logical(pairs$biological_outcomes_computed)) ||
        !setequal(pairs$homology_dimension, c("H0", "H1")) ||
        any(c("tissue", "approach") %in% names(pairs))) {
      stop("MV5-V pair-scope source violates the accepted boundary.",
           call. = FALSE)
    }
    groups <- split(pairs, pairs$group_id)
    if (length(groups) != 75L) {
      stop("MV5-V requires 75 accepted groups per representation.",
           call. = FALSE)
    }
    for (group in groups) {
      key <- paste(group$query_sample_id, group$training_sample_id, sep = "\r")
      counts <- table(group$homology_dimension)
      biological_pairs <- length(unique(key))
      if (any(counts[c("H0", "H1")] != biological_pairs) ||
          any(table(key) != 2L)) {
        stop("MV5-V pair scope must contain H0/H1 for every biological pair.",
             call. = FALSE)
      }
      axis <- sort(unique(paste(
        group$query_sample_id, group$training_sample_id, sep = "\r"
      )), method = "radix")
      cursor <- cursor + 1L
      rows[[cursor]] <- data.frame(
        contract_id = "mv05v_base_pair_scope_v1",
        source_group_id = group$group_id[[1L]],
        group_order = as.integer(group$group_order[[1L]]),
        fold_id = group$fold_id[[1L]], seed = as.integer(group$seed[[1L]]),
        representation = representation,
        heldout_samples = length(unique(group$query_sample_id)),
        training_samples = length(unique(group$training_sample_id)),
        biological_pairs = biological_pairs,
        landscape_request_rows = nrow(group),
        landscape_subchunks = sum(ceiling(
          as.numeric(counts[c("H0", "H1")]) / max_rows_per_subchunk
        )),
        max_rows_per_subchunk = max_rows_per_subchunk,
        base_pair_axis_sha256 = .mv05v_digest(axis),
        pair_scope = "held_out_query_to_training_reference",
        outcome_label_state = "closed", outcomes_computed = FALSE,
        stringsAsFactors = FALSE
      )
    }
  }
  result <- do.call(rbind, rows)
  result <- result[order(
    result$representation, result$group_order, method = "radix"
  ), , drop = FALSE]
  rownames(result) <- NULL
  if (nrow(result) != 150L ||
      any(table(result$representation) != 75L) ||
      sum(result$biological_pairs) != 70700L ||
      sum(result$landscape_request_rows) != 141400L ||
      sum(result$landscape_subchunks) != 720L ||
      !.mv05v_hash_ok(result$base_pair_axis_sha256)) {
    stop("MV5-V pair-scope totals do not reconcile.", call. = FALSE)
  }
  result
}

#' Build the exact 600-group streamed robustness execution queue.
#'
#' @param pair_scope Output of [mv05v_build_pair_scope_v1()].
#' @param configurations The four accepted MV5-T configurations.
#' @param private_inventory The accepted 150-file coordinate inventory.
#' @param source_freeze_sha256 Bound MV5-V source-freeze hash.
#' @param implementation_sha256 Bound full engine hash.
#' @param prospective_head Prospective engine Git commit.
#' @return The exact label-closed group queue.
#' @keywords internal
mv05v_build_group_queue_v1 <- function(
    pair_scope, configurations, private_inventory,
    source_freeze_sha256, implementation_sha256, prospective_head) {
  configurations <- mv05v_validate_configurations_v1(configurations)
  required_scope <- c(
    "fold_id", "seed", "representation", "heldout_samples",
    "training_samples", "biological_pairs", "landscape_request_rows",
    "landscape_subchunks", "base_pair_axis_sha256"
  )
  required_inventory <- c(
    "source_type", "fold_study", "seed", "private_locator", "sha256",
    "labels_opened", "outcomes_computed"
  )
  if (!is.data.frame(pair_scope) || nrow(pair_scope) != 150L ||
      !all(required_scope %in% names(pair_scope)) ||
      !is.data.frame(private_inventory) || nrow(private_inventory) != 150L ||
      !all(required_inventory %in% names(private_inventory)) ||
      !.mv05v_hash_ok(private_inventory$sha256) ||
      any(as.logical(private_inventory$labels_opened)) ||
      any(as.logical(private_inventory$outcomes_computed)) ||
      !.mv05v_hash_ok(source_freeze_sha256) ||
      !.mv05v_hash_ok(implementation_sha256) ||
      length(prospective_head) != 1L ||
      !grepl("^[0-9a-f]{40}$", prospective_head)) {
    stop("MV5-V group-queue inputs are invalid.", call. = FALSE)
  }
  inventory <- private_inventory
  inventory$representation <- ifelse(
    inventory$source_type == "sct", "sct_whole",
    ifelse(inventory$source_type == "integrated",
           "inductive_integrated", NA_character_)
  )
  inventory$fold_id <- paste0("large_loso_v1:", inventory$fold_study)
  joined <- merge(
    pair_scope, inventory,
    by = c("fold_id", "seed", "representation"), sort = FALSE,
    suffixes = c("", "_source")
  )
  if (nrow(joined) != 150L || any(is.na(joined$representation))) {
    stop("MV5-V coordinate sources do not join one-to-one to pair scope.",
         call. = FALSE)
  }
  rows <- lapply(seq_len(nrow(joined)), function(group_index) {
    group <- joined[group_index, , drop = FALSE]
    do.call(rbind, lapply(seq_len(nrow(configurations)), function(config_index) {
      config <- configurations[config_index, , drop = FALSE]
      identity <- list(
        contract_id = "mv05v_robustness_group_v1",
        fold_id = group$fold_id, seed = as.integer(group$seed),
        representation = group$representation,
        configuration_id = config$configuration_id,
        coordinate_source_sha256 = group$sha256,
        base_pair_axis_sha256 = group$base_pair_axis_sha256
      )
      data.frame(
        contract_id = "mv05v_full_group_queue_v1",
        robustness_group_id = paste0(
          "mv05v_group_v1:", .mv05v_digest(identity)
        ),
        fold_id = group$fold_id, seed = as.integer(group$seed),
        representation = group$representation,
        configuration_id = config$configuration_id,
        candidate_id = config$candidate_id,
        configuration_order = as.integer(config$configuration_order),
        cells = as.integer(config$cells),
        coordinates = as.integer(config$coordinates),
        point_metric = config$point_metric,
        heldout_samples = as.integer(group$heldout_samples),
        training_samples = as.integer(group$training_samples),
        view_count = 90L,
        biological_pairs = as.integer(group$biological_pairs),
        landscape_request_rows = as.integer(group$landscape_request_rows),
        landscape_subchunks = as.integer(group$landscape_subchunks),
        energy_request_rows = as.integer(group$biological_pairs),
        assembled_method_rows = as.integer(group$biological_pairs) * 4L,
        deterministic_repeat_required = FALSE,
        private_locator = group$private_locator,
        coordinate_source_sha256 = group$sha256,
        base_pair_axis_sha256 = group$base_pair_axis_sha256,
        source_freeze_sha256 = source_freeze_sha256,
        implementation_sha256 = implementation_sha256,
        prospective_head = prospective_head,
        outcome_label_state = "closed", outcomes_computed = FALSE,
        execution_authorized = FALSE, execution_completed = FALSE,
        stringsAsFactors = FALSE
      )
    }))
  })
  result <- do.call(rbind, rows)
  result <- result[order(
    result$configuration_order, result$representation, result$fold_id,
    result$seed, method = "radix"
  ), , drop = FALSE]
  result$execution_order <- seq_len(nrow(result))
  repeat_key <- paste(
    result$configuration_id, result$representation, sep = "\r"
  )
  maximum <- ave(
    result$training_samples, repeat_key,
    FUN = function(values) values == max(values)
  )
  lexical <- ave(
    seq_len(nrow(result)), repeat_key,
    FUN = function(indices) {
      chosen <- indices[order(
        -result$training_samples[indices], result$fold_id[indices],
        result$seed[indices], method = "radix"
      )][[1L]]
      indices == chosen
    }
  )
  result$deterministic_repeat_required <- as.logical(maximum) & as.logical(lexical)
  rownames(result) <- NULL
  mv05v_validate_group_queue_v1(result)
  result
}

#' Validate the exact MV5-V full group queue.
#'
#' @param queue A candidate MV5-V queue.
#' @return `invisible(queue)` on success.
#' @keywords internal
mv05v_validate_group_queue_v1 <- function(queue) {
  required <- c(
    "robustness_group_id", "execution_order", "fold_id", "seed",
    "representation", "configuration_id", "cells", "coordinates",
    "point_metric", "heldout_samples", "training_samples", "view_count",
    "biological_pairs", "landscape_request_rows", "landscape_subchunks",
    "energy_request_rows", "assembled_method_rows",
    "deterministic_repeat_required", "coordinate_source_sha256",
    "base_pair_axis_sha256", "source_freeze_sha256",
    "implementation_sha256", "prospective_head", "outcome_label_state",
    "outcomes_computed", "execution_authorized", "execution_completed"
  )
  if (!is.data.frame(queue) || !all(required %in% names(queue)) ||
      nrow(queue) != 600L || anyDuplicated(queue$robustness_group_id) ||
      !identical(as.integer(queue$execution_order), 1:600) ||
      length(unique(queue$fold_id)) != 15L ||
      !identical(sort(unique(as.integer(queue$seed))), 20260805:20260809) ||
      !setequal(queue$representation,
                c("sct_whole", "inductive_integrated")) ||
      length(unique(queue$configuration_id)) != 4L ||
      any(table(queue$configuration_id, queue$representation) != 75L) ||
      any(queue$view_count != 90L) ||
      sum(queue$view_count) != 54000L ||
      sum(queue$biological_pairs) != 282800L ||
      sum(queue$landscape_request_rows) != 565600L ||
      sum(queue$landscape_subchunks) != 2880L ||
      sum(queue$energy_request_rows) != 282800L ||
      sum(queue$assembled_method_rows) != 1131200L ||
      sum(as.logical(queue$deterministic_repeat_required)) != 8L ||
      any(queue$outcome_label_state != "closed") ||
      any(as.logical(queue$outcomes_computed)) ||
      any(as.logical(queue$execution_authorized)) ||
      any(as.logical(queue$execution_completed)) ||
      any(c("tissue", "approach") %in% names(queue)) ||
      !.mv05v_hash_ok(queue$coordinate_source_sha256) ||
      !.mv05v_hash_ok(queue$base_pair_axis_sha256) ||
      length(unique(queue$source_freeze_sha256)) != 1L ||
      length(unique(queue$implementation_sha256)) != 1L ||
      length(unique(queue$prospective_head)) != 1L) {
    stop("MV5-V full group queue violates its frozen contract.", call. = FALSE)
  }
  invisible(queue)
}

#' Construct the conservative MV5-V resource projection.
#'
#' @param historical Named stage seconds for the six accepted stage/representation
#'   measurements.
#' @param conservative_storage_bytes MV5-T four-setting storage projection.
#' @return Stage-level resource projection and total decision attributes.
#' @keywords internal
mv05v_resource_projection_v1 <- function(
    historical, conservative_storage_bytes = 10180000000) {
  required <- c(
    "sct_ph", "integrated_ph", "sct_landscape", "integrated_landscape",
    "sct_assembly", "integrated_assembly"
  )
  if (!is.numeric(historical) || !setequal(names(historical), required) ||
      any(!is.finite(historical) | historical <= 0)) {
    stop("MV5-V resource projection requires six positive stage measurements.",
         call. = FALSE)
  }
  stages <- data.frame(
    contract_id = "mv05v_resource_projection_v1",
    stage_id = c("ph", "landscape", "assembly"),
    sct_seconds_per_configuration = c(
      historical[["sct_ph"]], historical[["sct_landscape"]],
      historical[["sct_assembly"]]
    ),
    integrated_seconds_per_configuration = c(
      historical[["integrated_ph"]], historical[["integrated_landscape"]],
      historical[["integrated_assembly"]]
    ),
    configurations = 4L, workers = 1L,
    labels_opened = FALSE, outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
  stages$projected_seconds <- 4 * (
    stages$sct_seconds_per_configuration +
      stages$integrated_seconds_per_configuration
  )
  total <- sum(stages$projected_seconds)
  attr(stages, "projected_worker_seconds") <- total
  attr(stages, "projected_worker_hours") <- total / 3600
  attr(stages, "worker_hour_cap") <- 30
  attr(stages, "storage_projection_bytes") <-
    as.numeric(conservative_storage_bytes)
  attr(stages, "storage_cap_bytes") <- 16 * 1024^3
  stages
}

#' Issue the bounded MV5-V prefreeze decision.
#'
#' @param queue Valid MV5-V group queue.
#' @param projection Output of [mv05v_resource_projection_v1()].
#' @return One-row decision record.
#' @keywords internal
mv05v_prefreeze_decision_v1 <- function(queue, projection) {
  mv05v_validate_group_queue_v1(queue)
  hours <- attr(projection, "projected_worker_hours")
  cap <- attr(projection, "worker_hour_cap")
  storage <- attr(projection, "storage_projection_bytes")
  storage_cap <- attr(projection, "storage_cap_bytes")
  passed <- is.finite(hours) && hours <= cap && storage <= storage_cap
  data.frame(
    contract_id = "mv05v_prefreeze_decision_v1",
    gate_complete = passed,
    full_execution_authorized = FALSE,
    next_action = if (passed) {
      "commit_engine_and_bind_queue_before_configuration_stratified_execution"
    } else "revise_or_stop_resource_envelope",
    robustness_groups = nrow(queue), views = sum(queue$view_count),
    landscape_request_rows = sum(queue$landscape_request_rows),
    energy_request_rows = sum(queue$energy_request_rows),
    landscape_subchunks = sum(queue$landscape_subchunks),
    projected_worker_hours = hours, worker_hour_cap = cap,
    storage_projection_bytes = storage, storage_cap_bytes = storage_cap,
    labels_opened = FALSE, outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
}
