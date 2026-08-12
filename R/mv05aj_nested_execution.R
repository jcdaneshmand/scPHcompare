# MV5-AJ label-closed nested-256 configuration execution contracts.

.mv05aj_hash_ok <- function(values, width = 64L) {
  is.character(values) && length(values) > 0L &&
    all(!is.na(values) & grepl(paste0("^[0-9a-f]{", width, "}$"), values))
}

#' Validate the prospectively bound 150-group nested-256 execution queue.
mv05aj_validate_nested_queue_v1 <- function(queue) {
  required <- c(
    "contract_id", "robustness_group_id", "fold_id", "seed",
    "representation", "configuration_id", "configuration_order", "cells",
    "coordinates", "point_metric", "view_count", "biological_pairs",
    "landscape_request_rows", "landscape_subchunks", "energy_request_rows",
    "assembled_method_rows", "training_samples",
    "deterministic_repeat_required", "coordinate_source_sha256",
    "base_pair_axis_sha256", "source_freeze_sha256",
    "implementation_sha256", "prospective_head",
    "python_executable_sha256", "python_version", "persim_version",
    "numpy_version", "scipy_version", "configuration_execution_order",
    "execution_authorized", "execution_completed", "labels_opened",
    "rankings_computed", "outcomes_computed")
  id <- "nested_cells_256_pc30_euclidean_v1"
  if (!is.data.frame(queue) || !all(required %in% names(queue)) ||
      nrow(queue) != 150L || anyDuplicated(queue$robustness_group_id) ||
      any(queue$contract_id != "mv05aj_nested256_execution_queue_v1") ||
      any(queue$configuration_id != id) ||
      any(as.integer(queue$configuration_order) != 4L) ||
      length(unique(queue$fold_id)) != 15L ||
      length(unique(queue$seed)) != 5L ||
      !setequal(queue$representation,
                c("sct_whole", "inductive_integrated")) ||
      any(table(queue$fold_id, queue$representation) != 5L) ||
      any(as.integer(queue$cells) != 256L) ||
      any(as.integer(queue$coordinates) != 30L) ||
      any(queue$point_metric != "euclidean") ||
      any(as.integer(queue$view_count) != 90L) ||
      sum(as.integer(queue$biological_pairs)) != 70700L ||
      sum(as.integer(queue$landscape_request_rows)) != 141400L ||
      sum(as.integer(queue$landscape_subchunks)) != 720L ||
      sum(as.integer(queue$energy_request_rows)) != 70700L ||
      sum(as.integer(queue$assembled_method_rows)) != 282800L ||
      sum(as.logical(queue$deterministic_repeat_required)) != 2L ||
      !setequal(
        queue$training_samples[as.logical(queue$deterministic_repeat_required)],
        range(queue$training_samples)
      ) ||
      !setequal(
        queue$representation[as.logical(queue$deterministic_repeat_required)],
        c("sct_whole", "inductive_integrated")
      ) ||
      any(!.mv05aj_hash_ok(queue$coordinate_source_sha256)) ||
      any(!.mv05aj_hash_ok(queue$base_pair_axis_sha256)) ||
      length(unique(queue$source_freeze_sha256)) != 1L ||
      !.mv05aj_hash_ok(unique(queue$source_freeze_sha256)) ||
      length(unique(queue$implementation_sha256)) != 1L ||
      !.mv05aj_hash_ok(unique(queue$implementation_sha256)) ||
      length(unique(queue$prospective_head)) != 1L ||
      !.mv05aj_hash_ok(unique(queue$prospective_head), 40L) ||
      length(unique(queue$python_executable_sha256)) != 1L ||
      !.mv05aj_hash_ok(unique(queue$python_executable_sha256)) ||
      !identical(as.integer(queue$configuration_execution_order), 1:150) ||
      !all(as.logical(queue$execution_authorized)) ||
      any(as.logical(queue$execution_completed)) ||
      any(as.logical(queue$labels_opened)) ||
      any(as.logical(queue$rankings_computed)) ||
      any(as.logical(queue$outcomes_computed)) ||
      any(c("tissue", "approach", "estimate", "p_value") %in% names(queue))) {
    stop("MV5-AJ requires exactly its 150 unopened nested-256 groups.",
         call. = FALSE)
  }
  invisible(queue)
}

#' Summarize the bound nested-256 calculation.
mv05aj_nested_scope_v1 <- function(queue) {
  mv05aj_validate_nested_queue_v1(queue)
  data.frame(
    contract_id = "mv05aj_nested256_execution_scope_v1",
    configuration_id = unique(queue$configuration_id), groups = nrow(queue),
    representations = length(unique(queue$representation)),
    folds = length(unique(queue$fold_id)), seeds = length(unique(queue$seed)),
    views = sum(queue$view_count),
    biological_pairs = sum(queue$biological_pairs),
    landscape_rows = sum(queue$landscape_request_rows),
    landscape_subchunks = sum(queue$landscape_subchunks),
    energy_rows = sum(queue$energy_request_rows),
    method_rows = sum(queue$assembled_method_rows), maximum_workers = 1L,
    unit_elapsed_cap_seconds = 600,
    configuration_elapsed_cap_seconds = 6 * 60 * 60,
    unit_rss_cap_bytes = 4 * 1024^3,
    configuration_storage_cap_bytes = 4 * 1024^3,
    point_selection = "sha256_sample_seed_cell_nested_v1",
    other_configuration_execution_authorized = FALSE,
    labels_opened = FALSE, rankings_computed = FALSE,
    outcomes_computed = FALSE, stringsAsFactors = FALSE)
}
