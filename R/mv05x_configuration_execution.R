# MV5-X bounded PC20 configuration execution contracts.

.mv05x_hash_ok <- function(values, width = 64L) {
  is.character(values) && length(values) > 0L &&
    all(!is.na(values) & grepl(paste0("^[0-9a-f]{", width, "}$"), values))
}

#' Validate the prospectively bound MV5-X PC20 execution queue.
#'
#' @param queue A 600-row MV5-V queue rebound to the MV5-X runtime.
#' @return `invisible(queue)` on success.
#' @keywords internal
mv05x_validate_configuration_queue_v1 <- function(queue) {
  runtime <- c(
    "python_executable_sha256", "python_version", "persim_version",
    "numpy_version", "scipy_version", "configuration_execution_order"
  )
  authorized_id <- "cells384_pc20_euclidean_v1"
  if (!is.data.frame(queue) || !all(runtime %in% names(queue)) ||
      nrow(queue) != 600L ||
      any(queue$contract_id != "mv05x_pc20_execution_queue_v1") ||
      sum(as.logical(queue$execution_authorized)) != 150L ||
      any(as.logical(queue$execution_authorized) !=
            (queue$configuration_id == authorized_id)) ||
      any(as.logical(queue$execution_completed)) ||
      !identical(
        as.integer(queue$configuration_execution_order[
          as.logical(queue$execution_authorized)
        ]),
        1:150
      ) ||
      any(!is.na(queue$configuration_execution_order[
        !as.logical(queue$execution_authorized)
      ])) ||
      length(unique(queue$python_executable_sha256)) != 1L ||
      !.mv05x_hash_ok(unique(queue$python_executable_sha256)) ||
      length(unique(queue$implementation_sha256)) != 1L ||
      !.mv05x_hash_ok(unique(queue$implementation_sha256)) ||
      length(unique(queue$source_freeze_sha256)) != 1L ||
      !.mv05x_hash_ok(unique(queue$source_freeze_sha256)) ||
      length(unique(queue$prospective_head)) != 1L ||
      !.mv05x_hash_ok(unique(queue$prospective_head), 40L) ||
      any(queue$outcome_label_state != "closed") ||
      any(as.logical(queue$outcomes_computed)) ||
      any(c("tissue", "approach") %in% names(queue))) {
    stop("MV5-X requires exactly its 150 label-closed PC20 groups.",
         call. = FALSE)
  }
  frozen <- queue
  frozen$execution_authorized <- FALSE
  frozen$execution_completed <- FALSE
  mv05v_validate_group_queue_v1(frozen)
  invisible(queue)
}

#' Summarize a bound MV5-X PC20 execution queue.
#'
#' @param queue A valid MV5-X queue.
#' @return One label-closed configuration-scope record.
#' @keywords internal
mv05x_configuration_scope_v1 <- function(queue) {
  mv05x_validate_configuration_queue_v1(queue)
  authorized <- queue[as.logical(queue$execution_authorized), , drop = FALSE]
  data.frame(
    contract_id = "mv05x_pc20_configuration_scope_v1",
    configuration_id = unique(authorized$configuration_id),
    authorized_groups = nrow(authorized),
    representations = length(unique(authorized$representation)),
    folds = length(unique(authorized$fold_id)),
    seeds = length(unique(authorized$seed)),
    views = sum(authorized$view_count),
    biological_pairs = sum(authorized$biological_pairs),
    landscape_rows = sum(authorized$landscape_request_rows),
    landscape_subchunks = sum(authorized$landscape_subchunks),
    energy_rows = sum(authorized$energy_request_rows),
    method_rows = sum(authorized$assembled_method_rows),
    maximum_workers = 1L,
    unit_elapsed_cap_seconds = 600,
    program_elapsed_cap_seconds = 30 * 60 * 60,
    unit_rss_cap_bytes = 4 * 1024^3,
    program_storage_cap_bytes = 16 * 1024^3,
    labels_opened = FALSE,
    outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
}
