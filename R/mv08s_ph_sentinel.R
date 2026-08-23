# MV8-S same-axis external baseline and bounded PH-sentinel helpers.

.mv08s_sha_object <- function(value) {
  digest::digest(value, algo = "sha256", serialize = TRUE)
}

.mv08s_ph_payload <- function(record) {
  list(
    identity = record$identity,
    topology_result = record$topology_result,
    h0_mst_oracle = record$h0_mst_oracle,
    downstream_execution = record$downstream_execution
  )
}

mv08s_validate_residual_cache_v1 <- function(cache, binding = NULL) {
  required <- c(
    "contract_id", "identity", "panels", "representations",
    "payload_sha256", "cache_key", "downstream_execution"
  )
  if (!is.list(cache) || !all(required %in% names(cache)) ||
      !identical(cache$contract_id, "mv08o_residual_source_cache_v1") ||
      !is.list(cache$identity) || nrow(cache$panels) != 500L ||
      !identical(cache$identity$outcome_label_state, "closed") ||
      isTRUE(cache$identity$biological_outcomes_computed) ||
      any(unlist(cache$downstream_execution, use.names = FALSE) != 0)) {
    stop("MV8-S residual source cache violates the frozen source contract.",
         call. = FALSE)
  }
  expected <- .mv08s_sha_object(cache[c(
    "identity", "panels", "representations", "downstream_execution"
  )])
  if (!identical(cache$payload_sha256, expected) ||
      !identical(cache$cache_key,
                 paste0("mv08o_residual_source_cache_v1:", expected))) {
    stop("MV8-S residual source cache payload is stale.", call. = FALSE)
  }
  if (!is.null(binding)) {
    if (!is.data.frame(binding) || nrow(binding) != 1L ||
        cache$identity$unit_id != binding$unit_id ||
        cache$identity$dataset_scope != binding$dataset_scope ||
        cache$cache_key != binding$source_cache_key) {
      stop("MV8-S residual source cache identity differs from its binding.",
           call. = FALSE)
    }
  }
  invisible(TRUE)
}

mv08s_residual_gene_view_v1 <- function(cache, row, common_panel) {
  mv08s_validate_residual_cache_v1(cache)
  if (!is.data.frame(row) || nrow(row) != 1L ||
      row$view_kind != "gene_topology_v1" ||
      row$execution_role != "source_produced_gene_ph" ||
      !(row$panel_id %in% c("exact500", "common475")) ||
      !(row$representation_id %in% c(
        "sct_data_all_qc_fit_selected384",
        "sct_pearson_residual_all_qc_fit_selected384"
      ))) {
    stop("MV8-S residual view request is outside the frozen queue.",
         call. = FALSE)
  }
  seed_key <- as.character(row$seed)
  representations <- cache$representations[[seed_key]]
  if (is.null(representations) ||
      is.null(representations[[row$representation_id]])) {
    stop("MV8-S requested representation is absent from the source cache.",
         call. = FALSE)
  }
  value <- representations[[row$representation_id]]
  if (row$panel_id == "common475") {
    if (!is.data.frame(common_panel) || nrow(common_panel) != 475L ||
        !identical(as.integer(common_panel$panel_order_475), 1:475)) {
      stop("MV8-S ordered common475 panel is stale.", call. = FALSE)
    }
    index <- match(common_panel$feature_id, cache$panels$feature_id)
    if (anyNA(index) || anyDuplicated(index)) {
      stop("MV8-S common475 panel does not map uniquely into the source cache.",
           call. = FALSE)
    }
    value <- value[index, , drop = FALSE]
    rownames(value) <- common_panel$feature_id
    profile <- "scientific_common475"
  } else {
    rownames(value) <- cache$panels$feature_id
    profile <- "scientific"
  }
  value <- as.matrix(value)
  storage.mode(value) <- "double"
  expected_selection <- cache$identity$selection_sha256[[seed_key]]
  if (is.null(expected_selection)) {
    expected_selection <- cache$identity$selection_sha256[[
      match(row$seed, as.integer(names(cache$identity$selection_sha256)))
    ]]
  }
  if (!identical(as.character(expected_selection),
                 as.character(row$selected_cell_sha256)) ||
      nrow(value) != row$panel_genes || ncol(value) != 384L ||
      !identical(.mv08s_sha_object(value), row$matrix_sha256)) {
    stop("MV8-S residual matrix or selected-cell binding drifted.",
         call. = FALSE)
  }
  source <- new_dual_view_source(
    value,
    sample_id = row$unit_id,
    cohort = row$dataset_scope,
    representation = row$representation_id,
    fit_scope_id = cache$cache_key,
    subsample_seed = as.integer(row$seed),
    standardization_id = row$matrix_sha256,
    contract_profile = profile
  )
  view <- construct_gene_topology_view(source)
  validate_topology_view(view)
  view
}

.mv08s_baseline_payload <- function(record) {
  list(
    identity = record$identity,
    panel = record$panel,
    views = record$views,
    downstream_execution = record$downstream_execution
  )
}

mv08s_new_baseline_record_v1 <- function(identity, panel, views) {
  if (!is.list(identity) || !is.data.frame(panel) || nrow(panel) != 475L ||
      !is.list(views) ||
      !identical(sort(names(views)),
                 sort(c("cell_topology_v1", "gene_topology_v1")))) {
    stop("MV8-S baseline record inputs are incomplete.", call. = FALSE)
  }
  lapply(views, validate_topology_view)
  record <- list(
    contract_id = "mv08s_same_axis_external_baseline_record_v1",
    identity = identity,
    panel = panel,
    views = views,
    payload_sha256 = NULL,
    cache_key = NULL,
    downstream_execution = list(
      ph_jobs = 0L, landscape_jobs = 0L, comparison_jobs = 0L,
      clustering_jobs = 0L, fusion_jobs = 0L, label_jobs = 0L,
      biological_outcome_jobs = 0L
    )
  )
  record$payload_sha256 <- .mv08s_sha_object(.mv08s_baseline_payload(record))
  record$cache_key <- paste0(
    "mv08s_same_axis_external_baseline_record_v1:", record$payload_sha256
  )
  class(record) <- c("scph_mv08s_baseline_record_v1", "list")
  mv08s_validate_baseline_record_v1(record)
  record
}

mv08s_validate_baseline_record_v1 <- function(record, binding = NULL) {
  required <- c(
    "contract_id", "identity", "panel", "views", "payload_sha256",
    "cache_key", "downstream_execution"
  )
  if (!is.list(record) || !all(required %in% names(record)) ||
      !identical(record$contract_id,
                 "mv08s_same_axis_external_baseline_record_v1") ||
      !identical(record$identity$outcome_label_state, "closed") ||
      isTRUE(record$identity$biological_outcomes_computed) ||
      nrow(record$panel) != 475L ||
      !identical(sort(names(record$views)),
                 sort(c("cell_topology_v1", "gene_topology_v1"))) ||
      any(unlist(record$downstream_execution, use.names = FALSE) != 0)) {
    stop("MV8-S baseline record violates the frozen contract.",
         call. = FALSE)
  }
  lapply(record$views, validate_topology_view)
  if (any(vapply(record$views, `[[`, character(1L), "sample_id") !=
          record$identity$unit_id)) {
    stop("MV8-S baseline typed-view identity is stale.", call. = FALSE)
  }
  expected <- .mv08s_sha_object(.mv08s_baseline_payload(record))
  if (!identical(record$payload_sha256, expected) ||
      !identical(record$cache_key, paste0(
        "mv08s_same_axis_external_baseline_record_v1:", expected
      ))) {
    stop("MV8-S baseline payload is stale.", call. = FALSE)
  }
  if (!is.null(binding) &&
      (record$identity$unit_id != binding$unit_id ||
       record$identity$selected_cell_sha256 != binding$selected_cell_sha256 ||
       record$identity$panel_sha256 != binding$panel_sha256 ||
       record$identity$filtered_h5_sha256 != binding$filtered_h5_sha256 ||
       record$identity$raw_h5_sha256 != binding$raw_h5_sha256)) {
    stop("MV8-S baseline input identity differs from its binding.",
         call. = FALSE)
  }
  invisible(TRUE)
}

mv08s_new_ph_record_v1 <- function(row, source_cache_key, view, result) {
  validate_topology_view(view)
  oracle <- mv07g_validate_ph_against_view_v1(result, view)
  identity <- list(
    contract_id = "mv08s_ph_identity_v1",
    job_id = row$job_id,
    unit_id = row$unit_id,
    dataset_scope = row$dataset_scope,
    seed = as.integer(row$seed),
    representation_id = row$representation_id,
    panel_id = row$panel_id,
    panel_sha256 = row$panel_sha256,
    selected_cell_sha256 = row$selected_cell_sha256,
    execution_role = row$execution_role,
    source_cache_key = source_cache_key,
    view_id = view$view_id,
    view_cache_key = view$cache_key,
    point_axis_role = view$point_axis_role,
    coordinate_axis_role = view$coordinate_axis_role,
    point_metric = view$point_metric,
    point_count = length(view$point_ids),
    max_dim = 1L,
    threshold = -1,
    field = 2L,
    outcome_label_state = "closed",
    biological_outcomes_computed = FALSE
  )
  record <- list(
    contract_id = "mv08s_ph_record_v1",
    identity = identity,
    topology_result = result,
    h0_mst_oracle = oracle,
    payload_sha256 = NULL,
    cache_key = NULL,
    downstream_execution = list(
      landscape_jobs = 0L, distance_jobs = 0L, comparison_jobs = 0L,
      clustering_jobs = 0L, fusion_jobs = 0L, label_jobs = 0L,
      biological_outcome_jobs = 0L
    )
  )
  record$payload_sha256 <- .mv08s_sha_object(.mv08s_ph_payload(record))
  record$cache_key <- paste0("mv08s_ph_record_v1:", record$payload_sha256)
  class(record) <- c("scph_mv08s_ph_record_v1", "list")
  mv08s_validate_ph_record_v1(record, row, view)
  record
}

mv08s_validate_ph_record_v1 <- function(record, row = NULL, view = NULL) {
  required <- c(
    "contract_id", "identity", "topology_result", "h0_mst_oracle",
    "payload_sha256", "cache_key", "downstream_execution"
  )
  if (!is.list(record) || !all(required %in% names(record)) ||
      !identical(record$contract_id, "mv08s_ph_record_v1") ||
      !(record$identity$view_id %in%
          c("cell_topology_v1", "gene_topology_v1")) ||
      !identical(record$identity$outcome_label_state, "closed") ||
      isTRUE(record$identity$biological_outcomes_computed) ||
      !isTRUE(record$h0_mst_oracle$passed) ||
      any(unlist(record$downstream_execution, use.names = FALSE) != 0)) {
    stop("MV8-S PH record violates the frozen contract.", call. = FALSE)
  }
  expected <- .mv08s_sha_object(.mv08s_ph_payload(record))
  if (!identical(record$payload_sha256, expected) ||
      !identical(record$cache_key, paste0("mv08s_ph_record_v1:", expected))) {
    stop("MV8-S PH record payload is stale.", call. = FALSE)
  }
  if (!is.null(row) &&
      (record$identity$job_id != row$job_id ||
       record$identity$unit_id != row$unit_id ||
       record$identity$selected_cell_sha256 != row$selected_cell_sha256 ||
       record$identity$view_id != row$view_kind)) {
    stop("MV8-S PH record differs from its frozen queue row.", call. = FALSE)
  }
  if (!is.null(view)) {
    if (record$identity$view_cache_key != view$cache_key) {
      stop("MV8-S PH record belongs to another typed view.", call. = FALSE)
    }
    mv07g_validate_ph_against_view_v1(record$topology_result, view)
  }
  invisible(TRUE)
}
