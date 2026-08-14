# Internal MV5-H label-closed integrated cell-PH production contracts.

.mv05h_sha256 <- function(value) {
  digest::digest(value, algo = "sha256", serialize = TRUE)
}

.mv05h_hash <- function(value, name) {
  if (length(value) != 1L || is.na(value) ||
      !grepl("^[0-9a-f]{64}$", value)) {
    stop(name, " must be one lowercase SHA-256.", call. = FALSE)
  }
  value
}

.mv05h_string <- function(value, name) {
  if (length(value) != 1L || is.na(value) || !nzchar(value)) {
    stop(name, " must be one non-empty string.", call. = FALSE)
  }
  as.character(value)
}

mv05h_new_integrated_cell_view_v1 <- function(
    coordinates, sample_id, fold_id, fit_scope_id, seed,
    source_group_cache_key, coordinate_set_sha256) {
  coordinates <- .validate_named_numeric_matrix(coordinates, "coordinates")
  if (!identical(dim(coordinates), c(384L, 30L)) ||
      !identical(colnames(coordinates), paste0("PC", seq_len(30L)))) {
    stop("MV5-H coordinates must be named 384-by-30 PC matrices.",
         call. = FALSE)
  }
  sample_id <- .mv05h_string(sample_id, "sample_id")
  fold_id <- .mv05h_string(fold_id, "fold_id")
  fit_scope_id <- .mv05h_string(fit_scope_id, "fit_scope_id")
  source_group_cache_key <- .mv05h_string(
    source_group_cache_key, "source_group_cache_key"
  )
  coordinate_set_sha256 <- .mv05h_hash(
    coordinate_set_sha256, "coordinate_set_sha256"
  )
  seed <- as.integer(seed)
  if (length(seed) != 1L || is.na(seed)) {
    stop("seed must be one integer.", call. = FALSE)
  }
  source_identity <- list(
    contract_id = "mv05h_integrated_coordinate_source_v1",
    sample_id = sample_id, fold_id = fold_id, fit_scope_id = fit_scope_id,
    seed = seed, source_group_cache_key = source_group_cache_key,
    coordinate_set_sha256 = coordinate_set_sha256,
    coordinate_sha256 = .scientific_digest(coordinates)
  )
  source <- list(
    contract = list(profile = "scientific", scientific_eligible = TRUE),
    cache_key = paste0(
      "mv05h_integrated_coordinate_source_v1:",
      .scientific_digest(source_identity)
    ),
    sample_id = sample_id, cohort = fold_id,
    representation = "inductive_integrated",
    fit_scope_id = fit_scope_id, subsample_seed = seed
  )
  .new_topology_view(
    view_id = "cell_topology_v1", source = source,
    point_metric = "euclidean_frozen_shared_coordinates_v1",
    payload = coordinates, point_ids = rownames(coordinates),
    coordinate_ids = colnames(coordinates),
    transformations = list(
      coordinate_contract_id = "mv05f_label_closed_sct_transfer_v1",
      coordinate_fit_cache_key = source_group_cache_key,
      coordinate_set_sha256 = coordinate_set_sha256,
      n_components = 30L
    ),
    payload_sha256 = .scientific_digest(coordinates),
    diagnostics = list(duplicated_point_rows = sum(duplicated(coordinates)))
  )
}

mv05h_new_identity_v1 <- function(
    job, view, manifest_sha256, implementation_sha256,
    runtime = mv05d2_ph_runtime_v1()) {
  hashes <- c(
    source_group_sha256 = job$source_group_sha256,
    source_payload_sha256 = job$source_payload_sha256,
    coordinate_set_sha256 = job$coordinate_set_sha256,
    view_payload_sha256 = view$payload_sha256,
    manifest_sha256 = manifest_sha256,
    implementation_sha256 = implementation_sha256
  )
  invisible(lapply(names(hashes), function(name) {
    .mv05h_hash(hashes[[name]], name)
  }))
  identity <- list(
    contract_id = "mv05h_integrated_cell_ph_identity_v1",
    job_id = .mv05h_string(job$job_id, "job_id"),
    group_id = .mv05h_string(job$group_id, "group_id"),
    fold_id = .mv05h_string(job$fold_id, "fold_id"),
    fit_scope_id = .mv05h_string(job$fit_scope_id, "fit_scope_id"),
    held_out_study = .mv05h_string(job$held_out_study, "held_out_study"),
    seed = as.integer(job$seed), sample_id = .mv05h_string(
      job$sample_id, "sample_id"
    ),
    execution_role = match.arg(job$execution_role, c("held_out", "training")),
    missing_feature_count = as.integer(job$missing_feature_count),
    representation = "inductive_integrated", view_id = "cell_topology_v1",
    point_axis_role = "cells",
    coordinate_axis_role = "reference_fitted_inductive_integrated_coordinates",
    point_count = 384L, coordinate_count = 30L,
    point_metric = "euclidean_frozen_shared_coordinates_v1",
    max_dim = 1L, threshold = -1, field = 2L,
    source_group_cache_key = .mv05h_string(
      job$source_group_cache_key, "source_group_cache_key"
    ),
    source_group_sha256 = unname(hashes[["source_group_sha256"]]),
    source_payload_sha256 = unname(hashes[["source_payload_sha256"]]),
    coordinate_set_sha256 = unname(hashes[["coordinate_set_sha256"]]),
    view_cache_key = view$cache_key,
    view_payload_sha256 = unname(hashes[["view_payload_sha256"]]),
    manifest_sha256 = unname(hashes[["manifest_sha256"]]),
    implementation_sha256 = unname(hashes[["implementation_sha256"]]),
    runtime = runtime, outcome_label_state = "closed",
    biological_outcomes_computed = FALSE
  )
  identity$cache_key <- paste0(
    "mv05h_integrated_cell_ph_v1:", .mv05h_sha256(identity)
  )
  mv05h_validate_identity_v1(identity)
  identity
}

mv05h_validate_identity_v1 <- function(identity) {
  required <- c(
    "contract_id", "job_id", "group_id", "fold_id", "fit_scope_id",
    "held_out_study", "seed", "sample_id", "execution_role",
    "missing_feature_count", "representation", "view_id",
    "point_axis_role", "coordinate_axis_role", "point_count",
    "coordinate_count", "point_metric", "max_dim", "threshold", "field",
    "source_group_cache_key", "source_group_sha256",
    "source_payload_sha256", "coordinate_set_sha256", "view_cache_key",
    "view_payload_sha256", "manifest_sha256", "implementation_sha256",
    "runtime", "outcome_label_state", "biological_outcomes_computed",
    "cache_key"
  )
  if (!is.list(identity) || !all(required %in% names(identity))) {
    stop("MV5-H identity is incomplete.", call. = FALSE)
  }
  supplied <- identity$cache_key
  identity$cache_key <- NULL
  expected <- paste0(
    "mv05h_integrated_cell_ph_v1:", .mv05h_sha256(identity)
  )
  hash_fields <- c(
    "source_group_sha256", "source_payload_sha256", "coordinate_set_sha256",
    "view_payload_sha256", "manifest_sha256", "implementation_sha256"
  )
  if (!identical(supplied, expected) ||
      !identical(identity$contract_id,
                 "mv05h_integrated_cell_ph_identity_v1") ||
      !identity$execution_role %in% c("held_out", "training") ||
      !identical(identity$representation, "inductive_integrated") ||
      !identical(identity$view_id, "cell_topology_v1") ||
      !identical(identity$point_axis_role, "cells") ||
      !identical(identity$coordinate_axis_role,
                 "reference_fitted_inductive_integrated_coordinates") ||
      identity$point_count != 384L || identity$coordinate_count != 30L ||
      identity$max_dim != 1L || identity$threshold != -1 ||
      identity$field != 2L || identity$outcome_label_state != "closed" ||
      isTRUE(identity$biological_outcomes_computed) ||
      any(!grepl("^[0-9a-f]{64}$", unlist(identity[hash_fields])))) {
    stop("MV5-H identity violates its frozen contract.", call. = FALSE)
  }
  invisible(identity)
}

mv05h_new_record_v1 <- function(identity, view, result) {
  mv05h_validate_identity_v1(identity)
  validate_topology_view(view)
  if (!identical(identity$view_cache_key, view$cache_key) ||
      !identical(identity$view_payload_sha256, view$payload_sha256) ||
      !identical(identity$sample_id, view$sample_id) ||
      !identical(identity$fit_scope_id, view$fit_scope_id) ||
      !identical(identity$seed, view$subsample_seed)) {
    stop("MV5-H identity does not match its typed integrated view.",
         call. = FALSE)
  }
  oracle <- mv05d2_validate_ph_against_view_v1(result, view)
  payload <- list(identity = identity, topology_result = result,
                  h0_mst_oracle = oracle)
  payload_sha256 <- .mv05h_sha256(payload)
  record <- list(
    contract_id = "mv05h_integrated_cell_ph_record_v1",
    identity = identity, topology_result = result, h0_mst_oracle = oracle,
    payload_sha256 = payload_sha256,
    cache_key = paste0("mv05h_integrated_cell_ph_record_v1:", payload_sha256),
    downstream_execution = list(
      landscape_jobs = 0L, distance_jobs = 0L, retrieval_jobs = 0L,
      clustering_jobs = 0L, gene_view_jobs = 0L, fusion_jobs = 0L,
      new_data_jobs = 0L, biological_outcome_jobs = 0L
    )
  )
  class(record) <- c("scph_mv05h_integrated_cell_ph_record_v1", "list")
  mv05h_validate_record_static_v1(record)
  record
}

mv05h_validate_record_static_v1 <- function(record) {
  required <- c(
    "contract_id", "identity", "topology_result", "h0_mst_oracle",
    "payload_sha256", "cache_key", "downstream_execution"
  )
  if (!is.list(record) || !all(required %in% names(record)) ||
      record$contract_id != "mv05h_integrated_cell_ph_record_v1") {
    stop("MV5-H record is incomplete.", call. = FALSE)
  }
  mv05h_validate_identity_v1(record$identity)
  result <- record$topology_result
  diagram <- result$diagram
  h0 <- diagram[diagram[, "dimension"] == 0, , drop = FALSE]
  h1 <- diagram[diagram[, "dimension"] == 1, , drop = FALSE]
  oracle <- record$h0_mst_oracle
  payload <- list(identity = record$identity, topology_result = result,
                  h0_mst_oracle = oracle)
  zero <- unlist(record$downstream_execution, use.names = FALSE)
  if (!inherits(result, "scph_topology_result_v1") || !is.matrix(diagram) ||
      !identical(colnames(diagram), c("dimension", "birth", "death")) ||
      !identical(result$provenance$result_contract_id,
                 "corrected_topology_result_v1") ||
      !identical(result$provenance$view_cache_key,
                 record$identity$view_cache_key) ||
      !identical(result$provenance$sample_id, record$identity$sample_id) ||
      result$provenance$point_count != 384L ||
      result$provenance$max_dim != 1L || result$provenance$threshold != -1 ||
      result$provenance$field != 2L ||
      result$provenance$invalid_interval_count != 0L ||
      result$provenance$zero_persistence_count != 0L ||
      result$provenance$essential_h0_count != 1L ||
      nrow(h0) != 384L || sum(is.finite(h0[, "death"])) != 383L ||
      sum(is.infinite(h0[, "death"])) != 1L ||
      any(!is.finite(h1[, "death"])) || any(h1[, "death"] <= h1[, "birth"]) ||
      !identical(result$provenance$diagram_sha256,
                 .scientific_digest(diagram)) ||
      !is.list(oracle) || !isTRUE(oracle$passed) ||
      oracle$finite_h0_intervals != 383L ||
      oracle$finite_h1_intervals != nrow(h1) ||
      oracle$maximum_absolute_error > oracle$tolerance ||
      any(zero != 0L) || !identical(record$payload_sha256,
                                    .mv05h_sha256(payload)) ||
      !identical(record$cache_key, paste0(
        "mv05h_integrated_cell_ph_record_v1:", record$payload_sha256
      ))) {
    stop("MV5-H record violates diagram, oracle, identity, or scope rules.",
         call. = FALSE)
  }
  invisible(record)
}

mv05h_validate_manifest_v1 <- function(manifest) {
  required <- c(
    "contract_id", "job_id", "group_id", "group_order", "view_order",
    "fold_id", "fit_scope_id", "held_out_study", "seed", "sample_id",
    "execution_role", "missing_feature_count", "mapping_stratum",
    "representation", "view_id", "point_axis_role", "coordinate_axis_role",
    "point_count", "coordinate_count", "point_metric", "max_dim",
    "threshold", "field", "source_group_cache_key", "source_group_file",
    "source_group_sha256", "source_payload_sha256", "coordinate_set_sha256",
    "view_cache_key", "view_payload_sha256", "outcome_label_state",
    "biological_outcomes_computed"
  )
  if (!is.data.frame(manifest) || !all(required %in% names(manifest)) ||
      nrow(manifest) != 6750L ||
      any(manifest$contract_id != "mv05h_integrated_cell_ph_manifest_v1") ||
      anyDuplicated(manifest$job_id) ||
      length(unique(manifest$group_id)) != 75L ||
      any(table(manifest$group_id) != 90L) ||
      length(unique(manifest$fold_id)) != 15L ||
      !identical(sort(unique(as.integer(manifest$seed))), 20260805:20260809) ||
      any(table(manifest$seed) != 1350L) ||
      sum(manifest$execution_role == "held_out") != 450L ||
      sum(manifest$execution_role == "training") != 6300L ||
      any(manifest$representation != "inductive_integrated") ||
      any(manifest$point_count != 384L) ||
      any(manifest$coordinate_count != 30L) ||
      any(manifest$max_dim != 1L) || any(manifest$threshold != -1) ||
      any(manifest$field != 2L) ||
      any(manifest$outcome_label_state != "closed") ||
      any(as.logical(manifest$biological_outcomes_computed)) ||
      any(c("tissue", "approach") %in% names(manifest))) {
    stop("MV5-H manifest violates its frozen label-closed contract.",
         call. = FALSE)
  }
  groups <- unique(manifest[c(
    "group_id", "group_order", "fold_id", "seed", "source_group_cache_key",
    "source_group_file", "source_group_sha256"
  )])
  if (nrow(groups) != 75L || anyDuplicated(groups$group_order) ||
      !identical(sort(as.integer(groups$group_order)), 1:75) ||
      any(vapply(split(manifest$view_order, manifest$group_id), function(x) {
        !identical(sort(as.integer(x)), 1:90)
      }, logical(1L)))) {
    stop("MV5-H group/view order is incomplete.", call. = FALSE)
  }
  invisible(manifest)
}

mv05h_validate_view_metrics_v1 <- function(metrics, expected_jobs = 6750L,
                                             storage_cap_bytes = 1e9) {
  required <- c(
    "job_id", "group_id", "disposition", "operation_seconds",
    "h0_intervals", "h1_intervals", "h0_mst_oracle_passed",
    "result_size_bytes", "result_file_sha256", "record_cache_key",
    "landscape_jobs_executed", "distance_jobs_executed",
    "retrieval_jobs_executed", "clustering_jobs_executed",
    "gene_view_jobs_executed", "fusion_jobs_executed",
    "new_data_jobs_executed", "biological_outcomes_computed",
    "outcome_label_state"
  )
  zero <- c(
    "landscape_jobs_executed", "distance_jobs_executed",
    "retrieval_jobs_executed", "clustering_jobs_executed",
    "gene_view_jobs_executed", "fusion_jobs_executed",
    "new_data_jobs_executed"
  )
  if (!is.data.frame(metrics) || !all(required %in% names(metrics)) ||
      nrow(metrics) != as.integer(expected_jobs) ||
      anyDuplicated(metrics$job_id) ||
      any(!metrics$disposition %in% c("built_atomic", "reuse_validated")) ||
      any(!is.finite(metrics$operation_seconds)) ||
      any(metrics$operation_seconds < 0) ||
      any(metrics$h0_intervals != 384L) || any(metrics$h1_intervals < 0L) ||
      any(!as.logical(metrics$h0_mst_oracle_passed)) ||
      sum(metrics$result_size_bytes) > storage_cap_bytes ||
      any(as.matrix(metrics[zero]) != 0) ||
      any(metrics$outcome_label_state != "closed") ||
      any(as.logical(metrics$biological_outcomes_computed))) {
    stop("MV5-H view metrics violate completion, scope, or storage gates.",
         call. = FALSE)
  }
  invisible(metrics)
}

mv05h_validate_group_metrics_v1 <- function(
    metrics, expected_groups = 75L, elapsed_cap_seconds = 900,
    rss_cap_bytes = 4 * 1024^3, stage_cap_seconds = 12 * 3600) {
  required <- c(
    "group_id", "group_order", "fold_id", "seed", "disposition",
    "completed_views", "elapsed_seconds", "peak_process_tree_rss_bytes",
    "private_result_bytes", "outcome_label_state",
    "biological_outcomes_computed"
  )
  if (!is.data.frame(metrics) || !all(required %in% names(metrics)) ||
      nrow(metrics) != as.integer(expected_groups) ||
      anyDuplicated(metrics$group_id) ||
      any(metrics$disposition != "completed") ||
      any(metrics$completed_views != 90L) ||
      any(metrics$elapsed_seconds > elapsed_cap_seconds) ||
      sum(metrics$elapsed_seconds) > stage_cap_seconds ||
      any(metrics$peak_process_tree_rss_bytes > rss_cap_bytes) ||
      any(metrics$outcome_label_state != "closed") ||
      any(as.logical(metrics$biological_outcomes_computed))) {
    stop("MV5-H group metrics violate completion or resource gates.",
         call. = FALSE)
  }
  invisible(metrics)
}
