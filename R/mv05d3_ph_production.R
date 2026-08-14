# MV5-D3 full label-closed cell-PH production contracts.

mv05d3_validate_full_manifest_v1 <- function(manifest) {
  required <- c(
    "contract_id", "job_id", "group_id", "group_order", "view_order",
    "fold_id", "fit_scope_id", "held_out_study", "seed", "sample_id",
    "execution_role", "missing_feature_count", "mapping_stratum",
    "representation", "view_id", "point_axis_role",
    "coordinate_axis_role", "point_count", "coordinate_count",
    "point_metric", "max_dim", "threshold", "field", "fold_cache_key",
    "fold_cache_file", "fold_cache_sha256", "view_cache_key",
    "view_payload_sha256", "outcome_label_state",
    "biological_outcomes_computed"
  )
  if (!is.data.frame(manifest) || !all(required %in% names(manifest)) ||
      nrow(manifest) != 6750L ||
      any(manifest$contract_id != "mv05d3_cell_ph_full_manifest_v1") ||
      anyDuplicated(manifest$job_id) ||
      length(unique(manifest$group_id)) != 75L ||
      any(table(manifest$group_id) != 90L) ||
      length(unique(manifest$fold_id)) != 15L ||
      !identical(sort(unique(as.integer(manifest$seed))), 20260805:20260809) ||
      any(table(manifest$seed) != 1350L) ||
      any(manifest$point_count != 384L) ||
      any(manifest$coordinate_count != 30L) ||
      any(manifest$max_dim != 1L) || any(manifest$threshold != -1) ||
      any(manifest$field != 2L) ||
      sum(manifest$execution_role == "held_out") != 450L ||
      sum(manifest$execution_role == "training") != 6300L ||
      any(manifest$outcome_label_state != "closed") ||
      any(as.logical(manifest$biological_outcomes_computed)) ||
      any(c("tissue", "approach") %in% names(manifest))) {
    stop("MV5-D3 full manifest violates the frozen label-closed contract.",
         call. = FALSE)
  }
  groups <- unique(manifest[c(
    "group_id", "group_order", "fold_id", "seed", "fold_cache_key",
    "fold_cache_file", "fold_cache_sha256"
  )])
  if (nrow(groups) != 75L || anyDuplicated(groups$group_order) ||
      !identical(sort(groups$group_order), seq_len(75L)) ||
      any(vapply(split(manifest$view_order, manifest$group_id), function(x) {
        !identical(sort(as.integer(x)), seq_len(90L))
      }, logical(1L)))) {
    stop("MV5-D3 group or view order is not a complete deterministic axis.",
         call. = FALSE)
  }
  invisible(manifest)
}

mv05d3_validate_record_static_v1 <- function(record) {
  mv05d2_validate_ph_record_v1(record)
  identity <- record$identity
  result <- record$topology_result
  diagram <- result$diagram
  h0 <- diagram[diagram[, "dimension"] == 0, , drop = FALSE]
  h1 <- diagram[diagram[, "dimension"] == 1, , drop = FALSE]
  oracle <- record$h0_mst_oracle
  if (!inherits(result, "scph_topology_result_v1") ||
      !is.matrix(diagram) || !identical(colnames(diagram),
                                        c("dimension", "birth", "death")) ||
      !identical(result$provenance$result_contract_id,
                 "corrected_topology_result_v1") ||
      !identical(result$provenance$view_cache_key, identity$view_cache_key) ||
      !identical(result$provenance$sample_id, identity$sample_id) ||
      !identical(result$provenance$point_count, 384L) ||
      !identical(result$provenance$max_dim, 1L) ||
      !identical(result$provenance$threshold, -1) ||
      !identical(result$provenance$field, 2L) ||
      result$provenance$invalid_interval_count != 0L ||
      result$provenance$zero_persistence_count != 0L ||
      result$provenance$essential_h0_count != 1L ||
      nrow(h0) != 384L || sum(is.finite(h0[, "death"])) != 383L ||
      sum(is.infinite(h0[, "death"])) != 1L ||
      any(!is.finite(h1[, "death"])) ||
      any(h1[, "death"] <= h1[, "birth"]) ||
      !identical(result$provenance$diagram_sha256,
                 .scientific_digest(diagram)) ||
      !identical(result$cache_key, paste0(
        "corrected_topology_result_v1:",
        .scientific_digest(result$provenance)
      )) ||
      !is.list(oracle) || !isTRUE(oracle$passed) ||
      oracle$finite_h0_intervals != 383L ||
      oracle$finite_h1_intervals != nrow(h1) ||
      !is.finite(oracle$maximum_absolute_error) ||
      oracle$maximum_absolute_error > oracle$tolerance) {
    stop("MV5-D3 record violates diagram or stored-oracle invariants.",
         call. = FALSE)
  }
  invisible(record)
}

mv05d3_validate_view_metrics_v1 <- function(metrics, expected_jobs = 6750L,
                                             storage_cap_bytes = 1e9) {
  required <- c(
    "job_id", "group_id", "fold_id", "seed", "sample_id",
    "execution_role", "disposition", "operation_seconds", "h0_intervals",
    "h1_intervals", "h0_mst_oracle_passed", "result_size_bytes",
    "result_file_sha256", "record_cache_key", "outcome_label_state",
    "biological_outcomes_computed", "landscape_jobs_executed",
    "distance_jobs_executed", "clustering_jobs_executed",
    "integration_jobs_executed", "gene_view_jobs_executed"
  )
  zero_fields <- c(
    "landscape_jobs_executed", "distance_jobs_executed",
    "clustering_jobs_executed", "integration_jobs_executed",
    "gene_view_jobs_executed"
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
      any(metrics$outcome_label_state != "closed") ||
      any(as.logical(metrics$biological_outcomes_computed)) ||
      any(as.matrix(metrics[zero_fields]) != 0)) {
    stop("MV5-D3 view metrics violate completion, scope, or storage gates.",
         call. = FALSE)
  }
  invisible(metrics)
}

mv05d3_validate_group_metrics_v1 <- function(
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
    stop("MV5-D3 group metrics violate completion or resource gates.",
         call. = FALSE)
  }
  invisible(metrics)
}

mv05d3_group_summary_v1 <- function(view_metrics) {
  groups <- split(view_metrics, view_metrics$group_id)
  result <- do.call(rbind, lapply(groups, function(group) {
    data.frame(
      contract_id = "mv05d3_group_view_summary_v1",
      group_id = group$group_id[[1L]], fold_id = group$fold_id[[1L]],
      seed = group$seed[[1L]], views = nrow(group),
      held_out_views = sum(group$execution_role == "held_out"),
      training_views = sum(group$execution_role == "training"),
      h0_intervals = sum(group$h0_intervals),
      h1_intervals = sum(group$h1_intervals),
      operation_seconds_sum = sum(group$operation_seconds),
      private_result_bytes = sum(group$result_size_bytes),
      outcome_label_state = "closed", biological_outcomes_computed = FALSE,
      stringsAsFactors = FALSE
    )
  }))
  rownames(result) <- NULL
  result
}
