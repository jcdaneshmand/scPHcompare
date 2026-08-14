mv05d3_manifest_fixture <- function() {
  group_order <- rep(seq_len(75L), each = 90L)
  view_order <- rep(seq_len(90L), times = 75L)
  fold_number <- (group_order - 1L) %/% 5L + 1L
  seed <- 20260805L + (group_order - 1L) %% 5L
  group_id <- paste0("group_", sprintf("%02d", group_order))
  role <- ifelse(view_order <= 6L, "held_out", "training")
  data.frame(
    contract_id = "mv05d3_cell_ph_full_manifest_v1",
    job_id = paste(group_id, sprintf("view_%03d", view_order), sep = "__"),
    group_id = group_id, group_order = group_order, view_order = view_order,
    fold_id = paste0("fold_", sprintf("%02d", fold_number)),
    fit_scope_id = paste0("fit_", sprintf("%02d", fold_number)),
    held_out_study = paste0("study_", sprintf("%02d", fold_number)),
    seed = seed, sample_id = paste0("sample_", sprintf("%03d", view_order)),
    execution_role = role, missing_feature_count = 0L,
    mapping_stratum = "no_missing_training_features",
    representation = "sct_whole", view_id = "cell_topology_v1",
    point_axis_role = "cells",
    coordinate_axis_role = "shared_training_fitted_principal_components",
    point_count = 384L, coordinate_count = 30L,
    point_metric = "euclidean_frozen_shared_coordinates_v1",
    max_dim = 1L, threshold = -1, field = 2L,
    fold_cache_key = paste0("fold_key_", group_id),
    fold_cache_file = paste0(group_id, ".rds"),
    fold_cache_sha256 = paste(rep("a", 64L), collapse = ""),
    view_cache_key = paste0("view_key_", group_id, "_", view_order),
    view_payload_sha256 = paste(rep("b", 64L), collapse = ""),
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
}

test_that("MV5-D3 full manifest requires the complete label-closed axis", {
  manifest <- mv05d3_manifest_fixture()
  expect_invisible(mv05d3_validate_full_manifest_v1(manifest))
  expect_error(
    mv05d3_validate_full_manifest_v1(manifest[-1L, ]), "label-closed"
  )
  leaking <- manifest
  leaking$tissue <- "hidden"
  expect_error(mv05d3_validate_full_manifest_v1(leaking), "label-closed")
})

test_that("MV5-D3 view and group resource gates enforce the stop boundary", {
  views <- data.frame(
    job_id = paste0("j", 1:2), group_id = "g", fold_id = "f", seed = 1L,
    sample_id = paste0("s", 1:2), execution_role = "training",
    disposition = "built_atomic", operation_seconds = 1,
    h0_intervals = 384L, h1_intervals = 5L,
    h0_mst_oracle_passed = TRUE, result_size_bytes = 1000,
    result_file_sha256 = paste(rep("a", 64L), collapse = ""),
    record_cache_key = "key", outcome_label_state = "closed",
    biological_outcomes_computed = FALSE, landscape_jobs_executed = 0L,
    distance_jobs_executed = 0L, clustering_jobs_executed = 0L,
    integration_jobs_executed = 0L, gene_view_jobs_executed = 0L,
    stringsAsFactors = FALSE
  )
  expect_invisible(mv05d3_validate_view_metrics_v1(views, 2L, 10000))
  groups <- data.frame(
    group_id = "g", group_order = 1L, fold_id = "f", seed = 1L,
    disposition = "completed", completed_views = 90L,
    elapsed_seconds = 10, peak_process_tree_rss_bytes = 100,
    private_result_bytes = 1000, outcome_label_state = "closed",
    biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
  )
  expect_invisible(mv05d3_validate_group_metrics_v1(
    groups, 1L, 20, 1000, 100
  ))
  previous <- data.frame(
    scenario = "resource_safe_sct_cell_primary",
    normalization_worker_hours = 2.5,
    measured_cell_coordinate_worker_hours = 2.4,
    landscape_worker_hours = 3.6,
    known_components_lower_bound_worker_hours = 8.5,
    planning_cap_with_10_percent_reserve_hours = 21.6,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
  projected <- mv05d3_measured_primary_projection_v1(previous, groups)
  expect_equal(projected$measured_cell_ph_worker_hours, 10 / 3600)
  expect_true(projected$cap_passes)
  expect_false(projected$landscapes_complete)
  views$landscape_jobs_executed[[1L]] <- 1L
  expect_error(mv05d3_validate_view_metrics_v1(views, 2L), "scope")
})
