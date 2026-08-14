mv05h_hash_fixture <- function(letter = "a") {
  paste(rep(letter, 64L), collapse = "")
}

mv05h_manifest_fixture <- function() {
  groups <- expand.grid(
    held_out_study = sprintf("study_%02d", 1:15),
    seed = 20260805:20260809, stringsAsFactors = FALSE
  )
  groups <- groups[order(groups$held_out_study, groups$seed), , drop = FALSE]
  groups$group_order <- seq_len(nrow(groups))
  rows <- lapply(seq_len(nrow(groups)), function(index) {
    group <- groups[index, ]
    view_order <- 1:90
    data.frame(
      contract_id = "mv05h_integrated_cell_ph_manifest_v1",
      job_id = paste0("job_", index, "_", view_order),
      group_id = paste0("group_", index), group_order = index,
      view_order = view_order,
      fold_id = paste0("fold:", group$held_out_study),
      fit_scope_id = paste0("fold:", group$held_out_study, ":training"),
      held_out_study = group$held_out_study, seed = group$seed,
      sample_id = paste0("sample_", sprintf("%03d", view_order)),
      execution_role = ifelse(view_order <= 6L, "held_out", "training"),
      missing_feature_count = 0L,
      mapping_stratum = ifelse(
        view_order <= 6L, "held_out_full_panel", "training_reference"
      ),
      representation = "inductive_integrated", view_id = "cell_topology_v1",
      point_axis_role = "cells",
      coordinate_axis_role =
        "reference_fitted_inductive_integrated_coordinates",
      point_count = 384L, coordinate_count = 30L,
      point_metric = "euclidean_frozen_shared_coordinates_v1",
      max_dim = 1L, threshold = -1, field = 2L,
      source_group_cache_key = paste0("source:", index),
      source_group_file = paste0("source_", index, ".rds"),
      source_group_sha256 = mv05h_hash_fixture("a"),
      source_payload_sha256 = mv05h_hash_fixture("b"),
      coordinate_set_sha256 = mv05h_hash_fixture("c"),
      view_cache_key = paste0("view:", index, ":", view_order),
      view_payload_sha256 = mv05h_hash_fixture("d"),
      outcome_label_state = "closed", biological_outcomes_computed = FALSE,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

test_that("MV5-H constructs a typed cells-as-observations integrated view", {
  set.seed(20260809)
  coordinates <- matrix(rnorm(384L * 30L), nrow = 384L)
  rownames(coordinates) <- paste0("sample__cell_", seq_len(384L))
  colnames(coordinates) <- paste0("PC", seq_len(30L))
  view <- mv05h_new_integrated_cell_view_v1(
    coordinates, "sample", "fold:study", "fold:study:training", 20260805,
    "source:key", mv05h_hash_fixture("a")
  )
  expect_invisible(validate_topology_view(view))
  expect_identical(dim(view$payload), c(384L, 30L))
  expect_identical(view$point_axis_role, "cells")
  expect_identical(view$representation, "inductive_integrated")
  expect_error(
    mv05h_new_integrated_cell_view_v1(
      t(coordinates), "sample", "fold:study", "fold:study:training",
      20260805, "source:key", mv05h_hash_fixture("a")
    ),
    "384-by-30"
  )
})

test_that("MV5-H freezes complete label-closed manifest axes", {
  manifest <- mv05h_manifest_fixture()
  expect_invisible(mv05h_validate_manifest_v1(manifest))
  expect_equal(nrow(manifest), 6750L)
  expect_equal(sum(manifest$execution_role == "held_out"), 450L)
  expect_equal(as.integer(table(manifest$group_id)), rep(90L, 75L))
  leaking <- manifest
  leaking$tissue <- "hidden"
  expect_error(mv05h_validate_manifest_v1(leaking), "label-closed")
})

test_that("MV5-H view metrics reject downstream execution", {
  metrics <- data.frame(
    job_id = paste0("job_", 1:90), group_id = "group",
    disposition = "built_atomic", operation_seconds = 0.1,
    h0_intervals = 384L, h1_intervals = 10L,
    h0_mst_oracle_passed = TRUE, result_size_bytes = 1000,
    result_file_sha256 = mv05h_hash_fixture("a"),
    record_cache_key = paste0("record_", 1:90),
    landscape_jobs_executed = 0L, distance_jobs_executed = 0L,
    retrieval_jobs_executed = 0L, clustering_jobs_executed = 0L,
    gene_view_jobs_executed = 0L, fusion_jobs_executed = 0L,
    new_data_jobs_executed = 0L, biological_outcomes_computed = FALSE,
    outcome_label_state = "closed", stringsAsFactors = FALSE
  )
  expect_invisible(mv05h_validate_view_metrics_v1(metrics, expected_jobs = 90L))
  metrics$landscape_jobs_executed[[1L]] <- 1L
  expect_error(mv05h_validate_view_metrics_v1(metrics, 90L), "scope")
})
