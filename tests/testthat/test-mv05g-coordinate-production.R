mv05g_d1_resource_fixture <- function() {
  studies <- sprintf("study_%02d", 1:15)
  result <- expand.grid(
    held_out_study = studies, seed = 20260805:20260809,
    stringsAsFactors = FALSE
  )
  held_out <- rep(c(1L, 2L, 4L, 5L, 6L), length.out = 15L)
  result$held_out_samples <- held_out[match(result$held_out_study, studies)]
  result$training_samples <- 90L - result$held_out_samples
  result$fold_id <- paste0("large_loso_v1:", result$held_out_study)
  result$fit_scope_id <- paste0(result$fold_id, ":training")
  result$panel_genes <- 500L
  result$cells_per_view <- 384L
  result$pca_components <- 30L
  result$private_cache_file <- paste0(
    result$held_out_study, "__", result$seed, ".rds"
  )
  result$private_cache_sha256 <- paste(rep("a", 64L), collapse = "")
  result$fold_cache_key <- paste0(
    "mv05d1_sct_cell_fold_v1:", paste(rep("b", 64L), collapse = "")
  )
  result$outcome_label_state <- "closed"
  result$biological_outcomes_computed <- FALSE
  result
}

test_that("MV5-G freezes the exact 75-group coordinate-only axis", {
  manifest <- mv05g_build_full_manifest_v1(mv05g_d1_resource_fixture())
  expect_invisible(mv05g_validate_full_manifest_v1(manifest))
  expect_equal(nrow(manifest), 75L)
  expect_equal(length(unique(manifest$held_out_study)), 15L)
  expect_equal(as.integer(table(manifest$held_out_study)), rep(5L, 15L))
  expect_true(all(manifest$maximum_heavy_workers == 1L))
  expect_true(all(manifest$stage_cap_seconds == 43200))
  leaking <- manifest
  leaking$tissue <- "hidden"
  expect_error(mv05g_validate_full_manifest_v1(leaking), "label-closed")
})

test_that("MV5-G production metrics require every group and zero downstream work", {
  manifest <- mv05g_build_full_manifest_v1(mv05g_d1_resource_fixture())
  metrics <- data.frame(
    group_id = manifest$group_id, held_out_study = manifest$held_out_study,
    seed = manifest$seed, held_out_samples = manifest$held_out_samples,
    disposition = "completed", exit_status = 0L,
    completed_query_mappings = manifest$held_out_samples,
    completed_coordinate_views = 90L, elapsed_seconds = 100,
    peak_process_tree_rss_bytes = 1000, private_result_bytes = 1000,
    reference_immutable = TRUE, label_transfer_jobs_executed = 0L,
    ph_jobs_executed = 0L, landscape_jobs_executed = 0L,
    distance_jobs_executed = 0L, retrieval_jobs_executed = 0L,
    clustering_jobs_executed = 0L, gene_view_jobs_executed = 0L,
    fusion_jobs_executed = 0L, new_data_jobs_executed = 0L,
    biological_outcomes_computed = FALSE, outcome_label_state = "closed",
    stringsAsFactors = FALSE
  )
  expect_invisible(mv05g_validate_group_metrics_v1(metrics))
  expect_error(mv05g_validate_group_metrics_v1(metrics[-1L, ]), "completion")
  metrics$ph_jobs_executed[[1L]] <- 1L
  expect_error(mv05g_validate_group_metrics_v1(metrics), "scope")
})
