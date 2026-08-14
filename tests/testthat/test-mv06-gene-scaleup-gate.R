test_that("MV6-B SCT summaries separate incomplete and unresolved queries", {
  groups <- data.frame(
    fold_id = paste0("f", 1:75), held_out_study = paste0("s", 1:75),
    seed = rep(1:5, 15), training_views = 84L, query_views = 6L,
    exact_panel_variance_resolved_views = 84L,
    exact_panel_variance_unresolved_views = 5L,
    incomplete_panel_views = 1L, missing_feature_instances = 2L,
    maximum_missing_features = 2L, full_panel_group = FALSE,
    expression_payload_status = "reconstructable_from_accepted_d0_d1_sources",
    stringsAsFactors = FALSE
  )
  integrated <- data.frame(
    fold_id = groups$fold_id, held_out_study = groups$held_out_study,
    seed = groups$seed, training_views = 84L, query_views = 6L,
    full_panel_views = 89L, incomplete_panel_views = 1L,
    missing_feature_instances = 1L, maximum_missing_features = 1L,
    full_panel_group = FALSE, expression_undefined_views = 90L,
    expression_payload_status =
      "undefined_from_accepted_integrated_cell_coordinates",
    stringsAsFactors = FALSE
  )
  result <- mv06b_finalize_inventory_v1(groups, integrated)
  expect_equal(nrow(result$group), 75L)
  expect_equal(result$summary$view_instances, c(6750L, 6750L))
  expect_equal(result$summary$incomplete_panel_views, c(75L, 75L))
  expect_equal(result$summary$variance_unresolved_views, c(375L, 0L))
  expect_equal(result$summary$expression_undefined_views, c(0L, 6750L))
  expect_identical(result$decision$decision,
                   "stop_contract_revision_required")
  expect_equal(result$workload$dimension_pair_landscape_distances,
               c(70700L, 70700L))
  expect_true(all(result$decision$ph_jobs_executed == 0L))
  expect_false(result$decision$biological_outcomes_computed)
})

test_that("MV6-B rejects mismatched representation axes", {
  base <- data.frame(
    fold_id = paste0("f", 1:75), held_out_study = paste0("s", 1:75),
    seed = rep(1:5, 15), training_views = 84L, query_views = 6L,
    exact_panel_variance_resolved_views = 84L,
    exact_panel_variance_unresolved_views = 6L,
    incomplete_panel_views = 0L, missing_feature_instances = 0L,
    maximum_missing_features = 0L, full_panel_group = TRUE,
    stringsAsFactors = FALSE
  )
  other <- transform(
    base,
    full_panel_views = 90L,
    expression_undefined_views = 90L
  )
  other$fold_id[[1L]] <- "different"
  expect_error(mv06b_finalize_inventory_v1(base, other),
               "axes differ")
})
