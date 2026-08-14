test_that("MV6-C selects a deterministic technical-filtered global core", {
  feature_ids <- c(
    "A-ENSG1", "B-ENSG2", "C-ENSG3", "MT-D-ENSG4",
    "A-ENSG5", "E-ENSG6"
  )
  seeds <- rep(20260805:20260809, each = 90L)
  variances <- matrix(1, nrow = length(feature_ids), ncol = length(seeds))
  variances[1L, ] <- 6
  variances[2L, ] <- 5
  variances[3L, ] <- 4
  variances[4L, ] <- 100
  variances[5L, ] <- 3
  variances[6L, ] <- 2
  variances[3L, 1L] <- 0
  result <- mv06c_build_global_core_panel_v1(
    feature_ids, variances, seeds, panel_size = 2L
  )
  expect_identical(result$decision, "go_bounded_matched_sct_profile")
  expect_identical(result$panel$feature_id, c("A-ENSG1", "B-ENSG2"))
  expect_identical(result$panel$gene, c("A", "B"))
  expect_true(all(result$panel$positive_cache_count == 450L))
  expect_equal(result$eligibility$technical_category_features, 1L)
  expect_equal(result$eligibility$nonpositive_any_cache_features, 1L)
  expect_equal(result$eligibility$duplicate_canonical_features_removed, 1L)
  expect_equal(nrow(result$seed_stability), 5L)
  expect_true(all(result$seed_stability$cache_count == 90L))
})

test_that("MV6-C stops rather than shrinking an insufficient panel", {
  feature_ids <- c("A", "B", "C")
  seeds <- rep(20260805:20260809, each = 90L)
  variances <- matrix(1, nrow = 3L, ncol = length(seeds))
  variances[3L, ] <- 0
  result <- mv06c_build_global_core_panel_v1(
    feature_ids, variances, seeds, panel_size = 3L
  )
  expect_identical(result$decision, "stop_global_core_insufficient")
  expect_equal(nrow(result$panel), 2L)
  expect_equal(result$eligibility$eligibility_margin, -1L)
})

test_that("MV6-C future workload remains bounded and unauthorized", {
  workload <- mv06c_future_workload_v1()
  expect_equal(workload$matched_cell_views, 6750L)
  expect_equal(workload$matched_gene_views, 6750L)
  expect_equal(workload$four_component_landscape_distances, 141400L)
  expect_equal(workload$five_weight_fusion_pair_rows, 176750L)
  expect_false(workload$execution_authorized)
})
