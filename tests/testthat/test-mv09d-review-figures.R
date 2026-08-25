test_that("MV9-D review figures preserve replication and dimension strata", {
  root <- testthat::test_path("..", "..", "docs", "audits",
                              "mv09b-robustness-synthesis-v1")
  result <- mv09d_prepare_review_figure_data_v1(root)
  expect_equal(nrow(result$metric_contract), 4L)
  expect_identical(result$metric_contract$metric,
                   c("pearson", "spearman", "relative_stress",
                     "mean_neighbor_jaccard"))
  expect_equal(nrow(result$internal), 120L)
  expect_equal(length(unique(result$internal$seed)), 5L)
  expect_equal(nrow(result$internal_summary), 24L)
  expect_true(all(result$internal_summary$seeds == 5L))
  expect_identical(result$internal_summary$value,
                   result$internal_summary$median)
  expect_equal(nrow(result$external), 40L)
  expect_true(all(result$external$replication_units == 1L))
  expect_equal(nrow(result$delta), 80L)
  expect_equal(sum(result$delta$dataset_scope == "internal124"), 60L)
  expect_equal(sum(result$delta$dataset_scope == "external8"), 20L)
  expect_true(all(is.finite(result$delta$h1_minus_h0)))
})
