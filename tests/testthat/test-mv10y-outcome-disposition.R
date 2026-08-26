test_that("MV10-Y disposition reports contrasts without selection", {
  root <- testthat::test_path("..", "..")
  complete <- read.csv(file.path(
    root, "docs", "audits", "mv10t-outcome-review-synthesis-v1",
    "mv10t-complete-summary.csv"
  ), stringsAsFactors = FALSE, check.names = FALSE)
  endpoints <- read.csv(file.path(
    root, "docs", "audits", "mv10t-outcome-review-synthesis-v1",
    "mv10t-endpoint-coverage.csv"
  ), stringsAsFactors = FALSE, check.names = FALSE)
  result <- mv10y_build_outcome_disposition_v1(complete, endpoints)
  expect_equal(nrow(result$primary_envelope), 20L)
  expect_equal(nrow(result$method_envelope), 60L)
  expect_equal(nrow(result$context_contrast), 120L)
  expect_equal(nrow(result$approach_contrast), 120L)
  expect_equal(nrow(result$disposition), 1L)
  expect_true(all(!result$context_contrast$inference_performed))
  expect_true(all(!result$context_contrast$causal_interpretation_allowed))
  expect_true(all(!result$approach_contrast$inference_performed))
  expect_true(all(!result$approach_contrast$causal_interpretation_allowed))
  expect_false(result$disposition$magnitude_threshold_applied)
  expect_false(result$disposition$p_values_computed)
  expect_false(result$disposition$method_selection_executed)
  expect_false(result$disposition$representation_ranking_executed)
  expect_false(result$disposition$H0_H1_combined)
  expect_false(result$disposition$approach_causal_interpretation)
  expect_false(result$disposition$biological_claims)
  expect_false(result$disposition$manuscript_claims)
})
