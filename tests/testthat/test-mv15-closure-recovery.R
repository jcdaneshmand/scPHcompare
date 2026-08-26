test_that("MV15 closure recovery changes only the failed difference expression", {
  root <- testthat::test_path("..", "..")
  failure <- read.csv(file.path(
    root, "docs", "audits", "mv15-cell-distance-comparison-closure-failure-v1",
    "mv15-closure-failure.csv"
  ), stringsAsFactors = FALSE, check.names = FALSE)
  expect_equal(nrow(failure), 1L)
  expect_true(failure$production_complete)
  expect_equal(failure$production_comparisons, 36L)
  expect_equal(failure$published_closure_artifacts, 0L)
  expect_false(failure$scientific_values_changed)
  expect_false(failure$production_retry_authorized)

  closure <- readLines(file.path(
    root, "scripts", "build_mv15_cell_distance_comparison_closure.R"
  ), warn = FALSE)
  expression <- paste(closure, collapse = "\n")
  expect_match(expression, "neighbor_summary_delta <- max\\(abs\\(c\\(")
  expect_silent(parse(file.path(
    root, "scripts", "build_mv15_closure_recovery_prefreeze.R"
  )))
})
