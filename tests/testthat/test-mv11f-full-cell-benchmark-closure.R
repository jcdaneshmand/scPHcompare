test_that("MV11-F exactly closes the matched historical cell benchmark", {
  root <- testthat::test_path("..", "..")
  audit <- file.path(root, "docs", "audits",
                     "mv11f-cell-benchmark-closure-v1")
  read_audit <- function(name) read.csv(
    file.path(audit, name), stringsAsFactors = FALSE, check.names = FALSE
  )
  validation <- read_audit("mv11f-validation.csv")
  repeat_rows <- read_audit("mv11f-artifact-repeat.csv")
  summary <- read_audit("mv11f-closure-summary.csv")
  decision <- read_audit("mv11f-decision.csv")
  manifest <- read_audit("mv11f-artifact-manifest.csv")
  expect_equal(nrow(validation), 30L)
  expect_true(all(validation$passed))
  expect_equal(c(summary$matrices, summary$partition_fits,
                 summary$private_assignment_rows, summary$quality_rows,
                 summary$stability_rows, summary$primary_k_rows,
                 summary$agreement_rows),
               c(10, 450, 55800, 450, 90, 2, 900))
  expect_equal(c(summary$selected_H0_k, summary$selected_H1_k), c(2, 2))
  expect_true(summary$exact_repeat)
  expect_false(summary$labels_used)
  expect_false(summary$outcomes_used)
  expect_false(summary$cross_view_comparison_performed)
  expect_equal(summary$result_interpretation_state,
               "label_closed_cell_benchmark_only")
  expect_equal(nrow(repeat_rows), 5L)
  expect_true(all(repeat_rows$exact_repeat))
  expect_equal(repeat_rows$saved_sha256, repeat_rows$repeat_sha256)
  expect_true(decision$historical_cell_benchmark_closed)
  expect_true(decision$common_k_cross_view_prefreeze_eligible_next)
  expect_false(decision$cross_view_comparison_authorized_now)
  expect_false(decision$labels_authorized)
  expect_false(decision$outcomes_authorized)
  artifacts <- file.path(audit, manifest$artifact)
  expect_equal(as.numeric(file.info(artifacts)$size), manifest$bytes)
  expect_equal(unname(vapply(artifacts, function(path) digest::digest(
    file = path, algo = "sha256", serialize = FALSE
  ), character(1L))), manifest$sha256)
})
