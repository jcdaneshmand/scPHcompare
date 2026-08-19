library(testthat)

testthat::test_that("MV8-H approved matrix/QC review passes with firewalls", {
  root <- testthat::test_path("..", "..")
  evidence <- file.path(root, "docs", "audits", "mv08h-matrix-qc-review-v1")
  summary <- utils::read.csv(file.path(evidence,
    "mv08h-matrix-qc-review-summary.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  structure <- utils::read.csv(file.path(evidence,
    "mv08h-matrix-qc-review-structure.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  validation <- utils::read.csv(file.path(evidence,
    "mv08h-matrix-qc-review-validation.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  decision <- utils::read.csv(file.path(evidence,
    "mv08h-matrix-qc-review-decision.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  resource <- utils::read.csv(file.path(evidence,
    "mv08h-matrix-qc-review-resource.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  manifest <- utils::read.csv(file.path(evidence,
    "mv08h-matrix-qc-review-artifact-manifest.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  report <- paste(readLines(file.path(evidence,
    "MV08H_MATRIX_QC_REVIEW_2026-08-18.md"), warn = FALSE), collapse = "\n")

  testthat::expect_equal(nrow(summary), 1L)
  testthat::expect_equal(as.numeric(summary$filtered_cells), 5037)
  testthat::expect_equal(as.numeric(summary$post_qc_cells), 4614)
  testthat::expect_true(summary$depth_384_pass)
  testthat::expect_true(summary$expression_values_opened)
  testthat::expect_false(summary$barcode_identifiers_opened)
  testthat::expect_false(summary$labels_outcomes_opened)
  testthat::expect_false(summary$landscapes_computed)
  testthat::expect_equal(nrow(structure), 2L)
  testthat::expect_true(all(structure$matrix_values_opened))
  testthat::expect_true(all(!structure$barcode_identifiers_opened))
  testthat::expect_equal(nrow(validation), 8L)
  testthat::expect_true(all(validation$passed))
  testthat::expect_identical(decision$decision, "matrix_qc_review_pass")
  testthat::expect_false(decision$barcode_identifiers_opened)
  testthat::expect_false(decision$labels_outcomes_opened)
  testthat::expect_false(decision$pca_ph_landscapes_computed)
  testthat::expect_false(decision$remaining_units_authorized)
  testthat::expect_false(decision$deletion_authorized)
  testthat::expect_true(resource$rss_cap_passed)
  testthat::expect_equal(nrow(manifest), 6L)
  testthat::expect_true(all(!manifest$contains_private_path))
  testthat::expect_match(report, "No topology or biological interpretation", fixed = TRUE)
  testthat::expect_match(report, "Next gate", fixed = TRUE)
})
