library(testthat)

testthat::test_that("MV8-H matrix/QC review prefreeze is closed pending owner approval", {
  root <- testthat::test_path("..", "..")
  evidence <- file.path(root, "docs", "audits",
    "mv08h-matrix-qc-review-prefreeze-v1")
  decision <- utils::read.csv(file.path(evidence,
    "mv08h-matrix-qc-review-decision.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  input <- utils::read.csv(file.path(evidence,
    "mv08h-matrix-qc-review-input-binding.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  firewall <- utils::read.csv(file.path(evidence,
    "mv08h-matrix-qc-review-firewall.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  validation <- utils::read.csv(file.path(evidence,
    "mv08h-matrix-qc-review-validation.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  manifest <- utils::read.csv(file.path(evidence,
    "mv08h-matrix-qc-review-artifact-manifest.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  report <- paste(readLines(file.path(evidence,
    "MV08H_MATRIX_QC_REVIEW_PREFREEZE_2026-08-18.md"), warn = FALSE),
    collapse = "\n")

  testthat::expect_equal(nrow(decision), 1L)
  testthat::expect_true(decision$prefreeze_completed)
  testthat::expect_false(decision$matrix_content_review_authorized)
  testthat::expect_false(decision$qc_content_review_authorized)
  testthat::expect_false(decision$labels_outcomes_authorized)
  testthat::expect_false(decision$landscapes_authorized)
  testthat::expect_false(decision$remaining_units_authorized)
  testthat::expect_false(decision$deletion_authorized)
  testthat::expect_equal(nrow(input), 5L)
  testthat::expect_true(all(!input$opened_by_prefreeze))
  testthat::expect_false(input$permitted_after_owner_approval[
    input$input_name == "molecule_info.h5"])
  testthat::expect_false(input$permitted_after_owner_approval[
    input$input_name == "web_summary.html"])
  testthat::expect_equal(nrow(firewall), 5L)
  testthat::expect_true(all(firewall$stop_on_breach))
  testthat::expect_equal(nrow(validation), 8L)
  testthat::expect_true(all(validation$passed))
  testthat::expect_true(all(!manifest$contains_private_path))
  testthat::expect_equal(nrow(manifest), 7L)
  testthat::expect_match(report, "does not open matrix values", fixed = TRUE)
  testthat::expect_match(report, "explicit owner approval", fixed = TRUE)
})
