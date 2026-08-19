library(testthat)

testthat::test_that("MV8-H metadata-only structural/QC admission is closed", {
  root <- testthat::test_path("..", "..")
  evidence <- file.path(root, "docs", "audits", "mv08h-structural-qc-admission-v1")
  decision <- utils::read.csv(file.path(evidence,
    "mv08h-structural-qc-decision.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  validation <- utils::read.csv(file.path(evidence,
    "mv08h-structural-qc-validation.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  h5 <- utils::read.csv(file.path(evidence,
    "mv08h-structural-qc-h5-metadata.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  metadata <- utils::read.csv(file.path(evidence,
    "mv08h-structural-qc-metadata.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  outputs <- utils::read.csv(file.path(evidence,
    "mv08h-structural-qc-output-structure.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  manifest <- utils::read.csv(file.path(evidence,
    "mv08h-structural-qc-artifact-manifest.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  report <- paste(readLines(file.path(evidence,
    "MV08H_STRUCTURAL_QC_ADMISSION_2026-08-18.md"), warn = FALSE),
    collapse = "\n")

  testthat::expect_equal(nrow(decision), 1L)
  testthat::expect_identical(decision$decision,
    "metadata_only_structural_qc_admit")
  testthat::expect_true(decision$exact_output_structure_passed)
  testthat::expect_true(decision$h5_metadata_passed)
  testthat::expect_true(decision$aggregate_metrics_passed)
  testthat::expect_true(decision$resource_provenance_passed)
  testthat::expect_false(decision$expression_values_opened)
  testthat::expect_false(decision$labels_outcomes_opened)
  testthat::expect_false(decision$landscapes_computed)
  testthat::expect_false(decision$remaining_units_authorized)

  testthat::expect_equal(nrow(validation), 8L)
  testthat::expect_true(all(validation$passed))
  testthat::expect_equal(nrow(h5), 2L)
  testthat::expect_true(all(h5$metadata_schema_pass))
  testthat::expect_true(all(!h5$expression_values_opened))
  testthat::expect_true(all(!h5$barcode_identifiers_opened))
  testthat::expect_equal(as.numeric(metadata$filtered_estimated_cells), 5037)
  testthat::expect_equal(as.numeric(metadata$filtered_median_genes_per_cell), 801)
  testthat::expect_equal(nrow(outputs), 7L)
  testthat::expect_true(all(outputs$exists))
  testthat::expect_true(all(!outputs$expression_values_opened))
  testthat::expect_true(all(!outputs$labels_outcomes_opened))
  testthat::expect_true(all(!manifest$contains_private_path))
  testthat::expect_equal(nrow(manifest), 6L)
  testthat::expect_match(report, "aggregate Cell Ranger QC metadata", fixed = TRUE)
  testthat::expect_match(report, "Next gate", fixed = TRUE)
})
