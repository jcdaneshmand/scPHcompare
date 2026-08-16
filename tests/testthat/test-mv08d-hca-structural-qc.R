testthat::test_that("MV8-D freezes exact legacy-comparable QC before H5 inspection", {
  root <- testthat::test_path("..", "..")
  spec <- paste(readLines(file.path(root, "docs", "specifications",
    "MV08D_HCA_STRUCTURAL_QC_PREFREEZE_V1.md"), warn = FALSE),
    collapse = "\n")
  qc <- utils::read.csv(file.path(root, "docs", "specifications",
    "mv08d-hca-qc-contract-v1.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  testthat::expect_identical(qc$rule_order, 1:10)
  testthat::expect_identical(qc$threshold,
    c("Gene Expression", "200", "3", "500", "9000", "20", "5", "3",
      "384", "500"))
  testthat::expect_false(any(qc$post_observation_change_allowed))
  testthat::expect_match(spec, "not an EmptyDrops-style statistical droplet caller",
    fixed = TRUE)
  testthat::expect_match(spec, "prohibit symbol-only rescue", fixed = TRUE)
  testthat::expect_match(spec, "No HCA value contributes to a reference fit",
    fixed = TRUE)
})

testthat::test_that("MV8-D freezes all five immutable 124-sample transforms", {
  root <- testthat::test_path("..", "..")
  transforms <- utils::read.csv(file.path(root, "docs", "specifications",
    "mv08d-reference-transform-freeze-v1.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  testthat::expect_identical(transforms$seed, 20260805:20260809)
  testthat::expect_true(all(transforms$fit_samples == 124L))
  testthat::expect_true(all(transforms$center_features == 500L))
  testthat::expect_true(all(transforms$scale_features == 500L))
  testthat::expect_true(all(transforms$rotation_rows == 500L))
  testthat::expect_true(all(transforms$rotation_columns == 30L))
  testthat::expect_length(unique(transforms$source_file_sha256), 5L)
  testthat::expect_length(unique(transforms$pca_model_cache_key), 5L)
  testthat::expect_identical(unique(transforms$panel_sha256),
    "48e881ee753893bfaecd7101fc16fbcb552dd30a2a791fedfc7204d7b322a32e")
})
