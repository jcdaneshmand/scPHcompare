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

testthat::test_that("MV8-D implementation is fail-closed and aggregate-only", {
  root <- testthat::test_path("..", "..")
  fetch <- paste(readLines(file.path(root, "scripts",
    "fetch_mv08d_hca_h5.py"), warn = FALSE), collapse = "\n")
  audit <- paste(readLines(file.path(root, "scripts",
    "build_mv08d_hca_structural_qc.py"), warn = FALSE), collapse = "\n")
  transforms <- paste(readLines(file.path(root, "scripts",
    "validate_mv08d_reference_transforms.R"), warn = FALSE), collapse = "\n")
  testthat::expect_match(fetch, "existing cache file is incompatible and was not replaced",
    fixed = TRUE)
  testthat::expect_match(fetch, "Transient redirected URLs are intentionally never logged",
    fixed = TRUE)
  testthat::expect_match(audit, "symbol-only", fixed = TRUE)
  testthat::expect_match(audit, "panel_annotation_reconciliation_or_cohort_block",
    fixed = TRUE)
  testthat::expect_false(grepl("selected_cell", audit, fixed = TRUE))
  testthat::expect_match(transforms, "reference_refit_performed = FALSE",
    fixed = TRUE)
  testthat::expect_match(transforms, "hca_expression_accessed = FALSE",
    fixed = TRUE)
})

testthat::test_that("MV8-D publishes the blocked exact-500 result transparently", {
  evidence <- testthat::test_path("..", "..", "docs", "audits",
    "mv08d-hca-structural-qc-evidence")
  decision <- utils::read.csv(file.path(evidence, "mv08d-decision.csv"),
    stringsAsFactors = FALSE, check.names = FALSE)
  summary <- utils::read.csv(file.path(evidence, "mv08d-qc-summary.csv"),
    stringsAsFactors = FALSE, check.names = FALSE)
  missing <- utils::read.csv(file.path(evidence,
    "mv08d-missing-panel-summary.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  recovery <- utils::read.csv(file.path(evidence,
    "mv08d-recovery-options.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  validation <- utils::read.csv(file.path(evidence,
    "mv08d-independent-validation.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  testthat::expect_identical(decision$disposition,
    "blocked_exact500_annotation_incompatibility")
  testthat::expect_equal(summary$post_qc_cells_min, 3403L)
  testthat::expect_equal(summary$post_qc_cells_max, 4707L)
  testthat::expect_equal(summary$exact_stable_id_intersection, 475L)
  testthat::expect_equal(summary$exact_stable_id_intersection_fraction, 0.95)
  testthat::expect_equal(nrow(missing), 25L)
  testthat::expect_true(all(missing$donors_missing_stable_id == 8L))
  testthat::expect_false(any(missing$symbol_only_rescue_applied))
  recommended <- recovery[recovery$recommendation == "recommended", ]
  testthat::expect_identical(recommended$option_id,
    "common475_reference_only_refit")
  testthat::expect_equal(recommended$reference_ph_jobs, 1240L)
  testthat::expect_equal(recommended$hca_ph_jobs, 80L)
  testthat::expect_equal(nrow(validation), 17L)
  testthat::expect_true(all(validation$passed))
  testthat::expect_false(decision$hca_pca_coordinates_computed)
  testthat::expect_false(decision$ph_computed)
  testthat::expect_false(decision$landscapes_computed)
  testthat::expect_false(decision$biological_outcomes_computed)
})
