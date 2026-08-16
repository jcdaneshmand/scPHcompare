testthat::test_that("MV8-F recovery specification freezes exact cache identity", {
  root <- testthat::test_path("..", "..")
  spec <- paste(readLines(file.path(root, "docs", "specifications",
    "MV08F_REFERENCE_CACHE_RECOVERY_PREFREEZE_V1.md"), warn = FALSE),
    collapse = "\n")
  testthat::expect_match(spec, "90 samples, five seeds", fixed = TRUE)
  testthat::expect_match(spec, "selected\\s+cell SHA-256")
  testthat::expect_match(spec, "payload\\s+SHA-256")
  testthat::expect_match(spec, "file SHA-256", fixed = TRUE)
  testthat::expect_match(spec, "does not\\s+authorize PH")
  testthat::expect_match(spec, "HCA FASTQs", fixed = TRUE)
})

testthat::test_that("MV8-F builder preserves the label and landscape boundary", {
  root <- testthat::test_path("..", "..")
  builder <- paste(readLines(file.path(root, "scripts",
    "build_mv08f_cache_recovery_prefreeze.R"), warn = FALSE), collapse = "\n")
  testthat::expect_match(builder, "nrow(primary) != 450L", fixed = TRUE)
  testthat::expect_match(builder, "nrow(added) != 170L", fixed = TRUE)
  testthat::expect_match(builder, "runtime_equal <- identical", fixed = TRUE)
  testthat::expect_match(builder, "nrow(ph_metrics) != 1240L", fixed = TRUE)
  testthat::expect_match(builder, "nrow(landscape_inventory) != 20L",
    fixed = TRUE)
  testthat::expect_match(builder,
    "decision = \"authorize_primary90_cache_reconstruction_only\"",
    fixed = TRUE)
  testthat::expect_match(builder, "hca_fastq_download_authorized = FALSE",
    fixed = TRUE)
  testthat::expect_match(builder, "label_access_authorized = FALSE",
    fixed = TRUE)
  testthat::expect_match(builder,
    "validator = \"scripts/validate_mv08f_cache_recovery_prefreeze.R\"",
    fixed = TRUE)
})

testthat::test_that("MV8-F published prefreeze evidence is independently closed", {
  evidence <- testthat::test_path("..", "..", "docs", "audits",
    "mv08f-recovery-prefreeze-evidence")
  axis <- utils::read.csv(file.path(evidence, "mv08f-recovery-axis.csv"),
    stringsAsFactors = FALSE, check.names = FALSE)
  decision <- utils::read.csv(file.path(evidence, "mv08f-decision.csv"),
    stringsAsFactors = FALSE, check.names = FALSE)
  validation <- utils::read.csv(file.path(evidence,
    "mv08f-independent-validation.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  testthat::expect_equal(nrow(axis), 450L)
  testthat::expect_equal(length(unique(axis$sample_id)), 90L)
  testthat::expect_equal(as.integer(table(axis$seed)), rep(90L, 5L))
  testthat::expect_equal(decision$raw_jobs_authorized, 90L)
  testthat::expect_equal(decision$sct_jobs_authorized, 450L)
  testthat::expect_false(decision$panel475_source_authorized)
  testthat::expect_false(decision$ph_authorized)
  testthat::expect_false(decision$hca_fastq_download_authorized)
  testthat::expect_equal(nrow(validation), 10L)
  testthat::expect_true(all(validation$passed))
})
