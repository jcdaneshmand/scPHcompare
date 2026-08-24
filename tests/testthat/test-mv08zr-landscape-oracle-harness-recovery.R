test_that("MV8-ZR confines failure to the legacy engine predicate", {
  root <- testthat::test_path("..", "..", "docs", "audits",
    "mv08zr-landscape-oracle-harness-recovery-prefreeze-v1")
  checks <- read.csv(file.path(root, "mv08zr-validation.csv"), check.names = FALSE)
  evidence <- read.csv(file.path(root, "mv08zr-failure-evidence.csv"),
                       check.names = FALSE)
  decision <- read.csv(file.path(root, "mv08zr-decision.csv"), check.names = FALSE)
  expect_equal(nrow(checks), 20L)
  expect_true(all(checks$passed))
  expect_equal(evidence$scientific_passes, 28L)
  expect_equal(evidence$engine_valid_passes, 0L)
  expect_equal(decision$expected_engine_version, 2L)
  expect_true(decision$preserve_eight_argument_engine1_interface)
  expect_false(decision$scientific_kernel_change_authorized)
  expect_false(decision$candidate_hash_change_authorized)
  expect_false(decision$failed_run_rerun_authorized)
  expect_equal(decision$production_landscape_jobs_authorized, 0L)
})

test_that("MV8-ZR builder is read-only with respect to production", {
  path <- testthat::test_path("..", "..", "scripts",
    "build_mv08zr_landscape_oracle_harness_recovery_prefreeze.R")
  expect_silent(parse(path))
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")
  expect_false(grepl("unlink\\(|file.remove\\(", text, perl = TRUE))
  expect_match(text, "failed_run_rerun_authorized = FALSE", fixed = TRUE)
})
