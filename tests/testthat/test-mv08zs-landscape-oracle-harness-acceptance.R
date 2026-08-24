test_that("MV8-ZS accepts only the backward-compatible harness amendment", {
  root <- testthat::test_path("..", "..", "docs", "audits",
    "mv08zs-landscape-oracle-harness-recovery-acceptance-v1")
  checks <- read.csv(file.path(root, "mv08zs-validation.csv"), check.names = FALSE)
  binding <- read.csv(file.path(root, "mv08zs-harness-amendment.csv"),
                      check.names = FALSE)
  decision <- read.csv(file.path(root, "mv08zs-decision.csv"), check.names = FALSE)
  expect_equal(nrow(checks), 18L)
  expect_true(all(checks$passed))
  expect_false(binding$scientific_kernel_changed)
  expect_equal(binding$backward_compatible_eight_argument_default, 1L)
  expect_equal(decision$expected_engine_version, 2L)
  expect_equal(decision$fresh_oracle_runs_authorized, 2L)
  expect_false(decision$failed_run_rerun_authorized)
  expect_equal(decision$production_landscape_jobs, 0L)
})

test_that("oracle runner retains engine-1 default and explicit engine parameter", {
  path <- testthat::test_path("..", "..", "scripts",
    "run_mv08x_rust_landscape_oracles.R")
  expect_silent(parse(path))
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")
  expect_match(text, "!length(args) %in% 8:9", fixed = TRUE)
  expect_match(text, "else 1L", fixed = TRUE)
  expect_match(text, "candidate$engine_version == expected_engine_version", fixed = TRUE)
})
