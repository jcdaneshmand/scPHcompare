test_that("MV8-ZP freezes engine-v2 repair admission while production stays closed", {
  root <- testthat::test_path("..", "..", "docs", "audits",
                             "mv08zp-landscape-kernel-repair-prefreeze-v1")
  checks <- read.csv(file.path(root, "mv08zp-validation.csv"), check.names = FALSE)
  decision <- read.csv(file.path(root, "mv08zp-decision.csv"), check.names = FALSE)
  contract <- read.csv(file.path(root, "mv08zp-contract.csv"), check.names = FALSE)
  tests <- read.csv(file.path(root, "mv08zp-test-contract.csv"), check.names = FALSE)
  expect_equal(nrow(checks), 28L)
  expect_true(all(checks$passed))
  expect_equal(contract$scientific_engine_version, 2L)
  expect_equal(decision$oracle_runs_authorized, 2L)
  expect_equal(decision$private_failed_pair_jobs_authorized, 1L)
  expect_false(decision$old_root_resume_authorized)
  expect_false(decision$old_output_reuse_authorized)
  expect_false(decision$fresh_production_authorized)
  expect_equal(nrow(tests), 9L)
})

test_that("MV8-ZP diagnostic cannot invoke production or downstream work", {
  runner <- testthat::test_path("..", "..", "scripts",
    "run_mv08zp_landscape_kernel_repair_diagnostic.R")
  builder <- testthat::test_path("..", "..", "scripts",
    "build_mv08zp_landscape_kernel_repair_prefreeze.R")
  expect_silent(parse(runner))
  expect_silent(parse(builder))
  text <- paste(readLines(runner, warn = FALSE), collapse = "\n")
  expect_false(grepl("run_mv08zf_full_landscape_production", text, fixed = TRUE))
  expect_false(grepl("unlink\\(|file.remove\\(", text, perl = TRUE))
  expect_match(text, "production_landscape_jobs = 0L", fixed = TRUE)
})
