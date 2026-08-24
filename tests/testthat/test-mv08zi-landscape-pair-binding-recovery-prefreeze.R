test_that("MV8-ZI binds private pair reconstruction without retry", {
  root <- file.path("..", "..", "docs", "audits",
                    "mv08zi-landscape-pair-binding-recovery-prefreeze-v1")
  read <- function(name) utils::read.csv(file.path(root, name),
                                          stringsAsFactors = FALSE,
                                          check.names = FALSE)
  validation <- read("mv08zi-validation.csv")
  failure <- read("mv08zi-preserved-attempt.csv")
  binding <- read("mv08zi-private-pair-binding.csv")
  implementation <- read("mv08zi-implementation-bindings.csv")
  decision <- read("mv08zi-decision.csv")

  expect_equal(nrow(validation), 25L)
  expect_true(all(validation$passed))
  expect_equal(failure$public_rows_before, 163L)
  expect_equal(failure$public_rows_after, 163L)
  expect_equal(failure$public_receipts_written, 0L)
  expect_equal(failure$landscape_jobs_run, 0L)
  expect_equal(binding$expected_pair_rows, 250L)
  expect_identical(binding$expected_pair_subset_sha256,
                   binding$observed_pair_subset_sha256)
  expect_false(binding$public_private_identities_exposed)
  expect_true(implementation$private_pair_binding_fix[[1L]])
  expect_false(any(implementation$scientific_change))
  expect_true(decision$orphan_adoption_authorized)
  expect_true(decision$private_pair_binding_required)
  expect_false(decision$landscape_recomputation_authorized)
  expect_equal(decision$automatic_retries, 0L)
  expect_equal(decision$resume_at_production_order, 165L)
})

test_that("MV8-ZI amended recovery parses and reconstructs private pairs", {
  root <- testthat::test_path("..", "..")
  recovery <- file.path(root, "scripts",
                        "recover_mv08zh_landscape_launcher_interruption.R")
  builder <- file.path(root, "scripts",
                       "build_mv08zi_landscape_pair_binding_recovery_prefreeze.R")
  expect_silent(parse(recovery))
  expect_silent(parse(builder))
  text <- paste(readLines(recovery, warn = FALSE), collapse = "\n")
  expect_match(text, "private_unit_bindings", fixed = TRUE)
  expect_match(text, ".mv08z_group_pairs", fixed = TRUE)
  expect_match(text, "expected_pairs$pair_identity_sha256", fixed = TRUE)
  expect_false(grepl("mv08z-pair-queue.csv", text, fixed = TRUE))
  expect_false(grepl("run_mv08z_landscape_chunk", text, fixed = TRUE))
})
