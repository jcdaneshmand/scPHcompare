test_that("MV15 recovery freezes closure-only correction and immutable production", {
  root <- testthat::test_path("..", "..")
  audit <- file.path(root, "docs", "audits",
                     "mv15-cell-distance-comparison-closure-recovery-prefreeze-v1")
  read_at <- function(name) read.csv(
    file.path(audit, name), stringsAsFactors = FALSE, check.names = FALSE
  )
  contract <- read_at("mv15-recovery-contract.csv")
  production <- read_at("mv15-recovery-production-binding.csv")
  implementation <- read_at("mv15-recovery-implementation-binding.csv")
  validation <- read_at("mv15-recovery-validation.csv")
  decision <- read_at("mv15-recovery-decision.csv")
  manifest <- read_at("mv15-recovery-artifact-manifest.csv")

  expect_equal(nrow(validation), 14L)
  expect_true(all(validation$passed))
  expect_equal(contract$immutable_production_comparisons, 36L)
  expect_equal(contract$immutable_global_rows, 36L)
  expect_equal(contract$immutable_neighbor_rows, 42L)
  expect_false(contract$production_rerun_authorized)
  expect_true(contract$closure_rerun_authorized_after_commit)
  expect_equal(sum(implementation$changed_for_recovery), 1L)
  expect_equal(implementation$file[implementation$changed_for_recovery],
               "scripts/build_mv15_cell_distance_comparison_closure.R")
  expect_true(all(implementation$original_sha256[
    !implementation$changed_for_recovery] == implementation$recovery_sha256[
      !implementation$changed_for_recovery]))
  expect_equal(nrow(production), 5L)
  expect_true(all(production$immutable))
  expect_true(decision$immutable_production_reused)
  expect_false(decision$production_rerun_authorized)
  expect_true(decision$corrected_closure_authorized_after_commit)
  expect_false(decision$result_interpretation_authorized)

  files <- file.path(audit, manifest$artifact)
  expect_equal(as.numeric(file.info(files)$size), manifest$bytes)
  expect_equal(unname(vapply(files, function(path) digest::digest(
    file = path, algo = "sha256", serialize = FALSE
  ), character(1L))), manifest$sha256)
})
