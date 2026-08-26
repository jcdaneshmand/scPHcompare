test_that("MV11-B v2 prospectively binds the corrected sentinel recovery", {
  root <- testthat::test_path("..", "..")
  audit <- file.path(root, "docs", "audits",
                     "mv11b-cell-benchmark-prefreeze-v2")
  read_audit <- function(name) read.csv(
    file.path(audit, name), stringsAsFactors = FALSE, check.names = FALSE
  )
  contract <- read_audit("mv11b-contract.csv")
  queue <- read_audit("mv11b-workload-queue.csv")
  failure <- read_audit("mv11b-prior-failure-binding.csv")
  implementation <- read_audit("mv11b-implementation-bindings.csv")
  validation <- read_audit("mv11b-validation.csv")
  decision <- read_audit("mv11b-decision.csv")
  manifest <- read_audit("mv11b-artifact-manifest.csv")
  expect_equal(contract$execution_head,
               "57d8852df0f3744351119218fae96aa63cbd1f16")
  expect_equal(contract$execution_attempt, 2L)
  expect_equal(contract$recovery_change,
               "canonical_character_catalog_equality_only")
  expect_equal(nrow(queue), 10L)
  expect_equal(c(contract$partition_fits, contract$private_assignment_rows,
                 contract$public_quality_rows, contract$public_stability_rows,
                 contract$public_primary_k_rows,
                 contract$public_method_agreement_rows),
               c(450, 55800, 450, 90, 2, 900))
  expect_equal(nrow(failure), 3L)
  failure_root <- file.path(root, "docs", "audits",
    "mv11c-cell-benchmark-sentinel-attempt1-failure-v1")
  failure_paths <- file.path(failure_root, failure$artifact)
  expect_equal(as.numeric(file.info(failure_paths)$size), failure$bytes)
  expect_equal(unname(vapply(failure_paths, function(path) digest::digest(
    file = path, algo = "sha256", serialize = FALSE
  ), character(1L))), failure$sha256)
  paths <- file.path(root, implementation$file)
  expect_equal(as.numeric(file.info(paths)$size), implementation$bytes)
  expect_equal(unname(vapply(paths, function(path) digest::digest(
    file = path, algo = "sha256", serialize = FALSE
  ), character(1L))), implementation$sha256)
  expect_equal(nrow(validation), 28L)
  expect_true(all(validation$passed))
  expect_true(decision$sentinel_execution_authorized_after_commit)
  expect_false(decision$full_execution_authorized)
  expect_false(decision$cross_view_comparison_authorized)
  artifacts <- file.path(audit, manifest$artifact)
  expect_equal(as.numeric(file.info(artifacts)$size), manifest$bytes)
  expect_equal(unname(vapply(artifacts, function(path) digest::digest(
    file = path, algo = "sha256", serialize = FALSE
  ), character(1L))), manifest$sha256)
})
