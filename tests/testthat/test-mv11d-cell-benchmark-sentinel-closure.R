test_that("MV11-D exactly closes the corrected cell sentinel", {
  root <- testthat::test_path("..", "..")
  audit <- file.path(root, "docs", "audits",
                     "mv11d-cell-benchmark-sentinel-closure-v1")
  read_audit <- function(name) read.csv(
    file.path(audit, name), stringsAsFactors = FALSE, check.names = FALSE
  )
  validation <- read_audit("mv11d-validation.csv")
  summary <- read_audit("mv11d-sentinel-closure-summary.csv")
  decision <- read_audit("mv11d-decision.csv")
  manifest <- read_audit("mv11d-artifact-manifest.csv")
  expect_equal(nrow(validation), 20L)
  expect_true(all(validation$passed))
  expect_equal(summary$execution_head,
               "57d8852df0f3744351119218fae96aa63cbd1f16")
  expect_equal(summary$catalog_order, 10L)
  expect_equal(summary$seed, 20260809L)
  expect_equal(summary$homology_dimension, "H1")
  expect_equal(c(summary$partition_rows, summary$quality_rows), c(5580, 45))
  expect_true(summary$exact_repeat)
  expect_false(summary$labels_used)
  expect_false(summary$outcomes_used)
  expect_false(summary$cross_view_comparison_performed)
  expect_true(decision$sentinel_independently_closed)
  expect_true(decision$full_execution_authorized_after_commit)
  expect_false(decision$labels_authorized)
  expect_false(decision$outcomes_authorized)
  expect_false(decision$cross_view_comparison_authorized)
  artifacts <- file.path(audit, manifest$artifact)
  expect_equal(as.numeric(file.info(artifacts)$size), manifest$bytes)
  expect_equal(unname(vapply(artifacts, function(path) digest::digest(
    file = path, algo = "sha256", serialize = FALSE
  ), character(1L))), manifest$sha256)
})
