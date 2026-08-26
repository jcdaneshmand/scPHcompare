test_that("MV10-G freezes comprehensive review before opening values", {
  root <- testthat::test_path("..", "..")
  audit <- file.path(root, "docs", "audits",
                     "mv10g-clustering-review-prefreeze-v1")
  read_audit <- function(name) read.csv(
    file.path(audit, name), stringsAsFactors = FALSE, check.names = FALSE
  )
  manifest <- read_audit("mv10g-artifact-manifest.csv")
  files <- file.path(audit, manifest$artifact)
  expect_true(all(file.exists(files)))
  expect_equal(as.numeric(file.info(files)$size), as.numeric(manifest$bytes))
  expect_equal(unname(vapply(files, function(file) digest::digest(
    file = file, algo = "sha256", serialize = FALSE
  ), character(1L))), manifest$sha256)
  contract <- read_audit("mv10g-contract.csv")
  schemas <- read_audit("mv10g-source-schema.csv")
  metrics <- read_audit("mv10g-metric-contract.csv")
  aggregation <- read_audit("mv10g-aggregation-contract.csv")
  figures <- read_audit("mv10g-figure-contract.csv")
  review <- read_audit("mv10g-review-contract.csv")
  implementation <- read_audit("mv10g-implementation-bindings.csv")
  sources <- read_audit("mv10g-source-freeze.csv")
  decision <- read_audit("mv10g-decision.csv")
  validation <- read_audit("mv10g-validation.csv")

  expect_equal(contract$execution_head,
               "f3281e12574942935915fadda93f78268dc83cbe")
  expect_equal(contract$source_matrices, 30L)
  expect_equal(contract$source_partition_fits, 1350L)
  expect_equal(contract$figures, 7L)
  expect_false(contract$source_values_opened_before_metric_and_figure_freeze)
  expect_equal(nrow(schemas), 4L)
  expect_identical(schemas$source_id,
                   c("quality", "stability", "primary_k", "agreement"))
  expect_identical(schemas$expected_rows, c(1350L, 270L, 2L, 2700L))
  expect_true(all(!schemas$scientific_values_opened_during_prefreeze))
  expect_equal(nrow(metrics), 8L)
  expect_equal(nrow(aggregation), 6L)
  expect_identical(aggregation$expected_rows,
                   c(270L, 270L, 540L, 2L, 18L, 90L))
  expect_true(all(!aggregation$result_dependent_filter))
  expect_equal(nrow(figures), 7L)
  expect_true(all(figures$format == "png"))
  expect_true(all(figures$H0_H1_simultaneous))
  expect_true(all(figures$complete_k2_k10))
  expect_equal(nrow(review), 13L)
  expect_equal(nrow(implementation), 8L)
  expect_true(all(file.exists(file.path(root, implementation$file))))
  expect_true(all(grepl("^[0-9a-f]{64}$", implementation$sha256)))
  source_files <- file.path(root, sources$artifact)
  expect_true(all(file.exists(source_files)))
  expect_equal(unname(vapply(source_files, function(file) digest::digest(
    file = file, algo = "sha256", serialize = FALSE
  ), character(1L))), sources$sha256)
  expect_equal(nrow(validation), 40L)
  expect_true(all(validation$passed))
  expect_true(decision$synthesis_authorized_after_commit)
  expect_false(decision$figure_render_authorized)
  expect_false(decision$scientific_interpretation_authorized)
  expect_false(decision$biological_interpretation_authorized)
  expect_false(decision$manuscript_claims_authorized)
})
