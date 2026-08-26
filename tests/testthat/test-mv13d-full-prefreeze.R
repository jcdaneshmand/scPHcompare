test_that("MV13-D freezes full cell PH without opening landscapes", {
  root <- testthat::test_path("..", "..")
  audit <- file.path(root, "docs", "audits",
                     "mv13d-allqc-cell-full-prefreeze-v1")
  read_audit <- function(name) read.csv(
    file.path(audit, name), stringsAsFactors = FALSE, check.names = FALSE
  )
  contract <- read_audit("mv13d-contract.csv")
  groups <- read_audit("mv13d-group-queue.csv")
  validation <- read_audit("mv13d-validation.csv")
  decision <- read_audit("mv13d-decision.csv")
  implementation <- read_audit("mv13d-implementation-bindings.csv")
  manifest <- read_audit("mv13d-artifact-manifest.csv")
  expect_equal(contract$execution_head,
               "cf544ebf4faf0afbd554b977134c557898d3ed7c")
  expect_equal(c(contract$pca_models, contract$adopted_models,
                 contract$new_models, contract$cell_views,
                 contract$adopted_views, contract$new_views,
                 contract$dimension_records), c(7, 1, 6, 636, 1, 635, 1272))
  expect_equal(nrow(groups), 7L)
  expect_equal(sum(groups$unit_count), 636)
  expect_equal(sum(groups$adopt_closed_model), 1)
  expect_equal(sum(groups$new_model_count), 6)
  expect_equal(sum(groups$new_view_count), 635)
  expect_true(contract$full_execution_authorized)
  expect_true(contract$independent_closure_required)
  expect_false(any(c(contract$landscapes_authorized,
                     contract$comparisons_authorized,
                     contract$clustering_authorized,
                     contract$fusion_authorized,
                     contract$labels_authorized,
                     contract$outcomes_authorized,
                     contract$biological_claims_authorized,
                     contract$manuscript_claims_authorized)))
  expect_equal(nrow(validation), 25L)
  expect_true(all(validation$passed))
  expect_true(decision$full_execution_authorized_after_commit)
  expect_true(decision$independent_closure_required)
  expect_false(decision$landscapes_authorized)
  paths <- file.path(root, implementation$file)
  expect_equal(as.numeric(file.info(paths)$size), implementation$bytes)
  expect_equal(unname(vapply(paths, function(path) digest::digest(
    file = path, algo = "sha256", serialize = FALSE
  ), character(1L))), implementation$sha256)
  artifacts <- file.path(audit, manifest$artifact)
  expect_equal(as.numeric(file.info(artifacts)$size), manifest$bytes)
  expect_equal(unname(vapply(artifacts, function(path) digest::digest(
    file = path, algo = "sha256", serialize = FALSE
  ), character(1L))), manifest$sha256)
})
