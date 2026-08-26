test_that("MV12-A prospectively freezes historical fusion feasibility", {
  root <- testthat::test_path("..", "..")
  audit <- file.path(root, "docs", "audits",
                     "mv12a-historical-fusion-prefreeze-v1")
  read_audit <- function(name) read.csv(
    file.path(audit, name), stringsAsFactors = FALSE, check.names = FALSE
  )
  contract <- read_audit("mv12a-contract.csv")
  source <- read_audit("mv12a-source-binding.csv")
  closures <- read_audit("mv12a-closure-bindings.csv")
  implementation <- read_audit("mv12a-implementation-bindings.csv")
  trigger <- read_audit("mv12a-option2-trigger.csv")
  validation <- read_audit("mv12a-validation.csv")
  decision <- read_audit("mv12a-decision.csv")
  manifest <- read_audit("mv12a-artifact-manifest.csv")

  expect_equal(contract$execution_head,
               "f077455127ab4d3b14600be9358c9803cce6133b")
  expect_equal(c(contract$matrices, contract$partition_fits,
                 contract$private_assignment_rows), c(50, 500, 62000))
  expect_equal(c(contract$public_scale_rows, contract$public_catalog_rows,
                 contract$public_quality_rows, contract$public_stability_rows,
                 contract$public_consensus_rows,
                 contract$public_primary_detail_rows,
                 contract$public_disposition_rows),
               c(20, 50, 500, 100, 300, 2, 1))
  expect_equal(contract$homology_dimensions, "H0;H1_separate")
  expect_equal(contract$gene_weights, "0;0.25;0.5;0.75;1")
  expect_equal(contract$primary_gene_weight, 0.5)
  expect_equal(contract$primary_method, "pam_dissimilarity_v1")
  expect_equal(contract$common_k, "2;3")
  expect_equal(contract$workers, 1)
  expect_equal(contract$automatic_retries, 0)

  source_path <- file.path(root, source$private_path_not_published)
  expect_equal(as.numeric(file.info(source_path)$size), source$bytes)
  expect_equal(unname(digest::digest(file = source_path, algo = "sha256",
                                     serialize = FALSE)), source$sha256)
  closure_paths <- file.path(root, closures$file)
  expect_equal(as.numeric(file.info(closure_paths)$size), closures$bytes)
  expect_equal(unname(vapply(closure_paths, function(path) digest::digest(
    file = path, algo = "sha256", serialize = FALSE
  ), character(1L))), closures$sha256)
  implementation_paths <- file.path(root, implementation$file)
  expect_equal(as.numeric(file.info(implementation_paths)$size),
               implementation$bytes)
  expect_equal(unname(vapply(implementation_paths, function(path) digest::digest(
    file = path, algo = "sha256", serialize = FALSE
  ), character(1L))), implementation$sha256)

  expect_equal(nrow(trigger), 2L)
  expect_equal(trigger$k, c(2, 3))
  expect_true(all(trigger$primary_pass_requires_both_rules))
  expect_false(any(trigger$labels_used))
  expect_false(any(trigger$outcomes_used))
  expect_equal(nrow(validation), 25L)
  expect_true(all(validation$passed))
  expect_true(decision$historical_fusion_execution_authorized_after_commit)
  expect_true(decision$option2_authorized_only_if_triggered)
  expect_false(decision$labels_authorized)
  expect_false(decision$outcomes_authorized)
  expect_false(decision$method_or_weight_selection_authorized)
  expect_false(decision$biological_claims_authorized)
  expect_false(decision$manuscript_claims_authorized)

  artifact_paths <- file.path(audit, manifest$artifact)
  expect_equal(as.numeric(file.info(artifact_paths)$size), manifest$bytes)
  expect_equal(unname(vapply(artifact_paths, function(path) digest::digest(
    file = path, algo = "sha256", serialize = FALSE
  ), character(1L))), manifest$sha256)
})
