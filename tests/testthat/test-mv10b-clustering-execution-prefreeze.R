test_that("MV10-B freezes only one sentinel before full execution", {
  root <- testthat::test_path("..", "..")
  audit <- file.path(
    root, "docs", "audits", "mv10b-clustering-execution-prefreeze-v1"
  )
  read_audit <- function(name) read.csv(
    file.path(audit, name), stringsAsFactors = FALSE, check.names = FALSE
  )
  manifest <- read_audit("mv10b-artifact-manifest.csv")
  paths <- file.path(audit, manifest$artifact)
  expect_true(all(file.exists(paths)))
  expect_equal(as.numeric(file.info(paths)$size), as.numeric(manifest$bytes))
  expect_equal(unname(vapply(paths, function(path) digest::digest(
    file = path, algo = "sha256", serialize = FALSE
  ), character(1L))), manifest$sha256)

  contract <- read_audit("mv10b-contract.csv")
  queue <- read_audit("mv10b-workload-queue.csv")
  sentinel <- read_audit("mv10b-sentinel-queue.csv")
  resource <- read_audit("mv10b-resource-policy.csv")
  closure <- read_audit("mv10b-prospective-closure.csv")
  implementation <- read_audit("mv10b-implementation-bindings.csv")
  sources <- read_audit("mv10b-source-freeze.csv")
  decision <- read_audit("mv10b-decision.csv")
  validation <- read_audit("mv10b-validation.csv")

  expect_equal(contract$execution_head,
               "06c0dd6c129d8567651081156272cfe972e52c7d")
  expect_equal(contract$matrices, 30L)
  expect_equal(contract$sentinel_matrices, 1L)
  expect_equal(contract$sentinel_partition_fits, 45L)
  expect_equal(contract$full_partition_fits, 1350L)
  expect_equal(contract$full_private_assignment_rows, 167400L)
  expect_true(contract$H0_H1_separate)
  expect_false(contract$cell_gene_combined)

  expect_equal(nrow(queue), 30L)
  expect_identical(queue$execution_order, 1:30)
  expect_equal(nrow(sentinel), 1L)
  expect_equal(sentinel$stack_id, "allqc_residual_exact500")
  expect_equal(sentinel$homology_dimension, "H1")
  expect_equal(sentinel$seed, 20260809L)
  expect_match(sentinel$sentinel_rationale, "without_cluster_results",
               fixed = TRUE)
  expect_equal(resource$workers, 1L)
  expect_equal(resource$automatic_retries, 0L)
  expect_equal(resource$child_elapsed_cap_seconds, 1800L)
  expect_equal(resource$process_tree_rss_cap_bytes, 4 * 1024^3)
  expect_equal(resource$private_storage_cap_bytes, 512 * 1024^2)
  expect_true(all(closure$requires_source_reload))
  expect_true(all(closure$requires_private_partition_recomputation))

  implementation_paths <- file.path(root, implementation$file)
  expect_true(all(file.exists(implementation_paths)))
  expect_true(all(grepl("^[0-9a-f]{64}$", implementation$sha256)))
  source_paths <- file.path(root, sources$artifact)
  expect_true(all(file.exists(source_paths)))
  expect_equal(unname(vapply(source_paths, function(path) digest::digest(
    file = path, algo = "sha256", serialize = FALSE
  ), character(1L))), sources$sha256)

  expect_equal(nrow(validation), 42L)
  expect_true(all(validation$passed))
  expect_true(decision$sentinel_execution_authorized_after_commit)
  expect_false(decision$full_execution_authorized)
  expect_false(decision$labels_authorized)
  expect_false(decision$outcomes_authorized)
  expect_false(decision$biological_interpretation_authorized)
  expect_false(decision$manuscript_claims_authorized)
})
