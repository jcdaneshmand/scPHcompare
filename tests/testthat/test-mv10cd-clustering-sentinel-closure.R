test_that("MV10-C/D closes the sentinel and admits only label-closed execution", {
  root <- testthat::test_path("..", "..")
  sentinel <- file.path(root, "docs", "audits",
                        "mv10c-clustering-sentinel-v1")
  closure <- file.path(root, "docs", "audits",
                       "mv10d-clustering-sentinel-closure-v1")
  verify_manifest <- function(path, name) {
    manifest <- read.csv(file.path(path, name), stringsAsFactors = FALSE,
                         check.names = FALSE)
    files <- file.path(path, manifest$artifact)
    expect_true(all(file.exists(files)))
    expect_equal(as.numeric(file.info(files)$size), as.numeric(manifest$bytes))
    expect_equal(unname(vapply(files, function(file) digest::digest(
      file = file, algo = "sha256", serialize = FALSE
    ), character(1L))), manifest$sha256)
  }
  verify_manifest(sentinel, "mv10c-artifact-manifest.csv")
  verify_manifest(closure, "mv10d-artifact-manifest.csv")
  receipt <- read.csv(file.path(sentinel, "mv10c-resource-receipt.csv"),
                      stringsAsFactors = FALSE, check.names = FALSE)
  quality <- read.csv(file.path(sentinel, "mv10c-partition-quality.csv"),
                      stringsAsFactors = FALSE, check.names = FALSE)
  validation <- read.csv(file.path(closure, "mv10d-validation.csv"),
                         stringsAsFactors = FALSE, check.names = FALSE)
  decision <- read.csv(file.path(closure, "mv10d-decision.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
  rehash <- read.csv(file.path(closure, "mv10d-private-rehash.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)

  expect_equal(receipt$completion_state, "complete")
  expect_equal(receipt$partition_fits, 45L)
  expect_equal(receipt$private_assignment_rows, 5580L)
  expect_equal(receipt$quality_rows, 45L)
  expect_gt(receipt$elapsed_seconds, 0)
  expect_gt(receipt$peak_process_tree_rss_bytes, 0)
  expect_lte(receipt$elapsed_seconds, receipt$elapsed_cap_seconds)
  expect_lte(receipt$peak_process_tree_rss_bytes,
             receipt$process_tree_rss_cap_bytes)
  expect_lte(receipt$private_bytes, receipt$private_storage_cap_bytes)
  expect_equal(receipt$stderr_bytes, 0)
  expect_equal(receipt$workers, 1L)
  expect_equal(receipt$retries, 0L)
  expect_equal(nrow(quality), 45L)
  expect_equal(nrow(rehash), 3L)
  expect_equal(nrow(validation), 26L)
  expect_true(all(validation$passed))
  expect_true(decision$full_execution_authorized_after_commit)
  expect_equal(decision$maximum_workers, 1L)
  expect_equal(decision$automatic_retries, 0L)
  expect_false(decision$labels_authorized)
  expect_false(decision$outcomes_authorized)
  expect_false(decision$biological_interpretation_authorized)
  expect_false(decision$manuscript_claims_authorized)
})
