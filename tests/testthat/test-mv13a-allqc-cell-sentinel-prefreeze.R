test_that("MV13-A freezes only one all-QC cell sentinel", {
  root <- testthat::test_path("..", "..")
  audit <- file.path(root, "docs", "audits",
                     "mv13a-allqc-cell-sentinel-prefreeze-v1")
  read_audit <- function(name) read.csv(
    file.path(audit, name), stringsAsFactors = FALSE, check.names = FALSE
  )
  contract <- read_audit("mv13a-contract.csv")
  sources <- read_audit("mv13a-source-bindings.csv")
  locator <- read_audit("mv13a-private-locator-binding.csv")
  pca <- read_audit("mv13a-pca-queue.csv")
  views <- read_audit("mv13a-view-queue.csv")
  ph <- read_audit("mv13a-ph-queue.csv")
  sentinel <- read_audit("mv13a-sentinel.csv")
  validation <- read_audit("mv13a-validation.csv")
  manifest <- read_audit("mv13a-artifact-manifest.csv")

  expect_equal(contract$execution_head,
               "2b498e672c721f7f07dc2fa032b780ce979863ad")
  expect_equal(c(contract$source_caches, contract$pca_models,
                 contract$cell_views, contract$ph_records),
               c(132, 7, 636, 1272))
  expect_false(contract$full_execution_authorized)
  expect_true(contract$sentinel_authorized_after_commit)
  expect_false(any(c(contract$landscapes_authorized,
                     contract$comparisons_authorized,
                     contract$clustering_authorized,
                     contract$fusion_authorized,
                     contract$labels_authorized,
                     contract$outcomes_authorized,
                     contract$biological_claims_authorized,
                     contract$manuscript_claims_authorized)))
  expect_equal(nrow(sources), 132L)
  expect_false("private_cache_path" %in% names(sources))
  expect_true(all(sources$private_locator_state ==
                    "private_hash_bound_not_published"))
  expect_equal(locator$rows, 132L)
  expect_equal(locator$publication_state, "private_not_tracked")
  expect_equal(nrow(pca), 7L)
  expect_equal(nrow(views), 636L)
  expect_equal(nrow(ph), 1272L)
  expect_equal(sentinel$unit_id, "SRA742961_SRS3565197")
  expect_equal(sentinel$seed, 20260809)
  expect_equal(sentinel$pca_fit_cells, 47616)
  expect_equal(sentinel$rss_cap_bytes, 8 * 1024^3)
  expect_equal(sentinel$workers, 1)
  expect_equal(sentinel$retries, 0)
  expect_equal(nrow(validation), 25L)
  expect_true(all(validation$passed))
  paths <- file.path(audit, manifest$artifact)
  expect_equal(as.numeric(file.info(paths)$size), manifest$bytes)
  expect_equal(unname(vapply(paths, function(path) digest::digest(
    file = path, algo = "sha256", serialize = FALSE
  ), character(1L))), manifest$sha256)
})
