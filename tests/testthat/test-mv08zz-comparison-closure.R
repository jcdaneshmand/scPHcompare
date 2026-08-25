test_that("MV8-ZZ independently closes all label-closed comparisons", {
  root <- testthat::test_path(
    "..", "..", "docs", "audits",
    "mv08zz-distance-comparison-closure-v1"
  )
  manifest <- read.csv(file.path(root, "mv08zz-artifact-manifest.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
  paths <- file.path(root, manifest$artifact)
  expect_true(all(file.exists(paths)))
  expect_equal(as.numeric(file.info(paths)$size), as.numeric(manifest$bytes))
  expect_equal(unname(vapply(paths, function(path) digest::digest(
    file = path, algo = "sha256", serialize = FALSE
  ), character(1L))), manifest$sha256)
  validation <- read.csv(file.path(root, "mv08zz-validation.csv"),
                         stringsAsFactors = FALSE, check.names = FALSE)
  recompute <- read.csv(file.path(root, "mv08zz-recomputation-summary.csv"),
                        stringsAsFactors = FALSE, check.names = FALSE)
  rehash <- read.csv(file.path(root, "mv08zz-private-rehash.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)
  resource <- read.csv(file.path(root, "mv08zz-resource-summary.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
  decision <- read.csv(file.path(root, "mv08zz-decision.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
  expect_equal(nrow(validation), 12L)
  expect_true(all(validation$passed))
  expect_equal(nrow(recompute), 40L)
  expect_identical(as.integer(recompute$comparison_order), 1:40)
  expect_equal(sum(recompute$dataset_scope == "internal124"), 30L)
  expect_equal(sum(recompute$dataset_scope == "external8"), 10L)
  expect_equal(sum(recompute$homology_dimension == "H0"), 20L)
  expect_equal(sum(recompute$homology_dimension == "H1"), 20L)
  expect_true(all(recompute$independently_recomputed))
  expect_lte(max(recompute$maximum_summary_numeric_difference), 1e-12)
  expect_lte(max(recompute$maximum_neighbor_difference), 1e-12)
  expect_equal(nrow(rehash), 40L)
  expect_true(all(rehash$independently_rehashed))
  expect_true(resource$elapsed_cap_passed)
  expect_true(resource$rss_cap_passed)
  expect_true(resource$storage_cap_passed)
  expect_equal(resource$workers, 1L)
  expect_equal(resource$retries, 0L)
  expect_equal(decision$decision,
               "close_label_closed_distance_comparisons")
  expect_equal(decision$labels_outcomes_state, "closed")
  expect_equal(decision$clustering_state, "closed")
  expect_equal(decision$fusion_state, "closed")
})
