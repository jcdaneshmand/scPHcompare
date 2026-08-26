test_that("MV15 v2 independently closes all 36 comparisons", {
  root <- testthat::test_path("..", "..")
  audit <- file.path(root, "docs", "audits",
                     "mv15-cell-distance-comparison-closure-v2")
  read_at <- function(name) read.csv(
    file.path(audit, name), stringsAsFactors = FALSE, check.names = FALSE
  )
  validation <- read_at("mv15-validation.csv")
  recompute <- read_at("mv15-recomputation.csv")
  rehash <- read_at("mv15-private-rehash.csv")
  resource <- read_at("mv15-resource-closure.csv")
  decision <- read_at("mv15-decision.csv")
  manifest <- read_at("mv15-artifact-manifest.csv")

  expect_equal(nrow(validation), 22L)
  expect_true(all(validation$passed))
  expect_equal(nrow(recompute), 36L)
  expect_equal(sum(recompute$contrast_family == "cell_seed_stability"), 20L)
  expect_equal(sum(recompute$contrast_family == "cell_panel_sensitivity"), 2L)
  expect_equal(sum(recompute$contrast_family == "cell_gene_view_agreement"),
               14L)
  expect_equal(as.integer(table(recompute$homology_dimension)), c(18L, 18L))
  expect_true(all(recompute$pair_axis_identical))
  expect_true(all(recompute$independently_recomputed))
  expect_lte(max(recompute$maximum_global_difference,
                 recompute$maximum_neighbor_summary_difference,
                 recompute$maximum_unit_neighbor_difference), 1e-12)
  expect_equal(nrow(rehash), 36L)
  expect_true(all(rehash$independently_rehashed))
  expect_equal(resource$GNU_time_exit_status, 0L)
  expect_lte(resource$total_elapsed_seconds, resource$elapsed_cap_seconds)
  expect_lte(resource$peak_rss_bytes, resource$rss_cap_bytes)
  expect_lte(resource$private_bytes, resource$private_storage_cap_bytes)
  expect_lte(resource$public_bytes, resource$public_storage_cap_bytes)
  expect_equal(resource$workers, 1L)
  expect_equal(resource$retries, 0L)
  expect_true(decision$distance_comparisons_independently_closed)
  expect_true(decision$values_not_interpreted_by_closure)
  expect_true(decision$descriptive_synthesis_prefreeze_eligible_next)
  expect_false(any(decision[c(
    "clustering_authorized", "fusion_authorized", "labels_authorized",
    "outcomes_authorized", "biological_claims_authorized",
    "manuscript_claims_authorized"
  )]))

  files <- file.path(audit, manifest$artifact)
  expect_equal(as.numeric(file.info(files)$size), manifest$bytes)
  expect_equal(unname(vapply(files, function(path) digest::digest(
    file = path, algo = "sha256", serialize = FALSE
  ), character(1L))), manifest$sha256)
})
