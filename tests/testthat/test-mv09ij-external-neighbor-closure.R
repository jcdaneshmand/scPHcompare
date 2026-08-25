test_that("MV9-I/J publish and independently close informative neighborhoods", {
  production <- testthat::test_path(
    "..", "..", "docs", "audits",
    "mv09i-external-neighbor-sensitivity-v1"
  )
  closure <- testthat::test_path(
    "..", "..", "docs", "audits",
    "mv09j-external-neighbor-closure-v1"
  )
  verify_manifest <- function(root, filename) {
    manifest <- read.csv(file.path(root, filename), stringsAsFactors = FALSE,
                         check.names = FALSE)
    paths <- file.path(root, manifest$artifact)
    expect_true(all(file.exists(paths)))
    expect_equal(as.numeric(file.info(paths)$size), as.numeric(manifest$bytes))
    expect_equal(unname(vapply(paths, function(path) digest::digest(
      file = path, algo = "sha256", serialize = FALSE
    ), character(1L))), manifest$sha256)
  }
  verify_manifest(production, "mv09i-artifact-manifest.csv")
  verify_manifest(closure, "mv09j-artifact-manifest.csv")

  summary <- read.csv(file.path(
    production, "mv09i-external-neighbor-summary.csv"
  ), stringsAsFactors = FALSE, check.names = FALSE)
  proof <- read.csv(file.path(
    production, "mv09i-degeneracy-classification.csv"
  ), stringsAsFactors = FALSE, check.names = FALSE)
  receipt <- read.csv(file.path(production, "mv09i-terminal-receipt.csv"),
                      stringsAsFactors = FALSE, check.names = FALSE)
  rehash <- read.csv(file.path(closure, "mv09j-private-rehash.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)
  validation <- read.csv(file.path(closure, "mv09j-validation.csv"),
                         stringsAsFactors = FALSE, check.names = FALSE)
  decision <- read.csv(file.path(closure, "mv09j-decision.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)

  expect_equal(nrow(summary), 20L)
  expect_equal(length(unique(summary$comparison_id)), 10L)
  expect_identical(sort(unique(summary$k)), c(2L, 3L))
  expect_identical(sort(unique(summary$homology_dimension)), c("H0", "H1"))
  expect_true(all(summary$mean_neighbor_jaccard >= 0 &
                    summary$mean_neighbor_jaccard <= 1))
  expect_gt(length(unique(summary$mean_neighbor_jaccard)), 2L)
  expect_true(all(summary$replication_units == 1L))
  expect_false("unit_id" %in% names(summary))

  expect_true(proof$k_equals_all_other_units)
  expect_equal(proof$possible_neighbor_sets_per_unit, 1)
  expect_equal(proof$jaccard_for_any_two_complete_rankings, 1)
  expect_false(proof$informative_for_neighborhood_preservation)
  expect_equal(receipt$completion_state, "complete")
  expect_equal(receipt$comparisons, 10L)
  expect_equal(receipt$summary_rows, 20L)
  expect_equal(receipt$private_unit_rows, 160L)
  expect_equal(receipt$workers, 1L)
  expect_equal(receipt$retries, 0L)
  expect_false(receipt$labels_used)
  expect_false(receipt$outcomes_used)
  expect_false(receipt$biological_claims)
  expect_false(receipt$manuscript_claims)

  expect_equal(nrow(rehash), 10L)
  expect_true(all(rehash$unit_axis_identical))
  expect_lte(max(rehash$independent_summary_difference), 1e-12)
  expect_lte(max(rehash$independent_unit_difference), 1e-12)
  expect_equal(nrow(validation), 19L)
  expect_true(all(validation$passed))
  expect_lte(decision$maximum_numeric_difference, 1e-12)
  expect_equal(decision$decision,
               "close_k2_k3_sensitivity_and_exclude_k7")
  expect_equal(decision$biological_interpretation_state, "closed")
  expect_equal(decision$manuscript_claim_state, "closed")
})
