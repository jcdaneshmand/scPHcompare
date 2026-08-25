test_that("MV9-H prospectively freezes informative external neighborhoods", {
  root <- testthat::test_path("..", "..", "docs", "audits",
                              "mv09h-external-neighbor-prefreeze-v1")
  manifest <- read.csv(file.path(root, "mv09h-artifact-manifest.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
  paths <- file.path(root, manifest$artifact)
  expect_true(all(file.exists(paths)))
  expect_equal(as.numeric(file.info(paths)$size), as.numeric(manifest$bytes))
  expect_equal(unname(vapply(paths, function(path) digest::digest(
    file = path, algo = "sha256", serialize = FALSE
  ), character(1L))), manifest$sha256)

  contract <- read.csv(file.path(root, "mv09h-contract.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
  proof <- read.csv(file.path(root, "mv09h-degeneracy-proof.csv"),
                    stringsAsFactors = FALSE, check.names = FALSE)
  queue <- read.csv(file.path(root, "mv09h-external-comparison-queue.csv"),
                    stringsAsFactors = FALSE, check.names = FALSE)
  figures <- read.csv(file.path(root, "mv09h-figure-contract.csv"),
                      stringsAsFactors = FALSE, check.names = FALSE)
  validation <- read.csv(file.path(root, "mv09h-validation.csv"),
                         stringsAsFactors = FALSE, check.names = FALSE)
  implementation <- read.csv(file.path(
    root, "mv09h-implementation-bindings.csv"
  ), stringsAsFactors = FALSE, check.names = FALSE)

  expect_equal(contract$external_units, 8L)
  expect_equal(contract$comparisons, 10L)
  expect_equal(contract$frozen_sensitivity_k, "2;3")
  expect_equal(contract$structurally_degenerate_k, 7L)
  expect_equal(contract$k_selection_timing, "prospective_before_calculation")
  expect_false(contract$result_dependent_selection)
  expect_equal(contract$internal_k10_state, "preserved_unchanged")
  expect_false(contract$labels_used)
  expect_false(contract$outcomes_used)
  expect_false(contract$biological_claims)
  expect_false(contract$manuscript_claims)

  expect_true(proof$k_equals_all_other_units)
  expect_equal(proof$possible_neighbor_sets_per_unit, 1)
  expect_equal(proof$jaccard_for_any_two_complete_rankings, 1)
  expect_false(proof$informative_for_neighborhood_preservation)
  expect_equal(nrow(queue), 10L)
  expect_identical(as.integer(queue$comparison_order), 31:40)
  expect_equal(nrow(figures), 4L)
  expect_false(any(figures$external_k7_displayed))
  expect_equal(sum(figures$external_k2_k3_displayed), 1L)
  expect_true(all(figures$H0_H1_separate))
  expect_true(all(figures$format == "PNG_only_no_PDF"))
  expect_equal(nrow(validation), 20L)
  expect_true(all(validation$passed))

  implementation_paths <- testthat::test_path("..", "..",
                                               implementation$file)
  expect_equal(unname(vapply(implementation_paths, function(path) {
    digest::digest(file = path, algo = "sha256", serialize = FALSE)
  }, character(1L))), implementation$sha256)
})
