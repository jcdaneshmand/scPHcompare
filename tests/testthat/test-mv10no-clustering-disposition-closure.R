test_that("MV10-N/O closes the label-closed methodological disposition", {
  root <- testthat::test_path("..", "..")
  production <- file.path(root, "docs", "audits",
                          "mv10n-clustering-disposition-v1")
  closure <- file.path(root, "docs", "audits",
                       "mv10o-clustering-disposition-closure-v1")
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
  verify_manifest(production, "mv10n-artifact-manifest.csv")
  verify_manifest(closure, "mv10o-artifact-manifest.csv")
  summary <- read.csv(file.path(production, "mv10n-primary-summary.csv"),
                      stringsAsFactors = FALSE, check.names = FALSE)
  disposition <- read.csv(file.path(production, "mv10n-disposition.csv"),
                          stringsAsFactors = FALSE, check.names = FALSE)
  receipt <- read.csv(file.path(production, "mv10n-terminal-receipt.csv"),
                      stringsAsFactors = FALSE, check.names = FALSE)
  validation <- read.csv(file.path(closure, "mv10o-validation.csv"),
                         stringsAsFactors = FALSE, check.names = FALSE)
  decision <- read.csv(file.path(closure, "mv10o-decision.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
  rehash <- read.csv(file.path(closure, "mv10o-rehash.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)

  expect_equal(receipt$completion_state, "complete")
  expect_equal(receipt$output_tables, 5L)
  expect_equal(nrow(summary), 2L)
  expect_identical(summary$homology_dimension, c("H0", "H1"))
  expect_identical(summary$selected_k, c(2L, 3L))
  expect_true(all(summary$structurally_nondegenerate))
  expect_true(all(summary$minimum_cluster_size_across_seeds >= 2L))
  expect_true(all(summary$maximum_singleton_clusters_across_seeds == 0L))
  expect_false(any(summary$silhouette_threshold_applied))
  expect_equal(disposition$decision,
               "retain_frozen_PAM_for_separate_label_evaluation_prefreeze")
  expect_true(disposition$H0_H1_remain_separate)
  expect_true(disposition$all_selected_partitions_structurally_nondegenerate)
  expect_true(disposition$internal_only)
  expect_false(disposition$labels_used)
  expect_false(disposition$outcomes_used)
  expect_false(disposition$biological_interpretation)
  expect_false(disposition$manuscript_claims)
  expect_equal(nrow(rehash), 5L)
  expect_true(all(rehash$exact_column_names))
  expect_true(all(rehash$numeric_repeat))
  expect_equal(nrow(validation), 25L)
  expect_true(all(validation$passed))
  expect_equal(decision$next_stage, "separate_label_evaluation_prefreeze")
  expect_equal(decision$labels_outcomes_state, "closed")
  expect_equal(decision$biological_interpretation_state, "closed")
  expect_equal(decision$manuscript_claims_state, "closed")
})
