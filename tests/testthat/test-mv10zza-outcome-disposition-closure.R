test_that("MV10-Z/ZA closes transparent descriptive disposition", {
  root <- testthat::test_path("..", "..")
  production <- file.path(root, "docs", "audits",
                          "mv10z-outcome-disposition-v1")
  closure <- file.path(root, "docs", "audits",
                       "mv10za-outcome-disposition-closure-v1")
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
  verify_manifest(production, "mv10z-artifact-manifest.csv")
  verify_manifest(closure, "mv10za-artifact-manifest.csv")
  primary <- read.csv(file.path(
    production, "mv10z-primary-representation-envelope.csv"
  ), stringsAsFactors = FALSE, check.names = FALSE)
  methods <- read.csv(file.path(production, "mv10z-method-envelope.csv"),
                      stringsAsFactors = FALSE, check.names = FALSE)
  context <- read.csv(file.path(production, "mv10z-context-contrast.csv"),
                      stringsAsFactors = FALSE, check.names = FALSE)
  approach <- read.csv(file.path(production, "mv10z-approach-contrast.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
  disposition <- read.csv(file.path(production, "mv10z-disposition.csv"),
                          stringsAsFactors = FALSE, check.names = FALSE)
  receipt <- read.csv(file.path(production, "mv10z-terminal-receipt.csv"),
                      stringsAsFactors = FALSE, check.names = FALSE)
  rehash <- read.csv(file.path(closure, "mv10za-rehash.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)
  validation <- read.csv(file.path(closure, "mv10za-validation.csv"),
                         stringsAsFactors = FALSE, check.names = FALSE)
  decision <- read.csv(file.path(closure, "mv10za-decision.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
  expect_equal(receipt$completion_state, "complete")
  expect_equal(receipt$output_tables, 5L)
  expect_equal(nrow(primary), 20L)
  expect_equal(nrow(methods), 60L)
  expect_equal(nrow(context), 120L)
  expect_equal(nrow(approach), 120L)
  expect_equal(nrow(disposition), 1L)
  expect_true(all(!context$inference_performed))
  expect_true(all(!context$causal_interpretation_allowed))
  expect_true(all(!approach$inference_performed))
  expect_true(all(!approach$causal_interpretation_allowed))
  expect_false(receipt$magnitude_threshold_applied)
  expect_false(receipt$p_values_computed)
  expect_false(receipt$method_selection_executed)
  expect_false(receipt$ranking_executed)
  expect_false(receipt$H0_H1_combined)
  expect_false(receipt$biological_claims)
  expect_false(receipt$manuscript_claims)
  expect_equal(nrow(rehash), 5L)
  expect_true(all(rehash$exact_columns))
  expect_true(all(rehash$numeric_repeat))
  expect_true(all(rehash$byte_identical))
  expect_equal(nrow(validation), 26L)
  expect_true(all(validation$passed))
  expect_equal(decision$frozen_method_state,
               "retain_PAM_without_outcome_tuning")
  expect_equal(decision$next_stage,
               "cross_view_descriptive_synthesis_prefreeze_or_stop")
  expect_equal(decision$method_selection_state, "closed")
  expect_equal(decision$biological_claims_state, "closed")
  expect_equal(decision$manuscript_claims_state, "closed")
})
