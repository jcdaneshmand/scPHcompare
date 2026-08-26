test_that("MV10-M freezes a transparent value-aware disposition", {
  root <- testthat::test_path("..", "..")
  audit <- file.path(root, "docs", "audits",
                     "mv10m-clustering-disposition-prefreeze-v1")
  read_audit <- function(name) read.csv(
    file.path(audit, name), stringsAsFactors = FALSE, check.names = FALSE
  )
  manifest <- read_audit("mv10m-artifact-manifest.csv")
  files <- file.path(audit, manifest$artifact)
  expect_true(all(file.exists(files)))
  expect_equal(as.numeric(file.info(files)$size), as.numeric(manifest$bytes))
  expect_equal(unname(vapply(files, function(file) digest::digest(
    file = file, algo = "sha256", serialize = FALSE
  ), character(1L))), manifest$sha256)
  contract <- read_audit("mv10m-contract.csv")
  rules <- read_audit("mv10m-rule-contract.csv")
  outputs <- read_audit("mv10m-output-contract.csv")
  implementation <- read_audit("mv10m-implementation-bindings.csv")
  sources <- read_audit("mv10m-source-freeze.csv")
  validation <- read_audit("mv10m-validation.csv")
  decision <- read_audit("mv10m-decision.csv")

  expect_equal(contract$execution_head,
               "dfeb837c7bdd5f49419bcdce40de696931a3aef8")
  expect_true(contract$value_aware_prefreeze)
  expect_equal(contract$primary_method, "pam_dissimilarity_v1")
  expect_equal(contract$K_selection, "frozen_five_seed_one_se_rule")
  expect_equal(contract$silhouette_role, "descriptive_no_threshold")
  expect_equal(contract$representation_role, "sensitivity_no_ranking")
  expect_equal(contract$method_role, "sensitivity_no_ranking")
  expect_equal(contract$H0_H1, "separate")
  expect_equal(nrow(rules), 8L)
  expect_true(all(!rules$result_dependent_threshold))
  expect_equal(nrow(outputs), 5L)
  expect_identical(outputs$expected_rows, c(10L, 2L, 6L, 24L, 1L))
  expect_true(all(!outputs$labels_allowed))
  expect_true(all(!outputs$outcomes_allowed))
  expect_true(all(!outputs$ranking_allowed))
  expect_equal(nrow(implementation), 4L)
  expect_true(all(file.exists(file.path(root, implementation$file))))
  expect_true(all(grepl("^[0-9a-f]{64}$", implementation$sha256)))
  source_files <- ifelse(
    grepl("^(/|[A-Za-z]:[/\\\\])", sources$artifact),
    sources$artifact, file.path(root, sources$artifact)
  )
  expect_true(all(file.exists(source_files)))
  expect_equal(unname(vapply(source_files, function(file) digest::digest(
    file = file, algo = "sha256", serialize = FALSE
  ), character(1L))), sources$sha256)
  expect_equal(nrow(validation), 25L)
  expect_true(all(validation$passed))
  expect_true(decision$execution_authorized_after_commit)
  expect_false(decision$labels_authorized)
  expect_false(decision$outcomes_authorized)
  expect_false(decision$biological_interpretation_authorized)
  expect_false(decision$manuscript_claims_authorized)
})
