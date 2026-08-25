test_that("MV9-DB freezes only the specification storage-class recovery", {
  root <- testthat::test_path("..", "..", "docs", "audits",
                              "mv09db-specification-type-recovery-prefreeze-v1")
  manifest <- read.csv(file.path(root, "mv09d-artifact-manifest.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
  paths <- file.path(root, manifest$artifact)
  expect_true(all(file.exists(paths)))
  expect_equal(as.numeric(file.info(paths)$size), as.numeric(manifest$bytes))
  expect_equal(unname(vapply(paths, function(path) digest::digest(
    file = path, algo = "sha256", serialize = FALSE
  ), character(1L))), manifest$sha256)
  recovery <- read.csv(file.path(root, "mv09db-recovery-contract.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
  evidence <- read.csv(file.path(root, "mv09db-failed-repeat-evidence.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
  validation <- read.csv(file.path(root, "mv09db-validation.csv"),
                         stringsAsFactors = FALSE, check.names = FALSE)
  decision <- read.csv(file.path(root, "mv09d-decision.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
  implementation <- read.csv(file.path(root,
                                       "mv09d-implementation-bindings.csv"),
                             stringsAsFactors = FALSE, check.names = FALSE)
  expect_equal(recovery$failed_repeat_pngs, 3L)
  expect_true(recovery$repeat_pngs_byte_identical_to_primary)
  expect_false(recovery$source_data_affected)
  expect_false(recovery$metric_selection_affected)
  expect_false(recovery$figure_values_affected)
  expect_false(recovery$scientific_contract_affected)
  expect_match(recovery$rerun_scope, "independent_repeat_closure_only$")
  expect_equal(nrow(evidence), 3L)
  expect_true(all(evidence$byte_identical_to_primary))
  expect_equal(evidence$sha256, evidence$primary_sha256)
  expect_equal(nrow(validation), 12L)
  expect_true(all(validation$passed))
  expect_true(decision$render_authorized_after_commit)
  expect_false(decision$interpretation_authorized)
  expect_equal(nrow(implementation), 4L)
  expect_true(all(implementation$bytes > 0))
  expect_true(all(grepl("^[0-9a-f]{64}$", implementation$sha256)))
  expect_equal(anyDuplicated(implementation$file), 0L)
})
