test_that("MV9-DA freezes only the render-aesthetic recovery", {
  root <- testthat::test_path("..", "..", "docs", "audits",
                              "mv09da-render-aesthetic-recovery-prefreeze-v1")
  manifest <- read.csv(file.path(root, "mv09d-artifact-manifest.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
  paths <- file.path(root, manifest$artifact)
  expect_true(all(file.exists(paths)))
  expect_equal(as.numeric(file.info(paths)$size), as.numeric(manifest$bytes))
  expect_equal(unname(vapply(paths, function(path) digest::digest(
    file = path, algo = "sha256", serialize = FALSE
  ), character(1L))), manifest$sha256)
  recovery <- read.csv(file.path(root, "mv09da-recovery-contract.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
  validation <- read.csv(file.path(root, "mv09da-validation.csv"),
                         stringsAsFactors = FALSE, check.names = FALSE)
  decision <- read.csv(file.path(root, "mv09d-decision.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
  implementation <- read.csv(file.path(root,
                                       "mv09d-implementation-bindings.csv"),
                             stringsAsFactors = FALSE, check.names = FALSE)
  expect_equal(recovery$failed_artifacts, 1L)
  expect_gt(recovery$failed_artifact_bytes, 0)
  expect_false(recovery$source_data_affected)
  expect_false(recovery$metric_selection_affected)
  expect_false(recovery$scientific_contract_affected)
  expect_equal(recovery$rerun_scope,
               "clean_render_then_independent_repeat_only")
  expect_equal(nrow(validation), 10L)
  expect_true(all(validation$passed))
  expect_true(decision$render_authorized_after_commit)
  expect_false(decision$interpretation_authorized)
  implementation_paths <- testthat::test_path("..", "..",
                                               implementation$file)
  expect_equal(unname(vapply(implementation_paths, function(path) digest::digest(
    file = path, algo = "sha256", serialize = FALSE
  ), character(1L))), implementation$sha256)
})
