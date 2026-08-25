test_that("MV8-ZZA freezes a closure-only serialization recovery", {
  root <- testthat::test_path(
    "..", "..", "docs", "audits",
    "mv08zza-closure-serialization-recovery-prefreeze-v1"
  )
  manifest <- read.csv(file.path(root, "mv08zza-artifact-manifest.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
  paths <- file.path(root, manifest$artifact)
  expect_true(all(file.exists(paths)))
  expect_equal(as.numeric(file.info(paths)$size), as.numeric(manifest$bytes))
  expect_equal(unname(vapply(paths, function(path) digest::digest(
    file = path, algo = "sha256", serialize = FALSE
  ), character(1L))), manifest$sha256)
  contract <- read.csv(file.path(root, "mv08zza-contract.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
  production <- read.csv(file.path(root, "mv08zza-immutable-production.csv"),
                         stringsAsFactors = FALSE, check.names = FALSE)
  validation <- read.csv(file.path(root, "mv08zza-validation.csv"),
                         stringsAsFactors = FALSE, check.names = FALSE)
  decision <- read.csv(file.path(root, "mv08zza-decision.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
  implementation <- read.csv(file.path(
    root, "mv08zza-implementation-bindings.csv"
  ), stringsAsFactors = FALSE, check.names = FALSE)
  expect_equal(nrow(validation), 10L)
  expect_true(all(validation$passed))
  expect_equal(production$completed_jobs, 40L)
  expect_equal(production$workers, 1L)
  expect_equal(production$retries, 0L)
  expect_false(production$production_mutation_allowed)
  expect_false(contract$scientific_values_affected)
  expect_false(contract$source_hashes_affected)
  expect_false(contract$production_artifacts_affected)
  expect_equal(contract$rerun_scope,
               "independent_closure_only_no_production_retry")
  expect_false(decision$production_retry_authorized)
  expect_true(decision$closure_rerun_authorized_after_commit)
  implementation_paths <- testthat::test_path("..", "..",
                                               implementation$file)
  expect_equal(unname(vapply(implementation_paths, function(path) digest::digest(
    file = path, algo = "sha256", serialize = FALSE
  ), character(1L))), implementation$sha256)
})
