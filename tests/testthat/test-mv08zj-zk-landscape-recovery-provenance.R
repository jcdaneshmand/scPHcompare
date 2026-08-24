test_that("MV8-ZJ prospectively binds truthful recovery-aware closure", {
  root <- testthat::test_path(
    "..", "..", "docs", "audits",
    "mv08zj-landscape-closure-recovery-prefreeze-v1"
  )
  read <- function(name) utils::read.csv(
    file.path(root, name), stringsAsFactors = FALSE, check.names = FALSE
  )
  validation <- read("mv08zj-validation.csv")
  contract <- read("mv08zj-contract.csv")
  decision <- read("mv08zj-decision.csv")
  implementation <- read("mv08zj-implementation-bindings.csv")
  manifest <- read("mv08zj-artifact-manifest.csv")

  expect_equal(nrow(validation), 20L)
  expect_true(all(validation$passed))
  expect_equal(contract$adopted_production_order, 164L)
  expect_false(contract$upper_bounds_are_measurements)
  expect_true(contract$original_MV8_ZG_runs_unchanged)
  expect_false(contract$companion_recomputes_landscapes)
  expect_true(decision$companion_MV8_ZK_required)
  expect_true(decision$order_164_upper_bounds_not_measurements)
  expect_false(decision$scientific_contract_changed)
  expect_false(decision$resource_contract_changed)
  expect_equal(decision$comparison_jobs_authorized, 0L)
  expect_equal(nrow(implementation), 3L)
  expect_true(all(file.exists(testthat::test_path("..", "..", implementation$file))))

  observed <- unname(vapply(file.path(root, manifest$artifact), function(path) {
    digest::digest(file = path, algo = "sha256", serialize = FALSE)
  }, character(1L)))
  expect_identical(observed, manifest$sha256)
})

test_that("MV8-ZK implementation is bounded to provenance-only closure", {
  root <- testthat::test_path("..", "..")
  builder <- file.path(root, "scripts",
                       "build_mv08zk_landscape_recovery_provenance_closure.R")
  prefreeze <- file.path(root, "scripts",
                         "build_mv08zj_landscape_closure_recovery_prefreeze.R")
  expect_silent(parse(builder))
  expect_silent(parse(prefreeze))
  text <- paste(readLines(builder, warn = FALSE), collapse = "\n")
  expect_match(text, "upper_bounds_are_measurements = FALSE", fixed = TRUE)
  expect_match(text, "measured_chunks = nrow(measured_ledger)", fixed = TRUE)
  expect_match(text, "landscape_recomputation_records = 0L", fixed = TRUE)
  expect_false(grepl("run_mv08z_landscape_chunk", text, fixed = TRUE))
  expect_false(grepl("unlink\\(|file.remove\\(", text, perl = TRUE))
})

test_that("MV8-ZK closure publishes canonical recovery resource semantics", {
  root <- testthat::test_path(
    "..", "..", "docs", "audits",
    "mv08zk-landscape-recovery-provenance-closure-v1"
  )
  skip_if_not(dir.exists(root), "MV8-ZK closure has not been produced")
  read <- function(name) utils::read.csv(
    file.path(root, name), stringsAsFactors = FALSE, check.names = FALSE
  )
  validation <- read("mv08zk-validation.csv")
  resource <- read("mv08zk-resource-interpretation.csv")
  decision <- read("mv08zk-decision.csv")
  manifest <- read("mv08zk-artifact-manifest.csv")

  expect_equal(nrow(validation), 25L)
  expect_true(all(validation$passed))
  expect_equal(resource$chunks, 628L)
  expect_equal(resource$pairs, 152744L)
  expect_equal(resource$measured_chunks, 627L)
  expect_equal(resource$upper_bound_chunks, 1L)
  expect_false(resource$upper_bounds_are_measurements)
  expect_true(resource$all_caps_passed_conservatively)
  expect_true(decision$full_landscape_closure_bound)
  expect_equal(decision$landscape_recomputation_records, 0L)
  expect_equal(decision$retries, 0L)

  observed <- unname(vapply(file.path(root, manifest$artifact), function(path) {
    digest::digest(file = path, algo = "sha256", serialize = FALSE)
  }, character(1L)))
  expect_identical(observed, manifest$sha256)
})
