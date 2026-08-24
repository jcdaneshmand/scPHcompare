test_that("MV8-ZL freezes exact no-recompute order-280 receipt promotion", {
  root <- testthat::test_path(
    "..", "..", "docs", "audits",
    "mv08zl-landscape-receipt-recovery-prefreeze-v1"
  )
  read <- function(name) utils::read.csv(
    file.path(root, name), stringsAsFactors = FALSE, check.names = FALSE
  )
  validation <- read("mv08zl-validation.csv")
  snapshot <- read("mv08zl-stopped-snapshot.csv")
  binding <- read("mv08zl-order-280-binding.csv")
  policy <- read("mv08zl-recovery-policy.csv")
  decision <- read("mv08zl-decision.csv")
  manifest <- read("mv08zl-artifact-manifest.csv")

  expect_equal(nrow(validation), 30L)
  expect_true(all(validation$passed))
  expect_equal(snapshot$ledger_rows, 280L)
  expect_equal(snapshot$completion_rows, 279L)
  expect_equal(snapshot$completion_partial_rows, 280L)
  expect_equal(snapshot$public_partial_files, 1L)
  expect_equal(snapshot$private_partial_files, 0L)
  expect_equal(snapshot$active_runner_processes, 0L)
  expect_equal(binding$production_order, 280L)
  expect_true(binding$telemetry_is_measured)
  expect_false(policy$landscape_recomputation)
  expect_false(policy$automatic_retry)
  expect_false(policy$ledger_rewrite)
  expect_true(policy$preserve_prefix_copy_privately)
  expect_true(policy$WSL_only_publication_and_monitoring)
  expect_true(decision$receipt_promotion_authorized)
  expect_equal(decision$resume_at_production_order, 281L)
  expect_true(decision$companion_MV8_ZM_required)

  observed <- unname(vapply(file.path(root, manifest$artifact), function(path) {
    digest::digest(file = path, algo = "sha256", serialize = FALSE)
  }, character(1L)))
  expect_identical(observed, manifest$sha256)
})

test_that("MV8-ZL recovery implementation cannot compute landscapes or delete", {
  root <- testthat::test_path("..", "..")
  executor <- file.path(root, "scripts",
                        "recover_mv08zl_landscape_receipt_publication.R")
  builder <- file.path(root, "scripts",
                       "build_mv08zl_landscape_receipt_recovery_prefreeze.R")
  closure <- file.path(root, "scripts",
                       "build_mv08zm_landscape_receipt_recovery_closure.R")
  expect_silent(parse(executor))
  expect_silent(parse(builder))
  expect_silent(parse(closure))
  text <- paste(readLines(executor, warn = FALSE), collapse = "\n")
  expect_match(text, "resume_at=281", fixed = TRUE)
  expect_match(text, "preserved prefix drift", fixed = TRUE)
  expect_false(grepl("run_mv08z_landscape_chunk", text, fixed = TRUE))
  expect_false(grepl("unlink\\(|file.remove\\(", text, perl = TRUE))
})

test_that("MV8-ZM binds receipt recovery to full closure", {
  root <- testthat::test_path(
    "..", "..", "docs", "audits",
    "mv08zm-landscape-receipt-recovery-closure-v1"
  )
  skip_if_not(dir.exists(root), "MV8-ZM closure has not been produced")
  read <- function(name) utils::read.csv(
    file.path(root, name), stringsAsFactors = FALSE, check.names = FALSE
  )
  validation <- read("mv08zm-validation.csv")
  summary <- read("mv08zm-recovery-summary.csv")
  decision <- read("mv08zm-decision.csv")
  manifest <- read("mv08zm-artifact-manifest.csv")
  expect_equal(nrow(validation), 25L)
  expect_true(all(validation$passed))
  expect_equal(summary$recovery_order, 280L)
  expect_equal(summary$preserved_prefix_rows, 279L)
  expect_equal(summary$final_completion_rows, 628L)
  expect_false(summary$receipt_reconstructed)
  expect_false(summary$ledger_rewritten)
  expect_equal(summary$landscape_recomputations, 0L)
  expect_equal(summary$retries, 0L)
  expect_true(decision$full_landscape_closure_bound)
  expect_equal(decision$landscape_pairs, 152744L)

  observed <- unname(vapply(file.path(root, manifest$artifact), function(path) {
    digest::digest(file = path, algo = "sha256", serialize = FALSE)
  }, character(1L)))
  expect_identical(observed, manifest$sha256)
})
