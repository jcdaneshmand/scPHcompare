test_that("MV8-ZN accepts only the completed order-280 recovery state", {
  root <- testthat::test_path(
    "..", "..", "docs", "audits",
    "mv08zn-landscape-receipt-recovery-acceptance-v1"
  )
  read <- function(name) utils::read.csv(
    file.path(root, name), stringsAsFactors = FALSE, check.names = FALSE
  )
  validation <- read("mv08zn-validation.csv")
  stop_record <- read("mv08zn-validation-stop.csv")
  decision <- read("mv08zn-decision.csv")
  manifest <- read("mv08zn-artifact-manifest.csv")

  expect_equal(nrow(validation), 28L)
  expect_true(all(validation$passed))
  expect_equal(stop_record$executor_exit_status, 1L)
  expect_equal(stop_record$failure_class,
               "csv_roundtrip_exact_double_equality_only")
  expect_true(stop_record$completion_promotion_finished)
  expect_true(stop_record$progress_refresh_finished)
  expect_true(stop_record$tolerance_passed)
  expect_lt(stop_record$absolute_difference_seconds, 1e-9)
  expect_false(stop_record$scientific_failure)
  expect_false(stop_record$resource_failure)
  expect_false(stop_record$publication_failure)
  expect_true(decision$receipt_recovery_complete)
  expect_false(decision$executor_rerun_authorized)
  expect_equal(decision$resume_at_production_order, 281L)
  expect_equal(decision$landscape_recomputation_records, 0L)
  expect_equal(decision$retry_records, 0L)
  expect_true(decision$WSL_only_monitoring_required)

  observed <- unname(vapply(file.path(root, manifest$artifact), function(path) {
    digest::digest(file = path, algo = "sha256", serialize = FALSE)
  }, character(1L)))
  expect_identical(observed, manifest$sha256)
})

test_that("MV8-ZN acceptance builder is read-only with respect to production", {
  path <- testthat::test_path(
    "..", "..", "scripts",
    "build_mv08zn_landscape_receipt_recovery_acceptance.R"
  )
  expect_silent(parse(path))
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")
  expect_match(text, "aggregate_progress != aggregate_ledger", fixed = TRUE)
  expect_match(text, "abs(aggregate_difference) <= tolerance", fixed = TRUE)
  expect_match(text, "executor_rerun_authorized = FALSE", fixed = TRUE)
  expect_false(grepl("run_mv08z_landscape_chunk", text, fixed = TRUE))
  expect_false(grepl("file.copy\\(|unlink\\(|file.remove\\(", text, perl = TRUE))
  expect_false(grepl("atomic_csv\\((ledger|completed|progress)", text,
                     perl = TRUE))
})
