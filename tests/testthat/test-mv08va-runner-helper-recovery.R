test_that("MV8-VA preserves the helper omission and admits job 1 without retry", {
  root <- file.path("..", "..", "docs", "audits",
                    "mv08va-runner-helper-recovery-prefreeze-v1")
  read <- function(name) utils::read.csv(
    file.path(root, name), check.names = FALSE, stringsAsFactors = FALSE
  )
  stopped <- read("mv08va-stopped-evidence.csv")
  empty_launches <- read("mv08va-empty-launch-history.csv")
  failure <- read("mv08va-failure-receipt.csv")
  implementation <- read("mv08va-implementation-bindings.csv")
  checks <- read("mv08va-validation.csv")
  decision <- read("mv08va-decision.csv")
  manifest <- read("mv08va-artifact-manifest.csv")

  expect_equal(nrow(stopped), 6L)
  expect_false(any(stopped$private_absolute_path_published))
  expect_equal(nrow(empty_launches), 4L)
  expect_false(any(empty_launches$r_process_started |
                     empty_launches$output_root_created))
  expect_equal(sum(empty_launches$scientific_records), 0L)
  expect_identical(failure$failure_class, "missing_loaded_helper")
  expect_equal(failure$completed_child_attempts, 1L)
  expect_equal(failure$accepted_completed_records, 1L)
  expect_true(failure$child_stderr_empty)
  expect_true(failure$job1_independently_validated)
  expect_true(failure$h0_mst_passed)
  expect_equal(failure$later_ph_records, 0L)
  expect_equal(failure$landscape_records, 0L)
  expect_equal(failure$label_records, 0L)
  expect_false(failure$biological_outcomes_computed)
  expect_true(all(checks$passed))
  expect_equal(nrow(checks), 20L)
  expect_identical(decision$decision,
                   "authorize_no_retry_job1_bootstrap_and_resume_at_job2")
  expect_equal(decision$accepted_completed_records, 1L)
  expect_equal(decision$retry_records_authorized, 0L)
  expect_equal(decision$recomputed_records_authorized, 0L)
  expect_equal(decision$resume_at_production_order, 2L)
  expect_true(decision$fresh_private_public_roots_required)
  expect_true(decision$original_roots_immutable)
  expect_false(decision$scientific_contract_changed)
  expect_false(decision$resource_contract_changed)
  expect_equal(decision$landscape_groups_authorized, 0L)
  expect_equal(decision$outcome_jobs_authorized, 0L)

  observed_implementation <- vapply(implementation$file, function(path) {
    digest::digest(file = file.path("..", "..", path), algo = "sha256",
                   serialize = FALSE)
  }, character(1L))
  expect_identical(unname(observed_implementation), implementation$sha256)
  observed_manifest <- vapply(file.path(root, manifest$artifact), function(path) {
    digest::digest(file = path, algo = "sha256", serialize = FALSE)
  }, character(1L))
  expect_identical(unname(observed_manifest), manifest$sha256)
})

test_that("MV8-VA bootstrap and runner enforce no-retry fresh-root resume", {
  bootstrap <- paste(readLines(file.path(
    "..", "..", "scripts", "bootstrap_mv08va_full_ph_recovery.R"
  ), warn = FALSE), collapse = "\n")
  runner <- paste(readLines(file.path(
    "..", "..", "scripts", "run_mv08va_full_ph_production_recovery.R"
  ), warn = FALSE), collapse = "\n")
  expect_match(bootstrap, "file.copy", fixed = TRUE)
  expect_match(bootstrap, "recomputed_records = 0L", fixed = TRUE)
  expect_match(bootstrap, "retry_records = 0L", fixed = TRUE)
  expect_match(bootstrap, "resume_at_production_order = 2L", fixed = TRUE)
  expect_match(bootstrap,
               "mv08s_validate_ph_record_v1(record, row, view)", fixed = TRUE)
  expect_match(runner, 'source("R/mv08s_ph_sentinel.R")', fixed = TRUE)
  expect_match(runner, "MV08VA_RECOVERY_PREFREEZE", fixed = TRUE)
  expect_match(runner,
               "authorize_no_retry_job1_bootstrap_and_resume_at_job2",
               fixed = TRUE)
  expect_match(runner, "completed strict prefix", fixed = TRUE)
  expect_false(grepl("run_landscape", runner, fixed = TRUE))
})
