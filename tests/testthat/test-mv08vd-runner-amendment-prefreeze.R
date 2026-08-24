test_that("MV8-VD binds the recovery runner and retained strict prefix", {
  root <- file.path("..", "..", "docs", "audits",
                    "mv08vd-runner-amendment-prefreeze-v1")
  read <- function(name) utils::read.csv(
    file.path(root, name), check.names = FALSE, stringsAsFactors = FALSE
  )
  binding <- read("mv08vd-implementation-bindings.csv")
  checks <- read("mv08vd-validation.csv")
  decision <- read("mv08vd-decision.csv")
  manifest <- read("mv08vd-artifact-manifest.csv")

  expect_equal(nrow(binding), 2L)
  expect_setequal(binding$file, c(
    "scripts/bootstrap_mv08va_full_ph_recovery.R",
    "scripts/run_mv08va_full_ph_production_recovery.R"
  ))
  expect_true(all(binding$mv08va_sha256 != binding$sha256))
  expect_true(all(checks$passed))
  expect_equal(nrow(checks), 16L)
  expect_true(decision$runner_resume_authorized)
  expect_equal(decision$accepted_completed_records, 1L)
  expect_equal(decision$remaining_ph_records_authorized, 1256L)
  expect_equal(decision$resume_at_production_order, 2L)
  expect_equal(decision$workers, 1L)
  expect_equal(decision$automatic_retries, 0L)
  expect_false(decision$scientific_contract_changed)
  expect_false(decision$resource_contract_changed)
  expect_equal(decision$landscape_groups_authorized, 0L)
  expect_equal(decision$outcome_jobs_authorized, 0L)

  observed <- vapply(binding$file, function(path) {
    digest::digest(file = file.path("..", "..", path), algo = "sha256",
                   serialize = FALSE)
  }, character(1L))
  expect_identical(unname(observed), binding$sha256)
  observed_manifest <- vapply(file.path(root, manifest$artifact), function(path) {
    digest::digest(file = path, algo = "sha256", serialize = FALSE)
  }, character(1L))
  expect_identical(unname(observed_manifest), manifest$sha256)
})
