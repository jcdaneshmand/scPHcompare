test_that("MV8-VC binds the amendment chain without producing science", {
  root <- file.path("..", "..", "docs", "audits",
                    "mv08vc-bootstrap-amendment-recovery-prefreeze-v1")
  read <- function(name) utils::read.csv(
    file.path(root, name), check.names = FALSE, stringsAsFactors = FALSE
  )
  event <- read("mv08vc-bootstrap-stop-event.csv")
  binding <- read("mv08vc-implementation-binding.csv")
  checks <- read("mv08vc-validation.csv")
  decision <- read("mv08vc-decision.csv")
  manifest <- read("mv08vc-artifact-manifest.csv")

  expect_equal(nrow(event), 1L)
  expect_identical(event$raw_log_capture,
                   "interactive_tool_output_not_file_captured")
  expect_false(event$private_root_created)
  expect_false(event$public_root_created)
  expect_equal(event$copied_records, 0L)
  expect_equal(event$recomputed_records, 0L)
  expect_equal(event$retry_records, 0L)
  expect_equal(nrow(binding), 1L)
  expect_false(binding$mv08va_sha256 == binding$mv08vb_sha256)
  expect_false(binding$mv08vb_sha256 == binding$sha256)
  expect_identical(binding$exact_change,
                   "validate_hash_bound_bootstrap_amendment_chain")
  expect_true(all(checks$passed))
  expect_equal(nrow(checks), 14L)
  expect_true(decision$corrected_bootstrap_authorized)
  expect_equal(decision$copied_records_authorized, 1L)
  expect_equal(decision$recomputed_records_authorized, 0L)
  expect_equal(decision$retry_records_authorized, 0L)
  expect_equal(decision$resume_at_production_order, 2L)
  expect_false(decision$scientific_contract_changed)
  expect_false(decision$resource_contract_changed)
  expect_equal(decision$landscape_groups_authorized, 0L)
  expect_equal(decision$outcome_jobs_authorized, 0L)

  observed <- digest::digest(
    file = file.path("..", "..", binding$file), algo = "sha256",
    serialize = FALSE
  )
  expect_identical(observed, binding$sha256)
  observed_manifest <- vapply(file.path(root, manifest$artifact), function(path) {
    digest::digest(file = path, algo = "sha256", serialize = FALSE)
  }, character(1L))
  expect_identical(unname(observed_manifest), manifest$sha256)
})
