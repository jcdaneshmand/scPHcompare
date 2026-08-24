test_that("MV8-VB preserves zero-output bootstrap stops and binds value equality", {
  root <- file.path("..", "..", "docs", "audits",
                    "mv08vb-bootstrap-type-recovery-prefreeze-v1")
  read <- function(name) utils::read.csv(
    file.path(root, name), check.names = FALSE, stringsAsFactors = FALSE
  )
  history <- read("mv08vb-bootstrap-stop-history.csv")
  implementation <- read("mv08vb-implementation-binding.csv")
  checks <- read("mv08vb-validation.csv")
  decision <- read("mv08vb-decision.csv")
  manifest <- read("mv08vb-artifact-manifest.csv")

  expect_equal(nrow(history), 4L)
  expect_false(any(history$output_root_created))
  expect_equal(sum(history$copied_records), 0L)
  expect_equal(sum(history$recomputed_records), 0L)
  expect_equal(sum(history$retry_records), 0L)
  expect_false(any(history$supplied_head_matches_actual[1:2]))
  expect_true(all(history$supplied_head_matches_actual[3:4]))
  expect_equal(nrow(implementation), 1L)
  expect_false(implementation$prior_sha256 == implementation$sha256)
  expect_identical(implementation$exact_change,
                   "identical_numeric_types_to_numeric_value_equality")
  expect_true(all(checks$passed))
  expect_equal(nrow(checks), 12L)
  expect_true(decision$corrected_bootstrap_authorized)
  expect_equal(decision$copied_records_authorized, 1L)
  expect_equal(decision$recomputed_records_authorized, 0L)
  expect_equal(decision$retry_records_authorized, 0L)
  expect_equal(decision$resume_at_production_order, 2L)
  expect_false(decision$scientific_contract_changed)
  expect_false(decision$resource_contract_changed)
  expect_equal(decision$landscape_groups_authorized, 0L)
  expect_equal(decision$outcome_jobs_authorized, 0L)

  amendment <- utils::read.csv(file.path(
    "..", "..", "docs", "audits",
    "mv08vc-bootstrap-amendment-recovery-prefreeze-v1",
    "mv08vc-implementation-binding.csv"
  ), check.names = FALSE, stringsAsFactors = FALSE)
  expect_identical(amendment$mv08vb_sha256, implementation$sha256)
  observed <- digest::digest(
    file = file.path("..", "..", amendment$file), algo = "sha256",
    serialize = FALSE
  )
  expect_identical(observed, amendment$sha256)
  observed_manifest <- vapply(file.path(root, manifest$artifact), function(path) {
    digest::digest(file = path, algo = "sha256", serialize = FALSE)
  }, character(1L))
  expect_identical(unname(observed_manifest), manifest$sha256)
})
