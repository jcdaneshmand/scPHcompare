test_that("MV8-ZB preserves zero-output helper failure and binds fresh recovery", {
  root <- file.path("..", "..", "docs", "audits",
                    "mv08zb-landscape-helper-recovery-prefreeze-v1")
  read <- function(name) read.csv(file.path(root, name), check.names = FALSE,
                                  stringsAsFactors = FALSE)
  binding <- read("mv08zb-implementation-bindings.csv")
  failure <- read("mv08zb-failure.csv")
  checks <- read("mv08zb-validation.csv")
  decision <- read("mv08zb-decision.csv")
  manifest <- read("mv08zb-artifact-manifest.csv")
  expect_equal(nrow(binding), 4L)
  expect_true(all(!binding$scientific_change))
  expect_equal(failure$failure_class,
               "missing_finite_landscape_diagram_helper")
  expect_equal(failure$landscape_pair_outputs, 0L)
  expect_equal(failure$later_children_started, 0L)
  expect_equal(nrow(checks), 14L)
  expect_true(all(checks$passed))
  expect_equal(decision$replacement_children_authorized, 3L)
  expect_equal(decision$automatic_retries, 0L)
  expect_false(decision$scientific_contract_changed)
  expect_equal(decision$production_pairs_authorized, 0L)
  expect_equal(decision$downstream_jobs_authorized, 0L)
  observed <- unname(vapply(file.path(root, manifest$artifact), function(path) {
    digest::digest(file = path, algo = "sha256", serialize = FALSE)
  }, character(1L)))
  expect_identical(observed, manifest$sha256)
})
