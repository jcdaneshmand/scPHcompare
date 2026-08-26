test_that("MV13-D attempt-1 failure and attempt-2 recovery are exact", {
  root <- testthat::test_path("..", "..")
  failure <- file.path(root, "docs", "audits",
                       "mv13d-allqc-cell-full-failure-v1")
  recovery <- file.path(root, "docs", "audits",
                        "mv13d-allqc-cell-full-prefreeze-v2")
  read_at <- function(dir, name) read.csv(
    file.path(dir, name), stringsAsFactors = FALSE, check.names = FALSE
  )
  failure_validation <- read_at(failure, "mv13d-failure-validation.csv")
  failure_decision <- read_at(failure, "mv13d-failure-decision.csv")
  failure_manifest <- read_at(failure, "mv13d-failure-artifact-manifest.csv")
  contract <- read_at(recovery, "mv13d-contract.csv")
  binding <- read_at(recovery, "mv13d-failure-binding.csv")
  validation <- read_at(recovery, "mv13d-validation.csv")
  implementation <- read_at(recovery, "mv13d-implementation-bindings.csv")
  manifest <- read_at(recovery, "mv13d-artifact-manifest.csv")
  expect_equal(nrow(failure_validation), 12L)
  expect_true(all(failure_validation$passed))
  expect_equal(failure_decision$classification,
               "prewrite_vector_name_attribute_guard_failure")
  expect_false(failure_decision$scientific_failure)
  expect_false(failure_decision$resource_failure)
  expect_equal(failure_decision$completed_group_artifacts, 0)
  expect_equal(failure_decision$retries, 0)
  expect_equal(contract$execution_head,
               "f9e065c1c673db9145fa703373961c9acf01b385")
  expect_equal(contract$execution_attempt, 2)
  expect_equal(contract$recovery_change,
               "unit_model_axis_exact_values_ignore_vector_names_only")
  expect_equal(nrow(binding), 1L)
  expect_equal(nrow(validation), 28L)
  expect_true(all(validation$passed))
  paths <- file.path(root, implementation$file)
  expect_equal(as.numeric(file.info(paths)$size), implementation$bytes)
  expect_equal(unname(vapply(paths, function(path) digest::digest(
    file = path, algo = "sha256", serialize = FALSE
  ), character(1L))), implementation$sha256)
  for (pair in list(c(failure, failure_manifest$artifact),
                    c(recovery, manifest$artifact))) {
    base <- pair[[1L]]; names <- pair[-1L]
    files <- file.path(base, names)
    source_manifest <- if (identical(base, failure)) failure_manifest else manifest
    expect_equal(as.numeric(file.info(files)$size), source_manifest$bytes)
    expect_equal(unname(vapply(files, function(path) digest::digest(
      file = path, algo = "sha256", serialize = FALSE
    ), character(1L))), source_manifest$sha256)
  }
})
