test_that("MV9-A prospectively freezes descriptive robustness synthesis", {
  root <- testthat::test_path("..", "..", "docs", "audits",
                              "mv09a-robustness-synthesis-prefreeze-v1")
  manifest <- read.csv(file.path(root, "mv09a-artifact-manifest.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
  paths <- file.path(root, manifest$artifact)
  expect_true(all(file.exists(paths)))
  expect_equal(as.numeric(file.info(paths)$size), as.numeric(manifest$bytes))
  expect_equal(unname(vapply(paths, function(path) digest::digest(
    file = path, algo = "sha256", serialize = FALSE
  ), character(1L))), manifest$sha256)
  contract <- read.csv(file.path(root, "mv09a-contract.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
  ranges <- read.csv(file.path(root, "mv09a-metric-ranges.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)
  outputs <- read.csv(file.path(root, "mv09a-output-contract.csv"),
                      stringsAsFactors = FALSE, check.names = FALSE)
  figures <- read.csv(file.path(root, "mv09a-figure-contract.csv"),
                      stringsAsFactors = FALSE, check.names = FALSE)
  validation <- read.csv(file.path(root, "mv09a-validation.csv"),
                         stringsAsFactors = FALSE, check.names = FALSE)
  implementation <- read.csv(file.path(root,
                                       "mv09a-implementation-bindings.csv"),
                             stringsAsFactors = FALSE, check.names = FALSE)
  expect_equal(contract$source_rows, 40L)
  expect_equal(contract$metrics, 11L)
  expect_equal(contract$internal_rows, 30L)
  expect_equal(contract$external_rows, 10L)
  expect_equal(contract$H0_rows, 20L)
  expect_equal(contract$H1_rows, 20L)
  expect_equal(contract$inference, "none")
  expect_equal(contract$combined_robustness_score, "forbidden")
  expect_equal(nrow(ranges), 11L)
  expect_true(all(ranges$finite))
  expect_equal(nrow(outputs), 6L)
  expect_identical(as.integer(outputs$rows),
                   c(40L, 440L, 66L, 110L, 220L, 88L))
  expect_false(any(outputs$labels_allowed))
  expect_false(any(outputs$outcomes_allowed))
  expect_false(any(outputs$inference_allowed))
  expect_equal(nrow(figures), 3L)
  expect_equal(nrow(validation), 18L)
  expect_true(all(validation$passed))
  implementation_paths <- testthat::test_path("..", "..",
                                               implementation$file)
  expect_equal(unname(vapply(implementation_paths, function(path) digest::digest(
    file = path, algo = "sha256", serialize = FALSE
  ), character(1L))), implementation$sha256)
})
