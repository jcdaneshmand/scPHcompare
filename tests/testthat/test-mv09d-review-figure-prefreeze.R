test_that("MV9-D prospectively freezes claim-free review figures", {
  root <- testthat::test_path("..", "..", "docs", "audits",
                              "mv09d-review-figure-prefreeze-v1")
  manifest <- read.csv(file.path(root, "mv09d-artifact-manifest.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
  paths <- file.path(root, manifest$artifact)
  expect_true(all(file.exists(paths)))
  expect_equal(as.numeric(file.info(paths)$size), as.numeric(manifest$bytes))
  expect_equal(unname(vapply(paths, function(path) digest::digest(
    file = path, algo = "sha256", serialize = FALSE
  ), character(1L))), manifest$sha256)
  contract <- read.csv(file.path(root, "mv09d-contract.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
  metrics <- read.csv(file.path(root, "mv09d-metric-selection.csv"),
                      stringsAsFactors = FALSE, check.names = FALSE)
  figures <- read.csv(file.path(root, "mv09d-figure-contract.csv"),
                      stringsAsFactors = FALSE, check.names = FALSE)
  review <- read.csv(file.path(root, "mv09d-review-contract.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)
  validation <- read.csv(file.path(root, "mv09d-validation.csv"),
                         stringsAsFactors = FALSE, check.names = FALSE)
  implementation <- read.csv(file.path(root,
                                       "mv09d-implementation-bindings.csv"),
                             stringsAsFactors = FALSE, check.names = FALSE)
  expect_equal(contract$figures, 3L)
  expect_equal(contract$selected_metrics, 4L)
  expect_equal(contract$metric_selection_timing, "prospective_before_render")
  expect_equal(contract$thresholds, "none")
  expect_equal(contract$rankings, "none")
  expect_equal(contract$inference, "none")
  expect_equal(contract$combined_score, "forbidden")
  expect_true(contract$human_review_required)
  expect_equal(nrow(metrics), 4L)
  expect_equal(nrow(figures), 3L)
  expect_true(all(figures$format == "PNG_only_no_PDF"))
  expect_true(all(figures$dpi == 180L))
  expect_equal(nrow(review), 6L)
  expect_true(all(review$owner_review_state == "pending"))
  expect_equal(nrow(validation), 18L)
  expect_true(all(validation$passed))
  implementation_paths <- testthat::test_path("..", "..",
                                               implementation$file)
  expect_equal(unname(vapply(implementation_paths, function(path) digest::digest(
    file = path, algo = "sha256", serialize = FALSE
  ), character(1L))), implementation$sha256)
})
