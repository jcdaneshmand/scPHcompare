test_that("MV10-V/W/X closes exact and visual outcome-review figures", {
  root <- testthat::test_path("..", "..")
  figures <- file.path(root, "docs", "figures", "mv10-outcome-review-v1")
  exact <- file.path(root, "docs", "audits",
                     "mv10w-outcome-review-figure-closure-v1")
  visual <- file.path(root, "docs", "audits",
                      "mv10x-outcome-review-visual-closure-v1")
  verify_manifest <- function(path, name) {
    manifest <- read.csv(file.path(path, name), stringsAsFactors = FALSE,
                         check.names = FALSE)
    files <- file.path(path, manifest$artifact)
    expect_true(all(file.exists(files)))
    expect_equal(as.numeric(file.info(files)$size), as.numeric(manifest$bytes))
    expect_equal(unname(vapply(files, function(file) digest::digest(
      file = file, algo = "sha256", serialize = FALSE
    ), character(1L))), manifest$sha256)
  }
  verify_manifest(figures, "mv10v-artifact-manifest.csv")
  verify_manifest(exact, "mv10w-artifact-manifest.csv")
  image_validation <- read.csv(
    file.path(exact, "mv10w-image-validation.csv"),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  validation <- read.csv(file.path(exact, "mv10w-validation.csv"),
                         stringsAsFactors = FALSE, check.names = FALSE)
  review <- read.csv(file.path(visual, "mv10x-visual-review.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)
  decision <- read.csv(file.path(visual, "mv10x-decision.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
  expect_equal(nrow(image_validation), 6L)
  expect_true(all(image_validation$byte_identical_repeat))
  expect_true(all(image_validation$maximum_pixel_difference == 0))
  expect_true(all(image_validation$exact_pixel_repeat))
  expect_true(all(image_validation$dimensions_passed))
  expect_equal(nrow(validation), 25L)
  expect_true(all(validation$passed))
  expect_equal(nrow(review), 15L)
  expect_true(all(review$passed))
  expect_equal(decision$figures, 6L)
  expect_true(decision$visual_review_passed)
  expect_true(decision$complete_reporting_condition_satisfied)
  expect_equal(decision$next_stage,
               "prospective_descriptive_outcome_disposition_prefreeze")
  expect_equal(decision$method_selection_state, "closed")
  expect_equal(decision$biological_interpretation_state, "closed")
  expect_equal(decision$manuscript_claims_state, "closed")
})
