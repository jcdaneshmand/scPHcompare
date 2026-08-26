test_that("MV10-J/K/L closes exact and visual clustering-review figures", {
  root <- testthat::test_path("..", "..")
  figures <- file.path(root, "docs", "figures",
                       "mv10-clustering-review-v1")
  exact <- file.path(root, "docs", "audits",
                     "mv10k-clustering-review-figure-closure-v1")
  visual <- file.path(root, "docs", "audits",
                      "mv10l-clustering-review-visual-closure-v1")
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
  verify_manifest(figures, "mv10j-artifact-manifest.csv")
  verify_manifest(exact, "mv10k-artifact-manifest.csv")
  image_validation <- read.csv(
    file.path(exact, "mv10k-image-validation.csv"),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  validation <- read.csv(file.path(exact, "mv10k-validation.csv"),
                         stringsAsFactors = FALSE, check.names = FALSE)
  review <- read.csv(file.path(visual, "mv10l-visual-review.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)
  decision <- read.csv(file.path(visual, "mv10l-decision.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)

  expect_equal(nrow(image_validation), 7L)
  expect_true(all(image_validation$byte_identical_repeat))
  expect_true(all(image_validation$maximum_pixel_difference == 0))
  expect_true(all(image_validation$exact_pixel_repeat))
  expect_true(all(image_validation$dimensions_passed))
  expect_equal(nrow(validation), 27L)
  expect_true(all(validation$passed))
  expect_equal(nrow(review), 15L)
  expect_true(all(review$passed))
  expect_equal(decision$figures, 7L)
  expect_true(decision$visual_review_passed)
  expect_true(decision$owner_comprehensiveness_condition_satisfied)
  expect_equal(decision$next_stage,
               "prospective_label_closed_methodological_disposition_prefreeze")
  expect_equal(decision$labels_state, "closed")
  expect_equal(decision$outcomes_state, "closed")
  expect_equal(decision$biological_interpretation_state, "closed")
  expect_equal(decision$manuscript_claims_state, "closed")
})
