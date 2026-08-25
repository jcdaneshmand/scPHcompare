test_that("MV9-E/F/G closes deterministic claim-free review figures", {
  root <- testthat::test_path("..", "..")
  figure_root <- file.path(root, "docs", "figures", "mv09-review-v1")
  closure_root <- file.path(root, "docs", "audits",
                            "mv09f-review-figure-closure-v1")
  inspection_root <- file.path(root, "docs", "audits",
                               "mv09g-review-figure-visual-inspection-v1")
  verify_manifest <- function(path, manifest_name) {
    manifest <- read.csv(file.path(path, manifest_name),
                         stringsAsFactors = FALSE, check.names = FALSE)
    files <- file.path(path, manifest$artifact)
    expect_true(all(file.exists(files)))
    expect_equal(as.numeric(file.info(files)$size), as.numeric(manifest$bytes))
    expect_equal(unname(vapply(files, function(file) digest::digest(
      file = file, algo = "sha256", serialize = FALSE
    ), character(1L))), manifest$sha256)
  }
  verify_manifest(figure_root, "mv09e-artifact-manifest.csv")
  verify_manifest(closure_root, "mv09f-artifact-manifest.csv")
  verify_manifest(inspection_root, "mv09g-artifact-manifest.csv")
  specs <- read.csv(file.path(figure_root, "mv09e-figure-specifications.csv"),
                    stringsAsFactors = FALSE, check.names = FALSE)
  receipt <- read.csv(file.path(figure_root, "mv09e-terminal-receipt.csv"),
                      stringsAsFactors = FALSE, check.names = FALSE)
  image_validation <- read.csv(file.path(
    closure_root, "mv09f-image-validation.csv"
  ), stringsAsFactors = FALSE, check.names = FALSE)
  closure_checks <- read.csv(file.path(closure_root, "mv09f-validation.csv"),
                             stringsAsFactors = FALSE, check.names = FALSE)
  bindings <- read.csv(file.path(inspection_root, "mv09g-figure-bindings.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
  inspection <- read.csv(file.path(
    inspection_root, "mv09g-visual-inspection.csv"
  ), stringsAsFactors = FALSE, check.names = FALSE)
  decision <- read.csv(file.path(inspection_root, "mv09g-decision.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
  builder <- read.csv(file.path(inspection_root, "mv09g-builder-binding.csv"),
                      stringsAsFactors = FALSE, check.names = FALSE)
  expect_equal(nrow(specs), 3L)
  expect_equal(specs$width_inches, c(12, 14, 15))
  expect_equal(specs$height_inches, c(9, 8.5, 11))
  expect_true(all(specs$dpi == 180L))
  expect_equal(receipt$completion_state, "complete")
  expect_equal(receipt$figures, 3L)
  expect_equal(receipt$selected_metrics, 4L)
  firewall <- c("labels_used", "outcomes_used", "inference_performed",
                "combined_score", "ranking_performed", "biological_claims",
                "manuscript_claims")
  expect_true(all(!unlist(receipt[firewall], use.names = FALSE)))
  expect_equal(receipt$human_review_state, "pending")
  expect_equal(nrow(image_validation), 3L)
  expect_true(all(image_validation$dimensions_passed))
  expect_true(all(image_validation$byte_identical_repeat))
  expect_true(all(image_validation$exact_pixel_repeat))
  expect_true(all(image_validation$maximum_pixel_difference == 0))
  expect_equal(nrow(closure_checks), 16L)
  expect_true(all(closure_checks$passed))
  expect_equal(nrow(bindings), 3L)
  expect_true(all(bindings$opened_at_original_resolution))
  expect_true(all(!bindings$clipping_or_truncation_observed))
  expect_true(all(bindings$titles_axes_legible))
  expect_true(all(bindings$replication_disclosure_legible))
  bound_paths <- file.path(figure_root, bindings$filename)
  expect_equal(unname(vapply(bound_paths, function(file) digest::digest(
    file = file, algo = "sha256", serialize = FALSE
  ), character(1L))), bindings$sha256)
  expect_equal(nrow(inspection), 12L)
  expect_true(all(inspection$passed))
  expect_true(all(inspection$owner_review_state == "pending"))
  expect_equal(decision$owner_scientific_review_state, "pending")
  expect_equal(decision$biological_interpretation_state, "closed")
  expect_equal(decision$manuscript_claim_state, "closed")
  expect_length(list.files(figure_root, pattern = "\\.pdf$",
                           ignore.case = TRUE), 0L)
  builder_path <- file.path(root, builder$file)
  expect_equal(digest::digest(file = builder_path, algo = "sha256",
                              serialize = FALSE), builder$sha256)
})
