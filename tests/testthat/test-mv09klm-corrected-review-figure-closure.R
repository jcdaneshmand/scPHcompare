test_that("MV9-K/L/M closes corrected deterministic review figures", {
  root <- testthat::test_path("..", "..")
  figures <- file.path(root, "docs", "figures",
                       "mv09-corrected-review-v1")
  closure <- file.path(root, "docs", "audits",
                       "mv09l-corrected-review-figure-closure-v1")
  inspection <- file.path(root, "docs", "audits",
                          "mv09m-corrected-review-visual-inspection-v1")
  verify_manifest <- function(path, filename) {
    manifest <- read.csv(file.path(path, filename), stringsAsFactors = FALSE,
                         check.names = FALSE)
    paths <- file.path(path, manifest$artifact)
    expect_true(all(file.exists(paths)))
    expect_equal(as.numeric(file.info(paths)$size), as.numeric(manifest$bytes))
    expect_equal(unname(vapply(paths, function(file) digest::digest(
      file = file, algo = "sha256", serialize = FALSE
    ), character(1L))), manifest$sha256)
  }
  verify_manifest(figures, "mv09k-artifact-manifest.csv")
  verify_manifest(closure, "mv09l-artifact-manifest.csv")
  verify_manifest(inspection, "mv09m-artifact-manifest.csv")

  specs <- read.csv(file.path(figures, "mv09k-figure-specifications.csv"),
                    stringsAsFactors = FALSE, check.names = FALSE)
  metric_contract <- read.csv(file.path(figures, "mv09k-metric-contract.csv"),
                              stringsAsFactors = FALSE, check.names = FALSE)
  receipt <- read.csv(file.path(figures, "mv09k-terminal-receipt.csv"),
                      stringsAsFactors = FALSE, check.names = FALSE)
  images <- read.csv(file.path(closure, "mv09l-image-validation.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)
  closure_checks <- read.csv(file.path(closure, "mv09l-validation.csv"),
                             stringsAsFactors = FALSE, check.names = FALSE)
  bindings <- read.csv(file.path(inspection, "mv09m-figure-bindings.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
  visual_checks <- read.csv(file.path(
    inspection, "mv09m-visual-inspection.csv"
  ), stringsAsFactors = FALSE, check.names = FALSE)
  decision <- read.csv(file.path(inspection, "mv09m-decision.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
  builder <- read.csv(file.path(inspection, "mv09m-builder-binding.csv"),
                      stringsAsFactors = FALSE, check.names = FALSE)

  expect_equal(nrow(specs), 4L)
  expect_equal(specs$width_inches, c(12, 14, 14, 15))
  expect_equal(specs$height_inches, c(9, 7, 7.5, 9))
  expect_true(all(specs$dpi == 180L))
  expect_equal(receipt$completion_state, "complete")
  expect_equal(receipt$figures, 4L)
  expect_equal(receipt$internal_neighbor_k, 10L)
  expect_equal(receipt$external_neighbor_k, "2;3")
  expect_false(receipt$external_k7_displayed)
  expect_false(receipt$external_k7_interpretation)
  firewall <- c("labels_used", "outcomes_used", "inference_performed",
                "combined_score", "ranking_performed", "biological_claims",
                "manuscript_claims")
  expect_true(all(!unlist(receipt[firewall], use.names = FALSE)))
  expect_equal(receipt$human_review_state, "pending")
  expect_equal(metric_contract$disposition[which(metric_contract$k == 7L)],
               "structurally_noninformative_exclude")
  expect_identical(sort(metric_contract$k[
    metric_contract$disposition == "display_sensitivity"
  ]), c(2L, 3L))

  expect_equal(nrow(images), 4L)
  expect_true(all(images$dimensions_passed))
  expect_true(all(images$byte_identical_repeat))
  expect_true(all(images$exact_pixel_repeat))
  expect_true(all(images$maximum_pixel_difference == 0))
  expect_equal(nrow(closure_checks), 20L)
  expect_true(all(closure_checks$passed))
  expect_equal(nrow(bindings), 4L)
  expect_true(all(bindings$opened_at_original_resolution))
  expect_false(any(bindings$clipping_or_truncation_observed))
  expect_true(all(bindings$titles_axes_legible))
  expect_true(all(bindings$scientific_scope_disclosure_legible))
  expect_equal(nrow(visual_checks), 14L)
  expect_true(all(visual_checks$passed))
  expect_true(all(visual_checks$owner_review_state == "pending"))
  expect_equal(decision$computational_closure_checks, 20L)
  expect_equal(decision$visual_inspection_checks, 14L)
  expect_equal(decision$owner_scientific_review_state, "pending")
  expect_equal(decision$biological_interpretation_state, "closed")
  expect_equal(decision$manuscript_claim_state, "closed")
  expect_length(list.files(figures, pattern = "\\.pdf$", ignore.case = TRUE),
                0L)
  builder_path <- file.path(root, builder$file)
  expect_equal(digest::digest(file = builder_path, algo = "sha256",
                              serialize = FALSE), builder$sha256)
})
