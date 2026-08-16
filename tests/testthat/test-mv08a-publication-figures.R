testthat::test_that("MV8-A freezes publication figures before rendering", {
  spec_path <- testthat::test_path("..", "..", "docs", "specifications",
    "MV08A_CORRECTED_PUBLICATION_FIGURE_PREFREEZE_V1.md")
  renderer_path <- testthat::test_path("..", "..", "scripts",
    "render_mv08a_corrected_publication_figures.R")
  spec <- paste(readLines(spec_path, warn = FALSE), collapse = "\n")
  renderer <- paste(readLines(renderer_path, warn = FALSE), collapse = "\n")
  testthat::expect_match(spec, "all-active-level", fixed = TRUE)
  testthat::expect_match(spec, "H0 and H1 remain separate", fixed = TRUE)
  testthat::expect_match(spec, "Figure 8", fixed = TRUE)
  testthat::expect_match(renderer, "MV8-A prospective HEAD mismatch", fixed = TRUE)
  testthat::expect_match(renderer, "MV8-A H1 source hash mismatch", fixed = TRUE)
  testthat::expect_match(renderer, "MV8-A fixed PH artifact hash mismatch", fixed = TRUE)
})

testthat::test_that("MV8-A renders corrected landscape language", {
  path <- testthat::test_path("..", "..", "scripts",
    "render_mv08a_corrected_publication_figures.R")
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")
  testthat::expect_match(text, "all active levels", fixed = TRUE)
  testthat::expect_match(text, "No fixed scientific grid", fixed = TRUE)
  testthat::expect_match(text, "no level cap", fixed = TRUE)
  testthat::expect_match(text, "H0 and H1 remain separate", fixed = TRUE)
  testthat::expect_match(text, "Essential H0 class excluded", fixed = TRUE)
  testthat::expect_match(text, "visualization-only", fixed = TRUE)
})

testthat::test_that("MV8-A publishes complete figure families without selection", {
  path <- testthat::test_path("..", "..", "scripts",
    "render_mv08a_corrected_publication_figures.R")
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")
  testthat::expect_match(text, "All 7,626 median-across-seed sample pairs per view",
    fixed = TRUE)
  testthat::expect_match(text, "All 120 prespecified units", fixed = TRUE)
  testthat::expect_match(text, "No p-values or rankings", fixed = TRUE)
  testthat::expect_match(text, "favorable winner", fixed = TRUE)
  testthat::expect_match(text, "figure8_status = \"deferred_cross_stage_estimand\"",
    fixed = TRUE)
})

testthat::test_that("MV8-A requires independent structural and repeat validation", {
  path <- testthat::test_path("..", "..", "scripts",
    "validate_mv08a_corrected_publication_figures.R")
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")
  testthat::expect_match(text, "figure_hashes_dimensions_and_dpi", fixed = TRUE)
  testthat::expect_match(text, "corrected_landscape_contract_visible", fixed = TRUE)
  testthat::expect_match(text, "complete_unselected_families", fixed = TRUE)
  testthat::expect_match(text, "privacy_and_pdf_exclusion", fixed = TRUE)
  testthat::expect_match(text, "byte_identical_repeat", fixed = TRUE)
  testthat::expect_match(text, "authorize_author_review_of_figures", fixed = TRUE)
})
