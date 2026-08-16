testthat::test_that("MV7-J freezes the corrected landscape and claim hierarchy", {
  spec_path <- testthat::test_path(
    "..", "..", "docs", "specifications",
    "MV07J_CLAIM_FIGURE_LITERATURE_MAP_PREFREEZE_V1.md")
  builder_path <- testthat::test_path(
    "..", "..", "scripts", "build_mv07j_claim_figure_literature_map.R")
  spec <- paste(readLines(spec_path, warn = FALSE), collapse = "\n")
  builder <- paste(readLines(builder_path, warn = FALSE), collapse = "\n")
  testthat::expect_match(spec, "all consecutive active landscape levels",
                         fixed = TRUE)
  testthat::expect_match(spec, "no universal grid count", fixed = TRUE)
  testthat::expect_match(spec, "H0 and H1 landscapes", fixed = TRUE)
  testthat::expect_match(spec, "first-level, 100-point-grid", fixed = TRUE)
  testthat::expect_match(builder, "supported_method", fixed = TRUE)
  testthat::expect_match(builder, "supported_descriptive", fixed = TRUE)
  testthat::expect_match(builder, "hypothesis_only", fixed = TRUE)
  testthat::expect_match(builder, "retire", fixed = TRUE)
})

testthat::test_that("MV7-J reconstructs complete result sensitivity families", {
  path <- testthat::test_path(
    "..", "..", "scripts", "build_mv07j_claim_figure_literature_map.R")
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")
  testthat::expect_match(text, "nrow(h1) != 15252L", fixed = TRUE)
  testthat::expect_match(text, "nrow(stability) != 54L", fixed = TRUE)
  testthat::expect_match(text, "nrow(selected) != 7440L", fixed = TRUE)
  testthat::expect_match(text, "nrow(outcome_units) != 120L", fixed = TRUE)
  testthat::expect_match(text, "fraction_gt_0_50", fixed = TRUE)
  testthat::expect_match(text, "h0_composite_concordance", fixed = TRUE)
  testthat::expect_match(text, "algorithm_sensitivity", fixed = TRUE)
  testthat::expect_match(text, "favorable_algorithm_selected = FALSE",
                         fixed = TRUE)
})

testthat::test_that("MV7-J maps figures and current literature without overclaiming", {
  path <- testthat::test_path(
    "..", "..", "scripts", "build_mv07j_claim_figure_literature_map.R")
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")
  testthat::expect_match(text, "narrows_dual_view_novelty", fixed = TRUE)
  testthat::expect_match(text, "narrows_topology_integration_novelty",
                         fixed = TRUE)
  testthat::expect_match(text, "narrows_robustness_claim", fixed = TRUE)
  testthat::expect_match(text, "ready_for_corrected_render", fixed = TRUE)
  testthat::expect_match(text, "confidential_source_text_included = FALSE",
                         fixed = TRUE)
  testthat::expect_match(text,
    "read_only_dataset_admission_audit", fixed = TRUE)
  testthat::expect_match(text, "new_data_download_authorized = FALSE",
                         fixed = TRUE)
  testthat::expect_match(text, "new_ph_authorized = FALSE", fixed = TRUE)
})

testthat::test_that("MV7-J requires independent reconstruction and repeat", {
  path <- testthat::test_path(
    "..", "..", "scripts", "validate_mv07j_claim_figure_literature_map.R")
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")
  testthat::expect_match(text, "independent_ari", fixed = TRUE)
  testthat::expect_match(text, "all 15,252 pair rows independently summarized",
                         fixed = TRUE)
  testthat::expect_match(text, "all 120 outcome units", fixed = TRUE)
  testthat::expect_match(text, "byte-identical", fixed = TRUE)
  testthat::expect_match(text, "no confidential text locator quote",
                         fixed = TRUE)
  testthat::expect_match(text, "17/17 checks pass", fixed = TRUE)
})
