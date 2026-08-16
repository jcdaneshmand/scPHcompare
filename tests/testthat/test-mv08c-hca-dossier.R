testthat::test_that("MV8-C freezes the exact HCA whole-marrow cohort", {
  spec_path <- testthat::test_path("..", "..", "docs", "specifications",
    "MV08C_HCA_ADMISSION_COMPUTE_DOSSIER_PREFREEZE_V1.md")
  registry_path <- testthat::test_path("..", "..", "docs", "specifications",
    "mv08c-hca-query-contract-v1.csv")
  spec <- paste(readLines(spec_path, warn = FALSE), collapse = "\n")
  registry <- utils::read.csv(registry_path, stringsAsFactors = FALSE,
                              check.names = FALSE)
  testthat::expect_equal(nrow(registry), 4L)
  testthat::expect_identical(registry$query_order, 1:4)
  testthat::expect_identical(registry$expected_entities, c(1L, 24L, 63L, 8L))
  testthat::expect_false(any(registry$expression_download_authorized))
  testthat::expect_false(any(registry$outcome_access_authorized))
  primary <- registry[registry$query_id == "primary_h5_manifest", ]
  testthat::expect_match(primary$filters_json,
    '"donorCount":{"is":[1]}', fixed = TRUE)
  testthat::expect_match(primary$filters_json,
    '"fileFormat":{"is":["h5"]}', fixed = TRUE)
  testthat::expect_match(spec, "eight one-donor controls", fixed = TRUE)
  testthat::expect_match(spec, "sixteen eight-donor pooled", fixed = TRUE)
  testthat::expect_match(spec, "63 `bone marrow hematopoietic cell`", fixed = TRUE)
})

testthat::test_that("MV8-C preserves structural and scientific gates", {
  spec_path <- testthat::test_path("..", "..", "docs", "specifications",
    "MV08C_HCA_ADMISSION_COMPUTE_DOSSIER_PREFREEZE_V1.md")
  spec <- paste(readLines(spec_path, warn = FALSE), collapse = "\n")
  testthat::expect_match(spec, "at least 384 usable cells", fixed = TRUE)
  testthat::expect_match(spec, "all 500 ordered MV7-FP genes", fixed = TRUE)
  testthat::expect_match(spec, "80 PH jobs total", fixed = TRUE)
  testthat::expect_match(spec, "19,840 query-to-reference", fixed = TRUE)
  testthat::expect_match(spec, "every consecutive active landscape level",
    fixed = TRUE)
  testthat::expect_match(spec, "H0 and H1 calculated and reported separately",
    fixed = TRUE)
  testthat::expect_match(spec, "no universal fixed grid", fixed = TRUE)
  testthat::expect_match(spec, "may not issue a GET or Range request",
    fixed = TRUE)
  testthat::expect_match(spec, "No current-sprint outcome authorizes an expression download",
    fixed = TRUE)
})

testthat::test_that("MV8-C acquisition is metadata-only by construction", {
  fetch_path <- testthat::test_path("..", "..", "scripts",
    "fetch_mv08c_hca_metadata.ps1")
  builder_path <- testthat::test_path("..", "..", "scripts",
    "build_mv08c_hca_admission_dossier.R")
  validator_path <- testthat::test_path("..", "..", "scripts",
    "validate_mv08c_hca_admission_dossier.R")
  fetch <- paste(readLines(fetch_path, warn = FALSE), collapse = "\n")
  builder <- paste(readLines(builder_path, warn = FALSE), collapse = "\n")
  validator <- paste(readLines(validator_path, warn = FALSE), collapse = "\n")
  testthat::expect_match(fetch, "compact TSV contains metadata, not expression data",
    fixed = TRUE)
  testthat::expect_false(grepl("repository/files", fetch, fixed = TRUE))
  testthat::expect_match(fetch, "expression_download_authorized", fixed = TRUE)
  testthat::expect_match(validator, "202770089", fixed = TRUE)
  testthat::expect_match(builder, "distance_rows <- 8L * 124L * 5L * 2L * 2L",
    fixed = TRUE)
  testthat::expect_match(builder, "all_consecutive_active_levels", fixed = TRUE)
  testthat::expect_match(validator, "production and repeat fifteen-file dossiers",
    fixed = TRUE)
})
