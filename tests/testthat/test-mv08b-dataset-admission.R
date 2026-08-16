testthat::test_that("MV8-B freezes dataset selection before external outcomes", {
  spec_path <- testthat::test_path("..", "..", "docs", "specifications",
    "MV08B_READ_ONLY_DATASET_ADMISSION_AUDIT_PREFREEZE_V1.md")
  registry_path <- testthat::test_path("..", "..", "docs", "specifications",
    "mv08b-dataset-candidate-registry-v1.csv")
  spec <- paste(readLines(spec_path, warn = FALSE), collapse = "\n")
  registry <- utils::read.csv(registry_path, stringsAsFactors = FALSE,
                              check.names = FALSE)
  testthat::expect_match(spec, "technical_reprocessing_validation", fixed = TRUE)
  testthat::expect_match(spec, "same_tissue_external_replication", fixed = TRUE)
  testthat::expect_match(spec, "multi_tissue_external_generalization", fixed = TRUE)
  testthat::expect_equal(nrow(registry), 6L)
  testthat::expect_identical(registry$candidate_order, 1:6)
  testthat::expect_false(any(registry$selection_outcomes_opened))
  testthat::expect_false(any(registry$external_expression_download_authorized))
  testthat::expect_false(any(registry$new_ph_authorized))
})

testthat::test_that("MV8-B explicitly includes the separate all-bone-marrow cohort", {
  registry_path <- testthat::test_path("..", "..", "docs", "specifications",
    "mv08b-dataset-candidate-registry-v1.csv")
  registry <- utils::read.csv(registry_path, stringsAsFactors = FALSE,
                              check.names = FALSE)
  incumbent <- registry[registry$candidate_id == "existing_gse120221_25", ]
  testthat::expect_equal(nrow(incumbent), 1L)
  testthat::expect_identical(incumbent$repository_accession, "GSE120221")
  testthat::expect_identical(incumbent$planned_role,
    "technical_reprocessing_validation")
  testthat::expect_true(incumbent$whole_bone_marrow_required)
  testthat::expect_true(incumbent$cell_view_requested)
  testthat::expect_true(incumbent$gene_view_requested)
})

testthat::test_that("MV8-B classifies the incumbent and independent candidates honestly", {
  facts_path <- testthat::test_path("..", "..", "docs", "audits",
    "mv08b-dataset-source-facts-v1.csv")
  facts <- utils::read.csv(facts_path, stringsAsFactors = FALSE,
                           check.names = FALSE)
  incumbent <- facts[facts$candidate_id == "existing_gse120221_25", ]
  preferred <- facts[facts$candidate_id ==
    "hca_hematopoietic_immune_cell_atlas", ]
  testthat::expect_equal(nrow(facts), 6L)
  testthat::expect_false(anyDuplicated(facts$candidate_id) > 0L)
  testthat::expect_identical(incumbent$known_accession_overlap, 25L)
  testthat::expect_identical(incumbent$independent_donors, 0L)
  testthat::expect_identical(incumbent$disposition,
    "technical_reprocessing_only")
  testthat::expect_match(incumbent$cell_view_admission,
    "cell-view technical reprocessing", fixed = TRUE)
  testthat::expect_match(incumbent$gene_view_admission,
    "non-estimable under the frozen 500-gene", fixed = TRUE)
  testthat::expect_identical(preferred$disposition,
    "preferred_pending_download_authorization")
  testthat::expect_false(tolower(trimws(preferred$unresolved_fields)) == "none")
})

testthat::test_that("MV8-B preserves the corrected landscape and view contract", {
  spec_path <- testthat::test_path("..", "..", "docs", "specifications",
    "MV08B_READ_ONLY_DATASET_ADMISSION_AUDIT_PREFREEZE_V1.md")
  spec <- paste(readLines(spec_path, warn = FALSE), collapse = "\n")
  testthat::expect_match(spec, "every consecutive active landscape level",
    fixed = TRUE)
  testthat::expect_match(spec, "H0 and H1 calculated and reported separately",
    fixed = TRUE)
  testthat::expect_match(spec, "no universal grid", fixed = TRUE)
  testthat::expect_match(spec, "no universal landscape-level cap", fixed = TRUE)
  testthat::expect_match(spec, "Cell and gene topology are separate views",
    fixed = TRUE)
  testthat::expect_match(spec, "fixed 500-gene", fixed = TRUE)
})

testthat::test_that("MV8-B builder and validator stop before download and PH", {
  builder_path <- testthat::test_path("..", "..", "scripts",
    "build_mv08b_read_only_dataset_admission_audit.R")
  validator_path <- testthat::test_path("..", "..", "scripts",
    "validate_mv08b_read_only_dataset_admission_audit.R")
  builder <- paste(readLines(builder_path, warn = FALSE), collapse = "\n")
  validator <- paste(readLines(validator_path, warn = FALSE), collapse = "\n")
  testthat::expect_match(builder, "MV8-B prospective HEAD mismatch", fixed = TRUE)
  testthat::expect_match(builder, "external_expression_download_authorized = FALSE",
    fixed = TRUE)
  testthat::expect_match(builder, "new_ph_authorized = FALSE", fixed = TRUE)
  testthat::expect_match(builder, "all_consecutive_active_levels", fixed = TRUE)
  testthat::expect_match(validator, "exact_accession_overlap", fixed = TRUE)
  testthat::expect_match(validator, "local_source_compatibility", fixed = TRUE)
  testthat::expect_match(validator, "privacy_and_confidentiality", fixed = TRUE)
})
