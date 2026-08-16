testthat::test_that("MV8-E identifies the exact historical HCA reference", {
  evidence <- testthat::test_path("..", "..", "docs", "audits",
    "mv08e-reference-reconciliation-evidence")
  fingerprint <- utils::read.csv(file.path(evidence,
    "mv08e-reference-fingerprint.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  testthat::expect_equal(fingerprint$unfiltered_ensembl93_genes, 58395L)
  testthat::expect_equal(fingerprint$filtered_reference_genes, 33538L)
  testthat::expect_equal(fingerprint$hca_h5_feature_genes, 33538L)
  testthat::expect_true(fingerprint$id_axis_byte_exact)
  testthat::expect_true(fingerprint$name_axis_byte_exact)
  testthat::expect_identical(fingerprint$filtered_gtf_id_axis_sha256,
    fingerprint$hca_h5_id_axis_sha256)
  testthat::expect_equal(fingerprint$panel_in_unfiltered_ensembl93, 500L)
  testthat::expect_equal(fingerprint$panel_in_filtered_cellranger_reference,
    475L)
  testthat::expect_equal(fingerprint$panel_excluded_by_reference_filter, 25L)
})

testthat::test_that("MV8-E freezes the ordered common 475 without leakage", {
  evidence <- testthat::test_path("..", "..", "docs", "audits",
    "mv08e-reference-reconciliation-evidence")
  panel <- utils::read.csv(file.path(evidence, "mv08e-common475-panel.csv"),
    stringsAsFactors = FALSE, check.names = FALSE)
  missing <- utils::read.csv(file.path(evidence,
    "mv08e-missing-biotypes.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  testthat::expect_equal(nrow(panel), 475L)
  testthat::expect_identical(panel$panel_order_475, 1:475)
  testthat::expect_true(all(diff(panel$panel_order_500) > 0L))
  testthat::expect_length(unique(panel$common475_axis_sha256), 1L)
  testthat::expect_false(any(panel$hca_expression_used_for_selection))
  testthat::expect_false(any(panel$biological_outcomes_used_for_selection))
  testthat::expect_equal(nrow(missing), 25L)
  testthat::expect_true(all(missing$current_ensembl_lookup_present))
  testthat::expect_false(any(missing$crosswalk_can_restore_matrix_count))
  testthat::expect_equal(sum(missing$ensembl93_gene_biotype ==
    "processed_transcript"), 4L)
  testthat::expect_equal(sum(grepl("pseudogene",
    missing$ensembl93_gene_biotype, fixed = TRUE)), 21L)
})

testthat::test_that("MV8-E preserves landscapes and prospects raw reads", {
  root <- testthat::test_path("..", "..")
  spec <- paste(readLines(file.path(root, "docs", "specifications",
    "MV08E_REFERENCE_RECONCILIATION_AND_PANEL_SENSITIVITY_PREFREEZE_V1.md"),
    warn = FALSE), collapse = "\n")
  evidence <- file.path(root, "docs", "audits",
    "mv08e-reference-reconciliation-evidence")
  fastq <- utils::read.csv(file.path(evidence,
    "mv08e-hca-fastq-resource-summary.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  decision <- utils::read.csv(file.path(evidence, "mv08e-decision.csv"),
    stringsAsFactors = FALSE, check.names = FALSE)
  testthat::expect_equal(fastq$biological_units, 8L)
  testthat::expect_equal(fastq$fastq_files, 48L)
  testthat::expect_equal(fastq$fastq_bytes, 85034239918)
  testthat::expect_identical(fastq$download_status, "not_downloaded")
  testthat::expect_false(decision$hca_fastq_download_authorized_now)
  testthat::expect_true(decision$raw_read_reprocessing_preference_recorded)
  testthat::expect_match(spec, "every consecutive active landscape level",
    fixed = TRUE)
  testthat::expect_match(spec, "exact or error-controlled squared-L2",
    fixed = TRUE)
  testthat::expect_match(spec,
    "H0 and H1 calculated and reported separately", fixed = TRUE)
  testthat::expect_match(spec, "no primary fixed grid", fixed = TRUE)
  testthat::expect_match(spec, "1,240 typed views", fixed = TRUE)
  testthat::expect_match(spec, "152,520 component rows", fixed = TRUE)
})

testthat::test_that("MV8-E public evidence is deterministic and validated", {
  evidence <- testthat::test_path("..", "..", "docs", "audits",
    "mv08e-reference-reconciliation-evidence")
  repeat_validation <- utils::read.csv(file.path(evidence,
    "mv08e-repeat-validation.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  validation <- utils::read.csv(file.path(evidence,
    "mv08e-independent-validation.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  manifest <- utils::read.csv(file.path(evidence,
    "mv08e-artifact-manifest.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  testthat::expect_equal(nrow(repeat_validation), 9L)
  testthat::expect_true(all(repeat_validation$byte_identical))
  testthat::expect_equal(nrow(validation), 12L)
  testthat::expect_true(all(validation$passed))
  testthat::expect_equal(nrow(manifest), 8L)
  testthat::expect_false(any(manifest$contains_expression))
  testthat::expect_false(any(manifest$contains_cell_barcode))
  testthat::expect_false(any(manifest$contains_local_absolute_path))
})
