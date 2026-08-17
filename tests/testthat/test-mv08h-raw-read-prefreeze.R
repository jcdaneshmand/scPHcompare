testthat::test_that("MV8-H freezes the exact authorized FASTQ workload", {
  evidence <- testthat::test_path("..", "..", "docs", "audits",
    "mv08h-raw-read-prefreeze-v1")
  manifest <- utils::read.csv(file.path(evidence, "mv08h-fastq-manifest.csv"),
    stringsAsFactors = FALSE, check.names = FALSE)
  units <- utils::read.csv(file.path(evidence, "mv08h-unit-manifest.csv"),
    stringsAsFactors = FALSE, check.names = FALSE)
  decision <- utils::read.csv(file.path(evidence, "mv08h-decision.csv"),
    stringsAsFactors = FALSE, check.names = FALSE)
  testthat::expect_equal(nrow(manifest), 48L)
  testthat::expect_identical(manifest$file_order, 1:48)
  testthat::expect_equal(sum(manifest$file_size_bytes), 85034239918)
  testthat::expect_equal(length(unique(manifest$file_uuid)), 48L)
  testthat::expect_true(all(manifest$download_authorized))
  testthat::expect_equal(nrow(units), 8L)
  testthat::expect_true(all(units$fastq_files == 6L))
  testthat::expect_equal(sum(units$fastq_bytes), 85034239918)
  testthat::expect_true(decision$download_authorized)
  testthat::expect_false(decision$raw_reprocessing_authorized)
  testthat::expect_false(decision$label_access_authorized)
})

testthat::test_that("MV8-H preserves the exact-500 and landscape definitions", {
  root <- testthat::test_path("..", "..")
  spec <- paste(readLines(file.path(root, "docs", "specifications",
    "MV08H_EXACT500_HCA_RAW_READ_ACQUISITION_PREFREEZE_V1.md"), warn = FALSE),
    collapse = "\n")
  processing <- utils::read.csv(file.path(root, "docs", "audits",
    "mv08h-raw-read-prefreeze-v1", "mv08h-processing-contract.csv"),
    stringsAsFactors = FALSE, check.names = FALSE)
  text <- paste(processing$frozen_value, collapse = "\n")
  for (term in c("33,563-gene GTF", "all 500 ordered target stable IDs",
                 "Cell Ranger 3.0.0", "SC3Pv2", "exactly 384 cells",
                 "every consecutive active level", "no fixed grid",
                 "no level cap", "separate H0 and H1",
                 "exclusion of essential H0")) {
    testthat::expect_match(spec, term, fixed = TRUE)
  }
  for (term in c("every active landscape level",
                 "no fixed grid or level cap", "H0/H1 separate",
                 "essential H0 excluded")) {
    testthat::expect_match(text, term, fixed = TRUE)
  }
  testthat::expect_false(any(processing$execution_authorized))
})

testthat::test_that("MV8-H acquisition is bounded, resumable, and fail closed", {
  root <- testthat::test_path("..", "..")
  downloader <- paste(readLines(file.path(root, "scripts",
    "fetch_mv08h_hca_fastq.py"), warn = FALSE), collapse = "\n")
  for (term in c("EXPECTED_MANIFEST_SHA256", "MAX_ATTEMPTS = 5",
                 "MINIMUM_POST_DOWNLOAD_FREE_BYTES", "Range",
                 "Content-Range", "os.replace(partial, target)",
                 "partial retained", "--dry-run", "--max-files")) {
    testthat::expect_match(downloader, term, fixed = TRUE)
  }
  testthat::expect_match(downloader,
    "MV8-H exact 48-file HCA FASTQ download authorized", fixed = TRUE)
})

testthat::test_that("MV8-H public prefreeze evidence passes independently", {
  evidence <- testthat::test_path("..", "..", "docs", "audits",
    "mv08h-raw-read-prefreeze-v1")
  validation <- utils::read.csv(file.path(evidence,
    "mv08h-independent-validation.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  decision <- utils::read.csv(file.path(evidence,
    "mv08h-validation-decision.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  manifest <- utils::read.csv(file.path(evidence,
    "mv08h-artifact-manifest.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  repeat_validation <- utils::read.csv(file.path(evidence,
    "mv08h-repeat-validation.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  testthat::expect_equal(nrow(validation), 13L)
  testthat::expect_true(all(validation$passed))
  testthat::expect_equal(decision$checks_passed, 13L)
  testthat::expect_true(decision$fastq_download_authorized)
  testthat::expect_false(decision$raw_reprocessing_authorized)
  testthat::expect_false(any(manifest$contains_expression))
  testthat::expect_false(any(manifest$contains_cell_barcode))
  testthat::expect_false(any(manifest$contains_absolute_private_path))
  testthat::expect_false(any(manifest$contains_donor_attribute))
  testthat::expect_false(any(manifest$contains_outcome_label))
  testthat::expect_equal(nrow(repeat_validation), 9L)
  testthat::expect_true(all(repeat_validation$byte_identical))
  paths <- file.path(evidence, manifest$file)
  observed <- unname(vapply(paths, function(path) digest::digest(
    file = path, algo = "sha256", serialize = FALSE), character(1L)))
  testthat::expect_identical(observed, manifest$sha256)
})
