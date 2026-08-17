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

testthat::test_that("MV8-H sentinel independently opens only the remaining download", {
  root <- testthat::test_path("..", "..")
  evidence <- file.path(root, "docs", "audits", "mv08h-download-sentinel-v1")
  files <- utils::read.csv(file.path(evidence,
    "mv08h-sentinel-file-validation.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  summary <- utils::read.csv(file.path(evidence,
    "mv08h-sentinel-summary.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  decision <- utils::read.csv(file.path(evidence,
    "mv08h-sentinel-decision.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  artifacts <- utils::read.csv(file.path(evidence,
    "mv08h-sentinel-artifact-manifest.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  testthat::expect_equal(nrow(files), 1L)
  testthat::expect_equal(files$observed_bytes, 394373114)
  testthat::expect_identical(files$observed_sha256,
    "4464d4ea5000356ec74dd80f37dc69016ca13668bc8f33c8d2fcffdef31dd302")
  testthat::expect_true(files$gzip_magic)
  testthat::expect_true(files$receipt_exact)
  testthat::expect_true(files$passed)
  testthat::expect_equal(summary$verified_files, 1L)
  testthat::expect_equal(summary$partial_files, 0L)
  testthat::expect_true(summary$storage_gate_passed)
  testthat::expect_true(decision$remaining_fastq_download_authorized)
  testthat::expect_false(decision$mkref_authorized)
  testthat::expect_false(decision$raw_reprocessing_authorized)
  testthat::expect_false(decision$label_access_authorized)
  paths <- file.path(evidence, artifacts$file)
  observed <- unname(vapply(paths, function(path) digest::digest(
    file = path, algo = "sha256", serialize = FALSE), character(1L)))
  testthat::expect_identical(observed, artifacts$sha256)
  validator <- paste(readLines(file.path(root, "scripts",
    "validate_mv08h_fastq_cache.py"), warn = FALSE), collapse = "\n")
  for (term in c("EXPECTED_MANIFEST_SHA256", "gzip_magic",
                 "receipt_exact", "RESERVE_BYTES", "CACHE_CAP_BYTES",
                 "complete_download_exact_authorize_reference_input_validation_only")) {
    testthat::expect_match(validator, term, fixed = TRUE)
  }
})

testthat::test_that("MV8-H complete acquisition opens only reference input validation", {
  root <- testthat::test_path("..", "..")
  evidence <- file.path(root, "docs", "audits", "mv08h-download-closure-v1")
  files <- utils::read.csv(file.path(evidence,
    "mv08h-complete-file-validation.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  summary <- utils::read.csv(file.path(evidence,
    "mv08h-complete-summary.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  decision <- utils::read.csv(file.path(evidence,
    "mv08h-complete-decision.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  artifacts <- utils::read.csv(file.path(evidence,
    "mv08h-complete-artifact-manifest.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  testthat::expect_equal(nrow(files), 48L)
  testthat::expect_identical(files$file_order, 1:48)
  testthat::expect_true(all(files$gzip_magic))
  testthat::expect_true(all(files$receipt_exact))
  testthat::expect_true(all(files$passed))
  testthat::expect_equal(summary$verified_files, 48L)
  testthat::expect_equal(summary$expected_cache_bytes, 85034239918)
  testthat::expect_equal(summary$observed_cache_bytes, 85034239918)
  testthat::expect_equal(summary$remaining_manifest_bytes, 0)
  testthat::expect_equal(summary$partial_files, 0L)
  testthat::expect_true(summary$storage_gate_passed)
  testthat::expect_true(summary$all_file_gates_passed)
  testthat::expect_false(decision$remaining_fastq_download_authorized)
  testthat::expect_true(decision$reference_input_validation_authorized)
  testthat::expect_false(decision$mkref_authorized)
  testthat::expect_false(decision$raw_reprocessing_authorized)
  testthat::expect_false(decision$label_access_authorized)
  testthat::expect_false(decision$biological_outcomes_authorized)
  testthat::expect_true("MV08H_FASTQ_COMPLETE_CLOSURE_2026-08-17.md" %in%
    artifacts$file)
  testthat::expect_false(any(artifacts$contains_expression))
  testthat::expect_false(any(artifacts$contains_cell_barcode))
  testthat::expect_false(any(artifacts$contains_absolute_private_path))
  testthat::expect_false(any(artifacts$contains_donor_attribute))
  testthat::expect_false(any(artifacts$contains_outcome_label))
  paths <- file.path(evidence, artifacts$file)
  observed <- unname(vapply(paths, function(path) digest::digest(
    file = path, algo = "sha256", serialize = FALSE), character(1L)))
  testthat::expect_identical(observed, artifacts$sha256)
})

testthat::test_that("MV8-H closes the exact Ensembl-93 input before Cell Ranger", {
  root <- testthat::test_path("..", "..")
  evidence <- file.path(root, "docs", "audits", "mv08h-reference-input-v1")
  identity <- utils::read.csv(file.path(evidence,
    "mv08h-reference-input-identity.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  decision <- utils::read.csv(file.path(evidence,
    "mv08h-reference-input-decision.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  artifacts <- utils::read.csv(file.path(evidence,
    "mv08h-reference-input-artifact-manifest.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  testthat::expect_equal(identity$observed_bytes, 881214682)
  testthat::expect_identical(identity$observed_sha256,
    "2a27436d44f0d6350f86894fbe5edec56faa5467028879784508041562406aa0")
  testthat::expect_identical(sprintf("%05d", identity$observed_bsd_sum),
    "07119")
  testthat::expect_equal(identity$observed_bsd_blocks, 860562L)
  testthat::expect_true(identity$gzip_integrity)
  testthat::expect_true(identity$all_identity_gates_passed)
  testthat::expect_true(decision$reference_input_identity_closed)
  testthat::expect_false(decision$cellranger_acquisition_authorized)
  testthat::expect_false(decision$cellranger_runtime_validated)
  testthat::expect_false(decision$mkref_authorized)
  testthat::expect_false(decision$raw_reprocessing_authorized)
  testthat::expect_false(decision$label_access_authorized)
  testthat::expect_false(decision$biological_outcomes_authorized)
  testthat::expect_true("MV08H_ENSEMBL93_REFERENCE_INPUT_CLOSURE_2026-08-17.md" %in%
    artifacts$file)
  testthat::expect_false(any(artifacts$contains_expression))
  testthat::expect_false(any(artifacts$contains_cell_barcode))
  testthat::expect_false(any(artifacts$contains_absolute_private_path))
  testthat::expect_false(any(artifacts$contains_donor_attribute))
  testthat::expect_false(any(artifacts$contains_outcome_label))
  paths <- file.path(evidence, artifacts$file)
  observed <- unname(vapply(paths, function(path) digest::digest(
    file = path, algo = "sha256", serialize = FALSE), character(1L)))
  testthat::expect_identical(observed, artifacts$sha256)
  validator <- paste(readLines(file.path(root, "scripts",
    "validate_mv08h_reference_input.py"), warn = FALSE), collapse = "\n")
  for (term in c("EXPECTED_FASTA_SHA256", "EXPECTED_BSD_SUM",
                 "validate_gzip", "cellranger_acquisition_authorized",
                 "stop_before_runtime_or_mkref")) {
    testthat::expect_match(validator, term, fixed = TRUE)
  }
})

testthat::test_that("MV8-H prospectively binds the installed Cell Ranger 8.0.1 runtime", {
  root <- testthat::test_path("..", "..")
  evidence <- file.path(root, "docs", "audits",
    "mv08h-cellranger8-runtime-reference-prefreeze-v2")
  runtime <- utils::read.csv(file.path(evidence,
    "mv08h-runtime-identity.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  testthat::expect_equal(nrow(runtime), 1L)
  testthat::expect_identical(runtime$reported_version,
    "cellranger cellranger-8.0.1")
  testthat::expect_identical(runtime$launcher_sha256,
    "4ee3a1670b4f14c826004fe8e17b4759e1edc701b15ff2e9623753bf1b34d4d6")
  testthat::expect_identical(runtime$star_version, "2.7.2a")
  testthat::expect_identical(runtime$samtools_version, "samtools 1.16.1")
  testthat::expect_match(runtime$tree_sha256, "^[0-9a-f]{64}$")
  testthat::expect_gt(runtime$regular_files, 25000)
  testthat::expect_gt(runtime$regular_file_bytes, 1900000000)
  testthat::expect_true(runtime$required_cli_controls_verified)
  testthat::expect_true(runtime$all_runtime_gates_passed)
})

testthat::test_that("MV8-H runtime amendment rebinds exact reference and panel identities", {
  root <- testthat::test_path("..", "..")
  evidence <- file.path(root, "docs", "audits",
    "mv08h-cellranger8-runtime-reference-prefreeze-v2")
  reference <- utils::read.csv(file.path(evidence,
    "mv08h-reference-binding.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  testthat::expect_equal(nrow(reference), 5L)
  testthat::expect_true(all(reference$gate_passed))
  testthat::expect_identical(
    reference$sha256[reference$resource_id == "ensembl93_primary_assembly_fasta_gz"],
    "2a27436d44f0d6350f86894fbe5edec56faa5467028879784508041562406aa0")
  custom <- reference[reference$resource_id == "target_complete_33563_gtf", ]
  for (term in c("genes=33563", "feature_records=2565751",
                 "all_exact500_present=true",
                 "independent_repeat_byte_identical=true")) {
    testthat::expect_match(custom$structural_detail, term, fixed = TRUE)
  }
  testthat::expect_match(
    reference$structural_detail[reference$resource_id == "exact500_panel"],
    "features=500", fixed = TRUE)
  testthat::expect_match(
    reference$structural_detail[reference$resource_id == "common475_panel"],
    "features=475", fixed = TRUE)
})

testthat::test_that("MV8-H runtime amendment opens only mkref and preserves landscapes", {
  root <- testthat::test_path("..", "..")
  evidence <- file.path(root, "docs", "audits",
    "mv08h-cellranger8-runtime-reference-prefreeze-v2")
  processing <- utils::read.csv(file.path(evidence,
    "mv08h-processing-amendment.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  decision <- utils::read.csv(file.path(evidence,
    "mv08h-decision.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  count <- processing$frozen_value[processing$step_id == "count_sentinel"]
  for (term in c("--chemistry=SC3Pv2", "--include-introns=false",
                 "--create-bam=false", "--nosecondary",
                 "--localcores=16", "--localmem=64")) {
    testthat::expect_match(count, term, fixed = TRUE)
  }
  landscapes <- processing$frozen_value[
    processing$step_id == "topology_landscapes"]
  for (term in c("cell_topology_v1", "gene_topology_v1",
                 "complete VR H0/H1",
                 "every consecutive active landscape level",
                 "no fixed grid", "no level cap", "H0/H1 separate",
                 "essential H0 excluded")) {
    testthat::expect_match(landscapes, term, fixed = TRUE)
  }
  testthat::expect_equal(sum(processing$execution_authorized), 1L)
  testthat::expect_true(processing$execution_authorized[
    processing$step_id == "mkref"])
  testthat::expect_true(decision$mkref_authorized)
  testthat::expect_false(decision$count_sentinel_authorized)
  testthat::expect_false(decision$remaining_units_authorized)
  testthat::expect_false(decision$pca_ph_landscape_authorized)
  testthat::expect_false(decision$label_access_authorized)
  testthat::expect_false(decision$biological_outcomes_authorized)
})

testthat::test_that("MV8-H Cell Ranger 8.0.1 evidence is independently closed", {
  root <- testthat::test_path("..", "..")
  evidence <- file.path(root, "docs", "audits",
    "mv08h-cellranger8-runtime-reference-prefreeze-v2")
  validation <- utils::read.csv(file.path(evidence,
    "mv08h-independent-validation.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  repeat_validation <- utils::read.csv(file.path(evidence,
    "mv08h-repeat-validation.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  artifacts <- utils::read.csv(file.path(evidence,
    "mv08h-artifact-manifest.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  testthat::expect_equal(nrow(validation), 12L)
  testthat::expect_true(all(validation$passed))
  testthat::expect_equal(nrow(repeat_validation), 6L)
  testthat::expect_true(all(repeat_validation$byte_identical))
  testthat::expect_false(any(artifacts$contains_expression))
  testthat::expect_false(any(artifacts$contains_cell_barcode))
  testthat::expect_false(any(artifacts$contains_absolute_private_path))
  testthat::expect_false(any(artifacts$contains_donor_attribute))
  testthat::expect_false(any(artifacts$contains_outcome_label))
  paths <- file.path(evidence, artifacts$file)
  observed <- unname(vapply(paths, function(path) digest::digest(
    file = path, algo = "sha256", serialize = FALSE), character(1L)))
  testthat::expect_identical(observed, artifacts$sha256)
  builder <- paste(readLines(file.path(root, "scripts",
    "build_mv08h_cellranger8_runtime_prefreeze.py"), warn = FALSE),
    collapse = "\n")
  validator <- paste(readLines(file.path(root, "scripts",
    "validate_mv08h_cellranger8_runtime_prefreeze.py"), warn = FALSE),
    collapse = "\n")
  for (term in c("distribution_tree_identity", "--include-introns=false",
                 "--create-bam=false", "runtime_reference_exact_authorize_mkref_only")) {
    testthat::expect_match(builder, term, fixed = TRUE)
  }
  for (term in c("runtime_tree", "public_firewall",
                 "deterministic_builder_repeat", "12/12 checks")) {
    testthat::expect_match(validator, term, fixed = TRUE)
  }
})
