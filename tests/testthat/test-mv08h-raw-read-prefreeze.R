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

testthat::test_that("MV8-H mkref execution amendment is downward-only and fail closed", {
  root <- testthat::test_path("..", "..")
  evidence <- file.path(root, "docs", "audits",
    "mv08h-cellranger8-mkref-prefreeze-v1")
  freeze <- utils::read.csv(file.path(evidence,
    "mv08h-mkref-execution-prefreeze.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  spec <- paste(readLines(file.path(root, "docs", "specifications",
    "MV08H_CELLRANGER8_MKREF_EXECUTION_PREFREEZE_V1.md"), warn = FALSE),
    collapse = "\n")
  launcher <- paste(readLines(file.path(root, "scripts",
    "run_mv08h_cellranger8_mkref.py"), warn = FALSE), collapse = "\n")
  testthat::expect_equal(
    as.integer(freeze$selected_value[freeze$resource == "cores"]), 4L)
  testthat::expect_equal(
    as.integer(freeze$selected_value[freeze$resource == "memory_gib"]), 32L)
  testthat::expect_true(all(freeze$gate %in%
    c("exact", "downward_only", "unchanged", "non_destructive", "authorized")))
  for (term in c("--nthreads=4", "--memgb=32", "--localcores=4",
                 "--localmem=32", "50-GiB", "1-TiB", "does not kill",
                 "does not authorize", "persistence landscapes")) {
    testthat::expect_match(spec, term, fixed = TRUE)
  }
  for (term in c("EXPECTED_FASTA_SHA256", "EXPECTED_GTF_SHA256",
                 "REFERENCE_CAP_BYTES", "FREE_SPACE_FLOOR_BYTES",
                 "process_tree_rss_kib", "resource_breach_detected",
                 "memory_allocation_passed", "automatic_kill_used",
                 "deletion_used", "--dry-run")) {
    testthat::expect_match(launcher, term, fixed = TRUE)
  }
})

testthat::test_that("MV8-H custom reference closes exact inputs and feature axes", {
  root <- testthat::test_path("..", "..")
  evidence <- file.path(root, "docs", "audits",
    "mv08h-cellranger8-mkref-prefreeze-v1")
  tree <- utils::read.csv(file.path(evidence,
    "mv08h-reference-tree-identity.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  inputs <- utils::read.csv(file.path(evidence,
    "mv08h-reference-input-closure.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  features <- utils::read.csv(file.path(evidence,
    "mv08h-reference-feature-closure.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  components <- utils::read.csv(file.path(evidence,
    "mv08h-reference-components.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  testthat::expect_equal(tree$regular_files, 19L)
  testthat::expect_equal(tree$regular_file_bytes, 20765871518)
  testthat::expect_identical(tree$tree_sha256,
    "5e2aff9e7154e6b02f98552a4419bd48edce66e617e579ae562e714f79199f1c")
  testthat::expect_equal(tree$partial_files, 0L)
  testthat::expect_true(tree$all_tree_gates_passed)
  testthat::expect_equal(nrow(components), 19L)
  testthat::expect_true(all(!grepl("^[A-Za-z]:|^/mnt/|^/home/",
    components$relative_path)))
  testthat::expect_true(all(inputs$gate_passed))
  testthat::expect_identical(
    inputs$logical_sha256[inputs$resource == "embedded_fasta"],
    "78777b0886e8dfa5e14e4957fbbaa53736fcbaa5668d59e09b6b7945fca93d8c")
  testthat::expect_identical(
    inputs$logical_sha256[inputs$resource == "embedded_gtf_gz"],
    "e28e4c4faf0dd76884d5e94c481fce2db43ad303968067c1276092a234727182")
  testthat::expect_equal(features$unique_genes, 33563L)
  testthat::expect_equal(features$all_feature_records, 2565751L)
  testthat::expect_equal(features$exact500_present, 500L)
  testthat::expect_equal(features$common475_present, 475L)
  testthat::expect_true(features$all_feature_gates_passed)
})

testthat::test_that("MV8-H mkref resource and authorization closure stays narrow", {
  root <- testthat::test_path("..", "..")
  evidence <- file.path(root, "docs", "audits",
    "mv08h-cellranger8-mkref-prefreeze-v1")
  resource <- utils::read.csv(file.path(evidence,
    "mv08h-mkref-resource-closure.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  validation <- utils::read.csv(file.path(evidence,
    "mv08h-mkref-independent-validation.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  decision <- utils::read.csv(file.path(evidence,
    "mv08h-mkref-decision.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  testthat::expect_equal(resource$elapsed_seconds, 13603L)
  testthat::expect_equal(resource$monitor_samples, 450L)
  testthat::expect_equal(resource$selected_cores, 4L)
  testthat::expect_equal(resource$selected_memory_gib, 32L)
  testthat::expect_equal(resource$observed_peak_rss_kib, 31975732L)
  testthat::expect_lt(resource$observed_peak_rss_gib, 32)
  testthat::expect_lt(resource$observed_peak_run_tree_bytes,
    resource$reference_cap_bytes)
  testthat::expect_gt(resource$minimum_free_bytes,
    resource$free_space_floor_bytes)
  testthat::expect_false(resource$resource_breach_detected)
  testthat::expect_true(resource$all_resource_gates_passed)
  testthat::expect_equal(nrow(validation), 15L)
  testthat::expect_true(all(validation$passed))
  testthat::expect_identical(decision$decision,
    "reference_exact_authorize_count_sentinel_prefreeze_only")
  testthat::expect_true(decision$count_sentinel_prefreeze_authorized)
  testthat::expect_false(decision$count_sentinel_execution_authorized)
  testthat::expect_false(decision$remaining_units_authorized)
  testthat::expect_false(decision$qc_pca_ph_landscape_authorized)
  testthat::expect_false(decision$label_access_authorized)
  testthat::expect_false(decision$biological_outcomes_authorized)
  testthat::expect_false(decision$deletion_authorized)
})

testthat::test_that("MV8-H custom-reference evidence preserves firewall and landscapes", {
  root <- testthat::test_path("..", "..")
  evidence <- file.path(root, "docs", "audits",
    "mv08h-cellranger8-mkref-prefreeze-v1")
  artifacts <- utils::read.csv(file.path(evidence,
    "mv08h-reference-artifact-manifest.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  report <- paste(readLines(file.path(evidence,
    "MV08H_CELLRANGER8_REFERENCE_CLOSURE_2026-08-18.md"), warn = FALSE),
    collapse = "\n")
  validator <- paste(readLines(file.path(root, "scripts",
    "validate_mv08h_cellranger8_reference.py"), warn = FALSE),
    collapse = "\n")
  for (field in c("contains_expression", "contains_cell_barcode",
                  "contains_absolute_private_path", "contains_donor_attribute",
                  "contains_outcome_label")) {
    testthat::expect_false(any(artifacts[[field]]))
  }
  paths <- file.path(evidence, artifacts$file)
  observed <- unname(vapply(paths, function(path) digest::digest(
    file = path, algo = "sha256", serialize = FALSE), character(1L)))
  testthat::expect_identical(observed, artifacts$sha256)
  for (term in c("separate typed views", "H0 and H1 remain separate",
                 "every consecutive active level", "no fixed grid",
                 "no universal level cap", "does **not** execute",
                 "RSS in its composite", "forward-only")) {
    testthat::expect_match(report, term, fixed = TRUE)
  }
  for (term in c("tree_identity", "parse_embedded_gtf",
                 "MEMORY_ALLOCATION_KIB", "count_not_executed",
                 "reference_exact_authorize_count_sentinel_prefreeze_only")) {
    testthat::expect_match(validator, term, fixed = TRUE)
  }
})

testthat::test_that("MV8-H count sentinel is selected prospectively without labels", {
  root <- testthat::test_path("..", "..")
  evidence <- file.path(root, "docs", "audits",
    "mv08h-cellranger8-count-sentinel-prefreeze-v1")
  selection <- utils::read.csv(file.path(evidence,
    "mv08h-count-sentinel-selection.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  fastqs <- utils::read.csv(file.path(evidence,
    "mv08h-count-sentinel-fastqs.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  testthat::expect_equal(nrow(selection), 8L)
  testthat::expect_identical(selection$burden_rank, 1:8)
  testthat::expect_equal(sum(selection$selected), 1L)
  testthat::expect_identical(selection$unit_id[selection$selected],
    "HCA_BM_002")
  testthat::expect_equal(selection$fastq_bytes[selection$selected],
    11249623632)
  testthat::expect_false(any(selection$label_fields_consulted))
  testthat::expect_false(any(selection$expression_or_qc_consulted))
  testthat::expect_equal(nrow(fastqs), 6L)
  testthat::expect_equal(sum(fastqs$expected_bytes), 11249623632)
  testthat::expect_true(all(fastqs$cache_identity_passed))
  testthat::expect_false(any(fastqs$private_path_published))
  testthat::expect_true(all(grepl("^[0-9a-f]{64}$", fastqs$expected_sha256)))
})

testthat::test_that("MV8-H count command and conservative resources are exact", {
  root <- testthat::test_path("..", "..")
  evidence <- file.path(root, "docs", "audits",
    "mv08h-cellranger8-count-sentinel-prefreeze-v1")
  command <- utils::read.csv(file.path(evidence,
    "mv08h-count-sentinel-command.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  resources <- utils::read.csv(file.path(evidence,
    "mv08h-count-sentinel-resource-policy.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  reference <- utils::read.csv(file.path(evidence,
    "mv08h-count-sentinel-reference-binding.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  public <- command$frozen_value[command$parameter == "public_command"]
  for (term in c(
      "--id=mv08h_count_sentinel_hca_bm_002",
      "--transcriptome=<verified_custom_reference>",
      "--fastqs=<verified_hca_bm_002_fastq_directory>",
      "--sample=MantonBM2_HiSeq_9", "--chemistry=SC3Pv2",
      "--expect-cells=7000", "--include-introns=false",
      "--create-bam=false", "--nosecondary", "--localcores=4",
      "--localmem=32", "--disable-ui")) {
    testthat::expect_match(public, term, fixed = TRUE)
  }
  values <- stats::setNames(resources$selected_value, resources$resource)
  testthat::expect_equal(values[["local_cores"]], 4)
  testthat::expect_equal(values[["local_memory"]], 32)
  testthat::expect_equal(values[["process_tree_rss_absolute_cap"]],
    80 * 1024^3)
  testthat::expect_equal(values[["run_workspace_cap"]], 200 * 1024^3)
  testthat::expect_equal(values[["minimum_free_space"]], 1024^4)
  testthat::expect_equal(values[["elapsed_observation_cap"]], 96 * 60 * 60)
  testthat::expect_identical(values[["automatic_termination"]], "false")
  testthat::expect_identical(values[["automatic_deletion"]], "false")
  testthat::expect_identical(reference$tree_sha256,
    "5e2aff9e7154e6b02f98552a4419bd48edce66e617e579ae562e714f79199f1c")
  testthat::expect_true(reference$binding_passed)
})

testthat::test_that("MV8-H count prefreeze preserves firewall and stop boundary", {
  root <- testthat::test_path("..", "..")
  evidence <- file.path(root, "docs", "audits",
    "mv08h-cellranger8-count-sentinel-prefreeze-v1")
  validation <- utils::read.csv(file.path(evidence,
    "mv08h-count-sentinel-validation-contract.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  firewall <- utils::read.csv(file.path(evidence,
    "mv08h-count-sentinel-firewall.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  decision <- utils::read.csv(file.path(evidence,
    "mv08h-count-sentinel-decision.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  testthat::expect_equal(nrow(validation), 12L)
  testthat::expect_false(any(validation$count_execution_performed_by_prefreeze))
  validation_text <- paste(validation$frozen_requirement, collapse = " ")
  for (term in c("HDF5", "33,563", "exact500", "common475",
                 "barcodes unique and private", "persistence landscape",
                 "remaining-unit decision")) {
    testthat::expect_match(validation_text, term, fixed = TRUE)
  }
  public <- firewall$field_class[firewall$public_release_permitted]
  testthat::expect_identical(public, c(
    "FASTQ byte sizes and SHA-256",
    "reference/runtime relative identities",
    "aggregate resource samples and run status",
    "matrix dimensions and feature-axis identities"))
  for (field in c("expression or UMI values", "cell barcodes",
                  "donor attributes or identifiers",
                  "QC values or eligibility decisions",
                  "study/tissue/approach labels", "biological outcomes")) {
    row <- firewall[firewall$field_class == field, ]
    testthat::expect_equal(nrow(row), 1L)
    testthat::expect_false(row$public_release_permitted)
  }
  testthat::expect_identical(decision$decision,
    "count_sentinel_prefreeze_exact_await_execution_authorization")
  testthat::expect_true(decision$count_sentinel_prefreeze_completed)
  for (field in c("count_sentinel_execution_authorized",
                  "matrix_access_authorized", "qc_authorized",
                  "remaining_units_authorized",
                  "pca_ph_landscape_authorized",
                  "label_access_authorized",
                  "biological_outcomes_authorized",
                  "deletion_authorized")) {
    testthat::expect_false(decision[[field]])
  }
})

testthat::test_that("MV8-H count prefreeze is independently closed and runnable only explicitly", {
  root <- testthat::test_path("..", "..")
  evidence <- file.path(root, "docs", "audits",
    "mv08h-cellranger8-count-sentinel-prefreeze-v1")
  validation <- utils::read.csv(file.path(evidence,
    "mv08h-count-sentinel-independent-validation.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  artifacts <- utils::read.csv(file.path(evidence,
    "mv08h-count-sentinel-artifact-manifest.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  report <- paste(readLines(file.path(evidence,
    "MV08H_CELLRANGER8_COUNT_SENTINEL_PREFREEZE_2026-08-18.md"),
    warn = FALSE), collapse = "\n")
  runner <- paste(readLines(file.path(root, "scripts",
    "run_mv08h_cellranger8_count_sentinel.py"), warn = FALSE),
    collapse = "\n")
  testthat::expect_equal(nrow(validation), 14L)
  testthat::expect_true(all(validation$passed))
  for (field in c("contains_expression", "contains_cell_barcode",
                  "contains_absolute_private_path", "contains_donor_attribute",
                  "contains_qc_value", "contains_outcome_label")) {
    testthat::expect_false(any(artifacts[[field]]))
  }
  paths <- file.path(evidence, artifacts$file)
  observed <- unname(vapply(paths, function(path) digest::digest(
    file = path, algo = "sha256", serialize = FALSE), character(1L)))
  testthat::expect_identical(observed, artifacts$sha256)
  for (term in c("does **not** execute or authorize count", "4 cores",
                 "32 GiB", "96 hours", "does not kill or delete",
                 "separate typed observation views", "H0 and H1 remain separate",
                 "every consecutive active level", "no fixed grid",
                 "no universal level cap")) {
    testthat::expect_match(report, term, fixed = TRUE)
  }
  for (term in c("--dry-run", "EXECUTION_TOKEN",
                 "EXPECTED_REFERENCE_TREE_SHA256",
                 "competing_cellranger_processes", "process_tree_rss_kib",
                 "RSS_CAP_BYTES", "WORKSPACE_CAP_BYTES",
                 "ELAPSED_CAP_SECONDS", "FREE_SPACE_FLOOR_BYTES",
                 "resource_breach_detected", "automatic_kill_used",
                 "deletion_used", "--create-bam=false", "--nosecondary",
                 "--localcores=4", "--localmem=32")) {
    testthat::expect_match(runner, term, fixed = TRUE)
  }
})
