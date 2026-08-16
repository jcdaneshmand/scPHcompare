#!/usr/bin/env Rscript

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1")
options(warn = 2)
for (package in c("digest", "Matrix")) {
  if (!requireNamespace(package, quietly = TRUE)) stop(package, " required.")
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 11L) {
  stop(
    "usage: build_mv08f_cache_recovery_prefreeze.R RECONCILIATION ",
    "CACHE_MANIFEST ADDED_CACHE SOURCE_ROOT MV07H_SOURCE_ROOT MV07H_PH_ROOT ",
    "MV07H_LANDSCAPE_ROOT PRIVATE_INPUT_DIR PUBLIC_DIR HISTORICAL_SHA EXPECTED_HEAD"
  )
}
reconciliation_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
manifest_path <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
added_cache <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
source_root <- normalizePath(args[[4L]], winslash = "/", mustWork = TRUE)
source_artifact_root <- normalizePath(args[[5L]], winslash = "/", mustWork = TRUE)
ph_artifact_root <- normalizePath(args[[6L]], winslash = "/", mustWork = TRUE)
landscape_root <- normalizePath(args[[7L]], winslash = "/", mustWork = TRUE)
private_dir <- args[[8L]]
public_dir <- args[[9L]]
historical_sha <- tolower(trimws(args[[10L]]))
expected_head <- tolower(trimws(args[[11L]]))
if (!grepl("^[0-9a-f]{64}$", historical_sha) ||
    !grepl("^[0-9a-f]{40}$", expected_head)) stop("Invalid frozen digest.")
head <- tolower(trimws(system2("git", c("rev-parse", "HEAD"), stdout = TRUE)))
if (!identical(head, expected_head)) stop("MV8-F exact HEAD mismatch.")
for (path in c(private_dir, public_dir)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  if (length(list.files(path, all.files = TRUE, no.. = TRUE))) {
    stop("MV8-F output directory must be empty: ", path)
  }
}
source("R/provenance_utils.R")
source("R/mv05_resource_safe_execution.R")
sha <- function(path) digest::digest(file = path, algo = "sha256",
                                     serialize = FALSE)
truth <- function(value) if (is.logical(value)) !is.na(value) & value else
  tolower(trimws(as.character(value))) == "true"
write_once <- function(value, directory, name) {
  path <- file.path(directory, name)
  if (file.exists(path)) stop("Refusing overwrite: ", path)
  write_provenance_csv(value, path)
  path
}

reconciliation <- read.csv(reconciliation_path, stringsAsFactors = FALSE,
                           check.names = FALSE)
manifest <- read.csv(manifest_path, stringsAsFactors = FALSE,
                     check.names = FALSE)
primary_samples <- reconciliation[truth(reconciliation$corrected_primary_90),
                                  c("sample_id", "post_qc_cells"), drop = FALSE]
primary_samples <- primary_samples[order(primary_samples$sample_id,
                                         method = "radix"),, drop = FALSE]
primary <- manifest[manifest$source_tier == "primary90",, drop = FALSE]
added <- manifest[manifest$source_tier == "added34",, drop = FALSE]
primary <- primary[order(primary$seed, primary$sample_id, method = "radix"),,
                   drop = FALSE]
if (nrow(primary_samples) != 90L || nrow(primary) != 450L ||
    nrow(added) != 170L || length(unique(primary$sample_id)) != 90L ||
    !setequal(primary_samples$sample_id, unique(primary$sample_id)) ||
    !identical(sort(unique(as.integer(primary$seed))), 20260805:20260809) ||
    any(table(primary$seed) != 90L) || any(primary$selected_cells != 384L) ||
    any(primary$outcome_label_state != "closed") ||
    any(truth(primary$biological_outcomes_computed))) {
  stop("MV8-F primary recovery axis differs from the accepted 90 by five.")
}

all_sources <- list.files(source_root, pattern = "[.]rds$", recursive = TRUE,
                          full.names = TRUE, ignore.case = TRUE)
source_ids <- tools::file_path_sans_ext(basename(all_sources))
coverage_counts <- vapply(primary_samples$sample_id,
                          function(id) sum(source_ids == id), integer(1L))
if (any(coverage_counts != 1L)) stop("Primary source coverage is not unique.")

added_paths <- file.path(added_cache, added$private_cache_file)
if (any(!file.exists(added_paths))) stop("An accepted added34 cache is absent.")
added_bytes <- as.numeric(file.info(added_paths)$size)
added_hashes <- vapply(added_paths, sha, character(1L))
if (!identical(unname(added_hashes), unname(added$private_cache_sha256)) ||
    !identical(unname(added_bytes), as.numeric(added$private_cache_bytes))) {
  stop("An accepted added34 cache has drifted.")
}
runtime_record <- readRDS(added_paths[[1L]])
mv05d0_validate_normalization_cache_record_v2(runtime_record)
current_runtime <- mv05d0_normalization_runtime_v1()
runtime_equal <- identical(runtime_record$identity$runtime, current_runtime)
if (!runtime_equal) stop("Current normalization runtime differs from accepted.")

source_metrics_path <- "docs/audits/mv07h-full-ph-evidence/mv07h-source-metrics.csv"
ph_metrics_path <- "docs/audits/mv07h-full-ph-evidence/mv07h-ph-metrics.csv"
landscape_inventory_path <- paste0(
  "docs/audits/mv07h-landscape-complete-validation/",
  "mv07h-landscape-complete-group-inventory.csv")
source_metrics <- read.csv(source_metrics_path, stringsAsFactors = FALSE,
                           check.names = FALSE)
ph_metrics <- read.csv(ph_metrics_path, stringsAsFactors = FALSE,
                       check.names = FALSE)
landscape_inventory <- read.csv(landscape_inventory_path,
                                stringsAsFactors = FALSE, check.names = FALSE)
source_paths <- file.path(source_artifact_root,
  paste0("mv07h__", source_metrics$seed, "__source.rds"))
ph_paths <- file.path(ph_artifact_root, paste0(
  "mv07h__", ph_metrics$seed, "__", ph_metrics$sample_id, "__",
  ph_metrics$view_id, "__ph.rds"))
landscape_dirs <- file.path(landscape_root,
  gsub(":", "_", landscape_inventory$group_id, fixed = TRUE))
landscape_distance <- file.path(landscape_dirs, "distances.csv")
landscape_metrics <- file.path(landscape_dirs, "metrics.csv")
landscape_status <- file.path(landscape_dirs, "status.csv")
if (nrow(source_metrics) != 5L || nrow(ph_metrics) != 1240L ||
    nrow(landscape_inventory) != 20L ||
    any(!file.exists(c(source_paths, ph_paths, landscape_distance,
                       landscape_metrics, landscape_status)))) {
  stop("Accepted MV7-H comparator artifact axis is incomplete.")
}
source_hash_ok <- identical(unname(vapply(source_paths, sha, character(1L))),
                            unname(source_metrics$output_sha256))
ph_hash_ok <- identical(unname(vapply(ph_paths, sha, character(1L))),
                        unname(ph_metrics$output_sha256))
landscape_hash_ok <- all(
  vapply(landscape_distance, sha, character(1L)) ==
    landscape_inventory$distances_sha256,
  vapply(landscape_metrics, sha, character(1L)) ==
    landscape_inventory$metrics_sha256,
  vapply(landscape_status, sha, character(1L)) ==
    landscape_inventory$status_sha256)
if (!source_hash_ok || !ph_hash_ok || !landscape_hash_ok) {
  stop("Accepted MV7-H comparator artifact hash drift.")
}

candidates <- data.frame(
  contract_id = "mv08f_primary90_candidate_v1",
  sample_id = primary_samples$sample_id,
  filtered_cells = as.integer(primary_samples$post_qc_cells),
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE)
source_preflight <- data.frame(
  contract_id = "mv08f_historical_source_preflight_v1",
  source_sha256 = historical_sha, outcome_label_state = "closed",
  biological_outcomes_computed = FALSE, stringsAsFactors = FALSE)
write_once(candidates, private_dir, "mv08f-primary90-candidates.csv")
write_once(source_preflight, private_dir, "mv08f-source-preflight.csv")

recovery_axis <- primary
recovery_axis$contract_id <- "mv08f_primary90_recovery_axis_v1"
coverage <- data.frame(
  contract_id = "mv08f_primary90_source_coverage_v1",
  expected_samples = 90L, uniquely_located_sources = sum(coverage_counts == 1L),
  minimum_post_qc_cells = min(primary_samples$post_qc_cells),
  maximum_post_qc_cells = max(primary_samples$post_qc_cells),
  source_root_published = FALSE, source_paths_published = FALSE,
  labels_present = FALSE, stringsAsFactors = FALSE)
cache_status <- data.frame(
  contract_id = "mv08f_live_cache_status_v1",
  source_tier = c("primary90", "added34"), expected_cache_records = c(450L, 170L),
  live_cache_records = c(0L, 170L), exact_byte_matches = c(0L, 170L),
  exact_sha256_matches = c(0L, 170L), action = c("reconstruct_then_exact_validate",
    "reuse_immutable"), stringsAsFactors = FALSE)
runtime <- data.frame(
  contract_id = "mv08f_runtime_identity_v1",
  runtime_exact = runtime_equal, r_version = current_runtime$r_version,
  seurat_version = current_runtime$seurat_version,
  seurat_object_version = current_runtime$seurat_object_version,
  sctransform_version = current_runtime$sctransform_version,
  matrix_version = current_runtime$matrix_version,
  future_plan = current_runtime$future_plan,
  thread_identity = paste(current_runtime$omp_num_threads,
    current_runtime$openblas_num_threads, current_runtime$mkl_num_threads,
    sep = "/"), stringsAsFactors = FALSE)
artifacts <- data.frame(
  contract_id = "mv08f_existing_comparator_status_v1",
  artifact_kind = c("source_bundle", "ph_record", "landscape_group"),
  expected_count = c(5L, 1240L, 20L), live_count = c(length(source_paths),
    length(ph_paths), length(landscape_dirs)), exact_hash_identity = TRUE,
  action = "reuse_immutable", stringsAsFactors = FALSE)
caps <- data.frame(
  contract_id = "mv08f_recovery_resource_caps_v1",
  scope = c("raw_child", "sct_child", "aggregate_storage"),
  elapsed_cap_seconds = c(1800, 1800, NA), rss_cap_bytes = c(8, 8, NA) * 1024^3,
  storage_cap_bytes = c(NA, NA, 40 * 1024^3), maximum_workers = c(2L, 2L, 2L),
  automatic_retry = FALSE, stringsAsFactors = FALSE)
source_files <- c(
  reconciliation = reconciliation_path, accepted_cache_manifest = manifest_path,
  source_metrics = source_metrics_path, ph_metrics = ph_metrics_path,
  landscape_inventory = landscape_inventory_path,
  runtime_helper = "R/mv05_resource_safe_execution.R",
  raw_queue = "scripts/run_mv05d0_raw_shard_queue.R",
  raw_entry = "scripts/run_mv05d0_raw_shard_entry.R",
  selection_builder = "scripts/build_mv05d0_selection_summary_from_shards.R",
  sct_queue = "scripts/run_mv05d0_sct_cache_queue.R",
  sct_entry = "scripts/run_mv05d0_sct_cache_entry.R",
  builder = "scripts/build_mv08f_cache_recovery_prefreeze.R",
  validator = "scripts/validate_mv08f_cache_recovery_prefreeze.R",
  specification = "docs/specifications/MV08F_REFERENCE_CACHE_RECOVERY_PREFREEZE_V1.md",
  tests = "tests/testthat/test-mv08f-cache-recovery-prefreeze.R")
if (any(!file.exists(source_files))) stop("MV8-F source freeze incomplete.")
locators <- unname(source_files)
locators[1:2] <- c(
  "docs/audits/mv07d-prefreeze-evidence/mv07d-sample-reconciliation.csv",
  "docs/audits/mv07fp-prefreeze-evidence-v4/mv07fp-cache-manifest.csv")
freeze <- data.frame(
  contract_id = "mv08f_recovery_source_freeze_v1",
  source_id = names(source_files), artifact_locator = locators,
  sha256 = vapply(source_files, sha, character(1L)),
  bytes = as.numeric(file.info(source_files)$size), accepted_head = expected_head,
  private_source = FALSE, stringsAsFactors = FALSE)
decision <- data.frame(
  contract_id = "mv08f_recovery_prefreeze_decision_v1",
  decision = "authorize_primary90_cache_reconstruction_only",
  raw_jobs_authorized = 90L, sct_jobs_authorized = 450L,
  exact_manifest_validation_required = TRUE, panel475_source_authorized = FALSE,
  ph_authorized = FALSE, landscape_authorized = FALSE,
  clustering_authorized = FALSE, hca_fastq_download_authorized = FALSE,
  label_access_authorized = FALSE, outcome_jobs_authorized = 0L,
  next_gate = "MV8-F_cache_recovery_independent_validation",
  stringsAsFactors = FALSE)
outputs <- list(
  "mv08f-recovery-axis.csv" = recovery_axis,
  "mv08f-source-coverage.csv" = coverage,
  "mv08f-live-cache-status.csv" = cache_status,
  "mv08f-runtime-identity.csv" = runtime,
  "mv08f-existing-comparator-status.csv" = artifacts,
  "mv08f-resource-caps.csv" = caps,
  "mv08f-source-freeze.csv" = freeze,
  "mv08f-decision.csv" = decision)
paths <- vapply(names(outputs), function(name) write_once(outputs[[name]],
  public_dir, name), character(1L))
artifact_manifest <- data.frame(
  contract_id = "mv08f_recovery_prefreeze_artifact_manifest_v1",
  file = basename(paths), bytes = as.numeric(file.info(paths)$size),
  sha256 = vapply(paths, sha, character(1L)),
  contains_expression = FALSE, contains_cell_barcode = FALSE,
  contains_absolute_private_path = FALSE, stringsAsFactors = FALSE)
write_once(artifact_manifest, public_dir, "mv08f-artifact-manifest.csv")
message("MV8-F recovery prefreeze complete: 90 raw + 450 SCT jobs only")
