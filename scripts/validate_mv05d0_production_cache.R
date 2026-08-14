#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) {
  stop(
    "usage: validate_mv05d0_production_cache.R CACHE_DIR SELECTION_CSV ",
    "SCT_METRICS_CSV RAW_METRICS_CSV PREVIOUS_PROJECTION_CSV ",
    "SUMMARY_CSV PROJECTION_CSV", call. = FALSE
  )
}
cache_dir <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
selection_path <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
sct_metrics_path <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
raw_metrics_path <- normalizePath(args[[4L]], winslash = "/", mustWork = TRUE)
previous_path <- normalizePath(args[[5L]], winslash = "/", mustWork = TRUE)
summary_path <- args[[6L]]
projection_path <- args[[7L]]
source("R/provenance_utils.R")
source("R/mv05_resource_safe_execution.R")

selection <- utils::read.csv(
  selection_path, stringsAsFactors = FALSE, check.names = FALSE
)
sct <- utils::read.csv(
  sct_metrics_path, stringsAsFactors = FALSE, check.names = FALSE
)
raw <- utils::read.csv(
  raw_metrics_path, stringsAsFactors = FALSE, check.names = FALSE
)
previous <- utils::read.csv(
  previous_path, stringsAsFactors = FALSE, check.names = FALSE
)
if (nrow(selection) != 450L || nrow(sct) != 450L || nrow(raw) != 90L ||
    any(c("tissue", "approach") %in%
          tolower(c(names(selection), names(sct), names(raw)))) ||
    any(selection$outcome_label_state != "closed") ||
    any(sct$outcome_label_state != "closed") ||
    any(raw$outcome_label_state != "closed") ||
    any(as.logical(selection$biological_outcomes_computed)) ||
    any(as.logical(sct$biological_outcomes_computed)) ||
    any(as.logical(raw$biological_outcomes_computed)) ||
    any(!sct$disposition %in% c("built_atomic", "reuse_validated")) ||
    any(!raw$disposition %in% c("built_atomic", "reuse_validated"))) {
  stop("Production manifests violate completeness or label boundaries.",
       call. = FALSE)
}
key <- function(x) paste(x$sample_id, x$seed, sep = "\r")
if (anyDuplicated(key(selection)) || anyDuplicated(key(sct)) ||
    !setequal(key(selection), key(sct)) ||
    anyDuplicated(raw$sample_id) ||
    sum(raw$recovered_monolithic_comparison == "exact_content_identity") != 53L ||
    sum(raw$recovered_monolithic_comparison == "not_available") != 37L) {
  stop("Production identities do not reconcile to 90 samples × 5 seeds.",
       call. = FALSE)
}
cache_files <- sort(list.files(
  cache_dir, pattern = "__sct\\.rds$", full.names = FALSE
), method = "radix")
if (length(cache_files) != 450L ||
    !setequal(cache_files, sct$private_cache_file)) {
  stop("Private cache directory does not contain exactly the audited files.",
       call. = FALSE)
}
selection_lookup <- selection[match(key(sct), key(selection)), , drop = FALSE]
runtime_hashes <- character(nrow(sct))
resume_validated <- logical(nrow(sct))
for (index in seq_len(nrow(sct))) {
  path <- file.path(cache_dir, sct$private_cache_file[[index]])
  record <- readRDS(path)
  mv05d0_validate_normalization_cache_record_v2(record)
  if (!identical(record$identity$sample_id, sct$sample_id[[index]]) ||
      record$identity$seed != as.integer(sct$seed[[index]]) ||
      !identical(record$identity$selected_cell_sha256,
                 selection_lookup$selected_cell_sha256[[index]]) ||
      !identical(record$cache_key, sct$normalization_cache_key[[index]]) ||
      !identical(record$payload_sha256, sct$payload_sha256[[index]]) ||
      unname(file.info(path)$size) != sct$private_cache_size_bytes[[index]] ||
      !identical(digest::digest(
        file = path, algo = "sha256", serialize = FALSE
      ), sct$private_cache_sha256[[index]])) {
    stop("Production cache identity/hash failed at row ", index,
         call. = FALSE)
  }
  runtime_hashes[[index]] <- digest::digest(
    record$identity$runtime, algo = "sha256", serialize = TRUE
  )
  resume_validated[[index]] <- TRUE
  rm(record)
  if (index %% 25L == 0L) invisible(gc())
}
if (length(unique(runtime_hashes)) != 1L) {
  stop("Production caches do not share one frozen runtime identity.",
       call. = FALSE)
}
total_cache_bytes <- sum(sct$private_cache_size_bytes)
actual_normalization_hours <- sum(sct$elapsed_seconds) / 3600
projection <- mv05d0_reproject_scenarios_v1(
  previous, actual_normalization_hours
)
summary <- data.frame(
  contract_id = "mv05d0_production_cache_validation_v1",
  samples = length(unique(sct$sample_id)),
  seeds = length(unique(sct$seed)), entries = nrow(sct),
  built_entries = sum(sct$disposition == "built_atomic"),
  resume_validated_entries = sum(resume_validated),
  failed_entries = sum(!sct$disposition %in% c("built_atomic", "reuse_validated")),
  raw_shards = nrow(raw),
  exact_recovered_monolithic_comparisons =
    sum(raw$recovered_monolithic_comparison == "exact_content_identity"),
  newly_recovered_individual_sources =
    sum(raw$recovered_monolithic_comparison == "not_available"),
  raw_shard_worker_hours = sum(raw$elapsed_seconds) / 3600,
  raw_shard_peak_rss_bytes = max(raw$peak_process_tree_rss_bytes),
  raw_shard_total_bytes = sum(raw$private_raw_cache_size_bytes),
  normalization_worker_hours = actual_normalization_hours,
  normalization_operation_hours = sum(sct$operation_seconds) / 3600,
  normalization_seconds_median = stats::median(sct$elapsed_seconds),
  normalization_seconds_max = max(sct$elapsed_seconds),
  normalization_peak_rss_bytes = max(sct$peak_process_tree_rss_bytes),
  normalization_cache_total_bytes = total_cache_bytes,
  normalization_cache_max_bytes = max(sct$private_cache_size_bytes),
  runtime_sha256 = unique(runtime_hashes),
  selection_set_sha256 = digest::digest(
    selection[order(key(selection), method = "radix"), ],
    algo = "sha256", serialize = TRUE
  ),
  cache_set_sha256 = digest::digest(
    sct[order(sct$private_cache_file, method = "radix"),
        c("private_cache_file", "private_cache_sha256",
          "normalization_cache_key", "payload_sha256")],
    algo = "sha256", serialize = TRUE
  ),
  elapsed_cap_seconds = unique(sct$elapsed_cap_seconds),
  rss_cap_bytes = unique(sct$rss_cap_bytes),
  storage_cap_bytes = unique(sct$storage_cap_bytes),
  fold_jobs_executed = 0L, ph_jobs_executed = 0L,
  landscape_jobs_executed = 0L, distance_jobs_executed = 0L,
  clustering_jobs_executed = 0L, integration_jobs_executed = 0L,
  biological_outcomes_computed = FALSE,
  outcome_label_state = "closed",
  stringsAsFactors = FALSE
)
if (summary$samples != 90L || summary$seeds != 5L ||
    summary$entries != 450L || summary$resume_validated_entries != 450L ||
    summary$failed_entries != 0L ||
    summary$normalization_seconds_max > summary$elapsed_cap_seconds ||
    summary$normalization_peak_rss_bytes > summary$rss_cap_bytes ||
    summary$normalization_cache_total_bytes > summary$storage_cap_bytes) {
  stop("Production cache failed its final completion/resource gate.",
       call. = FALSE)
}
write_provenance_csv(summary, summary_path)
write_provenance_csv(projection, projection_path)
message("Validated 450/450 v2 caches and wrote the mandatory reprojection.")
