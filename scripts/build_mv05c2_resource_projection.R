#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 10L) {
  stop(
    "usage: build_mv05c2_resource_projection.R BUILD_AUDIT RESUME_AUDIT ",
    "BUILD_TIME RESUME_TIME CACHED_FOLD_AUDIT CACHED_FOLD_TIME ",
    "CHUNK_RESOURCES ARCHITECTURE OLD_PROJECTION OUTPUT_PREFIX",
    call. = FALSE
  )
}
paths <- vapply(args[1:9], normalizePath, character(1L), winslash = "/",
                mustWork = TRUE)
output_prefix <- args[[10L]]
source("R/provenance_utils.R")

read_elapsed <- function(path) {
  lines <- readLines(path, warn = FALSE)
  value <- grep(
    "Elapsed \\(wall clock\\) time", lines, value = TRUE
  )
  value <- trimws(sub("^.*: ", "", value[[1L]]))
  parts <- as.numeric(strsplit(value, ":", fixed = TRUE)[[1L]])
  if (length(parts) == 3L) {
    parts[[1L]] * 3600 + parts[[2L]] * 60 + parts[[3L]]
  } else if (length(parts) == 2L) {
    parts[[1L]] * 60 + parts[[2L]]
  } else {
    parts[[1L]]
  }
}
read_peak_kb <- function(path) {
  lines <- readLines(path, warn = FALSE)
  value <- grep("Maximum resident set size", lines, value = TRUE)
  as.numeric(trimws(sub("^.*: ", "", value[[1L]])))
}

build <- utils::read.csv(paths[[1L]], stringsAsFactors = FALSE)
resume <- utils::read.csv(paths[[2L]], stringsAsFactors = FALSE)
build_elapsed <- read_elapsed(paths[[3L]])
resume_elapsed <- read_elapsed(paths[[4L]])
fold <- utils::read.csv(paths[[5L]], stringsAsFactors = FALSE)
fold_elapsed <- read_elapsed(paths[[6L]])
chunks <- utils::read.csv(paths[[7L]], stringsAsFactors = FALSE)
architecture <- utils::read.csv(paths[[8L]], stringsAsFactors = FALSE)
old <- utils::read.csv(paths[[9L]], stringsAsFactors = FALSE)

if (nrow(build) != 6L || nrow(resume) != 6L ||
    any(build$disposition != "built_atomic") ||
    any(resume$disposition != "reuse_validated") ||
    any(!fold$exact_numeric_identity) ||
    any(chunks$status != "completed") ||
    nrow(architecture) != 1L || nrow(old) != 1L) {
  stop("MV5-C2 projection inputs are incomplete or invalid.", call. = FALSE)
}

normalization_hours <- build_elapsed * (450 / 6) / 3600
fold_times <- unique(fold[c("fold_id", "cached_fold_seconds")])
seconds_per_sample_fold <- sum(fold_times$cached_fold_seconds) / (2 * 6)
cached_sct_fold_hours <- seconds_per_sample_fold * (75 * 90) / 3600
pair_seconds <- old$projected_landscape_hours_pair_linear[[1L]] * 3600 /
  old$projected_landscape_pair_rows[[1L]]
all_pair_rows <- architecture$query_training_landscape_pair_rows[[1L]]
sct_cell_pair_rows <- all_pair_rows / 3
sct_cell_gene_pair_rows <- all_pair_rows * 2 / 3
cell_landscape_hours <- sct_cell_pair_rows * pair_seconds / 3600
sct_landscape_hours <- sct_cell_gene_pair_rows * pair_seconds / 3600
all_landscape_hours <- all_pair_rows * pair_seconds / 3600
planning_cap_hours <- 24 * 0.9

scenarios <- data.frame(
  contract_id = "mv05c2_resource_projection_v1",
  scenario = c(
    "naive_full_mv05d", "resource_safe_all_planned_views_lower_bound",
    "resource_safe_sct_cell_gene", "resource_safe_sct_cell_primary"
  ),
  normalization_worker_hours = c(
    NA_real_, normalization_hours, normalization_hours, normalization_hours
  ),
  cached_sct_fold_worker_hours = c(
    NA_real_, cached_sct_fold_hours, cached_sct_fold_hours,
    cached_sct_fold_hours
  ),
  landscape_worker_hours = c(
    old$projected_landscape_hours_pair_linear[[1L]], all_landscape_hours,
    sct_landscape_hours, cell_landscape_hours
  ),
  integrated_reference_mapping_worker_hours = c(
    NA_real_, NA_real_, 0, 0
  ),
  projected_lower_bound_worker_hours = c(
    old$projected_worker_hours_sample_linear[[1L]] +
      old$projected_landscape_hours_pair_linear[[1L]],
    normalization_hours + cached_sct_fold_hours + all_landscape_hours,
    normalization_hours + cached_sct_fold_hours + sct_landscape_hours,
    normalization_hours + cached_sct_fold_hours + cell_landscape_hours
  ),
  nominal_cap_hours = 24,
  planning_cap_with_10_percent_reserve_hours = planning_cap_hours,
  cap_passes = c(
    FALSE,
    normalization_hours + cached_sct_fold_hours + all_landscape_hours <=
      planning_cap_hours,
    normalization_hours + cached_sct_fold_hours + sct_landscape_hours <=
      planning_cap_hours,
    normalization_hours + cached_sct_fold_hours + cell_landscape_hours <=
      planning_cap_hours
  ),
  disposition = c(
    "prohibited_measured_projection",
    "prohibited_lower_bound_exceeds_cap_before_integrated_mapping",
    "prohibited_insufficient_reserve_for_orchestration",
    "conditional_go_label_closed_sct_cell_primary_only"
  ),
  outcome_label_state = "closed",
  biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)

resources <- data.frame(
  contract_id = "mv05c2_bounded_resource_evidence_v1",
  normalization_cache_build_seconds_6_samples = build_elapsed,
  normalization_cache_build_peak_rss_bytes = read_peak_kb(paths[[3L]]) * 1024,
  normalization_cache_resume_seconds_6_samples = resume_elapsed,
  normalization_cache_resume_peak_rss_bytes = read_peak_kb(paths[[4L]]) * 1024,
  normalization_cache_build_to_resume_speedup = build_elapsed / resume_elapsed,
  normalization_cache_bytes_6_samples = sum(build$private_cache_size_bytes),
  projected_normalization_cache_bytes_450_entries =
    sum(build$private_cache_size_bytes) * 450 / 6,
  cached_sct_two_fold_process_seconds = fold_elapsed,
  cached_sct_two_fold_peak_rss_bytes = read_peak_kb(paths[[6L]]) * 1024,
  cached_sct_fold_seconds_min = min(fold_times$cached_fold_seconds),
  cached_sct_fold_seconds_max = max(fold_times$cached_fold_seconds),
  chunk_count = nrow(chunks),
  chunk_pairs = sum(chunks$completed_count),
  chunk_seconds_sum = sum(chunks$elapsed_seconds),
  chunk_seconds_max = max(chunks$elapsed_seconds),
  chunk_peak_rss_bytes = max(chunks$peak_process_rss_bytes),
  exact_reference_pairs = 250L,
  exact_reference_maximum_absolute_difference = 0,
  outcome_label_state = "closed",
  biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)

write_provenance_csv(
  resources, paste0(output_prefix, "-bounded-resource-evidence.csv")
)
write_provenance_csv(
  scenarios, paste0(output_prefix, "-scenario-projection.csv")
)
message("Wrote MV5-C2 bounded resource evidence and scenario projection.")
