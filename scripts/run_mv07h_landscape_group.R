#!/usr/bin/env Rscript

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1", RAYON_NUM_THREADS = "1")
options(warn = 2)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required.")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop("usage: run_mv07h_landscape_group.R PREFREEZE PH_ROOT RUST_LIBRARY GROUP_ID OUTPUT_ROOT EXPECTED_IMPLEMENTATION_ROOT")
}
prefreeze <- args[[1L]]
ph_root <- args[[2L]]
rust_library <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
group_id <- args[[4L]]
output_root <- args[[5L]]
expected_implementation_root <- args[[6L]]
source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/mv07g_sentinel.R")
source("R/mv07h_full_topology.R")
source("R/landscape_rust_prototype.R")
queue <- read.csv(file.path(prefreeze, "mv07h-landscape-queue.csv"),
                  stringsAsFactors = FALSE, check.names = FALSE)
contract <- read.csv(file.path(prefreeze, "mv07h-contract.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)
row <- queue[queue$group_id == group_id, , drop = FALSE]
if (nrow(row) != 1L || contract$implementation_root_sha256 !=
      expected_implementation_root ||
    contract$rust_library_sha256 != .mv07h_sha256(rust_library)) {
  stop("MV7-H landscape group execution identities are stale.")
}
safe <- gsub(":", "_", group_id, fixed = TRUE)
final_dir <- file.path(output_root, safe)
if (dir.exists(final_dir)) {
  status <- read.csv(file.path(final_dir, "status.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)
  distance_path <- file.path(final_dir, "distances.csv")
  metric_path <- file.path(final_dir, "metrics.csv")
  if (nrow(status) != 1L || status$completion_state != "complete" ||
      status$group_id != group_id ||
      status$distances_sha256 != .mv07h_sha256(distance_path) ||
      status$metrics_sha256 != .mv07h_sha256(metric_path) ||
      status$component_rows != row$component_rows) {
    stop("Existing MV7-H landscape group is incomplete or stale.")
  }
  message("Reused MV7-H landscape group: ", group_id)
  quit(save = "no", status = 0L)
}
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
partial <- tempfile(pattern = paste0(safe, "__partial__"),
                    tmpdir = output_root)
dir.create(partial)
started <- proc.time()[["elapsed"]]
sample_ids <- sort(unique(read.csv(
  file.path(prefreeze, "mv07h-sample-seed-axis.csv"),
  stringsAsFactors = FALSE, check.names = FALSE)$sample_id), method = "radix")
if (length(sample_ids) != 124L) stop("MV7-H sample axis is incomplete.")
intervals <- vector("list", length(sample_ids))
names(intervals) <- sample_ids
ph_keys <- character(length(sample_ids))
for (index in seq_along(sample_ids)) {
  sample_id <- sample_ids[[index]]
  path <- file.path(ph_root, paste0(
    "mv07h__", row$seed, "__", sample_id, "__", row$view_id, "__ph.rds"))
  if (!file.exists(path)) stop("Missing MV7-H PH record: ", basename(path))
  record <- readRDS(path)
  mv07h_validate_ph_record_v1(record)
  if (record$identity$seed != row$seed ||
      record$identity$sample_id != sample_id ||
      record$identity$view_id != row$view_id) {
    stop("MV7-H PH record belongs to another landscape group.")
  }
  intervals[[index]] <- mv07h_finite_intervals_v1(
    record, row$homology_dimension)
  ph_keys[[index]] <- record$cache_key
}
pairs <- utils::combn(sample_ids, 2L)
count <- ncol(pairs)
if (count != 7626L) stop("MV7-H unordered pair count is not 7,626.")
result <- data.frame(
  contract_id = rep("mv07h_exact_landscape_distance_v1", count),
  engine_id = rep("rust_scph_landscape_kernel_v1", count),
  group_id = rep(group_id, count), pair_ordinal = seq_len(count),
  pair_id = character(count), seed = rep(as.integer(row$seed), count),
  first_sample_id = pairs[1L, ], second_sample_id = pairs[2L, ],
  view_id = rep(row$view_id, count),
  homology_dimension = rep(row$homology_dimension, count),
  first_ph_cache_key = character(count), second_ph_cache_key = character(count),
  first_finite_intervals = integer(count), second_finite_intervals = integer(count),
  squared_distance = numeric(count), distance = numeric(count),
  active_levels = integer(count), event_segments = integer(count),
  exact = TRUE, all_active_levels = TRUE, level_cap_applied = FALSE,
  rust_status = integer(count), rust_engine_version = integer(count),
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  clustering_jobs = 0L, label_jobs = 0L, outcome_jobs = 0L,
  stringsAsFactors = FALSE
)
dimension_number <- as.integer(sub("H", "", row$homology_dimension,
                                   fixed = TRUE))
for (index in seq_len(count)) {
  first_id <- pairs[1L, index]
  second_id <- pairs[2L, index]
  value <- landscape_rust_prototype_dimension(
    intervals[[first_id]], intervals[[second_id]], dimension_number,
    library = rust_library
  )
  if (!isTRUE(value$rust_used) || value$status != 0L ||
      !is.finite(value$squared_distance) || value$squared_distance < 0) {
    stop("MV7-H Rust landscape calculation failed closed at pair ", index)
  }
  result$pair_id[[index]] <- mv07h_pair_id_v1(
    row$seed, first_id, second_id, row$view_id, row$homology_dimension)
  result$first_ph_cache_key[[index]] <- ph_keys[[first_id]]
  result$second_ph_cache_key[[index]] <- ph_keys[[second_id]]
  result$first_finite_intervals[[index]] <- value$first_finite_intervals
  result$second_finite_intervals[[index]] <- value$second_finite_intervals
  result$squared_distance[[index]] <- value$squared_distance
  result$distance[[index]] <- sqrt(value$squared_distance)
  result$active_levels[[index]] <- as.integer(value$active_levels)
  result$event_segments[[index]] <- as.integer(value$event_segments)
  result$rust_status[[index]] <- value$status
  result$rust_engine_version[[index]] <- value$engine_version
}
if (anyDuplicated(result$pair_id) || any(!is.finite(result$distance)) ||
    any(result$squared_distance < 0)) {
  stop("MV7-H landscape output is incomplete or invalid.")
}
elapsed <- proc.time()[["elapsed"]] - started
status_lines <- if (file.exists("/proc/self/status")) {
  readLines("/proc/self/status", warn = FALSE)
} else character()
peak_line <- grep("^VmHWM:", status_lines, value = TRUE)
peak_rss <- if (length(peak_line) == 1L) {
  as.numeric(gsub("[^0-9]", "", peak_line)) * 1024
} else NA_real_
metrics <- data.frame(
  contract_id = "mv07h_landscape_group_metrics_v1", group_id = group_id,
  seed = row$seed, view_id = row$view_id,
  homology_dimension = row$homology_dimension,
  elapsed_seconds = elapsed, peak_self_rss_bytes = peak_rss,
  samples = 124L, component_rows = nrow(result),
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
distance_path <- file.path(partial, "distances.csv")
metric_path <- file.path(partial, "metrics.csv")
write.csv(result, distance_path, row.names = FALSE, na = "")
write.csv(metrics, metric_path, row.names = FALSE, na = "")
status <- data.frame(
  contract_id = "mv07h_landscape_group_status_v1", group_id = group_id,
  implementation_root_sha256 = expected_implementation_root,
  rust_library_sha256 = contract$rust_library_sha256,
  completion_state = "complete", distances_sha256 = .mv07h_sha256(distance_path),
  distances_bytes = as.numeric(file.info(distance_path)$size),
  metrics_sha256 = .mv07h_sha256(metric_path),
  metrics_bytes = as.numeric(file.info(metric_path)$size),
  component_rows = nrow(result), outcome_label_state = "closed",
  biological_outcomes_computed = FALSE, clustering_jobs = 0L,
  label_jobs = 0L, outcome_jobs = 0L, stringsAsFactors = FALSE
)
write.csv(status, file.path(partial, "status.csv"), row.names = FALSE, na = "")
if (!file.rename(partial, final_dir)) {
  unlink(partial, recursive = TRUE)
  stop("Failed to atomically publish MV7-H landscape group.")
}
message("Completed MV7-H landscape group: ", group_id, "; ", nrow(result),
        " component rows.")
