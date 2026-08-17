#!/usr/bin/env Rscript

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1", RAYON_NUM_THREADS = "1")
options(warn = 2)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required.")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) {
  stop("usage: run_mv08g_matched_shift_group.R PREFREEZE MV07H_PH_ROOT MV08G_PH_ROOT RUST_LIBRARY GROUP_ID OUTPUT_ROOT EXPECTED_IMPLEMENTATION_ROOT")
}
prefreeze <- args[[1L]]; ph500_root <- args[[2L]]; ph475_root <- args[[3L]]
rust_library <- normalizePath(args[[4L]], winslash = "/", mustWork = TRUE)
group_id <- args[[5L]]; output_root <- args[[6L]]
expected_implementation_root <- args[[7L]]
source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/mv07g_sentinel.R")
source("R/mv07h_full_topology.R")
source("R/mv08g_panel_sensitivity.R")
source("R/landscape_rust_prototype.R")
queue <- read.csv(file.path(prefreeze, "mv08g-matched-shift-queue.csv"),
                  stringsAsFactors = FALSE, check.names = FALSE)
contract <- read.csv(file.path(prefreeze, "mv08g-landscape-contract.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)
row <- queue[queue$group_id == group_id, , drop = FALSE]
sha <- function(path) digest::digest(file = path, algo = "sha256",
                                     serialize = FALSE)
if (nrow(row) != 1L || contract$implementation_root_sha256 !=
      expected_implementation_root || contract$rust_library_sha256 !=
      sha(rust_library)) stop("MV8-G matched-shift identities are stale.")
final_dir <- file.path(output_root, group_id)
if (dir.exists(final_dir)) {
  status <- read.csv(file.path(final_dir, "status.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(status) != 1L || status$completion_state != "complete" ||
      status$distances_sha256 != sha(file.path(final_dir, "distances.csv")) ||
      status$component_rows != 124L) stop("Existing matched-shift group is stale.")
  message("Reused MV8-G matched-shift group: ", group_id)
  quit(save = "no", status = 0L)
}
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
partial <- tempfile(pattern = paste0(group_id, "__partial__"),
                    tmpdir = output_root); dir.create(partial)
sample_ids <- sort(unique(read.csv(file.path(prefreeze,
  "mv08g-sample-seed-axis.csv"), stringsAsFactors = FALSE,
  check.names = FALSE)$sample_id), method = "radix")
started <- proc.time()[["elapsed"]]
result <- vector("list", length(sample_ids))
dimension_number <- as.integer(sub("H", "", row$homology_dimension,
                                   fixed = TRUE))
for (index in seq_along(sample_ids)) {
  sample_id <- sample_ids[[index]]
  path500 <- file.path(ph500_root, paste0("mv07h__", row$seed, "__",
    sample_id, "__", row$view_id, "__ph.rds"))
  path475 <- file.path(ph475_root, paste0("mv08g__", row$seed, "__",
    sample_id, "__", row$view_id, "__ph.rds"))
  record500 <- readRDS(path500); record475 <- readRDS(path475)
  mv07h_validate_ph_record_v1(record500); mv08g_validate_ph_record_v1(record475)
  if (record500$identity$sample_id != sample_id ||
      record475$identity$sample_id != sample_id ||
      record500$identity$seed != row$seed || record475$identity$seed != row$seed ||
      record500$identity$view_id != row$view_id ||
      record475$identity$view_id != row$view_id) {
    stop("MV8-G matched-shift PH identity mismatch.")
  }
  first <- mv07h_finite_intervals_v1(record500, row$homology_dimension)
  second <- mv08g_finite_intervals_v1(record475, row$homology_dimension)
  value <- landscape_rust_prototype_dimension(first, second, dimension_number,
                                               library = rust_library)
  if (!isTRUE(value$rust_used) || value$status != 0L ||
      !is.finite(value$squared_distance) || value$squared_distance < 0) {
    stop("MV8-G matched-shift Rust calculation failed at sample ", sample_id)
  }
  result[[index]] <- data.frame(
    contract_id = "mv08g_matched_panel_landscape_shift_v1", group_id = group_id,
    sample_id = sample_id, seed = as.integer(row$seed), view_id = row$view_id,
    homology_dimension = row$homology_dimension,
    panel500_ph_cache_key = record500$cache_key,
    panel475_ph_cache_key = record475$cache_key,
    panel500_finite_intervals = value$first_finite_intervals,
    panel475_finite_intervals = value$second_finite_intervals,
    squared_distance = value$squared_distance,
    distance = sqrt(value$squared_distance), active_levels = value$active_levels,
    event_segments = value$event_segments, exact = TRUE,
    all_active_levels = TRUE, level_cap_applied = FALSE,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE)
}
result <- do.call(rbind, result)
elapsed <- proc.time()[["elapsed"]] - started
distance_path <- file.path(partial, "distances.csv")
metric_path <- file.path(partial, "metrics.csv")
write.csv(result, distance_path, row.names = FALSE, na = "")
write.csv(data.frame(
  contract_id = "mv08g_matched_shift_metrics_v1", group_id = group_id,
  seed = row$seed, view_id = row$view_id,
  homology_dimension = row$homology_dimension, elapsed_seconds = elapsed,
  samples = nrow(result), component_rows = nrow(result),
  outcome_label_state = "closed", biological_outcomes_computed = FALSE),
  metric_path, row.names = FALSE, na = "")
status <- data.frame(
  contract_id = "mv08g_matched_shift_status_v1", group_id = group_id,
  implementation_root_sha256 = expected_implementation_root,
  rust_library_sha256 = contract$rust_library_sha256,
  completion_state = "complete", distances_sha256 = sha(distance_path),
  distances_bytes = as.numeric(file.info(distance_path)$size),
  metrics_sha256 = sha(metric_path),
  metrics_bytes = as.numeric(file.info(metric_path)$size),
  component_rows = nrow(result), outcome_label_state = "closed",
  biological_outcomes_computed = FALSE, stringsAsFactors = FALSE)
write.csv(status, file.path(partial, "status.csv"), row.names = FALSE, na = "")
if (!file.rename(partial, final_dir)) {
  unlink(partial, recursive = TRUE)
  stop("Failed to atomically publish MV8-G matched-shift group.")
}
message("Completed MV8-G matched-shift group: ", group_id, "; 124 rows.")
