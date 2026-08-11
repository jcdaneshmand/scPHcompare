#!/usr/bin/env Rscript

options(warn = 2, digits = 17)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("usage: snapshot_mv05p_units.R PRODUCTION_ROOT SNAPSHOT_OUTPUT",
       call. = FALSE)
}
if (!requireNamespace("digest", quietly = TRUE)) {
  stop("digest is required for MV5-P immutable snapshots.", call. = FALSE)
}
source("R/provenance_utils.R")
production_root <- normalizePath(args[[1L]], mustWork = TRUE)
snapshot_output <- args[[2L]]
safe_name <- function(value) gsub("[^A-Za-z0-9_.-]", "_", value)
read_one <- function(path) utils::read.csv(
  normalizePath(path, mustWork = TRUE), stringsAsFactors = FALSE,
  check.names = FALSE
)
file_sha <- function(path) digest::digest(file = path, algo = "sha256",
                                          serialize = FALSE)
landscape_queue <- read_one(
  "docs/audits/mv05o-landscape-chunk-queue-2026-08-10.csv")
baseline_queue <- read_one(
  "docs/audits/mv05o-baseline-group-queue-2026-08-10.csv")
units <- rbind(
  data.frame(unit_family = "landscape",
             production_group_id = landscape_queue$production_group_id,
             unit_id = landscape_queue$production_chunk_id,
             output_subdir = "landscape-output",
             status_subdir = "landscape-status", stringsAsFactors = FALSE),
  data.frame(unit_family = "baseline",
             production_group_id = baseline_queue$production_group_id,
             unit_id = baseline_queue$baseline_unit_id,
             output_subdir = "baseline-output",
             status_subdir = "baseline-status", stringsAsFactors = FALSE)
)
if (nrow(units) != 4565L || anyDuplicated(units$unit_id)) {
  stop("MV5-P immutable snapshot queue accounting failed.", call. = FALSE)
}
rows <- lapply(seq_len(nrow(units)), function(index) {
  unit <- units[index, , drop = FALSE]
  root <- file.path(production_root, "groups",
                    safe_name(unit$production_group_id))
  stem <- safe_name(unit$unit_id)
  output <- file.path(root, unit$output_subdir, paste0(stem, ".csv"))
  status <- file.path(root, unit$status_subdir,
                      paste0(stem, "__status.csv"))
  info_output <- file.info(normalizePath(output, mustWork = TRUE))
  info_status <- file.info(normalizePath(status, mustWork = TRUE))
  data.frame(
    contract_id = "mv05p_immutable_unit_snapshot_v1",
    unit_family = unit$unit_family,
    production_group_id = unit$production_group_id,
    unit_id = unit$unit_id,
    output_sha256 = file_sha(output), output_size_bytes = info_output$size,
    output_mtime_epoch_seconds = as.numeric(info_output$mtime),
    status_sha256 = file_sha(status), status_size_bytes = info_status$size,
    status_mtime_epoch_seconds = as.numeric(info_status$mtime),
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    clustering_jobs_executed = 0L, stringsAsFactors = FALSE
  )
})
snapshot <- do.call(rbind, rows)
write_provenance_csv(snapshot, snapshot_output)
message("Snapshotted 4,565 MV5-P output/status unit pairs.")

