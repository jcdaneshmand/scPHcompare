#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop(
    "usage: audit_mv05d0_raw_shards.R RAW_DIR SAMPLE_MANIFEST ",
    "SOURCE_PREFLIGHT SELECTION_CSV RAW_AUDIT_CSV SEEDS_CSV",
    call. = FALSE
  )
}
raw_dir <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
sample_path <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
preflight_path <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
selection_path <- args[[4L]]
audit_path <- args[[5L]]
seeds <- as.integer(strsplit(args[[6L]], ",", fixed = TRUE)[[1L]])
source("R/provenance_utils.R")
source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/mv05_resource_safe_execution.R")

samples <- utils::read.csv(
  sample_path, stringsAsFactors = FALSE, check.names = FALSE
)
preflight <- utils::read.csv(
  preflight_path, stringsAsFactors = FALSE, check.names = FALSE
)
if (!all(c("sample_id", "outcome_label_state",
           "biological_outcomes_computed") %in% names(samples)) ||
    any(c("tissue", "approach") %in% names(samples)) ||
    any(samples$outcome_label_state != "closed") ||
    any(as.logical(samples$biological_outcomes_computed)) ||
    anyNA(seeds) || !length(seeds) || nrow(preflight) != 1L) {
  stop("Raw-shard audit inputs violate the label-closed contract.",
       call. = FALSE)
}
source_sha <- as.character(preflight$source_sha256[[1L]])
raw_samples <- list()
rows <- list()
for (sample_id in sort(samples$sample_id, method = "radix")) {
  file <- paste0(gsub("[^A-Za-z0-9_.-]", "_", sample_id), "__raw.rds")
  path <- file.path(raw_dir, file)
  if (!file.exists(path)) stop("Missing raw shard: ", file, call. = FALSE)
  record <- readRDS(path)
  if (!identical(record$contract_id, "mv05d0_raw_sample_cache_v1") ||
      !identical(record$sample_id, sample_id) ||
      !identical(record$historical_source_sha256, source_sha) ||
      !identical(record$outcome_label_state, "closed") ||
      !identical(record$biological_outcomes_computed, FALSE)) {
    stop("Raw shard identity failed: ", sample_id, call. = FALSE)
  }
  raw_samples[[sample_id]] <- record$counts
  rows[[length(rows) + 1L]] <- data.frame(
    contract_id = "mv05d0_raw_sample_cache_audit_v1",
    sample_id = sample_id, genes = nrow(record$counts),
    cells = ncol(record$counts), historical_source_sha256 = source_sha,
    private_raw_cache_file = file,
    private_raw_cache_size_bytes = unname(file.info(path)$size),
    private_raw_cache_sha256 = digest::digest(
      file = path, algo = "sha256", serialize = FALSE
    ),
    disposition = "recovered_validated_after_guard_stop",
    outcome_label_state = "closed",
    biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
  rm(record)
  invisible(gc())
}
selection <- mv05d0_build_selection_summary_v1(
  raw_samples, samples$sample_id, seeds = seeds, n = 384L
)
write_provenance_csv(selection, selection_path)
write_provenance_csv(do.call(rbind, rows), audit_path)
message("Audited ", length(raw_samples), " recovered raw shards and ",
        nrow(selection), " sample-seed selections.")
