#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop(
    "usage: build_mv05d0_frozen_inputs.R RAW_LIST SOURCE_PREFLIGHT ",
    "CANDIDATES RAW_SAMPLE_DIR SELECTION_CSV RAW_AUDIT_CSV",
    call. = FALSE
  )
}
if (!identical(Sys.getenv("MV05D0_ALLOW_HIGH_MEMORY_SOURCE_MIGRATION"),
               "YES_OWNER_APPROVED")) {
  stop(
    "Monolithic source migration exceeds the normal 8-GiB guard. Set ",
    "MV05D0_ALLOW_HIGH_MEMORY_SOURCE_MIGRATION=YES_OWNER_APPROVED only ",
    "after explicit project-owner approval and separate monitoring.",
    call. = FALSE
  )
}

raw_list_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
source_preflight_path <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
candidate_path <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
raw_sample_dir <- args[[4L]]
selection_path <- args[[5L]]
raw_audit_path <- args[[6L]]
dir.create(raw_sample_dir, recursive = TRUE, showWarnings = FALSE)

source("R/provenance_utils.R")
source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/mv05_resource_safe_execution.R")

preflight <- utils::read.csv(
  source_preflight_path, stringsAsFactors = FALSE, check.names = FALSE
)
candidates <- utils::read.csv(
  candidate_path, stringsAsFactors = FALSE, check.names = FALSE
)
if (nrow(preflight) != 1L || nrow(candidates) != 90L ||
    !all(c("sample_id", "outcome_label_state",
           "biological_outcomes_computed") %in% names(candidates)) ||
    any(c("tissue", "approach") %in% names(candidates)) ||
    any(candidates$outcome_label_state != "closed") ||
    any(as.logical(candidates$biological_outcomes_computed)) ||
    anyDuplicated(candidates$sample_id)) {
  stop("Frozen candidates or source preflight violate MV5-D0 boundaries.",
       call. = FALSE)
}
source_size <- unname(file.info(raw_list_path)$size)
source_sha <- as.character(preflight$source_sha256[[1L]])
if (!identical(source_size, as.numeric(preflight$source_size_bytes[[1L]])) ||
    !grepl("^[0-9a-f]{64}$", source_sha)) {
  stop("Historical source does not match the frozen preflight.", call. = FALSE)
}

started <- proc.time()[["elapsed"]]
raw_list <- readRDS(raw_list_path)
if (!is.list(raw_list) || is.null(names(raw_list)) ||
    anyDuplicated(names(raw_list)) ||
    !all(candidates$sample_id %in% names(raw_list))) {
  stop("Historical raw list is missing frozen candidate samples.", call. = FALSE)
}
raw_samples <- raw_list[candidates$sample_id]
rm(raw_list)
invisible(gc())

rows <- list()
for (sample_id in sort(names(raw_samples), method = "radix")) {
  value <- raw_samples[[sample_id]]
  if ((!is.matrix(value) && !inherits(value, "Matrix")) ||
      is.null(rownames(value)) || is.null(colnames(value)) ||
      anyDuplicated(rownames(value)) || anyDuplicated(colnames(value)) ||
      ncol(value) < 384L) {
    stop("Raw candidate failed identity/shape validation: ", sample_id,
         call. = FALSE)
  }
  safe_id <- gsub("[^A-Za-z0-9_.-]", "_", sample_id)
  cache_file <- paste0(safe_id, "__raw.rds")
  cache_path <- file.path(raw_sample_dir, cache_file)
  record <- list(
    contract_id = "mv05d0_raw_sample_cache_v1",
    sample_id = sample_id,
    historical_source_sha256 = source_sha,
    counts = value,
    outcome_label_state = "closed",
    biological_outcomes_computed = FALSE
  )
  disposition <- "built_atomic"
  if (file.exists(cache_path)) {
    existing <- readRDS(cache_path)
    if (!identical(existing$contract_id, record$contract_id) ||
        !identical(existing$sample_id, record$sample_id) ||
        !identical(existing$historical_source_sha256,
                   record$historical_source_sha256) ||
        !identical(existing$counts, record$counts) ||
        !identical(existing$outcome_label_state, "closed") ||
        !identical(existing$biological_outcomes_computed, FALSE)) {
      stop("Existing raw sample cache is stale; refusing overwrite: ",
           sample_id, call. = FALSE)
    }
    disposition <- "reuse_validated"
  } else {
    partial <- tempfile(pattern = paste0(cache_file, "."),
                        tmpdir = raw_sample_dir)
    saveRDS(record, partial, compress = "gzip", version = 3)
    if (!file.rename(partial, cache_path)) {
      unlink(partial)
      stop("Failed to atomically publish raw sample cache: ", sample_id,
           call. = FALSE)
    }
  }
  rows[[length(rows) + 1L]] <- data.frame(
    contract_id = "mv05d0_raw_sample_cache_audit_v1",
    sample_id = sample_id,
    genes = nrow(value), cells = ncol(value),
    historical_source_sha256 = source_sha,
    private_raw_cache_file = cache_file,
    private_raw_cache_size_bytes = unname(file.info(cache_path)$size),
    private_raw_cache_sha256 = digest::digest(
      file = cache_path, algo = "sha256", serialize = FALSE
    ),
    disposition = disposition,
    outcome_label_state = "closed",
    biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
  message("Frozen raw sample ", length(rows), "/", length(raw_samples),
          ": ", sample_id, " (", disposition, ")")
}
selection <- mv05d0_build_selection_summary_v1(
  raw_samples, candidates$sample_id, seeds = 20260805:20260809, n = 384L
)
raw_audit <- do.call(rbind, rows)
raw_audit$source_read_seconds <- proc.time()[["elapsed"]] - started
write_provenance_csv(selection, selection_path)
write_provenance_csv(raw_audit, raw_audit_path)
message("Frozen 90 raw sample caches and 450 matched-cell identities.")
