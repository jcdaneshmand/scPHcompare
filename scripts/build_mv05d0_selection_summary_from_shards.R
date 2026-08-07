#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(
    "usage: build_mv05d0_selection_summary_from_shards.R RAW_DIR ",
    "CANDIDATES OUTPUT_CSV", call. = FALSE
  )
}
raw_dir <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
candidate_path <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
output_path <- args[[3L]]
source("R/provenance_utils.R")
source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/mv05_resource_safe_execution.R")
candidates <- utils::read.csv(
  candidate_path, stringsAsFactors = FALSE, check.names = FALSE
)
if (nrow(candidates) != 90L ||
    any(c("tissue", "approach") %in% names(candidates)) ||
    any(candidates$outcome_label_state != "closed") ||
    any(as.logical(candidates$biological_outcomes_computed))) {
  stop("Candidate manifest violates the frozen label boundary.",
       call. = FALSE)
}
rows <- list()
for (sample_id in sort(candidates$sample_id, method = "radix")) {
  path <- file.path(
    raw_dir, paste0(gsub("[^A-Za-z0-9_.-]", "_", sample_id), "__raw.rds")
  )
  record <- readRDS(path)
  mv05d0_validate_raw_sample_cache_v2(record)
  if (!identical(record$sample_id, sample_id)) {
    stop("Raw shard sample identity mismatch.", call. = FALSE)
  }
  for (seed in 20260805:20260809) {
    selected <- select_matched_cells(
      colnames(record$counts), n = 384L, seed = seed
    )
    rows[[length(rows) + 1L]] <- data.frame(
      contract_id = "mv05d0_matched_cell_selection_summary_v2",
      sample_id = sample_id, seed = seed,
      eligible_cells = ncol(record$counts), selected_cells = 384L,
      selected_cell_sha256 = attr(selected, "selected_cell_sha256"),
      raw_counts_sha256 = record$counts_sha256,
      outcome_label_state = "closed",
      biological_outcomes_computed = FALSE,
      stringsAsFactors = FALSE
    )
  }
  rm(record)
  invisible(gc())
}
result <- do.call(rbind, rows)
result <- result[order(result$seed, result$sample_id, method = "radix"),
                 , drop = FALSE]
rownames(result) <- NULL
if (nrow(result) != 450L ||
    anyDuplicated(paste(result$sample_id, result$seed, sep = "\r"))) {
  stop("Selection summary did not produce 450 unique entries.", call. = FALSE)
}
write_provenance_csv(result, output_path)
message("Frozen 450 deterministic selections from 90 v2 raw shards.")
