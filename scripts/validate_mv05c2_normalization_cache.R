#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(
    "usage: validate_mv05c2_normalization_cache.R CACHE_DIR MV05C_JOB OUTPUT_CSV",
    call. = FALSE
  )
}
cache_dir <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
job_path <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
output_path <- args[[3L]]

source("R/provenance_utils.R")
source("R/mv05_resource_safe_execution.R")

job <- readRDS(job_path)
if (!identical(job$contract_id, "mv05c_existing_data_job_v1") ||
    !identical(job$outcome_label_state, "closed") ||
    !identical(job$biological_outcomes_computed, FALSE) ||
    is.null(job$sct_fold$panel) || is.null(job$sct_fold$prepared$selected)) {
  stop("MV5-C reference job is invalid for cache equivalence.", call. = FALSE)
}
panel <- job$sct_fold$panel
seed <- as.integer(job$seed)
rows <- list()
for (sample_id in sort(names(job$sct_fold$prepared$selected), method = "radix")) {
  safe_id <- gsub("[^A-Za-z0-9_.-]", "_", sample_id)
  path <- file.path(cache_dir, paste0(safe_id, "__", seed, "__sct.rds"))
  record <- readRDS(path)
  mv05c2_validate_normalization_cache_record_v1(record)
  data <- Seurat::GetAssayData(
    record$payload$sct_object, assay = "SCT", layer = "data"
  )
  observed <- as.matrix(data[panel$feature_id, , drop = FALSE])
  rownames(observed) <- panel$gene
  expected <- job$sct_fold$prepared$selected[[sample_id]]
  if (!identical(dim(observed), dim(expected)) ||
      !identical(dimnames(observed), dimnames(expected))) {
    stop("Cache/reference axes differ for ", sample_id, call. = FALSE)
  }
  difference <- observed - expected
  rows[[length(rows) + 1L]] <- data.frame(
    contract_id = "mv05c2_normalization_equivalence_v1",
    reference_job_id = job$job_id,
    sample_id = sample_id, seed = seed,
    compared_genes = nrow(observed), compared_cells = ncol(observed),
    maximum_absolute_difference = max(abs(difference)),
    exact_numeric_identity = identical(observed, expected),
    observed_sha256 = digest::digest(
      observed, algo = "sha256", serialize = TRUE
    ),
    expected_sha256 = digest::digest(
      expected, algo = "sha256", serialize = TRUE
    ),
    normalization_cache_key = record$cache_key,
    outcome_label_state = "closed",
    biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
}
result <- do.call(rbind, rows)
if (any(result$maximum_absolute_difference != 0) ||
    !all(result$exact_numeric_identity)) {
  stop("One or more cached normalization results differ from MV5-C.",
       call. = FALSE)
}
write_provenance_csv(result, output_path)
message("MV5-C2 normalization cache is exactly equivalent to MV5-C.")
