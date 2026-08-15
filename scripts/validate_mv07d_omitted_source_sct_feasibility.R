#!/usr/bin/env Rscript

options(warn = 2)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required.")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop("usage: validate_mv07d_omitted_source_sct_feasibility.R SENTINEL_CSV ",
       "RESULT_CSV INDIVIDUAL_SOURCE_ROOT PRIVATE_ROOT PROJECTION_CSV OUTPUT_CSV",
       call. = FALSE)
}
sentinel_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
result_path <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
source_root <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
private_root <- normalizePath(args[[4L]], winslash = "/", mustWork = TRUE)
projection_path <- normalizePath(args[[5L]], winslash = "/", mustWork = TRUE)
output_path <- args[[6L]]
if (file.exists(output_path)) stop("Refusing to overwrite: ", output_path)
source("R/provenance_utils.R")
source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/mv05_resource_safe_execution.R")
truth <- function(x) if (is.logical(x)) !is.na(x) & x else
  tolower(trimws(as.character(x))) == "true"
sha <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)

sentinels <- read.csv(sentinel_path, stringsAsFactors = FALSE, check.names = FALSE)
result <- read.csv(result_path, stringsAsFactors = FALSE, check.names = FALSE)
projection <- read.csv(projection_path, stringsAsFactors = FALSE,
                       check.names = FALSE)
axis_ok <- nrow(sentinels) == 6L && nrow(result) == 6L &&
  !anyDuplicated(result$sample_id) && identical(result$sample_id, sentinels$sample_id) &&
  identical(result$tissue, sentinels$tissue) &&
  identical(result$selection_boundary, sentinels$selection_boundary)

paths <- list.files(source_root, pattern = "\\.rds$", recursive = TRUE,
                    full.names = TRUE, ignore.case = TRUE)
ids <- tools::file_path_sans_ext(basename(paths))
source_ok <- vapply(seq_len(nrow(result)), function(i) {
  hit <- paths[ids == result$sample_id[[i]]]
  length(hit) == 1L && basename(hit) == result$individual_source_file[[i]] &&
    identical(sha(hit), result$individual_source_sha256[[i]])
}, logical(1L))

raw_ok <- vapply(seq_len(nrow(result)), function(i) {
  id <- result$sample_id[[i]]
  audit_path <- file.path(private_root, "raw-audit", paste0(id, ".csv"))
  cache_path <- file.path(private_root, "raw", paste0(id, "__raw.rds"))
  if (!file.exists(audit_path) || !file.exists(cache_path)) return(FALSE)
  audit <- read.csv(audit_path, stringsAsFactors = FALSE, check.names = FALSE)
  cache <- readRDS(cache_path)
  valid <- tryCatch({ mv05d0_validate_raw_sample_cache_v2(cache); TRUE },
                    error = function(e) FALSE)
  valid && audit$sample_id == id && audit$cells == result$observed_source_cells[[i]] &&
    audit$genes == result$observed_source_genes[[i]] &&
    ncol(cache$counts) == result$expected_post_qc_cells[[i]]
}, logical(1L))

sct_ok <- vapply(seq_len(nrow(result)), function(i) {
  id <- result$sample_id[[i]]
  audit_path <- file.path(private_root, "sct-audit", paste0(id, ".csv"))
  if (!file.exists(audit_path)) return(FALSE)
  audit <- read.csv(audit_path, stringsAsFactors = FALSE, check.names = FALSE)
  cache_path <- file.path(private_root, "sct", audit$private_cache_file[[1L]])
  if (!file.exists(cache_path)) return(FALSE)
  cache <- readRDS(cache_path)
  matrix <- tryCatch(mv05d0_sct_matrix_from_cache_v1(cache),
                     error = function(e) NULL)
  if (is.null(matrix)) return(FALSE)
  values <- if (inherits(matrix, "sparseMatrix")) matrix@x else as.vector(matrix)
  audit$sample_id == id && audit$selected_cells == 384L && ncol(matrix) == 384L &&
    identical(audit$payload_sha256[[1L]], result$sct_payload_sha256[[i]]) &&
    identical(audit$selected_cell_sha256[[1L]], result$selected_cell_sha256[[i]]) &&
    all(is.finite(values))
}, logical(1L))

resource_ok <- all(is.finite(result$raw_elapsed_seconds)) &&
  all(is.finite(result$sct_elapsed_seconds)) &&
  all(result$raw_elapsed_seconds > 0 &
        result$raw_elapsed_seconds <= result$elapsed_cap_seconds) &&
  all(result$sct_elapsed_seconds > 0 &
        result$sct_elapsed_seconds <= result$elapsed_cap_seconds) &&
  all(result$raw_peak_rss_bytes > 0 &
        result$raw_peak_rss_bytes <= result$rss_cap_bytes) &&
  all(result$sct_peak_rss_bytes > 0 &
        result$sct_peak_rss_bytes <= result$rss_cap_bytes)
scope_ok <- all(truth(result$source_sct_feasible)) &&
  !any(truth(result$ph_computed)) && !any(truth(result$landscape_computed)) &&
  !any(truth(result$biological_outcomes_computed))
depth_ok <- identical(result$expected_post_qc_cells,
                      result$observed_source_cells) &&
  all(result$selected_cells == 384L) && all(result$sct_cells == 384L) &&
  all(result$seed == 20260805L)
hash_ok <- all(grepl("^[0-9a-f]{64}$", result$individual_source_sha256)) &&
  all(grepl("^[0-9a-f]{64}$", result$selected_cell_sha256)) &&
  all(grepl("^[0-9a-f]{64}$", result$sct_payload_sha256))
raw_bytes <- sum(file.info(list.files(file.path(private_root, "raw"),
  pattern = "\\.rds$", full.names = TRUE))$size)
sct_bytes <- sum(file.info(list.files(file.path(private_root, "sct"),
  pattern = "\\.rds$", full.names = TRUE))$size)
expected_projection <- data.frame(
  workload = c("individual_source_to_raw_shard",
    "one_seed_sct_to_five_seed_extension", "combined_upstream_extension"),
  observed_seconds = c(sum(result$raw_elapsed_seconds),
    sum(result$sct_elapsed_seconds),
    sum(result$raw_elapsed_seconds) + sum(result$sct_elapsed_seconds)),
  projected_seconds = c(sum(result$raw_elapsed_seconds) / 6 * 34,
    sum(result$sct_elapsed_seconds) / 6 * 170,
    sum(result$raw_elapsed_seconds) / 6 * 34 +
      sum(result$sct_elapsed_seconds) / 6 * 170),
  observed_storage_bytes = c(raw_bytes, sct_bytes, raw_bytes + sct_bytes),
  projected_storage_bytes = c(raw_bytes / 6 * 34, sct_bytes / 6 * 170,
    raw_bytes / 6 * 34 + sct_bytes / 6 * 170), stringsAsFactors = FALSE)
projection <- projection[match(expected_projection$workload, projection$workload),
                         , drop = FALSE]
projection_ok <- nrow(projection) == 3L && !anyNA(projection$workload) &&
  all(vapply(c("observed_seconds", "projected_seconds",
    "observed_storage_bytes", "projected_storage_bytes"), function(field) {
      isTRUE(all.equal(projection[[field]], expected_projection[[field]],
                       tolerance = 1e-12, check.attributes = FALSE))
    }, logical(1L))) && !any(truth(projection$authorized_by_mv07d))

checks <- data.frame(
  contract_id = "mv07d_omitted_source_sct_independent_validation_v1",
  category = c("sentinel_axis", "individual_source_identity", "raw_cache_identity",
    "sct_cache_identity", "depth_and_selection", "resource_caps",
    "resource_projection", "scientific_scope", "published_hashes"),
  passed = c(axis_ok, all(source_ok), all(raw_ok), all(sct_ok), depth_ok,
             resource_ok, projection_ok, scope_ok, hash_ok),
  detail = c("six prospectively frozen tissue-depth extremes",
    "six private Seurat sources independently rehashed",
    "six v2 raw caches reopened and validated",
    "six v2 SCT caches reopened with finite 384-cell payloads",
    "source counts exact and deterministic selections retained",
    "all child elapsed and peak RSS values below frozen caps",
    "three linear projections independently reconstructed and nonauthorizing",
    "PH landscapes labels and outcomes remain closed",
    "all source selection and payload hashes are well formed"),
  stringsAsFactors = FALSE)
if (!all(checks$passed)) stop("MV7-D feasibility validation failed: ",
  paste(checks$category[!checks$passed], collapse = ", "), call. = FALSE)
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
write.table(checks, output_path, sep = ",", row.names = FALSE,
            col.names = TRUE, quote = TRUE, na = "")
message("MV7-D source/SCT independent validation: 9/9 categories pass")
