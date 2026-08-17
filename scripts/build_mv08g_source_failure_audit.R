#!/usr/bin/env Rscript

options(warn = 2)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required.")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("usage: build_mv08g_source_failure_audit.R PREFREEZE PRIVATE_ROOT OUTPUT EXPECTED_HEAD")
}
prefreeze <- args[[1L]]; private_root <- args[[2L]]; output <- args[[3L]]
expected_head <- tolower(trimws(args[[4L]]))
head <- tolower(trimws(system2("git", c("rev-parse", "HEAD"), stdout = TRUE)))
if (head != expected_head) stop("MV8-G failure audit exact HEAD mismatch.")
if (dir.exists(output) && length(list.files(output, all.files = TRUE,
                                             no.. = TRUE))) {
  stop("MV8-G failure audit output must be empty.")
}
source("R/provenance_utils.R")
sha <- function(path) digest::digest(file = path, algo = "sha256",
                                     serialize = FALSE)
contract <- read.csv(file.path(prefreeze, "mv08g-contract.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)
metric_path <- file.path(private_root, "source-metrics.csv")
stderr_path <- file.path(private_root, "logs",
  "production__source__20260805__stderr.txt")
metric <- read.csv(metric_path, stringsAsFactors = FALSE, check.names = FALSE)
error_text <- paste(readLines(stderr_path, warn = FALSE), collapse = "\n")
expected_error <- paste(c(
  "Error: x must be genes_by_cells with shape 500 x 384; observed 475 x 384.",
  "Execution halted"), collapse = "\n")
if (nrow(contract) != 1L || contract$accepted_head != expected_head ||
    nrow(metric) != 1L || metric$job_id != "production__source__20260805" ||
    metric$seed != 20260805L || metric$disposition != "failed" ||
    metric$exit_status != 1L || metric$retries != 0L ||
    metric$elapsed_seconds >= metric$elapsed_cap_seconds ||
    metric$peak_process_tree_rss_bytes >= metric$rss_cap_bytes ||
    !is.na(metric$output_sha256) || !identical(error_text, expected_error) ||
    file.exists(file.path(private_root, "source", metric$output_file))) {
  stop("MV8-G source failure evidence differs from the stopped v1 attempt.")
}
dir.create(output, recursive = TRUE, showWarnings = FALSE)
attempt <- data.frame(
  contract_id = "mv08g_source_failure_audit_v1",
  accepted_head = expected_head, job_id = metric$job_id, seed = metric$seed,
  disposition = metric$disposition, exit_status = metric$exit_status,
  elapsed_seconds = metric$elapsed_seconds,
  peak_process_tree_rss_bytes = metric$peak_process_tree_rss_bytes,
  elapsed_cap_seconds = metric$elapsed_cap_seconds,
  rss_cap_bytes = metric$rss_cap_bytes, workers = metric$workers,
  retries = metric$retries, output_published = FALSE,
  failure_class = "scientific_profile_shape_contract_mismatch",
  required_shape = "475x384", rejected_against_shape = "500x384",
  private_metric_sha256 = sha(metric_path),
  private_stderr_sha256 = sha(stderr_path), outcome_label_state = "closed",
  biological_outcomes_computed = FALSE, stringsAsFactors = FALSE)
decision <- data.frame(
  contract_id = "mv08g_source_failure_decision_v1",
  decision = "stop_v1_no_retry_require_new_head_and_prefreeze",
  data_failure = FALSE, resource_failure = FALSE,
  implementation_contract_failure = TRUE, retry_in_place_authorized = FALSE,
  ph_jobs_authorized = 0L, landscape_jobs_authorized = 0L,
  comparison_jobs_authorized = 0L, hca_fastq_download_authorized = FALSE,
  raw_reprocessing_authorized = FALSE, label_access_authorized = FALSE,
  next_gate = "MV8-G_v2_scientific_common475_prefreeze",
  stringsAsFactors = FALSE)
paths <- c(
  attempt = file.path(output, "mv08g-source-failure.csv"),
  decision = file.path(output, "mv08g-source-failure-decision.csv"))
write_provenance_csv(attempt, paths[["attempt"]])
write_provenance_csv(decision, paths[["decision"]])
manifest <- data.frame(
  contract_id = "mv08g_source_failure_artifact_manifest_v1",
  file = basename(paths), bytes = as.numeric(file.info(paths)$size),
  sha256 = vapply(paths, sha, character(1L)), contains_expression = FALSE,
  contains_cell_barcode = FALSE, contains_absolute_private_path = FALSE,
  contains_biological_label = FALSE, stringsAsFactors = FALSE)
write_provenance_csv(manifest,
  file.path(output, "mv08g-source-failure-artifact-manifest.csv"))
message("MV8-G stopped source attempt audited; no retry authorized")
