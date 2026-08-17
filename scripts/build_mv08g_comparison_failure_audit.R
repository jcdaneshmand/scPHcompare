#!/usr/bin/env Rscript

options(warn = 2)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required.")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop(paste("usage: build_mv08g_comparison_failure_audit.R",
    "COMPARISON_PREFREEZE PRIVATE_ROOT FAILED_PUBLIC_OUTPUT OUTPUT",
    "FAILED_HEAD RESULT_SUBDIR"))
}
prefreeze <- args[[1L]]; private_root <- args[[2L]]
failed_public <- args[[3L]]; output <- args[[4L]]
failed_head <- tolower(trimws(args[[5L]])); result_subdir <- args[[6L]]
head <- tolower(trimws(system2("git", c("rev-parse", "HEAD"), stdout = TRUE)))
if (head != failed_head) stop("MV8-G comparison failure audit HEAD mismatch.")
if (dir.exists(output) && length(list.files(output, all.files = TRUE,
                                             no.. = TRUE))) {
  stop("MV8-G comparison failure audit output must be empty.")
}
source("R/provenance_utils.R")
sha <- function(path) digest::digest(file = path, algo = "sha256",
                                     serialize = FALSE)
contract_path <- file.path(prefreeze, "mv08g-comparison-contract.csv")
decision_path <- file.path(prefreeze, "mv08g-comparison-decision.csv")
metric_path <- file.path(private_root, "comparison-resource-metric.csv")
stderr_path <- file.path(private_root, "logs", "comparison__stderr.txt")
result_path <- file.path(private_root, result_subdir)
run_path <- "scripts/run_mv08g_comparison.R"
clustering_path <- "R/mv05n_clustering_gate.R"
required <- c(contract_path, decision_path, metric_path, stderr_path,
              run_path, clustering_path)
if (any(!file.exists(required))) stop("MV8-G comparison failure evidence incomplete.")
contract <- read.csv(contract_path, stringsAsFactors = FALSE, check.names = FALSE)
decision <- read.csv(decision_path, stringsAsFactors = FALSE, check.names = FALSE)
metric <- read.csv(metric_path, stringsAsFactors = FALSE, check.names = FALSE)
stderr <- paste(readLines(stderr_path, warn = FALSE), collapse = "\n")
run_text <- paste(readLines(run_path, warn = FALSE), collapse = "\n")
clustering_text <- paste(readLines(clustering_path, warn = FALSE), collapse = "\n")
public_empty <- dir.exists(failed_public) &&
  !length(list.files(failed_public, all.files = TRUE, no.. = TRUE))
result_empty <- dir.exists(result_path) &&
  !length(list.files(result_path, all.files = TRUE, no.. = TRUE))
if (nrow(contract) != 1L || contract$accepted_head != failed_head ||
    decision$decision != "authorize_one_label_closed_comparison_job" ||
    nrow(metric) != 1L || metric$disposition != "failed" ||
    metric$exit_status != 1L || metric$elapsed_seconds > metric$elapsed_cap_seconds ||
    metric$peak_process_tree_rss_bytes > metric$rss_cap_bytes ||
    metric$result_bytes != 0 || !is.na(metric$result_manifest_sha256) ||
    !grepl("numbers of columns of arguments do not match", stderr, fixed = TRUE) ||
    !grepl("fixed_partitions <- do.call(rbind, fixed_partition_rows)",
           run_text, fixed = TRUE) ||
    !grepl("is_medoid = names(clusters) %in% medoids", clustering_text,
           fixed = TRUE) || !public_empty || !result_empty) {
  stop("MV8-G comparison failure evidence differs from the stopped job.")
}
dir.create(output, recursive = TRUE, showWarnings = FALSE)
failure <- data.frame(
  contract_id = "mv08g_comparison_failure_audit_v1",
  failed_head = failed_head,
  failure_class = "heterogeneous_fixed_partition_schema",
  failed_phase = "fixed_k_partition_publication",
  pam_only_column = "is_medoid", average_linkage_column_present = FALSE,
  result_mismatch_observed = FALSE, topology_validity_changed = FALSE,
  clustering_fit_failure = FALSE, resource_cap_failure = FALSE,
  accepted_result_files = 0L, retry_in_place_authorized = FALSE,
  comparison_prefreeze_contract_sha256 = sha(contract_path),
  comparison_prefreeze_decision_sha256 = sha(decision_path),
  private_resource_metric_sha256 = sha(metric_path),
  private_stderr_sha256 = sha(stderr_path), failed_runner_sha256 = sha(run_path),
  clustering_helper_sha256 = sha(clustering_path),
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE)
failure_decision <- data.frame(
  contract_id = "mv08g_comparison_failure_decision_v1",
  decision = "require_common_fixed_partition_schema_and_new_exact_prefreeze",
  comparison_jobs_authorized = 0L, hca_expression_jobs_authorized = 0L,
  hca_fastq_download_authorized = FALSE, raw_reprocessing_authorized = FALSE,
  label_access_authorized = FALSE,
  next_gate = "MV8-G_comparison_schema_repair_and_prefreeze",
  stringsAsFactors = FALSE)
paths <- c(
  failure = file.path(output, "mv08g-comparison-failure.csv"),
  decision = file.path(output, "mv08g-comparison-failure-decision.csv"))
write_provenance_csv(failure, paths[["failure"]])
write_provenance_csv(failure_decision, paths[["decision"]])
manifest <- data.frame(
  contract_id = "mv08g_comparison_failure_artifact_manifest_v1",
  file = basename(paths), bytes = as.numeric(file.info(paths)$size),
  sha256 = vapply(paths, sha, character(1L)), contains_expression = FALSE,
  contains_cell_barcode = FALSE, contains_absolute_private_path = FALSE,
  contains_biological_label = FALSE, stringsAsFactors = FALSE)
write_provenance_csv(manifest,
  file.path(output, "mv08g-comparison-failure-artifact-manifest.csv"))
message("MV8-G comparison schema failure audited; no result accepted")
