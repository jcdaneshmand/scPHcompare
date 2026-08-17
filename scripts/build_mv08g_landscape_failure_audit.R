#!/usr/bin/env Rscript

options(warn = 2)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required.")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("usage: build_mv08g_landscape_failure_audit.R PREFREEZE PRIVATE_ROOT OUTPUT FAILED_HEAD")
}
prefreeze <- args[[1L]]; private_root <- args[[2L]]; output <- args[[3L]]
failed_head <- tolower(trimws(args[[4L]]))
head <- tolower(trimws(system2("git", c("rev-parse", "HEAD"), stdout = TRUE)))
if (head != failed_head) stop("MV8-G landscape failure audit HEAD mismatch.")
if (dir.exists(output) && length(list.files(output, all.files = TRUE,
                                             no.. = TRUE))) {
  stop("MV8-G landscape failure audit output must be empty.")
}
source("R/provenance_utils.R")
sha <- function(path) digest::digest(file = path, algo = "sha256",
                                     serialize = FALSE)
contract <- read.csv(file.path(prefreeze, "mv08g-landscape-contract.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)
metric_path <- file.path(private_root, "group-resource-metrics.csv")
job_id <- "production__within475__mv08g__20260807__gene_topology_v1__H1"
stderr_path <- file.path(private_root, "logs", paste0(job_id, "__stderr.txt"))
metric <- read.csv(metric_path, stringsAsFactors = FALSE, check.names = FALSE)
error_text <- paste(readLines(stderr_path, warn = FALSE), collapse = "\n")
expected_fragment <- paste0("cannot open file '", prefreeze,
  "/mv08g-contract.csv': No such file or directory")
group_root <- file.path(private_root, "landscape",
  "mv08g__20260807__gene_topology_v1__H1")
if (nrow(contract) != 1L || contract$accepted_head != failed_head ||
    nrow(metric) != 1L || metric$job_id != job_id ||
    metric$scope != "within475" || metric$seed != 20260807L ||
    metric$view_id != "gene_topology_v1" ||
    metric$homology_dimension != "H1" || metric$disposition != "failed" ||
    metric$exit_status != 1L || metric$component_rows != 7626L ||
    metric$elapsed_seconds >= metric$elapsed_cap_seconds ||
    metric$peak_process_tree_rss_bytes >= metric$rss_cap_bytes ||
    !grepl(expected_fragment, error_text, fixed = TRUE) || dir.exists(group_root)) {
  stop("MV8-G landscape failure evidence differs from the stopped attempt.")
}
dir.create(output, recursive = TRUE, showWarnings = FALSE)
attempt <- data.frame(
  contract_id = "mv08g_landscape_failure_audit_v1",
  failed_head = failed_head, job_id = metric$job_id, scope = metric$scope,
  seed = metric$seed, view_id = metric$view_id,
  homology_dimension = metric$homology_dimension,
  disposition = metric$disposition, exit_status = metric$exit_status,
  elapsed_seconds = metric$elapsed_seconds,
  peak_process_tree_rss_bytes = metric$peak_process_tree_rss_bytes,
  elapsed_cap_seconds = metric$elapsed_cap_seconds,
  rss_cap_bytes = metric$rss_cap_bytes, output_published = FALSE,
  failure_class = "landscape_prefreeze_entry_interface_mismatch",
  observed_missing_artifact = "mv08g-contract.csv",
  corrected_contract_artifact = "mv08g-landscape-contract.csv",
  additional_missing_axis_detected_by_static_audit = TRUE,
  corrected_axis_artifact = "mv08g-sample-seed-axis.csv",
  private_metric_sha256 = sha(metric_path), private_stderr_sha256 = sha(stderr_path),
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE)
decision <- data.frame(
  contract_id = "mv08g_landscape_failure_decision_v1",
  decision = "stop_failed_head_no_retry_require_new_head_and_prefreeze",
  data_failure = FALSE, ph_failure = FALSE, rust_failure = FALSE,
  landscape_definition_failure = FALSE, resource_failure = FALSE,
  implementation_contract_failure = TRUE, retry_in_place_authorized = FALSE,
  accepted_landscape_groups = 0L, comparison_jobs_authorized = 0L,
  hca_fastq_download_authorized = FALSE, raw_reprocessing_authorized = FALSE,
  label_access_authorized = FALSE,
  next_gate = "MV8-G_corrected_landscape_prefreeze",
  stringsAsFactors = FALSE)
paths <- c(
  attempt = file.path(output, "mv08g-landscape-failure.csv"),
  decision = file.path(output, "mv08g-landscape-failure-decision.csv"))
write_provenance_csv(attempt, paths[["attempt"]])
write_provenance_csv(decision, paths[["decision"]])
manifest <- data.frame(
  contract_id = "mv08g_landscape_failure_artifact_manifest_v1",
  file = basename(paths), bytes = as.numeric(file.info(paths)$size),
  sha256 = vapply(paths, sha, character(1L)), contains_expression = FALSE,
  contains_cell_barcode = FALSE, contains_absolute_private_path = FALSE,
  contains_biological_label = FALSE, stringsAsFactors = FALSE)
write_provenance_csv(manifest,
  file.path(output, "mv08g-landscape-failure-artifact-manifest.csv"))
message("MV8-G stopped landscape attempt audited; new head and prefreeze required")
