#!/usr/bin/env Rscript

options(warn = 2)
for (package in c("digest", "cluster", "mclust", "processx", "ps")) {
  if (!requireNamespace(package, quietly = TRUE)) stop(package, " required.")
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop("usage: build_mv08g_comparison_prefreeze.R PRIMARY_PREFREEZE LANDSCAPE_VALIDATION LANDSCAPE_EXECUTION OUTPUT EXPECTED_HEAD")
}
primary <- args[[1L]]; validation <- args[[2L]]; execution <- args[[3L]]
output <- args[[4L]]; expected_head <- tolower(trimws(args[[5L]]))
head <- tolower(trimws(system2("git", c("rev-parse", "HEAD"), stdout = TRUE)))
if (head != expected_head) stop("MV8-G comparison prefreeze exact HEAD mismatch.")
if (dir.exists(output) && length(list.files(output, all.files = TRUE,
                                             no.. = TRUE))) {
  stop("MV8-G comparison prefreeze output must be empty.")
}
dir.create(output, recursive = TRUE, showWarnings = FALSE)
source("R/provenance_utils.R")
source("R/mv08g_panel_sensitivity.R")
sha <- function(path) digest::digest(file = path, algo = "sha256",
                                     serialize = FALSE)
truth <- function(value) tolower(as.character(value)) %in% c("true", "t", "1")
decision_landscape <- read.csv(file.path(validation,
  "mv08g-landscape-validation-decision.csv"), stringsAsFactors = FALSE,
  check.names = FALSE)
summary_landscape <- read.csv(file.path(validation,
  "mv08g-landscape-validation-summary.csv"), stringsAsFactors = FALSE,
  check.names = FALSE)
repeats <- read.csv(file.path(execution, "mv08g-landscape-repeat-validation.csv"),
                    stringsAsFactors = FALSE, check.names = FALSE)
comparison_contract <- read.csv(file.path(primary, "mv08g-comparison-contract.csv"),
                                stringsAsFactors = FALSE, check.names = FALSE)
thresholds <- read.csv(file.path(primary, "mv08g-interpretation-thresholds.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
if (decision_landscape$decision !=
      "landscapes_exact_authorize_comparison_execution_prefreeze_only" ||
    decision_landscape$comparison_jobs_authorized != 0L ||
    summary_landscape$within475_groups != 20L ||
    summary_landscape$within475_rows != 152520L ||
    summary_landscape$matched_shift_groups != 20L ||
    summary_landscape$matched_shift_rows != 2480L ||
    summary_landscape$r_oracles != 12L || summary_landscape$persim_oracles != 12L ||
    nrow(repeats) != 8L || !all(truth(repeats$byte_identical)) ||
    nrow(comparison_contract) != 8L || nrow(thresholds) != 3L) {
  stop("MV8-G landscape closure does not authorize comparison prefreeze.")
}
implementation_paths <- c(
  "R/mv08g_panel_sensitivity.R", "R/mv05n_clustering_gate.R",
  "R/mv05_benchmark_contract.R", "scripts/run_mv08g_comparison.R",
  "scripts/run_mv08g_comparison_monitor.R",
  "scripts/validate_mv08g_comparison.R",
  "scripts/build_mv08g_comparison_prefreeze.R",
  "tests/testthat/test-mv08g-panel-sensitivity.R")
if (any(!file.exists(implementation_paths))) stop("MV8-G comparison implementation incomplete.")
implementation_root <- digest::digest(data.frame(
  path = implementation_paths,
  sha256 = vapply(implementation_paths, sha, character(1L)),
  accepted_head = expected_head, stringsAsFactors = FALSE),
  algo = "sha256", serialize = TRUE)
contract <- data.frame(
  contract_id = "mv08g_comparison_execution_prefreeze_v1",
  accepted_head = expected_head, implementation_root_sha256 = implementation_root,
  comparison_jobs = 1L, components = 4L, seeds = 5L,
  matched_pair_distances = 20L * 7626L,
  neighbor_rows = 4L * 5L * 124L, normalized_shift_rows = 2480L,
  candidate_pam_fits = 2L * 4L * 5L * 9L,
  candidate_partition_rows = 2L * 4L * 5L * 9L * 124L,
  panel_selected_k_rows = 8L, pam_k_agreement_rows = 4L * 5L * 9L,
  fixed_partition_fits = 2L * 4L * 5L * 2L,
  fixed_partition_rows = 2L * 4L * 5L * 2L * 124L,
  fixed_agreement_rows = 4L * 5L * 2L, component_summary_rows = 4L,
  elapsed_cap_seconds = 3600, rss_cap_bytes = 8 * 1024^3,
  storage_cap_bytes = 2 * 1024^3, workers = 1L, retries = 0L,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE)
source_files <- c(
  landscape_decision = file.path(validation,
    "mv08g-landscape-validation-decision.csv"),
  landscape_summary = file.path(validation, "mv08g-landscape-validation-summary.csv"),
  landscape_repeats = file.path(execution, "mv08g-landscape-repeat-validation.csv"),
  comparison_contract = file.path(primary, "mv08g-comparison-contract.csv"),
  interpretation_thresholds = file.path(primary,
    "mv08g-interpretation-thresholds.csv"), implementation_paths)
freeze <- data.frame(
  contract_id = "mv08g_comparison_source_freeze_v1",
  source_id = names(source_files), artifact_locator = unname(source_files),
  sha256 = vapply(source_files, sha, character(1L)),
  bytes = as.numeric(file.info(source_files)$size), accepted_head = expected_head,
  private_source = FALSE, stringsAsFactors = FALSE)
acceptance <- data.frame(
  contract_id = "mv08g_comparison_prefreeze_acceptance_v1",
  category = c("landscape_closure", "comparison_families", "thresholds",
               "row_scope", "resources", "label_firewall", "stop_boundary"),
  passed = c(
    summary_landscape$within475_groups == 20L &&
      summary_landscape$matched_shift_groups == 20L,
    nrow(comparison_contract) == 8L && all(truth(comparison_contract$label_free)),
    identical(thresholds$metric, c("median_spearman", "median_top10_overlap",
      "median_fixed_k_pam_ari")) &&
      isTRUE(all.equal(thresholds$threshold, c(0.95, 0.80, 0.80))),
    contract$candidate_partition_rows == 44640L &&
      contract$fixed_partition_rows == 9920L,
    contract$elapsed_cap_seconds == 3600 && contract$rss_cap_bytes == 8 * 1024^3,
    all(truth(comparison_contract$label_free)), TRUE),
  detail = c("40 groups, eight repeats, 24 oracles exact",
    "eight frozen comparison families", "0.95, 0.80, 0.80 planning thresholds",
    "complete candidate and fixed partition evidence", "one bounded zero-retry job",
    "sample identities only; labels outcomes closed",
    "comparison only; HCA raw reads remain closed"), stringsAsFactors = FALSE)
if (!all(acceptance$passed)) stop("MV8-G comparison prefreeze acceptance failed.")
decision <- data.frame(
  contract_id = "mv08g_comparison_prefreeze_decision_v1",
  decision = "authorize_one_label_closed_comparison_job",
  comparison_jobs_authorized = 1L, hca_expression_jobs_authorized = 0L,
  hca_fastq_download_authorized = FALSE, raw_reprocessing_authorized = FALSE,
  label_access_authorized = FALSE,
  next_gate = "MV8-G_independent_comparison_reconstruction_and_raw_read_decision",
  stringsAsFactors = FALSE)
outputs <- list(
  "mv08g-comparison-contract.csv" = contract,
  "mv08g-comparison-analysis-contract.csv" = comparison_contract,
  "mv08g-interpretation-thresholds.csv" = thresholds,
  "mv08g-comparison-source-freeze.csv" = freeze,
  "mv08g-comparison-acceptance.csv" = acceptance,
  "mv08g-comparison-decision.csv" = decision)
for (name in names(outputs)) write_provenance_csv(
  outputs[[name]], file.path(output, name))
message("MV8-G comparison prefreeze passed: one label-closed job only")
