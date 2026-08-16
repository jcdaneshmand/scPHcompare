#!/usr/bin/env Rscript

options(warn = 2)
for (package in c("digest", "ripserr", "TDA", "processx", "ps")) {
  if (!requireNamespace(package, quietly = TRUE)) stop(package, " required.")
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop("usage: build_mv08g_ph_prefreeze.R PRIMARY_PREFREEZE SOURCE_VALIDATION SOURCE_EXECUTION OUTPUT EXPECTED_HEAD")
}
primary <- args[[1L]]; source_validation <- args[[2L]]
source_execution <- args[[3L]]; output <- args[[4L]]
expected_head <- tolower(trimws(args[[5L]]))
head <- tolower(trimws(system2("git", c("rev-parse", "HEAD"), stdout = TRUE)))
if (head != expected_head) stop("MV8-G PH prefreeze exact HEAD mismatch.")
if (dir.exists(output) && length(list.files(output, all.files = TRUE,
                                             no.. = TRUE))) {
  stop("MV8-G PH prefreeze output must be empty.")
}
dir.create(output, recursive = TRUE, showWarnings = FALSE)
source("R/provenance_utils.R")
source("R/mv08g_panel_sensitivity.R")
sha <- function(path) digest::digest(file = path, algo = "sha256",
                                     serialize = FALSE)
truth <- function(value) tolower(as.character(value)) %in% c("true", "t", "1")
source_decision <- read.csv(file.path(source_validation,
  "mv08g-source-validation-decision.csv"), stringsAsFactors = FALSE,
  check.names = FALSE)
source_summary <- read.csv(file.path(source_validation,
  "mv08g-source-validation-summary.csv"), stringsAsFactors = FALSE,
  check.names = FALSE)
source_repeat <- read.csv(file.path(source_execution,
  "mv08g-source-repeat-validation.csv"), stringsAsFactors = FALSE,
  check.names = FALSE)
queue <- read.csv(file.path(primary, "mv08g-ph-queue.csv"),
                  stringsAsFactors = FALSE, check.names = FALSE)
axis <- read.csv(file.path(primary, "mv08g-sample-seed-axis.csv"),
                 stringsAsFactors = FALSE, check.names = FALSE)
if (source_decision$decision != "source_exact_authorize_PH_execution_prefreeze_only" ||
    source_decision$ph_jobs_authorized != 0L ||
    source_decision$prospective_ph_jobs != 1240L ||
    source_summary$source_jobs != 5L || source_summary$typed_views != 1240L ||
    !truth(source_summary$exact_repeat_passed) ||
    !truth(source_repeat$byte_identical) || nrow(queue) != 1240L ||
    any(queue$workers != 1L) || any(queue$retries != 0L)) {
  stop("MV8-G source closure does not authorize PH prefreeze.")
}
sample_ids <- sort(unique(axis$sample_id), method = "radix")
repeat_queue <- data.frame(
  contract_id = "mv08g_ph_repeat_queue_v1", repeat_order = 1:4,
  seed = c(.mv08g_seeds[[1L]], .mv08g_seeds[[1L]],
           .mv08g_seeds[[5L]], .mv08g_seeds[[5L]]),
  sample_id = c(sample_ids[[1L]], sample_ids[[1L]],
                sample_ids[[124L]], sample_ids[[124L]]),
  view_id = rep(.mv08g_views, 2L), repeat_basis =
    "prospective_lexicographic_axis_extremes_balanced_by_view",
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE)
fallback <- data.frame(
  contract_id = "mv08g_exact_ph_resource_fallback_v1",
  primary_engine = "ripserr", eligible_view_id = "gene_topology_v1",
  eligible_primary_disposition = "rss_cap_exceeded",
  fallback_engine = "TDA_ripsDiag_GUDHI",
  mathematical_estimand =
    "complete_vietoris_rips_H0_H1_field_2_threshold_minus_1",
  fallback_elapsed_cap_seconds = 1800,
  fallback_rss_cap_bytes = 12 * 1024^3, fallback_workers = 1L,
  fallback_attempts = 1L, selected_engine_repeat_required = TRUE,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE)
implementation_paths <- c(
  "R/mv08g_panel_sensitivity.R", "scripts/run_mv08g_ph_entry.R",
  "scripts/run_mv08g_ph_fallback_entry.R", "scripts/run_mv08g_ph_monitor.R",
  "scripts/validate_mv08g_ph.R", "scripts/build_mv08g_ph_prefreeze.R",
  "tests/testthat/test-mv08g-panel-sensitivity.R")
if (any(!file.exists(implementation_paths))) stop("MV8-G PH implementation incomplete.")
implementation_root <- digest::digest(data.frame(
  path = implementation_paths,
  sha256 = vapply(implementation_paths, sha, character(1L)),
  accepted_head = expected_head, stringsAsFactors = FALSE),
  algo = "sha256", serialize = TRUE)
contract <- data.frame(
  contract_id = "mv08g_ph_execution_prefreeze_v1",
  accepted_head = expected_head, implementation_root_sha256 = implementation_root,
  samples = 124L, seeds = 5L, typed_views = 1240L, ph_jobs = 1240L,
  cell_ph_jobs = 620L, gene_ph_jobs = 620L, ph_repeat_jobs = 4L,
  aggregate_elapsed_cap_seconds = 172800,
  aggregate_storage_cap_bytes = 4 * 1024^3,
  max_dim = 1L, threshold = -1, field = 2L,
  filtration = "complete_vietoris_rips",
  essential_h0_policy = "one_essential_H0_excluded_downstream",
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE)
source_files <- c(
  source_validation_decision = file.path(source_validation,
    "mv08g-source-validation-decision.csv"),
  source_validation_summary = file.path(source_validation,
    "mv08g-source-validation-summary.csv"),
  source_repeat = file.path(source_execution, "mv08g-source-repeat-validation.csv"),
  ph_queue = file.path(primary, "mv08g-ph-queue.csv"),
  sample_axis = file.path(primary, "mv08g-sample-seed-axis.csv"),
  implementation_paths)
freeze <- data.frame(
  contract_id = "mv08g_ph_source_freeze_v1", source_id = names(source_files),
  artifact_locator = unname(source_files),
  sha256 = vapply(source_files, sha, character(1L)),
  bytes = as.numeric(file.info(source_files)$size), accepted_head = expected_head,
  private_source = FALSE, stringsAsFactors = FALSE)
acceptance <- data.frame(
  contract_id = "mv08g_ph_prefreeze_acceptance_v1",
  category = c("source_closure", "queue", "view_balance", "resource_caps",
               "fallback", "repeat", "label_firewall", "stop_boundary"),
  passed = c(
    truth(source_summary$exact_repeat_passed), nrow(queue) == 1240L,
    all(table(queue$view_id) == 620L),
    max(queue$rss_cap_bytes[queue$view_id == "gene_topology_v1"]) == 8 * 1024^3 &&
      max(queue$rss_cap_bytes[queue$view_id == "cell_topology_v1"]) == 4 * 1024^3,
    fallback$fallback_rss_cap_bytes == 12 * 1024^3 &&
      fallback$fallback_attempts == 1L,
    nrow(repeat_queue) == 4L && all(table(repeat_queue$view_id) == 2L),
    all(queue$outcome_label_state == "closed") &&
      all(repeat_queue$outcome_label_state == "closed"), TRUE),
  detail = c("five source bundles and repeat exact", "1,240 serial zero-retry jobs",
    "620 cell plus 620 gene", "4-GiB cell; 8-GiB gene; 48-hour and 4-GiB aggregate caps",
    "gene RSS-only exact 12-GiB fallback", "four balanced selected-engine repeats",
    "labels and outcomes closed", "PH only; landscapes and HCA remain closed"),
  stringsAsFactors = FALSE)
if (!all(acceptance$passed)) stop("MV8-G PH prefreeze acceptance failed.")
decision <- data.frame(
  contract_id = "mv08g_ph_prefreeze_decision_v1",
  decision = "authorize_1240_PH_jobs_and_four_selected_engine_repeats",
  ph_jobs_authorized = 1240L, ph_repeat_jobs_authorized = 4L,
  gene_fallback_attempts_authorized_only_after_rss_cap = 1L,
  landscape_jobs_authorized = 0L, matched_shift_jobs_authorized = 0L,
  comparison_jobs_authorized = 0L, hca_fastq_download_authorized = FALSE,
  raw_reprocessing_authorized = FALSE, label_access_authorized = FALSE,
  next_gate = "MV8-G_full_PH_independent_validation",
  stringsAsFactors = FALSE)
outputs <- list(
  "mv08g-ph-contract.csv" = contract,
  "mv08g-ph-queue.csv" = queue,
  "mv08g-ph-repeat-queue.csv" = repeat_queue,
  "mv08g-ph-fallback-policy.csv" = fallback,
  "mv08g-ph-source-freeze.csv" = freeze,
  "mv08g-ph-acceptance.csv" = acceptance,
  "mv08g-ph-decision.csv" = decision)
for (name in names(outputs)) write_provenance_csv(
  outputs[[name]], file.path(output, name))
message("MV8-G PH prefreeze passed: 1,240 PH plus four selected-engine repeats")
