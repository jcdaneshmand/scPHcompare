#!/usr/bin/env Rscript

options(warn = 2)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required.")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("usage: validate_mv08g_prefreeze.R PREFREEZE RECOVERY_EVIDENCE OUTPUT")
}
prefreeze <- args[[1L]]; recovery_dir <- args[[2L]]; output <- args[[3L]]
if (dir.exists(output) && length(list.files(output, all.files = TRUE,
                                             no.. = TRUE))) {
  stop("MV8-G independent validation output must be empty.")
}
dir.create(output, recursive = TRUE, showWarnings = FALSE)
source("R/provenance_utils.R")
source("R/mv07h_full_topology.R")
source("R/mv08g_panel_sensitivity.R")
truth <- function(value) tolower(as.character(value)) %in% c("true", "t", "1")
read <- function(name) read.csv(file.path(prefreeze, name),
  stringsAsFactors = FALSE, check.names = FALSE)
panel <- read("mv08g-panel.csv"); manifest <- read("mv08g-cache-manifest.csv")
axis <- read("mv08g-sample-seed-axis.csv"); source_queue <- read("mv08g-source-queue.csv")
ph_queue <- read("mv08g-ph-queue.csv"); landscape <- read("mv08g-landscape-queue.csv")
shift <- read("mv08g-matched-shift-queue.csv"); contract <- read("mv08g-contract.csv")
inventory500 <- read("mv08g-accepted500-landscape-inventory.csv")
comparison <- read("mv08g-comparison-contract.csv")
thresholds <- read("mv08g-interpretation-thresholds.csv")
firewall <- read("mv08g-label-firewall.csv"); decision <- read("mv08g-decision.csv")
recovery <- read.csv(file.path(recovery_dir, "mv08f-recovery-decision.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)
checks <- data.frame(
  contract_id = "mv08g_independent_prefreeze_validation_v1",
  check_id = c("recovery_gate", "panel_axis", "cache_axis", "source_queue",
    "ph_balance", "landscape_balance", "matched_shift_balance",
    "landscape_definition", "comparison_freeze", "threshold_freeze",
    "label_firewall", "authorization_boundary"),
  passed = c(
    recovery$decision == "recovery_exact_authorize_475_source_prefreeze" &&
      recovery$sct_exact == 450L && recovery$unexpected_cache_files == 0L,
    tryCatch({mv08g_validate_common475_panel_v1(panel); TRUE},
             error = function(error) FALSE),
    tryCatch({mv07h_validate_cache_manifest_v1(manifest); TRUE},
             error = function(error) FALSE) && nrow(axis) == 620L &&
      all(axis$panel_genes == 475L),
    nrow(source_queue) == 5L && sum(source_queue$typed_view_count) == 1240L,
    nrow(ph_queue) == 1240L && all(table(ph_queue$view_id) == 620L),
    nrow(landscape) == 20L && sum(landscape$component_rows) == 152520L &&
      all(table(landscape$view_id) == 10L) &&
      all(table(landscape$homology_dimension) == 10L),
    nrow(shift) == 20L && sum(shift$component_rows) == 2480L &&
      nrow(inventory500) == 20L && sum(inventory500$component_rows) == 152520L,
    grepl("finite_positive", contract$landscape_definition, fixed = TRUE) &&
      grepl("all_active_levels", contract$landscape_definition, fixed = TRUE) &&
      grepl("no_grid", contract$landscape_definition, fixed = TRUE) &&
      grepl("no_level_cap", contract$landscape_definition, fixed = TRUE),
    nrow(comparison) == 8L && all(truth(comparison$label_free)),
    identical(thresholds$metric,
      c("median_spearman", "median_top10_overlap", "median_fixed_k_pam_ari")) &&
      isTRUE(all.equal(thresholds$threshold, c(0.95, 0.80, 0.80))),
    firewall$label_access_state == "closed" &&
      sum(firewall[c("hca_expression_jobs", "hca_fastq_download_jobs",
        "raw_reprocessing_jobs", "label_jobs", "outcome_jobs")]) == 0,
    decision$decision == "authorize_five_common475_source_bundles_and_one_repeat" &&
      decision$source_jobs_authorized == 5L &&
      decision$source_repeat_jobs_authorized == 1L &&
      decision$ph_jobs_authorized == 0L &&
      !truth(decision$hca_fastq_download_authorized) &&
      !truth(decision$raw_reprocessing_authorized)),
  stringsAsFactors = FALSE)
if (!all(checks$passed)) stop("MV8-G independent prefreeze validation failed: ",
  paste(checks$check_id[!checks$passed], collapse = ", "))
out_decision <- data.frame(
  contract_id = "mv08g_independent_prefreeze_decision_v1",
  decision = "prefreeze_pass_authorize_five_source_jobs_and_one_repeat",
  passed_checks = sum(checks$passed), total_checks = nrow(checks),
  source_jobs_authorized = 5L, source_repeat_jobs_authorized = 1L,
  ph_jobs_authorized = 0L,
  hca_fastq_download_authorized = FALSE, raw_reprocessing_authorized = FALSE,
  label_access_authorized = FALSE, stringsAsFactors = FALSE)
write_provenance_csv(checks, file.path(output, "mv08g-independent-validation.csv"))
write_provenance_csv(out_decision, file.path(output, "mv08g-validation-decision.csv"))
message("MV8-G independent prefreeze validation passed: ", nrow(checks),
        "/", nrow(checks), "; source-only gate")
