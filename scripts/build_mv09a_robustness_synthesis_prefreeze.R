#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop(paste(
  "usage: build_mv09a_robustness_synthesis_prefreeze.R",
  "<mv08zy-public-summary> <output-dir>"
), call. = FALSE)
source_path <- normalizePath(args[[1L]], mustWork = TRUE)
output <- args[[2L]]
if (dir.exists(output) && length(list.files(output, all.files = TRUE,
                                            no.. = TRUE))) {
  stop("refusing to overwrite MV9-A prefreeze")
}
if (!dir.exists(output)) dir.create(output, recursive = TRUE)
source("R/mv08z_landscape_production.R")
source("R/mv09_robustness_synthesis.R")
readc <- .mv08z_read_csv; sha <- .mv08z_sha256_file; atomic <- .mv08z_atomic_csv
execution_head <- tolower(trimws(Sys.getenv("MV09A_GIT_HEAD", unset = "")))
if (!grepl("^[0-9a-f]{40}$", execution_head)) {
  stop("MV09A_GIT_HEAD must bind one 40-character commit")
}
value <- readc(source_path)
# Structural validation is permitted; the synthesis function itself is not run
# until after this audit is committed.
metrics <- .mv09_metrics_v1()
required <- c("comparison_order", "dataset_scope", "contrast_id", "seed",
              "homology_dimension", metrics, "outcome_label_state",
              "biological_outcomes_computed", "clustering_jobs", "fusion_jobs",
              "label_jobs", "outcome_jobs")
if (nrow(value) != 40L || !all(required %in% names(value)) ||
    !identical(as.integer(value$comparison_order), 1:40) ||
    any(!is.finite(as.matrix(value[metrics]))) ||
    any(value$outcome_label_state != "closed") ||
    any(.mv08z_truth(value$biological_outcomes_computed)) ||
    any(value$clustering_jobs != 0L | value$fusion_jobs != 0L |
          value$label_jobs != 0L | value$outcome_jobs != 0L)) {
  stop("MV9-A aggregate source contract drift")
}
closure <- "docs/audits/mv08zz-distance-comparison-closure-v1"
.mv08z_verify_manifest(closure, "mv08zz-artifact-manifest.csv")
closure_validation <- readc(file.path(closure, "mv08zz-validation.csv"))
if (!all(closure_validation$passed)) stop("MV9-A closure prerequisite drift")

contract <- data.frame(
  contract_id = "mv09a_robustness_synthesis_prefreeze_v1",
  execution_head = execution_head, source_rows = 40L,
  source_summary_bytes = as.numeric(file.info(source_path)$size),
  source_summary_sha256 = sha(source_path), metrics = length(metrics),
  internal_rows = sum(value$dataset_scope == "internal124"),
  external_rows = sum(value$dataset_scope == "external8"),
  H0_rows = sum(value$homology_dimension == "H0"),
  H1_rows = sum(value$homology_dimension == "H1"),
  inference = "none", equivalence_threshold = "none",
  combined_robustness_score = "forbidden",
  external_generalization_claim = "forbidden_single_cohort",
  biological_interpretation = "closed",
  outcome_label_state = "closed", clustering_state = "closed",
  fusion_state = "closed", manuscript_claims_state = "closed",
  cleanup_deletion_state = "closed", stringsAsFactors = FALSE
)
metric_ranges <- data.frame(
  contract_id = "mv09a_metric_range_v1", metric_order = seq_along(metrics),
  metric = metrics, minimum = vapply(value[metrics], min, numeric(1L)),
  median = vapply(value[metrics], stats::median, numeric(1L)),
  maximum = vapply(value[metrics], max, numeric(1L)),
  finite = TRUE, interpretation = "descriptive_source_inventory_only",
  stringsAsFactors = FALSE
)
outputs <- data.frame(
  contract_id = "mv09a_output_contract_v1", output_order = 1:6,
  artifact = c("aggregate_comparisons", "plot_data_long",
               "internal_seed_summary", "external_singleton",
               "paired_dimension_delta", "dimension_delta_summary"),
  rows = c(40L, 440L, 66L, 110L, 220L, 88L),
  rule = c(
    "exact public aggregate copy with provenance",
    "eleven metrics x forty strata; H0/H1 explicit",
    "three contrasts x H0/H1 x eleven metrics; five seeds; min/q25/median/q75/max",
    "five contrasts x H0/H1 x eleven metrics; replication_units=1",
    "twenty matched scope/contrast/seed pairs x eleven metrics; H1-H0",
    "scope/contrast/metric summaries; singleton external remains n=1"
  ),
  labels_allowed = FALSE, outcomes_allowed = FALSE,
  inference_allowed = FALSE, stringsAsFactors = FALSE
)
figures <- data.frame(
  contract_id = "mv09a_figure_contract_v1", figure_order = 1:3,
  figure_id = c("internal_seed_small_multiples",
                "external_singleton_small_multiples",
                "paired_H1_minus_H0_small_multiples"),
  geometry = c("seed points plus median and IQR", "one point per contrast",
               "seed points plus median and IQR; external singleton point"),
  facet = c("metric x homology_dimension", "metric x homology_dimension",
            "metric x dataset_scope"),
  forbidden = c("thresholds;rankings;combined score", rep(
    "error bars implying replication;generalization;combined score", 2L)),
  generation_state = "closed_pending_mv09a_commit",
  stringsAsFactors = FALSE
)
implementation_files <- c(
  "R/mv09_robustness_synthesis.R",
  "scripts/build_mv09a_robustness_synthesis_prefreeze.R",
  "scripts/run_mv09b_robustness_synthesis.R",
  "scripts/build_mv09c_robustness_synthesis_closure.R"
)
implementation <- data.frame(
  contract_id = "mv09a_implementation_binding_v1",
  implementation_order = seq_along(implementation_files),
  file = implementation_files,
  bytes = as.numeric(file.info(implementation_files)$size),
  sha256 = vapply(implementation_files, sha, character(1L)),
  stringsAsFactors = FALSE
)
source_files <- c(
  file.path(closure, "mv08zz-artifact-manifest.csv"),
  file.path(closure, "mv08zz-validation.csv"),
  file.path(closure, "mv08zz-decision.csv")
)
source_freeze <- data.frame(
  contract_id = "mv09a_source_freeze_v1", source_order = seq_along(source_files),
  artifact = source_files, bytes = as.numeric(file.info(source_files)$size),
  sha256 = vapply(source_files, sha, character(1L)), stringsAsFactors = FALSE
)
validation <- data.frame(
  contract_id = "mv09a_validation_v1",
  check_id = c("mv08zz_closure", "source_hash", "row_cardinality",
               "metric_cardinality", "finite_complete", "internal_30",
               "external_10", "H0_H1_separate", "five_seed_internal_design",
               "external_singleton_disclosed", "six_outputs_frozen",
               "three_figures_frozen", "implementation_bound", "no_inference",
               "no_combined_score", "label_outcome_firewall",
               "clustering_fusion_firewall", "claim_cleanup_firewall"),
  passed = c(all(closure_validation$passed), nzchar(contract$source_summary_sha256),
             nrow(value) == 40L, length(metrics) == 11L,
             all(metric_ranges$finite), contract$internal_rows == 30L,
             contract$external_rows == 10L,
             contract$H0_rows == 20L && contract$H1_rows == 20L,
             all(table(value$seed[value$dataset_scope == "internal124"]) == 6L),
             all(outputs$rows[outputs$artifact == "external_singleton"] == 110L),
             nrow(outputs) == 6L, nrow(figures) == 3L,
             nrow(implementation) == 4L, contract$inference == "none",
             contract$combined_robustness_score == "forbidden",
             contract$outcome_label_state == "closed",
             contract$clustering_state == "closed" && contract$fusion_state == "closed",
             contract$manuscript_claims_state == "closed" &&
               contract$cleanup_deletion_state == "closed"),
  stringsAsFactors = FALSE
)
if (!all(validation$passed)) stop("MV9-A validation failed")
decision <- data.frame(
  contract_id = "mv09a_decision_v1",
  decision = "authorize_aggregate_label_closed_synthesis_after_commit",
  execution_authorized_after_commit = TRUE,
  next_if_closed = "separate_scientific_review_and_figure_prefreeze",
  stringsAsFactors = FALSE
)
artifacts <- list(
  "mv09a-contract.csv" = contract, "mv09a-metric-ranges.csv" = metric_ranges,
  "mv09a-output-contract.csv" = outputs, "mv09a-figure-contract.csv" = figures,
  "mv09a-implementation-bindings.csv" = implementation,
  "mv09a-source-freeze.csv" = source_freeze,
  "mv09a-validation.csv" = validation, "mv09a-decision.csv" = decision
)
for (name in names(artifacts)) atomic(artifacts[[name]], file.path(output, name))
writeLines(c(
  "# MV9-A aggregate robustness-synthesis prefreeze", "",
  "This gate publishes and summarizes only the already-closed 40-row aggregate",
  "comparison corpus. Internal five-seed variation and external singleton",
  "evidence remain visibly distinct. H0/H1 are never combined into one score.",
  "No synthesis or figure was generated while building this audit."
), file.path(output, "MV09A_ROBUSTNESS_SYNTHESIS_PREFREEZE_2026-08-25.md"))
files <- sort(setdiff(list.files(output), "mv09a-artifact-manifest.csv"))
manifest <- data.frame(
  contract_id = "mv09a_artifact_manifest_v1", artifact = files,
  bytes = as.numeric(file.info(file.path(output, files))$size),
  sha256 = vapply(file.path(output, files), sha, character(1L)),
  stringsAsFactors = FALSE
)
atomic(manifest, file.path(output, "mv09a-artifact-manifest.csv"))
message("Built MV9-A robustness-synthesis prefreeze; checks=18")
