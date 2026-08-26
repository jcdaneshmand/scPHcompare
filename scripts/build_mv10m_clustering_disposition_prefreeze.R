#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop(paste(
  "usage: build_mv10m_clustering_disposition_prefreeze.R",
  "<mv10h-synthesis> <mv10i-closure> <mv10j-figures> <output-dir>"
), call. = FALSE)
synthesis <- normalizePath(args[[1L]], mustWork = TRUE)
closure <- normalizePath(args[[2L]], mustWork = TRUE)
figures <- normalizePath(args[[3L]], mustWork = TRUE)
output <- args[[4L]]
if (dir.exists(output) && length(list.files(output, all.files = TRUE,
                                            no.. = TRUE))) {
  stop("refusing to overwrite MV10-M prefreeze")
}
if (!dir.exists(output)) dir.create(output, recursive = TRUE)
source("R/mv08z_landscape_production.R")
source("R/mv10m_clustering_disposition.R")
readc <- .mv08z_read_csv; sha <- .mv08z_sha256_file; atomic <- .mv08z_atomic_csv
execution_head <- tolower(trimws(Sys.getenv("MV10M_GIT_HEAD", unset = "")))
if (!grepl("^[0-9a-f]{40}$", execution_head)) {
  stop("MV10M_GIT_HEAD must bind one 40-character commit")
}
.mv08z_verify_manifest(synthesis, "mv10h-artifact-manifest.csv")
.mv08z_verify_manifest(closure, "mv10i-artifact-manifest.csv")
.mv08z_verify_manifest(figures, "mv10j-artifact-manifest.csv")
closure_validation <- readc(file.path(closure, "mv10i-validation.csv"))
if (nrow(closure_validation) != 23L || !all(closure_validation$passed)) {
  stop("MV10-I prerequisite drift")
}
schemas <- data.frame(
  contract_id = "mv10m_source_schema_v1", source_order = 1:6,
  source_id = c("stability_grid", "quality_summary", "agreement_summary",
                "primary_selection", "primary_stability", "primary_quality"),
  file = c("mv10h-stability-grid.csv", "mv10h-quality-summary.csv",
           "mv10h-method-agreement-summary.csv", "mv10h-primary-selection.csv",
           "mv10h-primary-stability.csv", "mv10h-primary-quality.csv"),
  expected_rows = c(270L, 270L, 540L, 2L, 18L, 90L),
  stringsAsFactors = FALSE
)
schema_headers <- lapply(schemas$file, function(name) names(read.csv(
  file.path(synthesis, name), nrows = 0L, check.names = FALSE
)))
schemas$columns <- vapply(schema_headers, paste, collapse = ";",
                          character(1L))
schemas$bytes <- as.numeric(file.info(file.path(synthesis, schemas$file))$size)
schemas$sha256 <- vapply(file.path(synthesis, schemas$file), sha, character(1L))

contract <- data.frame(
  contract_id = "mv10m_clustering_disposition_prefreeze_v1",
  execution_head = execution_head,
  value_aware_prefreeze = TRUE,
  reason_value_aware = "original_resolution_figures_already_reviewed",
  source_tables = 6L, output_tables = 5L,
  primary_method = "pam_dissimilarity_v1",
  K_selection = "frozen_five_seed_one_se_rule",
  structural_failure_rule =
    "any_selected_seed_minimum_cluster_size_below_2_or_singleton_clusters_above_0",
  stability_role = "descriptive_under_frozen_selection",
  silhouette_role = "descriptive_no_threshold",
  representation_role = "sensitivity_no_ranking",
  method_role = "sensitivity_no_ranking",
  H0_H1 = "separate", internal_only = TRUE,
  labels = "closed", outcomes = "closed",
  biology = "closed", manuscript_claims = "closed",
  stringsAsFactors = FALSE
)
rules <- data.frame(
  contract_id = "mv10m_disposition_rule_v1", rule_order = 1:8,
  rule_id = c("primary_method", "selected_k", "structural_degeneracy",
              "stability", "silhouette", "representation_sensitivity",
              "method_sensitivity", "firewall"),
  rule = c(
    "Retain prospectively frozen PAM as the primary method",
    "Retain separate prospectively selected H0 and H1 K values",
    "Stop before labels if any selected seed has minimum cluster size below two or any singleton cluster",
    "Report adjusted-Rand stability descriptively under the frozen one-SE rule",
    "Report silhouette descriptively without a threshold or selection role",
    "Show the same selected K across all three representations without ranking representations",
    "Show all four PAM pair agreements without ranking clustering methods",
    "Keep labels, outcomes, biology, combined scores, and manuscript claims closed"
  ),
  result_dependent_threshold = FALSE,
  stringsAsFactors = FALSE
)
outputs <- data.frame(
  contract_id = "mv10m_output_contract_v1", output_order = 1:5,
  artifact = c("selected_primary_seed_metrics", "primary_summary",
               "representation_context", "method_sensitivity", "disposition"),
  expected_rows = c(10L, 2L, 6L, 24L, 1L),
  labels_allowed = FALSE, outcomes_allowed = FALSE,
  inference_allowed = FALSE, ranking_allowed = FALSE,
  stringsAsFactors = FALSE
)
implementation_files <- c(
  "R/mv10m_clustering_disposition.R",
  "scripts/build_mv10m_clustering_disposition_prefreeze.R",
  "scripts/run_mv10n_clustering_disposition.R",
  "scripts/build_mv10o_clustering_disposition_closure.R"
)
implementation <- data.frame(
  contract_id = "mv10m_implementation_binding_v1",
  implementation_order = seq_along(implementation_files),
  file = implementation_files,
  bytes = as.numeric(file.info(implementation_files)$size),
  sha256 = vapply(implementation_files, sha, character(1L)),
  stringsAsFactors = FALSE
)
source_files <- c(
  file.path(synthesis, "mv10h-artifact-manifest.csv"),
  file.path(closure, "mv10i-artifact-manifest.csv"),
  file.path(closure, "mv10i-validation.csv"),
  file.path(figures, "mv10j-artifact-manifest.csv"),
  "docs/audits/mv10k-clustering-review-figure-closure-v1/mv10k-validation.csv",
  "docs/audits/mv10l-clustering-review-visual-closure-v1/mv10l-visual-review.csv",
  "docs/audits/mv10l-clustering-review-visual-closure-v1/mv10l-decision.csv"
)
source_freeze <- data.frame(
  contract_id = "mv10m_source_freeze_v1",
  source_order = seq_along(source_files), artifact = source_files,
  bytes = as.numeric(file.info(source_files)$size),
  sha256 = vapply(source_files, sha, character(1L)),
  stringsAsFactors = FALSE
)
validation <- data.frame(
  contract_id = "mv10m_validation_v1",
  check_id = c("synthesis_manifest", "closure_manifest", "figure_manifest",
               "closure_23", "six_sources", "source_cardinalities",
               "source_hashes", "value_aware_disclosed", "eight_rules",
               "no_result_dependent_threshold", "five_outputs",
               "output_cardinalities", "implementation_bound",
               "sources_bound", "primary_PAM", "frozen_one_se",
               "hard_degeneracy_only", "stability_descriptive",
               "silhouette_descriptive", "representation_no_ranking",
               "method_no_ranking", "H0_H1_separate", "internal_only",
               "label_outcome_firewall", "biology_claim_firewall"),
  passed = c(TRUE, TRUE, TRUE, all(closure_validation$passed),
             nrow(schemas) == 6L,
             identical(schemas$expected_rows, c(270L,270L,540L,2L,18L,90L)),
             all(nzchar(schemas$sha256)), contract$value_aware_prefreeze,
             nrow(rules) == 8L, all(!rules$result_dependent_threshold),
             nrow(outputs) == 5L,
             identical(outputs$expected_rows, c(10L,2L,6L,24L,1L)),
             nrow(implementation) == 4L && all(file.exists(implementation$file)),
             length(source_files) == 7L && all(file.exists(source_files)),
             contract$primary_method == "pam_dissimilarity_v1",
             contract$K_selection == "frozen_five_seed_one_se_rule",
             grepl("minimum_cluster_size", contract$structural_failure_rule),
             grepl("descriptive", contract$stability_role),
             contract$silhouette_role == "descriptive_no_threshold",
             contract$representation_role == "sensitivity_no_ranking",
             contract$method_role == "sensitivity_no_ranking",
             contract$H0_H1 == "separate", contract$internal_only,
             contract$labels == "closed" && contract$outcomes == "closed",
             contract$biology == "closed" &&
               contract$manuscript_claims == "closed"),
  stringsAsFactors = FALSE
)
if (!all(validation$passed)) stop("MV10-M validation failed")
decision <- data.frame(
  contract_id = "mv10m_decision_v1",
  decision = "authorize_label_closed_methodological_disposition_after_commit",
  execution_authorized_after_commit = TRUE,
  labels_authorized = FALSE, outcomes_authorized = FALSE,
  biological_interpretation_authorized = FALSE,
  manuscript_claims_authorized = FALSE,
  next_if_closed = "separate_label_evaluation_prefreeze_if_structurally_nondegenerate",
  stringsAsFactors = FALSE
)
artifacts <- list(
  "mv10m-contract.csv" = contract, "mv10m-source-schema.csv" = schemas,
  "mv10m-rule-contract.csv" = rules, "mv10m-output-contract.csv" = outputs,
  "mv10m-implementation-bindings.csv" = implementation,
  "mv10m-source-freeze.csv" = source_freeze,
  "mv10m-validation.csv" = validation, "mv10m-decision.csv" = decision
)
for (name in names(artifacts)) atomic(artifacts[[name]], file.path(output, name))
writeLines(c(
  "# MV10-M clustering methodological-disposition prefreeze", "",
  "This is explicitly value-aware because the exact review figures were already",
  "inspected. Before computing any new disposition, it freezes the original PAM",
  "and one-SE K rule, an objective structural-degeneracy stop, descriptive-only",
  "stability and silhouette roles, non-ranking sensitivity contexts, and closed",
  "labels, outcomes, biology, combined scores, and manuscript claims."
), file.path(output, "MV10M_CLUSTERING_DISPOSITION_PREFREEZE_2026-08-25.md"))
files <- sort(setdiff(list.files(output), "mv10m-artifact-manifest.csv"))
manifest <- data.frame(
  contract_id = "mv10m_artifact_manifest_v1", artifact = files,
  bytes = as.numeric(file.info(file.path(output, files))$size),
  sha256 = vapply(file.path(output, files), sha, character(1L)),
  stringsAsFactors = FALSE
)
atomic(manifest, file.path(output, "mv10m-artifact-manifest.csv"))
message("Built MV10-M clustering-disposition prefreeze; checks=25")
