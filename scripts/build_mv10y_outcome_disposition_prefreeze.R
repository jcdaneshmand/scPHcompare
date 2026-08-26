#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) stop(paste(
  "usage: build_mv10y_outcome_disposition_prefreeze.R <mv10t-synthesis>",
  "<mv10u-closure> <mv10w-closure> <output> <execution-head>"
), call. = FALSE)
synthesis <- normalizePath(args[[1L]], mustWork = TRUE)
synthesis_closure <- normalizePath(args[[2L]], mustWork = TRUE)
figure_closure <- normalizePath(args[[3L]], mustWork = TRUE)
output <- args[[4L]]
execution_head <- tolower(trimws(args[[5L]]))
if (!grepl("^[0-9a-f]{40}$", execution_head)) stop("invalid execution head")
if (dir.exists(output) && length(list.files(output, all.files = TRUE,
                                            no.. = TRUE))) {
  stop("refusing to overwrite MV10-Y prefreeze")
}
source("R/mv08z_landscape_production.R")
readc <- .mv08z_read_csv; sha <- .mv08z_sha256_file; atomic <- .mv08z_atomic_csv
.mv08z_verify_manifest(synthesis, "mv10t-artifact-manifest.csv")
.mv08z_verify_manifest(synthesis_closure, "mv10u-artifact-manifest.csv")
.mv08z_verify_manifest(figure_closure, "mv10w-artifact-manifest.csv")
u_validation <- readc(file.path(synthesis_closure, "mv10u-validation.csv"))
w_validation <- readc(file.path(figure_closure, "mv10w-validation.csv"))
implementation_files <- c(
  "R/mv08z_landscape_production.R", "R/mv10y_outcome_disposition.R",
  "scripts/build_mv10y_outcome_disposition_prefreeze.R",
  "scripts/run_mv10z_outcome_disposition.R",
  "scripts/build_mv10za_outcome_disposition_closure.R"
)
implementation <- data.frame(
  contract_id = "mv10y_implementation_binding_v1",
  implementation_order = seq_along(implementation_files),
  file = implementation_files,
  bytes = as.numeric(file.info(implementation_files)$size),
  sha256 = vapply(implementation_files, sha, character(1L)),
  stringsAsFactors = FALSE
)
source_files <- c(
  file.path(synthesis, "mv10t-artifact-manifest.csv"),
  file.path(synthesis, "mv10t-complete-summary.csv"),
  file.path(synthesis, "mv10t-endpoint-coverage.csv"),
  file.path(synthesis_closure, "mv10u-artifact-manifest.csv"),
  file.path(figure_closure, "mv10w-artifact-manifest.csv")
)
source_freeze <- data.frame(
  contract_id = "mv10y_source_freeze_v1", source_order = seq_along(source_files),
  artifact = source_files, bytes = as.numeric(file.info(source_files)$size),
  sha256 = vapply(source_files, sha, character(1L)), stringsAsFactors = FALSE
)
rules <- data.frame(
  contract_id = "mv10y_rule_contract_v1", rule_order = 1:10,
  rule_id = c(
    "primary_PAM_only_for_primary_envelope", "all_methods_complete_sensitivity",
    "representations_not_ranked", "H0_H1_separate", "two_metrics_separate",
    "primary90_minus_full124_context", "approach_minus_tissue_study_diagnostic",
    "no_magnitude_threshold", "no_inference_or_selection", "no_claims"
  ), result_dependent_threshold = FALSE, p_value = FALSE,
  method_selection = FALSE, ranking = FALSE, biological_claim = FALSE,
  stringsAsFactors = FALSE
)
outputs <- data.frame(
  contract_id = "mv10y_output_contract_v1",
  output_id = c(
    "primary_representation_envelope", "method_envelope", "context_contrast",
    "approach_contrast", "disposition"
  ), expected_rows = c(20L, 60L, 120L, 120L, 1L),
  inference = FALSE, ranking = FALSE, stringsAsFactors = FALSE
)
contract <- data.frame(
  contract_id = "mv10y_outcome_disposition_prefreeze_v1",
  execution_head = execution_head, value_aware_prefreeze = TRUE,
  figures_already_reviewed = TRUE, source_summary_rows = 300L,
  primary_rows = 60L, output_tables = 5L,
  magnitude_threshold = "none", p_values = FALSE,
  method_selection = FALSE, representation_ranking = FALSE,
  H0_H1_combined = FALSE, approach_causal_interpretation = FALSE,
  biological_claims = FALSE, manuscript_claims = FALSE,
  stringsAsFactors = FALSE
)
decision <- data.frame(
  contract_id = "mv10y_decision_v1",
  decision = "authorize_transparent_descriptive_disposition_after_commit",
  execution_authorized_after_commit = TRUE, method_selection_authorized = FALSE,
  biological_claims_authorized = FALSE, manuscript_claims_authorized = FALSE,
  stringsAsFactors = FALSE
)
validation <- data.frame(
  contract_id = "mv10y_validation_v1",
  check_id = c(
    "synthesis_manifest", "synthesis_closure_manifest", "figure_closure_manifest",
    "synthesis_25_of_25", "figures_25_of_25", "five_implementation_files",
    "implementation_hashes", "five_source_freezes", "source_hashes",
    "ten_rules", "no_result_thresholds", "five_outputs",
    "output_cardinalities", "value_awareness_disclosed", "figures_disclosed",
    "three_hundred_source_rows", "sixty_primary_rows", "no_magnitude_threshold",
    "no_p_values", "no_method_selection", "no_representation_ranking",
    "H0_H1_separate", "approach_noncausal", "claim_firewall",
    "execution_only_authorized"
  ),
  passed = c(
    TRUE, TRUE, TRUE, nrow(u_validation) == 25L && all(u_validation$passed),
    nrow(w_validation) == 25L && all(w_validation$passed),
    nrow(implementation) == 5L, all(file.exists(implementation$file)),
    nrow(source_freeze) == 5L, all(file.exists(source_freeze$artifact)),
    nrow(rules) == 10L, all(!rules$result_dependent_threshold),
    nrow(outputs) == 5L,
    identical(outputs$expected_rows, c(20L, 60L, 120L, 120L, 1L)),
    contract$value_aware_prefreeze, contract$figures_already_reviewed,
    contract$source_summary_rows == 300L, contract$primary_rows == 60L,
    contract$magnitude_threshold == "none", !contract$p_values,
    !contract$method_selection, !contract$representation_ranking,
    !contract$H0_H1_combined, !contract$approach_causal_interpretation,
    !contract$biological_claims && !contract$manuscript_claims,
    decision$execution_authorized_after_commit &&
      !decision$method_selection_authorized
  ), stringsAsFactors = FALSE
)
if (!all(validation$passed)) stop("MV10-Y prefreeze validation failed")
if (!dir.exists(output)) dir.create(output, recursive = TRUE)
artifacts <- list(
  "mv10y-contract.csv" = contract, "mv10y-rule-contract.csv" = rules,
  "mv10y-output-contract.csv" = outputs,
  "mv10y-implementation-bindings.csv" = implementation,
  "mv10y-source-freeze.csv" = source_freeze,
  "mv10y-decision.csv" = decision, "mv10y-validation.csv" = validation
)
for (name in names(artifacts)) atomic(artifacts[[name]], file.path(output, name))
writeLines(c(
  "# MV10-Y descriptive outcome-disposition prefreeze", "",
  "This is transparently value-aware because the exact figures were already",
  "reviewed. It adds no magnitude threshold, p-value, selection rule, ranking,",
  "causal approach interpretation, biological claim, or manuscript claim."
), file.path(output, "MV10Y_OUTCOME_DISPOSITION_PREFREEZE_2026-08-25.md"))
files <- sort(setdiff(list.files(output), "mv10y-artifact-manifest.csv"))
manifest <- data.frame(
  contract_id = "mv10y_artifact_manifest_v1", artifact = files,
  bytes = as.numeric(file.info(file.path(output, files))$size),
  sha256 = vapply(file.path(output, files), sha, character(1L)),
  stringsAsFactors = FALSE
)
atomic(manifest, file.path(output, "mv10y-artifact-manifest.csv"))
cat("Built MV10-Y descriptive outcome prefreeze; checks=25\n")
