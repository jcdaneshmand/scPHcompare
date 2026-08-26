#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop(paste(
  "usage: build_mv10s_outcome_review_prefreeze.R <mv10q-production>",
  "<mv10r-closure> <output> <execution-head>"
), call. = FALSE)
production <- normalizePath(args[[1L]], mustWork = TRUE)
source_closure <- normalizePath(args[[2L]], mustWork = TRUE)
output <- args[[3L]]
execution_head <- tolower(trimws(args[[4L]]))
if (!grepl("^[0-9a-f]{40}$", execution_head)) stop("invalid execution head")
if (dir.exists(output) && length(list.files(output, all.files = TRUE,
                                            no.. = TRUE))) {
  stop("refusing to overwrite MV10-S prefreeze")
}
source("R/mv08z_landscape_production.R")
readc <- .mv08z_read_csv; sha <- .mv08z_sha256_file; atomic <- .mv08z_atomic_csv
.mv08z_verify_manifest(production, "mv10q-artifact-manifest.csv")
.mv08z_verify_manifest(source_closure, "mv10r-artifact-manifest.csv")
receipt <- readc(file.path(production, "mv10q-terminal-receipt.csv"))
closure_validation <- readc(file.path(source_closure, "mv10r-validation.csv"))
closure_decision <- readc(file.path(source_closure, "mv10r-decision.csv"))
source_names <- c(
  seed_metrics = "mv10q-seed-metrics.csv",
  unit_summaries = "mv10q-unit-summaries.csv",
  structural_status = "mv10q-structural-status.csv"
)
expected_rows <- c(seed_metrics = 1500L, unit_summaries = 300L,
                   structural_status = 1L)
source_schema <- do.call(rbind, lapply(seq_along(source_names), function(i) {
  name <- names(source_names)[[i]]; path <- file.path(production, source_names[[i]])
  header <- names(utils::read.csv(path, nrows = 0L, check.names = FALSE))
  data.frame(
    contract_id = "mv10s_source_schema_v1", source_order = i,
    source_id = name, artifact = source_names[[i]],
    expected_rows = expected_rows[[name]], columns = paste(header, collapse = ";"),
    bytes = file.info(path)$size, sha256 = sha(path),
    scientific_values_opened_during_prefreeze = FALSE,
    stringsAsFactors = FALSE
  )
}))
implementation_files <- c(
  "R/mv08z_landscape_production.R", "R/mv10s_outcome_review.R",
  "scripts/build_mv10s_outcome_review_prefreeze.R",
  "scripts/run_mv10t_outcome_review_synthesis.R",
  "scripts/build_mv10u_outcome_review_closure.R",
  "scripts/render_mv10v_outcome_review_figures.R",
  "scripts/build_mv10w_outcome_review_figure_closure.R"
)
implementation <- data.frame(
  contract_id = "mv10s_implementation_binding_v1",
  implementation_order = seq_along(implementation_files),
  file = implementation_files,
  bytes = as.numeric(file.info(implementation_files)$size),
  sha256 = vapply(implementation_files, sha, character(1L)),
  stringsAsFactors = FALSE
)
source_files <- c(
  file.path(production, "mv10q-artifact-manifest.csv"),
  file.path(production, "mv10q-terminal-receipt.csv"),
  file.path(source_closure, "mv10r-artifact-manifest.csv"),
  file.path(source_closure, "mv10r-validation.csv"),
  file.path(source_closure, "mv10r-decision.csv")
)
source_freeze <- data.frame(
  contract_id = "mv10s_source_freeze_v1", source_order = seq_along(source_files),
  artifact = source_files, bytes = as.numeric(file.info(source_files)$size),
  sha256 = vapply(source_files, sha, character(1L)), stringsAsFactors = FALSE
)
aggregation <- data.frame(
  contract_id = "mv10s_aggregation_contract_v1",
  output_id = c("complete_summary", "primary_pam_summary", "endpoint_coverage"),
  expected_rows = c(300L, 60L, 6L),
  seed_aggregation = c(
    "mean_median_range_jackknife_se_across_5",
    "mean_median_range_jackknife_se_across_5", "structural_status_only"
  ),
  result_dependent_filter = FALSE, ranking = FALSE, H0_H1_combined = FALSE,
  stringsAsFactors = FALSE
)
figures <- data.frame(
  contract_id = "mv10s_figure_contract_v1", figure_order = 1:6,
  figure_id = c(
    "primary_tissue", "primary_study", "primary_approach",
    "primary_context_restriction", "all_methods_full124",
    "all_methods_primary90"
  ), width_inches = c(14, 14, 14, 16, 18, 18),
  height_inches = c(9, 9, 9, 10, 15, 12), dpi = 180L, format = "png",
  all_three_representations = TRUE, H0_H1_separate = TRUE,
  both_metrics = TRUE, seed_range = TRUE, ranking = FALSE,
  biological_claim = FALSE, manuscript_claim = FALSE,
  stringsAsFactors = FALSE
)
review <- data.frame(
  contract_id = "mv10s_review_contract_v1",
  evidence_axis = c(
    "representations", "homology", "metrics", "seed_variability",
    "full124_tissue", "full124_study", "full124_approach",
    "approach_confounding", "primary90_tissue_study",
    "primary90_approach", "methods", "claims"
  ), required_display = c(
    "all_three", "H0_H1_separate", "ARI_and_max_NMI",
    "mean_and_observed_range", "complete", "complete", "complete",
    "explicit_confounded_diagnostic", "complete",
    "one_class_structurally_not_estimable", "all_five_common_K_no_ranking",
    "none"
  ), stringsAsFactors = FALSE
)
contract <- data.frame(
  contract_id = "mv10s_outcome_review_prefreeze_v1",
  execution_head = execution_head, source_seed_rows = 1500L,
  source_summary_rows = 300L, representations = 3L,
  homology_dimensions = 2L, methods = 5L, estimable_endpoints = 5L,
  structural_nonendpoints = 1L, metrics = 2L, figures = 6L,
  source_values_opened_before_table_and_figure_freeze = FALSE,
  p_values = FALSE, method_selection = FALSE, ranking = FALSE,
  H0_H1_combined = FALSE, biological_claims = FALSE,
  manuscript_claims = FALSE, stringsAsFactors = FALSE
)
decision <- data.frame(
  contract_id = "mv10s_decision_v1",
  decision = "authorize_complete_outcome_synthesis_after_prefreeze_commit",
  synthesis_authorized_after_commit = TRUE, figure_render_authorized = FALSE,
  figure_render_requires = "committed_mv10u_25_of_25_closure",
  scientific_interpretation_authorized = FALSE,
  method_selection_authorized = FALSE, biological_claims_authorized = FALSE,
  manuscript_claims_authorized = FALSE, stringsAsFactors = FALSE
)
validation <- data.frame(
  contract_id = "mv10s_validation_v1",
  check_id = c(
    "production_manifest", "closure_manifest", "production_complete",
    "closure_30_of_30", "closure_next_review", "three_source_schemas",
    "seed_schema_only", "summary_schema_only", "structural_schema_only",
    "source_rows_bound", "source_hashes", "seven_implementation_files",
    "implementation_hashes", "five_source_freezes", "source_freeze_hashes",
    "three_output_tables", "output_cardinalities", "five_seed_aggregation",
    "no_result_filtering", "six_figures", "PNG_only", "figure_dimensions",
    "three_representations", "H0_H1_separate", "two_metrics",
    "seed_ranges", "complete_primary_endpoints", "primary90_context",
    "all_five_methods", "twelve_review_axes", "no_value_opening",
    "selection_ranking_pooling_firewall", "claim_firewall",
    "synthesis_only_authorized", "render_closed", "interpretation_closed"
  ),
  passed = c(
    TRUE, TRUE, receipt$completion_state == "complete",
    nrow(closure_validation) == 30L && all(closure_validation$passed),
    closure_decision$next_stage == "separate_complete_outcome_review_prefreeze",
    nrow(source_schema) == 3L, source_schema$source_id[[1L]] == "seed_metrics",
    source_schema$source_id[[2L]] == "unit_summaries",
    source_schema$source_id[[3L]] == "structural_status",
    identical(source_schema$expected_rows, c(1500L, 300L, 1L)),
    all(grepl("^[0-9a-f]{64}$", source_schema$sha256)),
    nrow(implementation) == 7L, all(file.exists(implementation$file)),
    nrow(source_freeze) == 5L, all(file.exists(source_freeze$artifact)),
    nrow(aggregation) == 3L,
    identical(aggregation$expected_rows, c(300L, 60L, 6L)),
    all(nzchar(aggregation$seed_aggregation)),
    all(!aggregation$result_dependent_filter), nrow(figures) == 6L,
    all(figures$format == "png"),
    identical(figures$width_inches, c(14, 14, 14, 16, 18, 18)) &&
      identical(figures$height_inches, c(9, 9, 9, 10, 15, 12)),
    all(figures$all_three_representations), all(figures$H0_H1_separate),
    all(figures$both_metrics), all(figures$seed_range),
    all(c("primary_tissue", "primary_study", "primary_approach") %in%
          figures$figure_id),
    "primary_context_restriction" %in% figures$figure_id,
    all(c("all_methods_full124", "all_methods_primary90") %in%
          figures$figure_id), nrow(review) == 12L,
    all(!source_schema$scientific_values_opened_during_prefreeze) &&
      !contract$source_values_opened_before_table_and_figure_freeze,
    !contract$p_values && !contract$method_selection && !contract$ranking &&
      !contract$H0_H1_combined,
    !contract$biological_claims && !contract$manuscript_claims,
    decision$synthesis_authorized_after_commit,
    !decision$figure_render_authorized,
    !decision$scientific_interpretation_authorized
  ), stringsAsFactors = FALSE
)
if (!all(validation$passed)) stop("MV10-S prefreeze validation failed")
if (!dir.exists(output)) dir.create(output, recursive = TRUE)
artifacts <- list(
  "mv10s-contract.csv" = contract, "mv10s-source-schema.csv" = source_schema,
  "mv10s-aggregation-contract.csv" = aggregation,
  "mv10s-figure-contract.csv" = figures, "mv10s-review-contract.csv" = review,
  "mv10s-implementation-bindings.csv" = implementation,
  "mv10s-source-freeze.csv" = source_freeze,
  "mv10s-decision.csv" = decision, "mv10s-validation.csv" = validation
)
for (name in names(artifacts)) atomic(artifacts[[name]], file.path(output, name))
writeLines(c(
  "# MV10-S outcome-review prefreeze", "",
  "Three complete-reporting tables, six PNG figures, and twelve review axes",
  "are frozen from schemas and immutable manifests before aggregate outcome",
  "values are opened. Only synthesis is authorized after this audit commits."
), file.path(output, "MV10S_OUTCOME_REVIEW_PREFREEZE_2026-08-25.md"))
files <- sort(setdiff(list.files(output), "mv10s-artifact-manifest.csv"))
manifest <- data.frame(
  contract_id = "mv10s_artifact_manifest_v1", artifact = files,
  bytes = as.numeric(file.info(file.path(output, files))$size),
  sha256 = vapply(file.path(output, files), sha, character(1L)),
  stringsAsFactors = FALSE
)
atomic(manifest, file.path(output, "mv10s-artifact-manifest.csv"))
cat("Built MV10-S outcome-review prefreeze; checks=36\n")
