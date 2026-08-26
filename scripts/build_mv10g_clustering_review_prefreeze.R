#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop(paste(
  "usage: build_mv10g_clustering_review_prefreeze.R <output>",
  "<execution-head>"
), call. = FALSE)
output <- args[[1L]]
execution_head <- tolower(trimws(args[[2L]]))
if (!grepl("^[0-9a-f]{40}$", execution_head)) stop("invalid execution head")
if (dir.exists(output) && length(list.files(output, all.files = TRUE,
                                            no.. = TRUE))) {
  stop("refusing to overwrite MV10-G prefreeze")
}
source("R/mv08z_landscape_production.R")
readc <- .mv08z_read_csv; sha <- .mv08z_sha256_file; atomic <- .mv08z_atomic_csv
production <- "docs/audits/mv10e-full-clustering-benchmark-v1"
source_closure <- "docs/audits/mv10f-full-clustering-closure-v1"
.mv08z_verify_manifest(production, "mv10e-artifact-manifest.csv")
.mv08z_verify_manifest(source_closure, "mv10f-artifact-manifest.csv")
receipt <- readc(file.path(production, "mv10e-terminal-receipt.csv"))
closure_validation <- readc(file.path(source_closure, "mv10f-validation.csv"))
closure_decision <- readc(file.path(source_closure, "mv10f-decision.csv"))
source_names <- c(
  quality = "mv10e-partition-quality.csv",
  stability = "mv10e-seed-stability.csv",
  primary_k = "mv10e-primary-k-selection.csv",
  agreement = "mv10e-method-agreement.csv"
)
expected_rows <- c(quality = 1350L, stability = 270L, primary_k = 2L,
                   agreement = 2700L)
source_schema <- do.call(rbind, lapply(seq_along(source_names), function(i) {
  name <- names(source_names)[[i]]
  path <- file.path(production, source_names[[i]])
  header <- names(utils::read.csv(path, nrows = 0L, check.names = FALSE))
  data.frame(
    contract_id = "mv10g_source_schema_v1", source_order = i,
    source_id = name, artifact = source_names[[i]], expected_rows =
      expected_rows[[name]], columns = paste(header, collapse = ";"),
    bytes = file.info(path)$size, sha256 = sha(path),
    scientific_values_opened_during_prefreeze = FALSE,
    stringsAsFactors = FALSE
  )
}))
implementation_files <- c(
  "R/mv08z_landscape_production.R",
  "R/mv10_clustering_benchmark.R",
  "R/mv10g_clustering_review.R",
  "scripts/build_mv10g_clustering_review_prefreeze.R",
  "scripts/run_mv10h_clustering_review_synthesis.R",
  "scripts/build_mv10i_clustering_review_closure.R",
  "scripts/render_mv10j_clustering_review_figures.R",
  "scripts/build_mv10k_clustering_review_figure_closure.R"
)
implementation <- data.frame(
  contract_id = "mv10g_implementation_binding_v1",
  implementation_order = seq_along(implementation_files),
  file = implementation_files,
  bytes = as.numeric(file.info(implementation_files)$size),
  sha256 = vapply(implementation_files, sha, character(1L)),
  stringsAsFactors = FALSE
)
source_files <- c(
  file.path(production, "mv10e-artifact-manifest.csv"),
  file.path(production, "mv10e-terminal-receipt.csv"),
  file.path(source_closure, "mv10f-artifact-manifest.csv"),
  file.path(source_closure, "mv10f-validation.csv"),
  file.path(source_closure, "mv10f-decision.csv")
)
source_freeze <- data.frame(
  contract_id = "mv10g_source_freeze_v1", source_order = seq_along(source_files),
  artifact = source_files, bytes = as.numeric(file.info(source_files)$size),
  sha256 = vapply(source_files, sha, character(1L)), stringsAsFactors = FALSE
)
metric_contract <- data.frame(
  contract_id = "mv10g_metric_contract_v1", metric_order = 1:8,
  metric_id = c(
    "mean_adjusted_rand", "mean_silhouette", "negative_silhouette_fraction",
    "singleton_clusters", "minimum_cluster_size", "method_pair_adjusted_rand",
    "primary_one_se_threshold", "primary_selected_k"
  ),
  role = c(
    "five_seed_stability_primary", "descriptive_internal_quality",
    "discordant_assignment_diagnostic", "fragmentation_diagnostic",
    "small_cluster_diagnostic", "method_concordance_diagnostic",
    "frozen_PAM_selection_threshold", "frozen_PAM_selection_result"
  ),
  selection_metric = c(TRUE, rep(FALSE, 5L), TRUE, TRUE),
  labels_allowed = FALSE, outcomes_allowed = FALSE, stringsAsFactors = FALSE
)
aggregation <- data.frame(
  contract_id = "mv10g_aggregation_contract_v1",
  output_id = c(
    "stability_grid", "quality_summary", "agreement_summary",
    "primary_selection", "primary_stability", "primary_quality"
  ),
  expected_rows = c(270L, 270L, 540L, 2L, 18L, 90L),
  seed_aggregation = c(
    "already_all_10_seed_pairs", "median_and_range_across_5_seeds",
    "median_and_range_across_5_seeds", "none", "already_all_10_seed_pairs",
    "none_show_all_5_seed_trajectories"
  ),
  result_dependent_filter = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
  H0_H1_combined = FALSE, cell_gene_combined = FALSE,
  stringsAsFactors = FALSE
)
figures <- data.frame(
  contract_id = "mv10g_figure_contract_v1", figure_order = 1:7,
  figure_id = c(
    "stability_grid", "silhouette_grid", "negative_silhouette_heatmap",
    "singleton_heatmap", "minimum_cluster_size_heatmap",
    "method_agreement_heatmap", "primary_pam_selection_dossier"
  ),
  width_inches = c(16, 16, 15, 15, 15, 17, 13),
  height_inches = c(10, 10, 10, 10, 10, 12, 10), dpi = 180L,
  format = "png", all_three_representations = c(rep(TRUE, 6L), FALSE),
  H0_H1_simultaneous = TRUE, complete_k2_k10 = TRUE,
  biological_interpretation = FALSE, manuscript_claim = FALSE,
  stringsAsFactors = FALSE
)
review <- data.frame(
  contract_id = "mv10g_review_contract_v1",
  evidence_axis = c(
    "representations", "homology", "methods", "K_grid", "seed_stability",
    "silhouette", "negative_silhouette", "singletons", "minimum_size",
    "method_agreement", "primary_selection", "privacy", "claims"
  ),
  required_display = c(
    "all_three", "H0_H1_separate", "all_five", "2_through_10",
    "mean_and_observed_range", "descriptive_median_and_range",
    "separate_heatmap", "separate_heatmap", "separate_heatmap",
    "all_ten_pairs", "threshold_selected_K_and_silhouette_context",
    "no_sample_ids_labels_or_outcomes", "none"
  ),
  stringsAsFactors = FALSE
)
contract <- data.frame(
  contract_id = "mv10g_clustering_review_prefreeze_v1",
  execution_head = execution_head, source_matrices = receipt$matrices,
  source_partition_fits = receipt$partition_fits, representations = 3L,
  homology_dimensions = 2L, methods = 5L, k_grid = "2:10", figures = 7L,
  source_values_opened_before_metric_and_figure_freeze = FALSE,
  labels_used = FALSE, outcomes_used = FALSE, inference_performed = FALSE,
  ranking_performed = FALSE, combined_score = FALSE,
  H0_H1_combined = FALSE, cell_gene_combined = FALSE,
  biological_claims = FALSE, manuscript_claims = FALSE,
  stringsAsFactors = FALSE
)
decision <- data.frame(
  contract_id = "mv10g_decision_v1",
  decision = "authorize_label_closed_synthesis_after_prefreeze_commit",
  synthesis_authorized_after_commit = TRUE,
  figure_render_authorized = FALSE,
  figure_render_requires = "committed_mv10i_23_of_23_closure",
  scientific_interpretation_authorized = FALSE,
  biological_interpretation_authorized = FALSE,
  manuscript_claims_authorized = FALSE, stringsAsFactors = FALSE
)
validation <- data.frame(
  contract_id = "mv10g_validation_v1",
  check_id = c(
    "production_manifest", "closure_manifest", "production_complete",
    "closure_40_of_40", "closure_interpretation_closed", "four_source_schemas",
    "quality_schema_only", "stability_schema_only", "primary_k_schema_only",
    "agreement_schema_only", "source_rows_bound", "source_hashes",
    "eight_implementation_files", "implementation_hashes",
    "five_source_freezes", "source_freeze_hashes", "eight_metrics",
    "ARI_stability_primary", "silhouette_descriptive", "three_pathology_metrics",
    "method_agreement_metric", "primary_selection_metrics", "six_output_tables",
    "output_row_contracts", "seed_aggregation_frozen", "no_result_filtering",
    "seven_figures", "PNG_only", "figure_dimensions", "all_representations",
    "H0_H1_simultaneous", "complete_k2_k10", "thirteen_review_axes",
    "no_source_value_opening", "label_outcome_firewall",
    "selection_inference_ranking_firewall", "combination_firewall",
    "claim_firewall", "synthesis_only_authorized", "render_closed"
  ),
  passed = c(
    TRUE, TRUE, receipt$completion_state == "complete",
    nrow(closure_validation) == 40L && all(closure_validation$passed),
    closure_decision$result_interpretation_state ==
      "closed_pending_separate_review",
    nrow(source_schema) == 4L,
    source_schema$source_id[[1L]] == "quality",
    source_schema$source_id[[2L]] == "stability",
    source_schema$source_id[[3L]] == "primary_k",
    source_schema$source_id[[4L]] == "agreement",
    identical(as.integer(source_schema$expected_rows),
              c(1350L, 270L, 2L, 2700L)),
    all(grepl("^[0-9a-f]{64}$", source_schema$sha256)),
    nrow(implementation) == 8L, all(file.exists(implementation$file)),
    nrow(source_freeze) == 5L, all(file.exists(source_freeze$artifact)),
    nrow(metric_contract) == 8L,
    metric_contract$role[[1L]] == "five_seed_stability_primary",
    metric_contract$role[[2L]] == "descriptive_internal_quality",
    all(grepl("diagnostic", metric_contract$role[3:5])),
    metric_contract$role[[6L]] == "method_concordance_diagnostic",
    all(metric_contract$selection_metric[c(1L, 7L, 8L)]),
    nrow(aggregation) == 6L,
    identical(aggregation$expected_rows, c(270L, 270L, 540L, 2L, 18L, 90L)),
    all(nzchar(aggregation$seed_aggregation)),
    all(!aggregation$result_dependent_filter), nrow(figures) == 7L,
    all(figures$format == "png"),
    identical(figures$width_inches, c(16, 16, 15, 15, 15, 17, 13)) &&
      identical(figures$height_inches, c(10, 10, 10, 10, 10, 12, 10)),
    all(figures$all_three_representations[1:6]),
    all(figures$H0_H1_simultaneous), all(figures$complete_k2_k10),
    nrow(review) == 13L,
    all(!source_schema$scientific_values_opened_during_prefreeze) &&
      !contract$source_values_opened_before_metric_and_figure_freeze,
    !contract$labels_used && !contract$outcomes_used,
    !contract$inference_performed && !contract$ranking_performed &&
      !contract$combined_score,
    !contract$H0_H1_combined && !contract$cell_gene_combined,
    !contract$biological_claims && !contract$manuscript_claims,
    decision$synthesis_authorized_after_commit,
    !decision$figure_render_authorized
  ),
  stringsAsFactors = FALSE
)
if (!all(validation$passed)) stop("MV10-G prefreeze validation failed")
if (!dir.exists(output)) dir.create(output, recursive = TRUE)
artifacts <- list(
  "mv10g-contract.csv" = contract,
  "mv10g-source-schema.csv" = source_schema,
  "mv10g-metric-contract.csv" = metric_contract,
  "mv10g-aggregation-contract.csv" = aggregation,
  "mv10g-figure-contract.csv" = figures,
  "mv10g-review-contract.csv" = review,
  "mv10g-implementation-bindings.csv" = implementation,
  "mv10g-source-freeze.csv" = source_freeze,
  "mv10g-decision.csv" = decision,
  "mv10g-validation.csv" = validation
)
for (name in names(artifacts)) atomic(artifacts[[name]], file.path(output, name))
writeLines(c(
  "# MV10-G clustering-review prefreeze", "",
  "Metrics, aggregations, seven PNG figures, and review questions are frozen",
  "from schemas and immutable manifests before scientific values are opened.",
  "Only label-closed synthesis is authorized after commit; rendering requires",
  "a separately committed independent MV10-I closure."
), file.path(output, "MV10G_CLUSTERING_REVIEW_PREFREEZE_2026-08-25.md"))
files <- sort(setdiff(list.files(output), "mv10g-artifact-manifest.csv"))
manifest <- data.frame(
  contract_id = "mv10g_artifact_manifest_v1", artifact = files,
  bytes = as.numeric(file.info(file.path(output, files))$size),
  sha256 = vapply(file.path(output, files), sha, character(1L)),
  stringsAsFactors = FALSE
)
atomic(manifest, file.path(output, "mv10g-artifact-manifest.csv"))
cat("Built MV10-G clustering review prefreeze; checks=40\n")
