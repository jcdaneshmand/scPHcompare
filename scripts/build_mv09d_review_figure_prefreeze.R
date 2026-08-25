#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop(paste(
  "usage: build_mv09d_review_figure_prefreeze.R <mv09b-root> <output-dir>"
), call. = FALSE)
production <- normalizePath(args[[1L]], mustWork = TRUE)
output <- args[[2L]]
if (dir.exists(output) && length(list.files(output, all.files = TRUE,
                                            no.. = TRUE))) {
  stop("refusing to overwrite MV9-D prefreeze")
}
if (!dir.exists(output)) dir.create(output, recursive = TRUE)
source("R/mv08z_landscape_production.R")
source("R/mv09d_review_figures.R")
readc <- .mv08z_read_csv; sha <- .mv08z_sha256_file; atomic <- .mv08z_atomic_csv
execution_head <- tolower(trimws(Sys.getenv("MV09D_GIT_HEAD", unset = "")))
if (!grepl("^[0-9a-f]{40}$", execution_head)) {
  stop("MV09D_GIT_HEAD must bind one 40-character commit")
}
.mv08z_verify_manifest(production, "mv09b-artifact-manifest.csv")
closure <- "docs/audits/mv09c-robustness-synthesis-closure-v1"
.mv08z_verify_manifest(closure, "mv09c-artifact-manifest.csv")
closure_validation <- readc(file.path(closure, "mv09c-validation.csv"))
if (!all(closure_validation$passed)) stop("MV9-D closure prerequisite drift")
data <- mv09d_prepare_review_figure_data_v1(production)

contract <- data.frame(
  contract_id = "mv09d_review_figure_prefreeze_v1",
  execution_head = execution_head,
  mv09b_manifest_sha256 = sha(file.path(production,
                                        "mv09b-artifact-manifest.csv")),
  mv09c_manifest_sha256 = sha(file.path(closure,
                                        "mv09c-artifact-manifest.csv")),
  figures = 3L, selected_metrics = 4L,
  metric_selection_timing = "prospective_before_render",
  metric_selection_policy =
    "one linear_global;one_rank_global;one_scaled_residual;one_local_neighbor",
  H0_H1_policy = "simultaneous_explicit_never_combined",
  internal_policy = "show_five_seed_points_plus_median_IQR",
  external_policy = "show_one_cohort_point_no_error_bar_no_generalization",
  dimension_delta_policy = "H1_minus_H0_zero_reference_not_threshold",
  thresholds = "none", rankings = "none", inference = "none",
  combined_score = "forbidden", biological_claims = "closed",
  manuscript_claims = "closed", labels_outcomes = "closed",
  clustering_fusion = "closed", human_review_required = TRUE,
  stringsAsFactors = FALSE
)
metrics <- mv09d_figure_metrics_v1()
metrics$contract_id <- "mv09d_metric_selection_v1"
figures <- data.frame(
  contract_id = "mv09d_figure_contract_v1", figure_order = 1:3,
  figure_id = c("internal_seed_sensitivity", "external_singleton_sensitivity",
                "paired_dimension_shift"),
  data_rows = c(120L, 40L, 80L),
  geometry = c("five seed points plus median and IQR",
               "one cohort point per contrast and dimension; no error bars",
               "paired H1-minus-H0 points with zero reference line"),
  facets = c("four metrics", "four metrics", "four metrics x evidence stratum"),
  width_inches = c(12, 14, 15), height_inches = c(9, 8.5, 11), dpi = 180L,
  format = "PNG_only_no_PDF", visual_review_required = TRUE,
  stringsAsFactors = FALSE
)
review <- data.frame(
  contract_id = "mv09d_review_contract_v1", review_order = 1:6,
  review_item = c("metric_legibility", "H0_H1_visibility",
                  "internal_seed_visibility", "external_singleton_disclosure",
                  "zero_line_not_threshold", "claim_firewall"),
  acceptance = c(
    "four quantities legible under free metric scales",
    "both dimensions shown simultaneously and never pooled",
    "all five seed marks visible with median and IQR",
    "one cohort and absence of error bars stated on figure",
    "subtitle explicitly identifies descriptive reference",
    "no thresholds rankings inference biology or manuscript claims"
  ), owner_review_state = "pending", stringsAsFactors = FALSE
)
implementation_files <- c(
  "R/mv09d_review_figures.R", "scripts/build_mv09d_review_figure_prefreeze.R",
  "scripts/render_mv09e_review_figures.R",
  "scripts/build_mv09f_review_figure_closure.R"
)
implementation <- data.frame(
  contract_id = "mv09d_implementation_binding_v1",
  implementation_order = seq_along(implementation_files),
  file = implementation_files,
  bytes = as.numeric(file.info(implementation_files)$size),
  sha256 = vapply(implementation_files, sha, character(1L)),
  stringsAsFactors = FALSE
)
source_files <- c(
  file.path(production, "mv09b-artifact-manifest.csv"),
  file.path(production, "mv09b-plot-data.csv"),
  file.path(production, "mv09b-internal-seed-summary.csv"),
  file.path(production, "mv09b-external-singleton.csv"),
  file.path(production, "mv09b-dimension-delta.csv"),
  file.path(closure, "mv09c-artifact-manifest.csv"),
  file.path(closure, "mv09c-validation.csv")
)
source_freeze <- data.frame(
  contract_id = "mv09d_source_freeze_v1", source_order = seq_along(source_files),
  artifact = source_files, bytes = as.numeric(file.info(source_files)$size),
  sha256 = vapply(source_files, sha, character(1L)), stringsAsFactors = FALSE
)
validation <- data.frame(
  contract_id = "mv09d_validation_v1",
  check_id = c("mv09b_manifest", "mv09c_closure", "four_metrics",
               "selection_prospective", "internal_120", "internal_summary_24",
               "external_40", "external_singleton", "dimension_delta_80",
               "H0_H1_explicit", "three_figures", "PNG_no_PDF",
               "fixed_dimensions", "implementation_bound", "no_threshold_rank",
               "no_inference_combined_score", "label_outcome_downstream_firewall",
               "human_review_gate"),
  passed = c(TRUE, all(closure_validation$passed), nrow(metrics) == 4L,
             contract$metric_selection_timing == "prospective_before_render",
             nrow(data$internal) == 120L, nrow(data$internal_summary) == 24L,
             nrow(data$external) == 40L,
             all(data$external$replication_units == 1L), nrow(data$delta) == 80L,
             contract$H0_H1_policy == "simultaneous_explicit_never_combined",
             nrow(figures) == 3L, all(figures$format == "PNG_only_no_PDF"),
             all(figures$width_inches > 0 & figures$height_inches > 0 &
                   figures$dpi == 180L), nrow(implementation) == 4L,
             contract$thresholds == "none" && contract$rankings == "none",
             contract$inference == "none" && contract$combined_score == "forbidden",
             contract$labels_outcomes == "closed" &&
               contract$clustering_fusion == "closed",
             contract$human_review_required),
  stringsAsFactors = FALSE
)
if (!all(validation$passed)) stop("MV9-D validation failed")
decision <- data.frame(
  contract_id = "mv09d_decision_v1",
  decision = "authorize_claim_free_review_figure_render_after_commit",
  render_authorized_after_commit = TRUE,
  interpretation_authorized = FALSE,
  next_if_closed = "owner_scientific_review",
  stringsAsFactors = FALSE
)
artifacts <- list(
  "mv09d-contract.csv" = contract, "mv09d-metric-selection.csv" = metrics,
  "mv09d-figure-contract.csv" = figures, "mv09d-review-contract.csv" = review,
  "mv09d-implementation-bindings.csv" = implementation,
  "mv09d-source-freeze.csv" = source_freeze,
  "mv09d-validation.csv" = validation, "mv09d-decision.csv" = decision
)
for (name in names(artifacts)) atomic(artifacts[[name]], file.path(output, name))
writeLines(c(
  "# MV9-D claim-free scientific-review figure prefreeze", "",
  "Four metrics are selected prospectively by measurement role, not observed",
  "magnitude. Three PNG figures preserve five-seed internal and singleton",
  "external evidence while keeping H0/H1 explicit. No PDF is produced.",
  "Rendering and all scientific interpretation remain closed until commit."
), file.path(output, "MV09D_REVIEW_FIGURE_PREFREEZE_2026-08-25.md"))
files <- sort(setdiff(list.files(output), "mv09d-artifact-manifest.csv"))
manifest <- data.frame(
  contract_id = "mv09d_artifact_manifest_v1", artifact = files,
  bytes = as.numeric(file.info(file.path(output, files))$size),
  sha256 = vapply(file.path(output, files), sha, character(1L)),
  stringsAsFactors = FALSE
)
atomic(manifest, file.path(output, "mv09d-artifact-manifest.csv"))
message("Built MV9-D figure prefreeze; checks=18")
