#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop(paste(
  "usage: render_mv09e_review_figures.R <prefreeze> <mv09b-root>",
  "<output-dir> <execution-head>"
), call. = FALSE)
prefreeze <- normalizePath(args[[1L]], mustWork = TRUE)
production <- normalizePath(args[[2L]], mustWork = TRUE)
output <- args[[3L]]
execution_head <- tolower(trimws(args[[4L]]))
if (dir.exists(output) && length(list.files(output, all.files = TRUE,
                                            no.. = TRUE))) {
  stop("refusing to overwrite MV9-E figures")
}
if (!dir.exists(output)) dir.create(output, recursive = TRUE)
source("R/mv08z_landscape_production.R")
source("R/mv09d_review_figures.R")
.mv08z_verify_manifest(prefreeze, "mv09d-artifact-manifest.csv")
.mv08z_verify_manifest(production, "mv09b-artifact-manifest.csv")
readc <- .mv08z_read_csv; sha <- .mv08z_sha256_file; atomic <- .mv08z_atomic_csv
contract <- readc(file.path(prefreeze, "mv09d-contract.csv"))
implementation <- readc(file.path(prefreeze, "mv09d-implementation-bindings.csv"))
decision <- readc(file.path(prefreeze, "mv09d-decision.csv"))
if (execution_head != contract$execution_head ||
    sha(file.path(production, "mv09b-artifact-manifest.csv")) !=
      contract$mv09b_manifest_sha256 ||
    !all(file.exists(implementation$file)) ||
    !all(vapply(implementation$file, sha, character(1L)) ==
           implementation$sha256) ||
    !.mv08z_truth(decision$render_authorized_after_commit)) {
  stop("MV9-E render binding drift")
}
data <- mv09d_prepare_review_figure_data_v1(production)
specifications <- mv09d_render_review_figures_v1(data, output, "mv09e")
atomic(specifications, file.path(output, "mv09e-figure-specifications.csv"))
atomic(data$metric_contract, file.path(output, "mv09e-metric-contract.csv"))
review <- data.frame(
  contract_id = "mv09e_review_dossier_v1", review_order = 1:6,
  review_question = c(
    "Are all four metrics legible without sharing incompatible y scales?",
    "Are H0 and H1 simultaneously visible and uncombined?",
    "Are all five internal seeds visible behind median and IQR summaries?",
    "Is external evidence visibly labeled as one cohort with no error bars?",
    "Is the H1-minus-H0 zero line described as a reference, not a threshold?",
    "Are thresholds, rankings, inference, biology, and manuscript claims absent?"
  ),
  required_answer = "yes_before_scientific_interpretation",
  owner_review_state = "pending", stringsAsFactors = FALSE
)
atomic(review, file.path(output, "mv09e-human-review-checklist.csv"))
receipt <- data.frame(
  contract_id = "mv09e_terminal_receipt_v1", execution_head = execution_head,
  completion_state = "complete", figures = 3L, selected_metrics = 4L,
  labels_used = FALSE, outcomes_used = FALSE, inference_performed = FALSE,
  combined_score = FALSE, ranking_performed = FALSE,
  biological_claims = FALSE, manuscript_claims = FALSE,
  human_review_state = "pending", stringsAsFactors = FALSE
)
atomic(receipt, file.path(output, "mv09e-terminal-receipt.csv"))
files <- sort(setdiff(list.files(output), "mv09e-artifact-manifest.csv"))
manifest <- data.frame(
  contract_id = "mv09e_artifact_manifest_v1", artifact = files,
  bytes = as.numeric(file.info(file.path(output, files))$size),
  sha256 = vapply(file.path(output, files), sha, character(1L)),
  stringsAsFactors = FALSE
)
atomic(manifest, file.path(output, "mv09e-artifact-manifest.csv"))
message("Rendered MV9-E claim-free review figures; figures=3")
