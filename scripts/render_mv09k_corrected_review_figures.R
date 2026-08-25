#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) stop(paste(
  "usage: render_mv09k_corrected_review_figures.R <prefreeze> <mv09b-root>",
  "<mv09i-public> <mv09j-closure> <output-dir> <execution-head>"
), call. = FALSE)
prefreeze <- normalizePath(args[[1L]], mustWork = TRUE)
mv09b <- normalizePath(args[[2L]], mustWork = TRUE)
neighbor <- normalizePath(args[[3L]], mustWork = TRUE)
neighbor_closure <- normalizePath(args[[4L]], mustWork = TRUE)
output <- args[[5L]]
execution_head <- tolower(trimws(args[[6L]]))
if (dir.exists(output) && length(list.files(output, all.files = TRUE,
                                            no.. = TRUE))) {
  stop("refusing to overwrite MV9-K figures")
}
if (!dir.exists(output)) dir.create(output, recursive = TRUE)
source("R/mv08z_landscape_production.R")
source("R/mv09d_review_figures.R")
source("R/mv09h_corrected_review_figures.R")
readc <- .mv08z_read_csv; sha <- .mv08z_sha256_file; atomic <- .mv08z_atomic_csv
.mv08z_verify_manifest(prefreeze, "mv09h-artifact-manifest.csv")
.mv08z_verify_manifest(mv09b, "mv09b-artifact-manifest.csv")
.mv08z_verify_manifest(neighbor, "mv09i-artifact-manifest.csv")
.mv08z_verify_manifest(neighbor_closure, "mv09j-artifact-manifest.csv")
contract <- readc(file.path(prefreeze, "mv09h-contract.csv"))
implementation <- readc(file.path(prefreeze, "mv09h-implementation-bindings.csv"))
decision <- readc(file.path(prefreeze, "mv09h-decision.csv"))
closure_decision <- readc(file.path(neighbor_closure, "mv09j-decision.csv"))
if (execution_head != contract$execution_head ||
    !all(file.exists(implementation$file)) ||
    !all(vapply(implementation$file, sha, character(1L)) ==
           implementation$sha256) ||
    !.mv08z_truth(decision$corrected_figure_render_authorized_after_numeric_closure) ||
    closure_decision$corrected_figure_render_state !=
      "authorized_after_closure_commit") {
  stop("MV9-K render binding drift")
}
data <- mv09h_prepare_corrected_review_data_v1(mv09b, neighbor)
specifications <- mv09h_render_corrected_review_figures_v1(
  data, output, "mv09k"
)
atomic(specifications, file.path(output, "mv09k-figure-specifications.csv"))
atomic(data$metric_contract, file.path(output, "mv09k-metric-contract.csv"))
atomic(data$degeneracy, file.path(output, "mv09k-degeneracy-disclosure.csv"))
review <- data.frame(
  contract_id = "mv09k_review_dossier_v1", review_order = 1:8,
  review_question = c(
    "Is external k=7 clearly identified as structurally non-informative?",
    "Are both prospectively frozen external k=2 and k=3 displayed?",
    "Is internal k=10 evidence preserved and explicitly labeled?",
    "Are global and local metrics visually separated?",
    "Are H0 and H1 simultaneously visible and uncombined?",
    "Is external evidence labeled as one cohort with no inference?",
    "Is zero described as a reference rather than a threshold?",
    "Are rankings, biology, and manuscript claims absent?"
  ),
  required_answer = "yes_before_scientific_interpretation",
  owner_review_state = "pending", stringsAsFactors = FALSE
)
atomic(review, file.path(output, "mv09k-human-review-checklist.csv"))
receipt <- data.frame(
  contract_id = "mv09k_terminal_receipt_v1", execution_head = execution_head,
  completion_state = "complete", figures = 4L, global_metrics = 3L,
  internal_neighbor_k = 10L, external_neighbor_k = "2;3",
  external_k7_displayed = FALSE, external_k7_interpretation = FALSE,
  labels_used = FALSE, outcomes_used = FALSE, inference_performed = FALSE,
  combined_score = FALSE, ranking_performed = FALSE,
  biological_claims = FALSE, manuscript_claims = FALSE,
  human_review_state = "pending", stringsAsFactors = FALSE
)
atomic(receipt, file.path(output, "mv09k-terminal-receipt.csv"))
files <- sort(setdiff(list.files(output), "mv09k-artifact-manifest.csv"))
manifest <- data.frame(
  contract_id = "mv09k_artifact_manifest_v1", artifact = files,
  bytes = as.numeric(file.info(file.path(output, files))$size),
  sha256 = vapply(file.path(output, files), sha, character(1L)),
  stringsAsFactors = FALSE
)
atomic(manifest, file.path(output, "mv09k-artifact-manifest.csv"))
message("Rendered MV9-K corrected review figures; figures=4")
