#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) stop(paste(
  "usage: render_mv10v_outcome_review_figures.R <prefreeze>",
  "<mv10q-production> <mv10r-closure> <mv10t-synthesis>",
  "<mv10u-closure> <output> <execution-head>"
), call. = FALSE)
prefreeze <- normalizePath(args[[1L]], mustWork = TRUE)
production <- normalizePath(args[[2L]], mustWork = TRUE)
source_closure <- normalizePath(args[[3L]], mustWork = TRUE)
synthesis <- normalizePath(args[[4L]], mustWork = TRUE)
synthesis_closure <- normalizePath(args[[5L]], mustWork = TRUE)
output <- args[[6L]]
execution_head <- tolower(trimws(args[[7L]]))
if (dir.exists(output)) stop("MV10-V output already exists")
for (package in c("ggplot2", "png")) {
  if (!requireNamespace(package, quietly = TRUE)) stop(package, " is required")
}
source("R/mv08z_landscape_production.R")
source("R/mv10s_outcome_review.R")
readc <- .mv08z_read_csv; sha <- .mv08z_sha256_file; atomic <- .mv08z_atomic_csv
.mv08z_verify_manifest(prefreeze, "mv10s-artifact-manifest.csv")
.mv08z_verify_manifest(production, "mv10q-artifact-manifest.csv")
.mv08z_verify_manifest(source_closure, "mv10r-artifact-manifest.csv")
.mv08z_verify_manifest(synthesis, "mv10t-artifact-manifest.csv")
.mv08z_verify_manifest(synthesis_closure, "mv10u-artifact-manifest.csv")
contract <- readc(file.path(prefreeze, "mv10s-contract.csv"))
implementation <- readc(file.path(prefreeze, "mv10s-implementation-bindings.csv"))
decision <- readc(file.path(synthesis_closure, "mv10u-decision.csv"))
if (execution_head != contract$execution_head ||
    !.mv08z_truth(decision$figure_render_authorized_after_commit) ||
    !all(file.exists(implementation$file)) ||
    !all(vapply(implementation$file, sha, character(1L)) ==
           implementation$sha256)) stop("MV10-V render binding drift")
data <- mv10s_build_outcome_review_v1(
  readc(file.path(production, "mv10q-seed-metrics.csv")),
  readc(file.path(production, "mv10q-unit-summaries.csv")),
  readc(file.path(production, "mv10q-structural-status.csv"))
)
dir.create(output, recursive = TRUE)
specifications <- mv10s_render_outcome_figures_v1(data, output, "mv10v")
atomic(specifications, file.path(output, "mv10v-figure-specifications.csv"))
review <- data.frame(
  contract_id = "mv10v_human_review_checklist_v1", review_order = 1:12,
  review_question = c(
    "Are all three representations present?",
    "Are H0 and H1 displayed separately?",
    "Are both ARI and max-normalized NMI present?",
    "Are five-seed ranges visible?",
    "Are tissue, study, and approach shown for full 124?",
    "Is approach explicitly marked as confounded diagnostic?",
    "Are tissue and study shown for primary 90?",
    "Is primary-90 approach explicitly non-estimable?",
    "Are all five frozen methods included in sensitivity panels?",
    "Are method-specific retuning and ranking absent?",
    "Are sample identities and label values absent?",
    "Are biological and manuscript claims absent?"
  ), required_answer = "yes_before_scientific_interpretation",
  review_state = "pending", stringsAsFactors = FALSE
)
atomic(review, file.path(output, "mv10v-human-review-checklist.csv"))
receipt <- data.frame(
  contract_id = "mv10v_terminal_receipt_v1", execution_head = execution_head,
  completion_state = "complete", figures = 6L, representations = 3L,
  homology_dimensions = 2L, methods = 5L, estimable_endpoints = 5L,
  metrics = 2L, p_values = FALSE, method_selection = FALSE,
  ranking = FALSE, H0_H1_combined = FALSE, biological_claims = FALSE,
  manuscript_claims = FALSE, human_review_state = "pending",
  stringsAsFactors = FALSE
)
atomic(receipt, file.path(output, "mv10v-terminal-receipt.csv"))
files <- sort(setdiff(list.files(output), "mv10v-artifact-manifest.csv"))
manifest <- data.frame(
  contract_id = "mv10v_artifact_manifest_v1", artifact = files,
  bytes = as.numeric(file.info(file.path(output, files))$size),
  sha256 = vapply(file.path(output, files), sha, character(1L)),
  stringsAsFactors = FALSE
)
atomic(manifest, file.path(output, "mv10v-artifact-manifest.csv"))
cat("Rendered MV10-V outcome-review figures; figures=6\n")
