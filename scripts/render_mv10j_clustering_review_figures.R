#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) stop(paste(
  "usage: render_mv10j_clustering_review_figures.R <prefreeze>",
  "<mv10e-production> <mv10f-closure> <mv10h-synthesis>",
  "<mv10i-closure> <output> <execution-head>"
), call. = FALSE)
prefreeze <- normalizePath(args[[1L]], mustWork = TRUE)
production <- normalizePath(args[[2L]], mustWork = TRUE)
source_closure <- normalizePath(args[[3L]], mustWork = TRUE)
synthesis <- normalizePath(args[[4L]], mustWork = TRUE)
synthesis_closure <- normalizePath(args[[5L]], mustWork = TRUE)
output <- args[[6L]]
execution_head <- tolower(trimws(args[[7L]]))
if (dir.exists(output)) stop("MV10-J output already exists")
for (package in c("ggplot2", "patchwork", "png")) {
  if (!requireNamespace(package, quietly = TRUE)) stop(package, " is required")
}
source("R/mv08z_landscape_production.R")
source("R/mv10_clustering_benchmark.R")
source("R/mv10g_clustering_review.R")
readc <- .mv08z_read_csv; sha <- .mv08z_sha256_file; atomic <- .mv08z_atomic_csv
.mv08z_verify_manifest(prefreeze, "mv10g-artifact-manifest.csv")
.mv08z_verify_manifest(production, "mv10e-artifact-manifest.csv")
.mv08z_verify_manifest(source_closure, "mv10f-artifact-manifest.csv")
.mv08z_verify_manifest(synthesis, "mv10h-artifact-manifest.csv")
.mv08z_verify_manifest(synthesis_closure, "mv10i-artifact-manifest.csv")
contract <- readc(file.path(prefreeze, "mv10g-contract.csv"))
implementation <- readc(file.path(prefreeze, "mv10g-implementation-bindings.csv"))
decision <- readc(file.path(synthesis_closure, "mv10i-decision.csv"))
if (execution_head != contract$execution_head ||
    !.mv08z_truth(decision$figure_render_authorized_after_commit) ||
    !all(file.exists(implementation$file)) ||
    !all(vapply(implementation$file, sha, character(1L)) ==
           implementation$sha256)) stop("MV10-J render binding drift")
data <- mv10g_build_review_data_v1(
  readc(file.path(production, "mv10e-partition-quality.csv")),
  readc(file.path(production, "mv10e-seed-stability.csv")),
  readc(file.path(production, "mv10e-primary-k-selection.csv")),
  readc(file.path(production, "mv10e-method-agreement.csv"))
)
dir.create(output, recursive = TRUE)
specifications <- mv10g_render_review_figures_v1(data, output, "mv10j")
atomic(specifications, file.path(output, "mv10j-figure-specifications.csv"))
review <- data.frame(
  contract_id = "mv10j_human_review_checklist_v1", review_order = 1:12,
  review_question = c(
    "Are all three representations displayed without omission?",
    "Are H0 and H1 simultaneously visible and never combined?",
    "Are all five authorized clustering methods displayed?",
    "Is the complete K=2:10 grid retained?",
    "Are five-seed stability ranges visible?",
    "Is silhouette explicitly descriptive rather than a K-selection rule?",
    "Are negative-silhouette, singleton, and minimum-size diagnostics separate?",
    "Are all ten pairwise method agreements displayed?",
    "Does the primary dossier show the frozen threshold and selected K?",
    "Are sample identities, labels, and outcomes absent?",
    "Are H0/H1 and cell/gene combined scores absent?",
    "Are biological and manuscript claims absent?"
  ),
  required_answer = "yes_before_scientific_interpretation",
  owner_review_state = "pending", stringsAsFactors = FALSE
)
atomic(review, file.path(output, "mv10j-human-review-checklist.csv"))
receipt <- data.frame(
  contract_id = "mv10j_terminal_receipt_v1", execution_head = execution_head,
  completion_state = "complete", figures = 7L, representations = 3L,
  homology_dimensions = 2L, methods = 5L, k_values = 9L,
  labels_used = FALSE, outcomes_used = FALSE, inference_performed = FALSE,
  ranking_performed = FALSE, combined_score = FALSE,
  H0_H1_combined = FALSE, cell_gene_combined = FALSE,
  biological_claims = FALSE, manuscript_claims = FALSE,
  human_review_state = "pending", stringsAsFactors = FALSE
)
atomic(receipt, file.path(output, "mv10j-terminal-receipt.csv"))
files <- sort(setdiff(list.files(output), "mv10j-artifact-manifest.csv"))
manifest <- data.frame(
  contract_id = "mv10j_artifact_manifest_v1", artifact = files,
  bytes = as.numeric(file.info(file.path(output, files))$size),
  sha256 = vapply(file.path(output, files), sha, character(1L)),
  stringsAsFactors = FALSE
)
atomic(manifest, file.path(output, "mv10j-artifact-manifest.csv"))
cat("Rendered MV10-J clustering review figures; figures=7\n")
