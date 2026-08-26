#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) stop(paste(
  "usage: build_mv10i_clustering_review_closure.R <prefreeze>",
  "<mv10e-production> <mv10f-closure> <mv10h-synthesis> <output>",
  "<execution-head>"
), call. = FALSE)
prefreeze <- normalizePath(args[[1L]], mustWork = TRUE)
production <- normalizePath(args[[2L]], mustWork = TRUE)
source_closure <- normalizePath(args[[3L]], mustWork = TRUE)
synthesis <- normalizePath(args[[4L]], mustWork = TRUE)
output <- args[[5L]]
execution_head <- tolower(trimws(args[[6L]]))
if (dir.exists(output)) stop("MV10-I output already exists")
source("R/mv08z_landscape_production.R")
source("R/mv10_clustering_benchmark.R")
source("R/mv10g_clustering_review.R")
readc <- .mv08z_read_csv; sha <- .mv08z_sha256_file; atomic <- .mv08z_atomic_csv
.mv08z_verify_manifest(prefreeze, "mv10g-artifact-manifest.csv")
.mv08z_verify_manifest(production, "mv10e-artifact-manifest.csv")
.mv08z_verify_manifest(source_closure, "mv10f-artifact-manifest.csv")
.mv08z_verify_manifest(synthesis, "mv10h-artifact-manifest.csv")
contract <- readc(file.path(prefreeze, "mv10g-contract.csv"))
receipt <- readc(file.path(synthesis, "mv10h-terminal-receipt.csv"))
if (execution_head != contract$execution_head ||
    receipt$execution_head != execution_head ||
    receipt$completion_state != "complete") stop("MV10-I binding drift")
fresh <- mv10g_build_review_data_v1(
  readc(file.path(production, "mv10e-partition-quality.csv")),
  readc(file.path(production, "mv10e-seed-stability.csv")),
  readc(file.path(production, "mv10e-primary-k-selection.csv")),
  readc(file.path(production, "mv10e-method-agreement.csv"))
)
mapping <- c(
  stability_grid = "mv10h-stability-grid.csv",
  quality_summary = "mv10h-quality-summary.csv",
  agreement_summary = "mv10h-method-agreement-summary.csv",
  primary_selection = "mv10h-primary-selection.csv",
  primary_stability = "mv10h-primary-stability.csv",
  primary_quality = "mv10h-primary-quality.csv"
)
rehash <- do.call(rbind, lapply(seq_along(mapping), function(i) {
  name <- names(mapping)[[i]]; path <- file.path(synthesis, mapping[[i]])
  saved <- readc(path)
  equal <- isTRUE(all.equal(saved, fresh[[name]], tolerance = 1e-14,
                            check.attributes = FALSE))
  data.frame(
    contract_id = "mv10i_synthesis_rehash_v1", artifact_order = i,
    artifact = mapping[[i]], rows = nrow(saved), bytes = file.info(path)$size,
    sha256 = sha(path), independent_numeric_repeat = equal,
    stringsAsFactors = FALSE
  )
}))
validation <- data.frame(
  contract_id = "mv10i_validation_v1",
  check_id = c(
    "prefreeze_manifest", "production_manifest", "source_closure_manifest",
    "synthesis_manifest", "execution_head", "terminal_complete",
    "six_output_tables", "two_hundred_seventy_stability_rows",
    "two_hundred_seventy_quality_rows", "five_hundred_forty_agreement_rows",
    "two_primary_selection_rows", "eighteen_primary_stability_rows",
    "ninety_primary_quality_rows", "six_independent_rehashes",
    "all_numeric_repeats", "complete_three_representations",
    "complete_five_methods", "complete_k2_k10", "H0_H1_separate",
    "label_outcome_firewall", "selection_inference_ranking_firewall",
    "combination_firewall", "claim_firewall"
  ),
  passed = c(
    TRUE, TRUE, TRUE, TRUE, receipt$execution_head == execution_head,
    receipt$completion_state == "complete", receipt$output_tables == 6L,
    nrow(fresh$stability_grid) == 270L,
    nrow(fresh$quality_summary) == 270L,
    nrow(fresh$agreement_summary) == 540L,
    nrow(fresh$primary_selection) == 2L,
    nrow(fresh$primary_stability) == 18L,
    nrow(fresh$primary_quality) == 90L,
    nrow(rehash) == 6L, all(rehash$independent_numeric_repeat),
    length(unique(fresh$stability_grid$stack_id)) == 3L,
    length(unique(fresh$stability_grid$method_id)) == 5L,
    identical(sort(unique(fresh$stability_grid$k)), 2:10),
    setequal(fresh$stability_grid$homology_dimension, c("H0", "H1")),
    !receipt$labels_used && !receipt$outcomes_used,
    !receipt$inference_performed && !receipt$ranking_performed &&
      !receipt$combined_score,
    !receipt$H0_H1_combined && !receipt$cell_gene_combined,
    !receipt$biological_claims && !receipt$manuscript_claims
  ),
  stringsAsFactors = FALSE
)
if (!all(validation$passed)) stop("MV10-I closure failed")
decision <- data.frame(
  contract_id = "mv10i_decision_v1",
  decision = "close_label_closed_review_synthesis",
  figure_render_authorized_after_commit = TRUE,
  interpretation_state = "closed_pending_exact_figures_and_visual_review",
  biological_claims_state = "closed", manuscript_claims_state = "closed",
  stringsAsFactors = FALSE
)
dir.create(output, recursive = TRUE)
atomic(rehash, file.path(output, "mv10i-synthesis-rehash.csv"))
atomic(validation, file.path(output, "mv10i-validation.csv"))
atomic(decision, file.path(output, "mv10i-decision.csv"))
writeLines(c(
  "# MV10-I clustering-review synthesis closure", "",
  "All six label-closed review tables independently reproduce from the",
  "immutable MV10-E aggregates. Seven prospectively frozen figures may render",
  "only after this closure is committed. Interpretation remains closed."
), file.path(output, "MV10I_CLUSTERING_REVIEW_CLOSURE_2026-08-25.md"))
files <- sort(setdiff(list.files(output), "mv10i-artifact-manifest.csv"))
manifest <- data.frame(
  contract_id = "mv10i_artifact_manifest_v1", artifact = files,
  bytes = as.numeric(file.info(file.path(output, files))$size),
  sha256 = vapply(file.path(output, files), sha, character(1L)),
  stringsAsFactors = FALSE
)
atomic(manifest, file.path(output, "mv10i-artifact-manifest.csv"))
cat("Closed MV10-I clustering review synthesis; checks=23\n")
