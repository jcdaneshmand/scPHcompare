#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) stop(paste(
  "usage: build_mv10u_outcome_review_closure.R <prefreeze>",
  "<mv10q-production> <mv10r-closure> <mv10t-synthesis> <output>",
  "<execution-head>"
), call. = FALSE)
prefreeze <- normalizePath(args[[1L]], mustWork = TRUE)
production <- normalizePath(args[[2L]], mustWork = TRUE)
source_closure <- normalizePath(args[[3L]], mustWork = TRUE)
synthesis <- normalizePath(args[[4L]], mustWork = TRUE)
output <- args[[5L]]
execution_head <- tolower(trimws(args[[6L]]))
if (dir.exists(output)) stop("MV10-U output already exists")
source("R/mv08z_landscape_production.R")
source("R/mv10s_outcome_review.R")
readc <- .mv08z_read_csv; sha <- .mv08z_sha256_file; atomic <- .mv08z_atomic_csv
.mv08z_verify_manifest(prefreeze, "mv10s-artifact-manifest.csv")
.mv08z_verify_manifest(production, "mv10q-artifact-manifest.csv")
.mv08z_verify_manifest(source_closure, "mv10r-artifact-manifest.csv")
.mv08z_verify_manifest(synthesis, "mv10t-artifact-manifest.csv")
contract <- readc(file.path(prefreeze, "mv10s-contract.csv"))
receipt <- readc(file.path(synthesis, "mv10t-terminal-receipt.csv"))
if (execution_head != contract$execution_head ||
    receipt$execution_head != execution_head ||
    receipt$completion_state != "complete") stop("MV10-U binding drift")
fresh <- mv10s_build_outcome_review_v1(
  readc(file.path(production, "mv10q-seed-metrics.csv")),
  readc(file.path(production, "mv10q-unit-summaries.csv")),
  readc(file.path(production, "mv10q-structural-status.csv"))
)
mapping <- c(
  complete_summary = "mv10t-complete-summary.csv",
  primary_summary = "mv10t-primary-pam-summary.csv",
  endpoint_coverage = "mv10t-endpoint-coverage.csv"
)
rehash <- do.call(rbind, lapply(seq_along(mapping), function(i) {
  name <- names(mapping)[[i]]; path <- file.path(synthesis, mapping[[i]])
  saved <- readc(path)
  data.frame(
    contract_id = "mv10u_synthesis_rehash_v1", artifact_order = i,
    artifact = mapping[[i]], rows = nrow(saved), bytes = file.info(path)$size,
    sha256 = sha(path), exact_columns = identical(names(saved), names(fresh[[name]])),
    independent_numeric_repeat = isTRUE(all.equal(
      saved, fresh[[name]], tolerance = 1e-14, check.attributes = FALSE
    )), stringsAsFactors = FALSE
  )
}))
validation <- data.frame(
  contract_id = "mv10u_validation_v1",
  check_id = c(
    "prefreeze_manifest", "production_manifest", "source_closure_manifest",
    "synthesis_manifest", "execution_head", "terminal_complete",
    "three_output_tables", "three_hundred_complete_rows",
    "sixty_primary_rows", "six_endpoint_rows", "three_rehashes",
    "exact_columns", "numeric_repeat", "three_representations",
    "two_homology_dimensions", "five_methods", "five_estimable_endpoints",
    "two_metrics", "five_seeds_represented", "selected_K_2_3",
    "one_structural_nonendpoint", "no_p_values", "no_method_selection",
    "no_ranking_or_pooling", "claim_firewall"
  ),
  passed = c(
    TRUE, TRUE, TRUE, TRUE, receipt$execution_head == execution_head,
    receipt$completion_state == "complete", receipt$output_tables == 3L,
    nrow(fresh$complete_summary) == 300L, nrow(fresh$primary_summary) == 60L,
    nrow(fresh$endpoint_coverage) == 6L, nrow(rehash) == 3L,
    all(rehash$exact_columns), all(rehash$independent_numeric_repeat),
    length(unique(fresh$complete_summary$stack_id)) == 3L,
    setequal(fresh$complete_summary$homology_dimension, c("H0", "H1")),
    length(unique(fresh$complete_summary$method_id)) == 5L,
    length(unique(fresh$complete_summary$endpoint_id)) == 5L,
    length(unique(fresh$complete_summary$metric_id)) == 2L,
    all(fresh$complete_summary$completed_seeds == 5L),
    setequal(fresh$complete_summary$selected_k, c(2L, 3L)),
    sum(grepl("structurally_not_estimable",
              fresh$endpoint_coverage$execution_status)) == 1L,
    !receipt$p_values, !receipt$method_selection,
    !receipt$ranking && !receipt$H0_H1_combined,
    !receipt$biological_claims && !receipt$manuscript_claims
  ), stringsAsFactors = FALSE
)
if (!all(validation$passed)) stop("MV10-U closure failed")
decision <- data.frame(
  contract_id = "mv10u_decision_v1",
  decision = "close_complete_outcome_review_synthesis",
  figure_render_authorized_after_commit = TRUE,
  interpretation_state = "closed_pending_exact_figures_and_visual_review",
  method_selection_state = "closed", biological_claims_state = "closed",
  manuscript_claims_state = "closed", stringsAsFactors = FALSE
)
dir.create(output, recursive = TRUE)
atomic(rehash, file.path(output, "mv10u-synthesis-rehash.csv"))
atomic(validation, file.path(output, "mv10u-validation.csv"))
atomic(decision, file.path(output, "mv10u-decision.csv"))
writeLines(c(
  "# MV10-U outcome-review synthesis closure", "",
  "All three prospectively frozen review tables independently reproduce.",
  "Six figures may render only after this closure commits; interpretation,",
  "selection, biological claims, and manuscript claims remain closed."
), file.path(output, "MV10U_OUTCOME_REVIEW_CLOSURE_2026-08-25.md"))
files <- sort(setdiff(list.files(output), "mv10u-artifact-manifest.csv"))
manifest <- data.frame(
  contract_id = "mv10u_artifact_manifest_v1", artifact = files,
  bytes = as.numeric(file.info(file.path(output, files))$size),
  sha256 = vapply(file.path(output, files), sha, character(1L)),
  stringsAsFactors = FALSE
)
atomic(manifest, file.path(output, "mv10u-artifact-manifest.csv"))
cat("Closed MV10-U outcome-review synthesis; checks=25\n")
