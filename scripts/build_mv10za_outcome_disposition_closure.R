#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) stop(paste(
  "usage: build_mv10za_outcome_disposition_closure.R <prefreeze>",
  "<mv10t-synthesis> <mv10z-production> <repeat-root> <output>",
  "<execution-head>"
), call. = FALSE)
prefreeze <- normalizePath(args[[1L]], mustWork = TRUE)
synthesis <- normalizePath(args[[2L]], mustWork = TRUE)
production <- normalizePath(args[[3L]], mustWork = TRUE)
repeat_root <- args[[4L]]; output <- args[[5L]]
execution_head <- tolower(trimws(args[[6L]]))
for (path in c(repeat_root, output)) {
  if (dir.exists(path)) stop("MV10-ZA output root already exists")
  dir.create(path, recursive = TRUE)
}
source("R/mv08z_landscape_production.R")
source("R/mv10y_outcome_disposition.R")
readc <- .mv08z_read_csv; sha <- .mv08z_sha256_file; atomic <- .mv08z_atomic_csv
.mv08z_verify_manifest(prefreeze, "mv10y-artifact-manifest.csv")
.mv08z_verify_manifest(synthesis, "mv10t-artifact-manifest.csv")
.mv08z_verify_manifest(production, "mv10z-artifact-manifest.csv")
contract <- readc(file.path(prefreeze, "mv10y-contract.csv"))
receipt <- readc(file.path(production, "mv10z-terminal-receipt.csv"))
if (execution_head != contract$execution_head ||
    receipt$execution_head != execution_head ||
    receipt$completion_state != "complete") stop("MV10-ZA binding drift")
fresh <- mv10y_build_outcome_disposition_v1(
  readc(file.path(synthesis, "mv10t-complete-summary.csv")),
  readc(file.path(synthesis, "mv10t-endpoint-coverage.csv"))
)
mapping <- c(
  primary_envelope = "mv10z-primary-representation-envelope.csv",
  method_envelope = "mv10z-method-envelope.csv",
  context_contrast = "mv10z-context-contrast.csv",
  approach_contrast = "mv10z-approach-contrast.csv",
  disposition = "mv10z-disposition.csv"
)
rehash <- do.call(rbind, lapply(seq_along(mapping), function(i) {
  name <- names(mapping)[[i]]; saved_path <- file.path(production, mapping[[i]])
  repeat_path <- file.path(repeat_root, mapping[[i]])
  atomic(fresh[[name]], repeat_path)
  saved <- readc(saved_path); repeated <- readc(repeat_path)
  data.frame(
    contract_id = "mv10za_rehash_v1", artifact_order = i,
    artifact = mapping[[i]], rows = nrow(saved),
    exact_columns = identical(names(saved), names(repeated)),
    numeric_repeat = isTRUE(all.equal(saved, repeated, tolerance = 1e-14,
                                      check.attributes = FALSE)),
    production_sha256 = sha(saved_path), repeat_sha256 = sha(repeat_path),
    byte_identical = sha(saved_path) == sha(repeat_path),
    stringsAsFactors = FALSE
  )
}))
validation <- data.frame(
  contract_id = "mv10za_validation_v1",
  check_id = c(
    "prefreeze_manifest", "synthesis_manifest", "production_manifest",
    "execution_head", "terminal_complete", "five_outputs",
    "twenty_primary_envelopes", "sixty_method_envelopes",
    "one_hundred_twenty_context_contrasts",
    "one_hundred_twenty_approach_contrasts", "one_disposition",
    "five_rehashes", "exact_columns", "numeric_repeat", "byte_repeat",
    "three_representations", "two_homology_dimensions", "five_methods",
    "five_endpoints", "two_metrics", "no_magnitude_threshold",
    "no_p_values", "no_method_selection", "no_ranking_or_pooling",
    "no_causal_approach_interpretation", "claim_firewall"
  ),
  passed = c(
    TRUE, TRUE, TRUE, receipt$execution_head == execution_head,
    receipt$completion_state == "complete", receipt$output_tables == 5L,
    nrow(fresh$primary_envelope) == 20L,
    nrow(fresh$method_envelope) == 60L,
    nrow(fresh$context_contrast) == 120L,
    nrow(fresh$approach_contrast) == 120L,
    nrow(fresh$disposition) == 1L, nrow(rehash) == 5L,
    all(rehash$exact_columns), all(rehash$numeric_repeat),
    all(rehash$byte_identical),
    length(unique(fresh$context_contrast$stack_id)) == 3L,
    length(unique(fresh$context_contrast$homology_dimension)) == 2L,
    length(unique(fresh$context_contrast$method_id)) == 5L,
    nrow(fresh$primary_envelope) / 4L == 5L,
    length(unique(fresh$context_contrast$metric_id)) == 2L,
    !receipt$magnitude_threshold_applied, !receipt$p_values_computed,
    !receipt$method_selection_executed,
    !receipt$ranking_executed && !receipt$H0_H1_combined,
    !fresh$disposition$approach_causal_interpretation,
    !receipt$biological_claims && !receipt$manuscript_claims
  ), stringsAsFactors = FALSE
)
if (!all(validation$passed)) stop("MV10-ZA closure failed")
decision <- data.frame(
  contract_id = "mv10za_decision_v1",
  decision = "close_descriptive_outcome_disposition",
  frozen_method_state = "retain_PAM_without_outcome_tuning",
  next_stage = "cross_view_descriptive_synthesis_prefreeze_or_stop",
  method_selection_state = "closed", biological_claims_state = "closed",
  manuscript_claims_state = "closed", stringsAsFactors = FALSE
)
atomic(rehash, file.path(output, "mv10za-rehash.csv"))
atomic(validation, file.path(output, "mv10za-validation.csv"))
atomic(decision, file.path(output, "mv10za-decision.csv"))
writeLines(c(
  "# MV10-ZA descriptive outcome-disposition closure", "",
  "All five transparent, value-aware descriptive outputs repeat exactly.",
  "The closure retains frozen PAM without label-driven tuning and authorizes",
  "no threshold, inference, ranking, causal approach interpretation, biology,",
  "or manuscript claim."
), file.path(output, "MV10ZA_OUTCOME_DISPOSITION_CLOSURE_2026-08-25.md"))
files <- sort(setdiff(list.files(output), "mv10za-artifact-manifest.csv"))
manifest <- data.frame(
  contract_id = "mv10za_artifact_manifest_v1", artifact = files,
  bytes = as.numeric(file.info(file.path(output, files))$size),
  sha256 = vapply(file.path(output, files), sha, character(1L)),
  stringsAsFactors = FALSE
)
atomic(manifest, file.path(output, "mv10za-artifact-manifest.csv"))
cat("Closed MV10-ZA descriptive outcome disposition; checks=26\n")
