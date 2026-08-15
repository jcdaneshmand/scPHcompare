#!/usr/bin/env Rscript

options(warn = 2)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required.")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("usage: validate_mv07d_full_corpus_reconciliation.R AUDIT_DIR ",
       "RETAINED_CSV EXPECTED_HEAD OUTPUT", call. = FALSE)
}
audit_dir <- args[[1L]]
retained_path <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
expected_head <- tolower(trimws(args[[3L]])); output <- args[[4L]]
readc <- function(name) read.csv(file.path(audit_dir, name), stringsAsFactors = FALSE,
                                 check.names = FALSE)
truth <- function(x) if (is.logical(x)) !is.na(x) & x else
  tolower(trimws(as.character(x))) == "true"
sha <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
required <- c("mv07d-source-freeze.csv", "mv07d-sample-reconciliation.csv",
  "mv07d-tissue-study-summary.csv", "mv07d-estimand-populations.csv",
  "mv07d-source-coverage.csv", "mv07d-artifact-coverage.csv",
  "mv07d-sentinel-prefreeze.csv", "mv07d-landscape-contract.csv",
  "mv07d-code-rule-trace.csv", "mv07d-gap-register.csv",
  "mv07d-acceptance-criteria.csv", "mv07d-decision.csv")
if (any(!file.exists(file.path(audit_dir, required)))) stop("MV7-D outputs incomplete.")
s <- readc(required[[1L]]); x <- readc(required[[2L]]); t <- readc(required[[3L]])
p <- readc(required[[4L]]); c <- readc(required[[5L]]); a <- readc(required[[6L]])
z <- readc(required[[7L]]); l <- readc(required[[8L]]); r <- readc(required[[9L]])
g <- readc(required[[10L]]); q <- readc(required[[11L]]); d <- readc(required[[12L]])

public_rows <- s[!truth(s$private_source), , drop = FALSE]
source_ok <- nrow(s) == 12L && all(s$accepted_head == expected_head) &&
  all(file.exists(public_rows$artifact_locator)) &&
  all(vapply(seq_len(nrow(public_rows)), function(i)
    identical(sha(public_rows$artifact_locator[[i]]), public_rows$sha256[[i]]),
    logical(1L))) &&
  identical(sha(retained_path), s$sha256[s$source_id == "retained_metadata"])
axis_ok <- nrow(x) == 127L && !anyDuplicated(x$sample_id) &&
  sum(truth(x$historical_retained_124)) == 124L &&
  sum(truth(x$corrected_primary_90)) == 90L &&
  sum(truth(x$corrected_descriptive_124)) == 124L &&
  sum(truth(x$threshold_sensitivity_only)) == 3L &&
  sum(x$corpus_class == "single_study_tissue_descriptive_only") == 34L &&
  sum(truth(x$approach_metadata_conflict)) == 16L
tissue_expected <- c("bone marrow" = 31L, "colon" = 13L, "liver" = 6L,
  "pancreatic islets" = 12L, "pbmc" = 12L, "prostate" = 16L,
  "substantia nigra" = 9L, "testis" = 28L)
tissue_ok <- nrow(t) == 8L && identical(stats::setNames(t$candidate_samples, t$tissue),
                                       tissue_expected) &&
  sum(truth(t$primary_cross_study_eligible)) == 5L &&
  sum(t$primary_samples) == 90L
population_ok <- nrow(p) == 5L && identical(p$samples, c(127L,124L,90L,124L,3L)) &&
  identical(which(truth(p$primary_claim_eligible)), 3L)
coverage_ok <- nrow(c) == 3L && all(c$expected_samples == c$observed_samples) &&
  identical(c$observed_samples, c(127L, 127L, 90L)) && nrow(a) == 127L
sentinel_ok <- nrow(z) == 6L && !anyDuplicated(z$sample_id) &&
  all(z$selected_cells == 384L) && all(z$seed == 20260805L) &&
  setequal(z$selection_boundary, c("minimum_post_qc_cells", "maximum_post_qc_cells")) &&
  length(unique(z$tissue)) == 3L && !any(truth(z$ph_authorized))
landscape_ok <- nrow(l) == 8L && all(truth(l$applies_to_full_corpus_expansion)) &&
  !any(truth(l$changed_by_mv07d)) &&
  all(c("all_consecutive_active_levels", "h0_h1_separate",
        "no_universal_fixed_grid", "no_universal_level_cap") %in% l$required_state)
trace_ok <- nrow(r) == 8L && !any(truth(r$changed_by_mv07d)) &&
  "bounded_knn" %in% r$rule_id
gap_ok <- nrow(g) == 7L && identical(g$gap_order, 1:7) &&
  !any(truth(g$can_change_primary_90_result))
criteria_ok <- nrow(q) == 11L && all(truth(q$passed))
decision_ok <- nrow(d) == 1L &&
  d$decision == "authorize_six_sample_source_sct_feasibility_only" &&
  !truth(d$primary_90_recalculation_authorized) &&
  !truth(d$omitted_34_ph_authorized) && !truth(d$omitted_34_outcome_authorized) &&
  !truth(d$below_250_sensitivity_authorized) && !truth(d$new_data_authorized)
checks <- data.frame(
  contract_id = "mv07d_independent_validation_v1",
  category = c("source_freeze", "sample_axis", "tissue_study_axis",
    "estimand_populations", "source_and_artifact_coverage", "sentinel_prefreeze",
    "landscape_contract", "code_rule_trace", "gap_register",
    "acceptance_criteria", "decision_gate"),
  passed = c(source_ok, axis_ok, tissue_ok, population_ok, coverage_ok,
    sentinel_ok, landscape_ok, trace_ok, gap_ok, criteria_ok, decision_ok),
  detail = c("12 identities including private retained metadata",
    "127 candidates 124 retained 90 primary 34 descriptive-only 3 threshold-only",
    "eight tissues and five multi-study primary tissues",
    "five non-interchangeable populations", "127 existing source pairs and 90 corrected",
    "six depth-extreme sentinels at 384 cells", "eight fixed landscape requirements",
    "eight code and design rules", "seven prospectively ordered gaps",
    "11 of 11 criteria pass", "source/SCT sentinel feasibility only"),
  stringsAsFactors = FALSE)
if (!all(checks$passed)) stop("MV7-D independent validation failed: ",
  paste(checks$category[!checks$passed], collapse = ", "), call. = FALSE)
if (file.exists(output)) stop("Refusing to overwrite: ", output, call. = FALSE)
write.table(checks, output, sep = ",", row.names = FALSE, col.names = TRUE,
            quote = TRUE, na = "")
message("MV7-D independent validation: 11/11 categories pass")
