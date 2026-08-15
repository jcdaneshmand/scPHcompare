#!/usr/bin/env Rscript

options(warn = 2)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required.")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("usage: validate_mv07a_robustness_confounding_synthesis.R AUDIT_DIR EXPECTED_HEAD OUTPUT",
    call. = FALSE)
}
audit_dir <- args[[1L]]; expected_head <- tolower(trimws(args[[2L]])); output <- args[[3L]]
readc <- function(name) read.csv(file.path(audit_dir, name), stringsAsFactors = FALSE,
  check.names = FALSE)
truth <- function(x) if (is.logical(x)) !is.na(x) & x else tolower(trimws(x)) == "true"
sha <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
required <- c("mv07a-source-freeze.csv", "mv07a-landscape-contract.csv",
  "mv07a-robustness-coverage.csv", "mv07a-confounding-coverage.csv",
  "mv07a-computation-budget.csv", "mv07a-evidence-gap-registry.csv",
  "mv07a-claim-boundaries.csv", "mv07a-selection-firewall.csv",
  "mv07a-acceptance-criteria.csv", "mv07a-decision.csv")
if (any(!file.exists(file.path(audit_dir, required)))) stop("MV7-A outputs incomplete.")
s <- readc(required[[1L]]); l <- readc(required[[2L]]); r <- readc(required[[3L]])
c <- readc(required[[4L]]); b <- readc(required[[5L]]); g <- readc(required[[6L]])
cl <- readc(required[[7L]]); f <- readc(required[[8L]]); a <- readc(required[[9L]])
d <- readc(required[[10L]])
source_ok <- nrow(s) == 23L && !anyDuplicated(s$source_id) &&
  all(file.exists(s$artifact_locator)) &&
  all(vapply(seq_len(nrow(s)), function(i) identical(sha(s$artifact_locator[[i]]),
    s$sha256[[i]]), logical(1L))) && all(s$accepted_head == expected_head) &&
  !any(truth(s$consumed_by_continuation_decision))
landscape_expected <- c("all_finite_positive_persistence_intervals",
  "exclude_infinite_interval", "all_consecutive_active_levels",
  "exact_or_error_controlled_squared_l2_on_dimension_support", "h0_h1_separate",
  "no_universal_fixed_grid", "no_universal_level_cap",
  "unweighted_h0_h1_euclidean_composite_descriptive_only")
landscape_ok <- nrow(l) == 8L && identical(l$required_state, landscape_expected) &&
  all(truth(l$preserved)) && !any(truth(l$changed_by_mv07a))
robustness_ok <- nrow(r) == 14L && !anyDuplicated(r$axis_id) &&
  identical(r$axis_order, 1:14) &&
  r$coverage[r$axis_id == "gene_panel_size"] == "gap" &&
  r$coverage[r$axis_id == "homology_dimension"] == "complete"
confounding_ok <- nrow(c) == 10L && !anyDuplicated(c$axis_id) &&
  identical(c$axis_order, 1:10) &&
  c$coverage[c$axis_id == "library_size"] == "unavailable" &&
  !any(truth(c$causal_adjustment_authorized))
budget_ok <- nrow(b) == 5L && all(is.finite(b$worker_or_wall_hours)) &&
  all(b$worker_or_wall_hours > 0) &&
  !truth(b$new_ph[b$workload == "mv07b_no_new_ph_diagnostics"])
gap_ok <- nrow(g) == 6L && identical(g$prerequisite_order, 1:6) &&
  identical(which(truth(g$next_sprint_eligible)), 1:3) &&
  !any(truth(g$requires_new_ph[1:3])) &&
  !any(truth(g$result_value_used_for_priority))
claim_ok <- nrow(cl) == 9L && !anyDuplicated(cl$claim_family) &&
  cl$current_status[cl$claim_family == "fusion"] == "negative" &&
  cl$current_status[cl$claim_family == "external_generalization"] == "prohibited"
firewall_ok <- nrow(f) == 8L && !any(truth(f$consumed_by_continuation_decision))
criteria_ok <- nrow(a) == 10L && all(truth(a$passed))
decision_ok <- nrow(d) == 1L && d$authorized_next_sprint == "MV7-B" &&
  d$numerical_result_rows_consumed == 0L && truth(d$prefreeze_only) &&
  !truth(d$new_ph_authorized) && !truth(d$new_data_authorized) &&
  !truth(d$method_or_weight_selection_authorized) &&
  !truth(d$advanced_fusion_authorized) && !truth(d$default_change_authorized) &&
  !truth(d$manuscript_claim_promotion_authorized)
checks <- data.frame(
  contract_id = "mv07a_independent_validation_v1",
  category = c("source_freeze", "landscape_contract", "robustness_registry",
    "confounding_registry", "computation_budget", "gap_order", "claim_boundaries",
    "selection_firewall", "acceptance_criteria", "continuation_decision"),
  passed = c(source_ok, landscape_ok, robustness_ok, confounding_ok, budget_ok,
    gap_ok, claim_ok, firewall_ok, criteria_ok, decision_ok),
  detail = c("23 source identities independently rehashed", "8 fixed items preserved",
    "14 unique robustness axes", "10 unique confounding axes",
    "four measured workloads plus bounded no-PH diagnostic", "first three gaps reuse existing artifacts",
    "nine separate claim families", "zero result inputs reach decision",
    "10 of 10 criteria pass", "MV7-B prefreeze only; all expansion false"),
  stringsAsFactors = FALSE)
if (!all(checks$passed)) stop("MV7-A independent validation failed: ",
  paste(checks$category[!checks$passed], collapse = ", "), call. = FALSE)
if (file.exists(output)) stop("Refusing to overwrite: ", output, call. = FALSE)
write.table(checks, output, sep = ",", row.names = FALSE, col.names = TRUE,
  quote = TRUE, na = "")
message("MV7-A independent validation: 10/10 categories pass")
