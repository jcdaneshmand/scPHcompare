#!/usr/bin/env Rscript

options(warn = 2)
for (package in c("digest", "pkgload")) {
  if (!requireNamespace(package, quietly = TRUE)) stop(package, " required.")
}
pkgload::load_all(".", quiet = TRUE, export_all = TRUE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("usage: build_mv05am_four_panel_synthesis_gate.R AUDIT_DIR EXPECTED_HEAD",
       call. = FALSE)
}
audit_dir <- args[[1L]]; expected_head <- tolower(trimws(args[[2L]]))
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
readc <- function(path) read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
sha <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
write_once <- function(value, name) {
  path <- file.path(audit_dir, name)
  if (file.exists(path)) stop("Refusing to overwrite: ", path, call. = FALSE)
  write_provenance_csv(value, path)
}
if (!grepl("^[0-9a-f]{40}$", expected_head)) stop("Full EXPECTED_HEAD required.")
head <- trimws(system2("git", c("rev-parse", "HEAD"), stdout = TRUE))
if (!identical(head, expected_head)) {
  stop("MV5-AM requires prospective engine HEAD ", expected_head,
       "; observed ", head, ".", call. = FALSE)
}

panel_prefix <- c(pc20_vs_pc30 = "mv05z",
  cosine_chord_vs_euclidean = "mv05ad",
  nested192_vs_384_cells = "mv05ah", nested256_vs_384_cells = "mv05al")
report_names <- c(
  pc20_vs_pc30 = "MV05Z_PC20_RETRIEVAL_ROBUSTNESS_EXECUTION_2026-08-11.md",
  cosine_chord_vs_euclidean = "MV05AD_COSINE_RETRIEVAL_ROBUSTNESS_EXECUTION_2026-08-11.md",
  nested192_vs_384_cells = "MV05AH_NESTED192_RETRIEVAL_ROBUSTNESS_EXECUTION_2026-08-11.md",
  nested256_vs_384_cells = "MV05AL_NESTED256_RETRIEVAL_ROBUSTNESS_EXECUTION_2026-08-11.md")
paths <- c(
  am_code = "R/mv05am_four_panel_synthesis_gate.R",
  am_builder = "scripts/build_mv05am_four_panel_synthesis_gate.R",
  am_validator = "scripts/validate_mv05am_four_panel_synthesis_gate.R",
  am_tests = "tests/testthat/test-mv05am-four-panel-synthesis-gate.R",
  am_spec = "docs/specifications/MV05AM_SELECTION_RESISTANT_FOUR_PANEL_SYNTHESIS_GATE_SPECIFICATION_V1.md",
  v_queue = "docs/audits/mv05v-full-group-queue-2026-08-10.csv",
  ai_decision = "docs/audits/mv05ai-continuation-decision-2026-08-11.csv",
  landscape_exact = "scripts/mv05w_exact_landscape_group.py",
  landscape_reference = "R/landscape_reference.R",
  landscape_contract = "R/landscape_contract.R")
for (panel_id in names(panel_prefix)) {
  prefix <- panel_prefix[[panel_id]]
  paths[paste0(panel_id, "__production")] <-
    paste0("docs/audits/", prefix, "-production-summary-2026-08-11.csv")
  paths[paste0(panel_id, "__macro")] <-
    paste0("docs/audits/", prefix, "-macro-estimands-2026-08-11.csv")
  paths[paste0(panel_id, "__intervals")] <-
    paste0("docs/audits/", prefix, "-estimand-intervals-2026-08-11.csv")
  paths[paste0(panel_id, "__primary")] <-
    paste0("docs/audits/", prefix, "-primary-contrasts-2026-08-11.csv")
  paths[paste0(panel_id, "__validation")] <-
    paste0("docs/audits/", prefix, "-outcome-independent-validation-2026-08-11.csv")
  paths[paste0(panel_id, "__prediction_lock")] <-
    paste0("docs/audits/", prefix, "-prediction-lock-2026-08-11.csv")
  paths[paste0(panel_id, "__report")] <-
    file.path("docs/audits", report_names[[panel_id]])
}
if (any(!file.exists(paths))) {
  stop("MV5-AM source set incomplete: ",
       paste(names(paths)[!file.exists(paths)], collapse = ", "), call. = FALSE)
}
source <- data.frame(
  contract_id = "mv05am_source_freeze_v1", source_id = names(paths),
  artifact_locator = unname(paths), sha256 = vapply(paths, sha, character(1L)),
  bytes = as.numeric(file.info(paths)$size), accepted_head = expected_head,
  numerical_result_source = grepl("__(macro|intervals|primary)$", names(paths)),
  consumed_by_continuation_helper = FALSE, stringsAsFactors = FALSE)

registry <- mv05am_panel_registry_v1()
bindings <- list(); macros <- list(); intervals <- list(); primaries <- list()
for (panel_id in registry$panel_id) {
  get_path <- function(kind) paths[[paste0(panel_id, "__", kind)]]
  production <- readc(get_path("production")); macro <- readc(get_path("macro"))
  interval <- readc(get_path("intervals")); primary <- readc(get_path("primary"))
  validation <- readc(get_path("validation")); lock <- readc(get_path("prediction_lock"))
  if (!all(.mv05am_true(validation$passed)) || nrow(lock) != 1L ||
      !.mv05am_true(production$outcomes_computed) ||
      !.mv05am_true(production$evaluation_executed)) {
    stop("MV5-AM accepted panel validation or execution identity failed.",
         call. = FALSE)
  }
  bindings[[panel_id]] <- mv05am_bind_complete_panel_v1(
    panel_id, production, macro, interval, primary)
  macros[[panel_id]] <- mv05am_stack_panel_v1(panel_id, macro, "macro")
  intervals[[panel_id]] <- mv05am_stack_panel_v1(panel_id, interval, "interval")
  primaries[[panel_id]] <- mv05am_stack_panel_v1(panel_id, primary, "primary")
}
bindings <- do.call(rbind, bindings); macros <- do.call(rbind, macros)
intervals <- do.call(rbind, intervals); primaries <- do.call(rbind, primaries)
rownames(bindings) <- rownames(macros) <- rownames(intervals) <- rownames(primaries) <- NULL
if (nrow(bindings) != 4L || nrow(macros) != 96L ||
    nrow(intervals) != 96L || nrow(primaries) != 16L) {
  stop("MV5-AM complete synthesis cardinality failed.", call. = FALSE)
}

comparability <- mv05am_comparability_v1()
landscape <- mv05am_landscape_contract_v1()
gaps <- mv05am_evidence_gaps_v1()
clustering <- data.frame(
  contract_id = "mv05am_clustering_identifiability_v1",
  panel_order = registry$panel_order, panel_id = registry$panel_id,
  directed_retrieval_pairs_complete = TRUE,
  complete_within_training_pair_matrix = FALSE,
  clustering_identifiable_from_accepted_rows = FALSE,
  clustering_authorized = FALSE,
  required_future_action =
    "separate_prospective_complete_pair_distance_contract_after_public_landscape_contract",
  stringsAsFactors = FALSE)
criteria <- mv05am_continuation_criteria_v1()
criteria$observed_evidence <- c(
  "4_of_4_complete_24_estimand_4_test_panels_bound",
  "MV5-V_positions_1_through_4_complete_no_position_5",
  "decision_helper_schema_rejects_estimate_interval_pvalue_tissue_winner_rank",
  "8_of_8_landscape_contract_items_preserved",
  "four_distinct_interventions_reported_separately",
  "4_of_4_panels_lack_complete_within_training_pair_matrices",
  "exact_engine_and_general_public_reference_default_state_not_yet_reconciled",
  "new_calculation_ranking_label_outcome_clustering_default_change_all_zero")
criteria$passed <- TRUE
firewall <- data.frame(
  contract_id = "mv05am_selection_firewall_v1",
  prohibited_decision_input = c("representation", "homology_dimension",
    "endpoint", "tissue", "seed", "estimate", "confidence_interval",
    "raw_p_value", "adjusted_p_value", "method_rank_or_winner"),
  present_in_public_synthesis = TRUE,
  consumed_by_continuation_helper = FALSE,
  privileged_or_filtered = FALSE, stringsAsFactors = FALSE)
decision <- mv05am_decide_v1(registry, bindings, comparability, gaps, criteria)
execution <- data.frame(
  contract_id = "mv05am_prohibited_execution_counters_v1",
  operation = c("ph", "landscape", "distance", "ranking", "label_access",
    "outcome", "clustering", "gene_view", "fusion", "new_data", "rust",
    "optimization", "public_default_change", "manuscript_claim"),
  executed = FALSE, authorized = FALSE, count = 0L,
  stringsAsFactors = FALSE)

outputs <- list(
  "mv05am-source-freeze-2026-08-12.csv" = source,
  "mv05am-panel-registry-2026-08-12.csv" = registry,
  "mv05am-complete-panel-binding-2026-08-12.csv" = bindings,
  "mv05am-complete-macro-synthesis-2026-08-12.csv" = macros,
  "mv05am-complete-interval-synthesis-2026-08-12.csv" = intervals,
  "mv05am-complete-primary-synthesis-2026-08-12.csv" = primaries,
  "mv05am-cross-panel-comparability-2026-08-12.csv" = comparability,
  "mv05am-landscape-definition-2026-08-12.csv" = landscape,
  "mv05am-clustering-identifiability-2026-08-12.csv" = clustering,
  "mv05am-evidence-gap-registry-2026-08-12.csv" = gaps,
  "mv05am-continuation-criteria-2026-08-12.csv" = criteria,
  "mv05am-selection-firewall-2026-08-12.csv" = firewall,
  "mv05am-continuation-decision-2026-08-12.csv" = decision,
  "mv05am-prohibited-execution-counters-2026-08-12.csv" = execution)
for (name in names(outputs)) write_once(outputs[[name]], name)
message("MV5-AM bound 4 panels: macro=96 intervals=96 primary=16; new execution=0")
