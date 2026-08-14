#!/usr/bin/env Rscript

options(warn = 2)
for (package in c("digest", "pkgload")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(package, " is required for MV5-AI.", call. = FALSE)
  }
}
pkgload::load_all(".", quiet = TRUE, export_all = TRUE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("usage: build_mv05ai_robustness_continuation_gate.R AUDIT_DIR EXPECTED_HEAD",
       call. = FALSE)
}
audit_dir <- args[[1L]]
expected_head <- tolower(trimws(args[[2L]]))
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
sha <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
readc <- function(path) read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
write_once <- function(value, name) {
  path <- file.path(audit_dir, name)
  if (file.exists(path)) stop("Refusing to overwrite: ", path, call. = FALSE)
  write_provenance_csv(value, path)
}

if (!grepl("^[0-9a-f]{40}$", expected_head)) {
  stop("EXPECTED_HEAD must be a full 40-character Git commit identity.",
       call. = FALSE)
}
head <- trimws(system2("git", c("rev-parse", "HEAD"), stdout = TRUE))
if (!identical(head, expected_head)) {
  stop("MV5-AI requires its prospectively committed engine HEAD ", expected_head,
       "; observed ", head, ".", call. = FALSE)
}

paths <- c(
  t_code = "R/mv05t_robustness_gate.R",
  t_configuration = "docs/audits/mv05t-configuration-registry-2026-08-10.csv",
  u_admission_resources = "docs/audits/mv05u-admission-resources-2026-08-10.csv",
  u_validation = "docs/audits/mv05u-independent-validation-2026-08-10.csv",
  v_spec = "docs/specifications/MV05V_STREAMED_FULL_ROBUSTNESS_PREFREEZE_SPECIFICATION_V1.md",
  v_queue = "docs/audits/mv05v-full-group-queue-2026-08-10.csv",
  v_source_freeze = "docs/audits/mv05v-source-freeze-2026-08-10.csv",
  aa_spec = "docs/specifications/MV05AA_SELECTION_RESISTANT_ROBUSTNESS_CONTINUATION_GATE_SPECIFICATION_V1.md",
  aa_order = "docs/audits/mv05aa-configuration-order-2026-08-11.csv",
  aa_decision = "docs/audits/mv05aa-continuation-decision-2026-08-11.csv",
  ae_order = "docs/audits/mv05ae-configuration-order-2026-08-11.csv",
  ae_decision = "docs/audits/mv05ae-continuation-decision-2026-08-11.csv",
  x_resources = "docs/audits/mv05x-pc20-production-resources-2026-08-11.csv",
  z_production = "docs/audits/mv05z-production-summary-2026-08-11.csv",
  z_primary = "docs/audits/mv05z-primary-contrasts-2026-08-11.csv",
  z_macro = "docs/audits/mv05z-macro-estimands-2026-08-11.csv",
  z_intervals = "docs/audits/mv05z-estimand-intervals-2026-08-11.csv",
  z_manifest = "docs/audits/mv05z-artifact-manifest-2026-08-11.csv",
  ab_resources = "docs/audits/mv05ab-cosine-production-resources-2026-08-11.csv",
  ad_spec = "docs/specifications/MV05AD_COSINE_RETRIEVAL_ROBUSTNESS_EXECUTION_SPECIFICATION_V1.md",
  ad_report = "docs/audits/MV05AD_COSINE_RETRIEVAL_ROBUSTNESS_EXECUTION_2026-08-11.md",
  ad_production = "docs/audits/mv05ad-production-summary-2026-08-11.csv",
  ad_primary = "docs/audits/mv05ad-primary-contrasts-2026-08-11.csv",
  ad_macro = "docs/audits/mv05ad-macro-estimands-2026-08-11.csv",
  ad_intervals = "docs/audits/mv05ad-estimand-intervals-2026-08-11.csv",
  ad_manifest = "docs/audits/mv05ad-artifact-manifest-2026-08-11.csv",
  ad_validation = "docs/audits/mv05ad-outcome-independent-validation-2026-08-11.csv",
  af_completion = "docs/audits/mv05af-nested192-unit-completion-2026-08-11.csv",
  af_resources = "docs/audits/mv05af-nested192-production-resources-2026-08-11.csv",
  af_validation = "docs/audits/mv05af-independent-validation-2026-08-11.csv",
  ah_spec = "docs/specifications/MV05AH_NESTED192_RETRIEVAL_ROBUSTNESS_EXECUTION_SPECIFICATION_V1.md",
  ah_report = "docs/audits/MV05AH_NESTED192_RETRIEVAL_ROBUSTNESS_EXECUTION_2026-08-11.md",
  ah_production = "docs/audits/mv05ah-production-summary-2026-08-11.csv",
  ah_primary = "docs/audits/mv05ah-primary-contrasts-2026-08-11.csv",
  ah_macro = "docs/audits/mv05ah-macro-estimands-2026-08-11.csv",
  ah_intervals = "docs/audits/mv05ah-estimand-intervals-2026-08-11.csv",
  ah_manifest = "docs/audits/mv05ah-outcome-clean-repeat-2026-08-11.csv",
  ah_validation = "docs/audits/mv05ah-outcome-independent-validation-2026-08-11.csv",
  ai_code = "R/mv05ai_robustness_continuation_gate.R",
  ai_tests = "tests/testthat/test-mv05ai-robustness-continuation-gate.R",
  ai_spec = "docs/specifications/MV05AI_SELECTION_RESISTANT_NESTED256_CONTINUATION_GATE_SPECIFICATION_V1.md",
  ai_builder = "scripts/build_mv05ai_robustness_continuation_gate.R",
  ai_validator = "scripts/validate_mv05ai_robustness_continuation_gate.R")
if (any(!file.exists(paths))) {
  stop("MV5-AI source set incomplete: ",
       paste(names(paths)[!file.exists(paths)], collapse = ", "), call. = FALSE)
}
source_freeze <- data.frame(
  contract_id = "mv05ai_source_freeze_v1", source_id = names(paths),
  artifact_locator = unname(paths), sha256 = vapply(paths, sha, character(1L)),
  bytes = as.numeric(file.info(paths)$size), accepted_head = expected_head,
  outcomes_read_for_decision = names(paths) %in% c(
    "z_production", "z_primary", "z_macro", "z_intervals",
    "ad_production", "ad_primary", "ad_macro", "ad_intervals",
    "ah_production", "ah_primary", "ah_macro", "ah_intervals"),
  stringsAsFactors = FALSE)

v_queue <- readc(paths[["v_queue"]]); v_sources <- readc(paths[["v_source_freeze"]])
if (nrow(v_queue) != 600L || nrow(v_sources) != 176L ||
    !all(grepl("^[0-9a-f]{64}$", v_sources$sha256))) {
  stop("MV5-V scope or source identity drifted.", call. = FALSE)
}
order <- mv05ai_configuration_order_v1(v_queue)
t_registry <- readc(paths[["t_configuration"]])
aa_order <- readc(paths[["aa_order"]]); aa_decision <- readc(paths[["aa_decision"]])
ae_order <- readc(paths[["ae_order"]]); ae_decision <- readc(paths[["ae_decision"]])
if (!setequal(t_registry$configuration_id, order$configuration_id) ||
    !identical(aa_order$configuration_id, order$configuration_id) ||
    aa_decision$authorized_configuration_id != order$configuration_id[[2L]] ||
    !identical(ae_order$configuration_id, order$configuration_id) ||
    ae_decision$authorized_configuration_id != order$configuration_id[[3L]]) {
  stop("Earlier configuration-order evidence drifted.", call. = FALSE)
}

bind <- function(prefix, analysis_id) {
  production <- readc(paths[[paste0(prefix, "_production")]])
  primary <- readc(paths[[paste0(prefix, "_primary")]])
  macro <- readc(paths[[paste0(prefix, "_macro")]])
  intervals <- readc(paths[[paste0(prefix, "_intervals")]])
  result <- mv05ai_bind_complete_evidence_v1(
    analysis_id, production, primary, macro, intervals)
  result$production_sha256 <- sha(paths[[paste0(prefix, "_production")]])
  result$primary_sha256 <- sha(paths[[paste0(prefix, "_primary")]])
  result$macro_sha256 <- sha(paths[[paste0(prefix, "_macro")]])
  result$interval_sha256 <- sha(paths[[paste0(prefix, "_intervals")]])
  result$completion_evidence_sha256 <- sha(paths[[paste0(prefix, "_manifest")]])
  result
}
evidence <- rbind(bind("z", "pc20_vs_pc30"),
                  bind("ad", "cosine_chord_vs_euclidean"),
                  bind("ah", "nested192_vs_384_cells"))

criteria <- mv05ai_continuation_criteria_v1()
criterion_pass <- criteria
criterion_pass$observed_evidence <- c(
  "MV5-V positions nested 256 fourth after PC20, cosine, and nested 192 complete",
  "all three analyses bind 24/24 estimands, intervals, and 4/4 tests without slicing",
  "nested 256 is the prespecified intermediate cell depth and was not selected by outcomes",
  "nested 256 retains 30 coordinates and Euclidean geometry while changing only 384 to 256 cells",
  "same SHA-256 ordering makes every accepted 192 set a strict subset of its closed 256 set",
  "six real nested-256 admissions and observed nested-192 full execution support bounded caps",
  "decision helper accepts no estimates, signs, intervals, p-values, tissues, or subgroups",
  "later outcome contract must reuse the complete 24-estimand registry",
  "authorization leaves calculation, ranking, labels, outcomes, and clustering false")
criterion_pass$passed <- TRUE
decision <- mv05ai_decide_v1(order, evidence, criterion_pass)
queue <- mv05ai_nested_256_queue_v1(v_queue, decision)

scope <- data.frame(
  contract_id = "mv05ai_nested_256_execution_scope_v1",
  configuration_id = decision$authorized_configuration_id,
  groups = nrow(queue), folds = length(unique(queue$fold_id)),
  seeds = length(unique(queue$seed)),
  representations = length(unique(queue$representation)),
  cells = unique(queue$cells), coordinates = unique(queue$coordinates),
  point_metric = unique(queue$point_metric),
  views = sum(as.integer(queue$view_count)),
  biological_pairs = sum(as.integer(queue$biological_pairs)),
  landscape_request_rows = sum(as.integer(queue$landscape_request_rows)),
  landscape_subchunks = sum(as.integer(queue$landscape_subchunks)),
  energy_request_rows = sum(as.integer(queue$energy_request_rows)),
  assembled_method_rows = sum(as.integer(queue$assembled_method_rows)),
  labels_opened = FALSE, rankings_computed = FALSE,
  outcomes_computed = FALSE, execution_completed = FALSE,
  stringsAsFactors = FALSE)
expected <- c(groups = 150L, folds = 15L, seeds = 5L,
  representations = 2L, cells = 256L, coordinates = 30L, views = 13500L,
  biological_pairs = 70700L, landscape_request_rows = 141400L,
  landscape_subchunks = 720L, energy_request_rows = 70700L,
  assembled_method_rows = 282800L)
if (any(vapply(names(expected), function(name) {
  as.integer(scope[[name]]) != expected[[name]]
}, logical(1L))) || scope$point_metric != "euclidean") {
  stop("MV5-AI nested-256 scope count drifted.", call. = FALSE)
}

admission <- readc(paths[["u_admission_resources"]])
admission <- admission[admission$configuration_id ==
  "nested_cells_256_pc30_euclidean_v1", , drop = FALSE]
x_resource <- readc(paths[["x_resources"]]); ab_resource <- readc(paths[["ab_resources"]])
af_resource <- readc(paths[["af_resources"]])
valid_resource <- function(x, n) nrow(x) == n &&
  all(x$disposition == "completed") && all(as.integer(x$exit_status) == 0L) &&
  !any(.mv05ai_true(x$labels_opened)) && !any(.mv05ai_true(x$outcomes_computed))
if (!valid_resource(admission, 6L) || !valid_resource(x_resource, 150L) ||
    !valid_resource(ab_resource, 150L) || !valid_resource(af_resource, 150L) ||
    !all(af_resource$configuration_id ==
         "nested_cells_192_pc30_euclidean_v1")) {
  stop("MV5-AI resource precedents are incomplete.", call. = FALSE)
}
projection_factor <- (256 / 192)^2
observed_nested192_seconds <- sum(as.numeric(af_resource$elapsed_seconds))
observed_nested192_max_seconds <- max(as.numeric(af_resource$elapsed_seconds))
observed_nested192_peak_rss <- max(as.numeric(
  af_resource$peak_process_tree_rss_bytes))
observed_nested192_private_bytes <- max(as.numeric(
  af_resource$cumulative_private_bytes))
resource <- data.frame(
  contract_id = "mv05ai_nested_256_resource_envelope_v1",
  admission_units = nrow(admission),
  admission_seconds = sum(as.numeric(admission$elapsed_seconds)),
  admission_max_group_seconds = max(as.numeric(admission$elapsed_seconds)),
  admission_peak_rss_bytes = max(as.numeric(admission$peak_process_tree_rss_bytes)),
  pc20_full_worker_hours = sum(as.numeric(x_resource$elapsed_seconds)) / 3600,
  cosine_full_worker_hours = sum(as.numeric(ab_resource$elapsed_seconds)) / 3600,
  nested192_full_groups = nrow(af_resource),
  nested192_full_worker_hours = observed_nested192_seconds / 3600,
  nested192_max_group_seconds = observed_nested192_max_seconds,
  nested192_peak_rss_bytes = observed_nested192_peak_rss,
  nested192_private_bytes = observed_nested192_private_bytes,
  conservative_projection_factor = projection_factor,
  projected_nested256_worker_hours =
    observed_nested192_seconds * projection_factor / 3600,
  projected_nested256_max_group_seconds =
    observed_nested192_max_seconds * projection_factor,
  projected_nested256_peak_rss_bytes =
    observed_nested192_peak_rss * projection_factor,
  projected_nested256_private_bytes =
    observed_nested192_private_bytes * projection_factor,
  later_max_workers = 1L, later_group_cap_seconds = 600,
  later_group_rss_cap_bytes = 4294967296,
  later_configuration_cap_worker_hours = 6,
  later_storage_cap_bytes = 4294967296,
  precedent_is_feasibility_evidence_not_runtime_guarantee = TRUE,
  resource_gate_pass = TRUE, stringsAsFactors = FALSE)

validation <- data.frame(
  contract_id = "mv05ai_nested_256_validation_plan_v1",
  validation_order = seq_len(12L),
  validation_id = c("source_hashes", "queue_axes", "nested_inclusion",
    "point_identity", "ph_and_mst", "exact_landscapes", "matched_energy",
    "method_rows", "atomicity", "clean_repeat", "immutable_resume",
    "independent_reconstruction"),
  required = TRUE, status = "required_before_any_nested_256_outcome_prefreeze",
  stringsAsFactors = FALSE)
aborts <- data.frame(
  contract_id = "mv05ai_nested_256_abort_rules_v1",
  abort_order = seq_len(10L),
  abort_id = c("source_or_commit_drift", "queue_or_axis_drift",
    "nested_cell_inclusion_failure", "shape_or_point_identity_failure",
    "ph_mst_or_landscape_failure", "energy_or_method_row_failure",
    "atomic_repeat_or_resume_failure", "resource_cap_breach",
    "label_outcome_or_ranking_access", "other_configuration_or_clustering_access"),
  required_action = "abort_without_substitution_or_automatic_repair",
  stringsAsFactors = FALSE)
firewall <- data.frame(
  contract_id = "mv05ai_selection_firewall_v1",
  prohibited_selection_input = c(
    "representation", "homology_dimension", "tissue", "endpoint", "seed",
    "estimate_sign_or_magnitude", "interval_exclusion", "p_value",
    "method_ranking"),
  consumed_by_decision_helper = FALSE,
  complete_prior_panels_bound = TRUE,
  stringsAsFactors = FALSE)

outputs <- list(
  "mv05ai-source-freeze-2026-08-11.csv" = source_freeze,
  "mv05ai-configuration-order-2026-08-11.csv" = order,
  "mv05ai-complete-evidence-binding-2026-08-11.csv" = evidence,
  "mv05ai-continuation-criteria-2026-08-11.csv" = criterion_pass,
  "mv05ai-selection-firewall-2026-08-11.csv" = firewall,
  "mv05ai-resource-envelope-2026-08-11.csv" = resource,
  "mv05ai-continuation-decision-2026-08-11.csv" = decision,
  "mv05ai-execution-scope-2026-08-11.csv" = scope,
  "mv05ai-execution-queue-2026-08-11.csv" = queue,
  "mv05ai-validation-plan-2026-08-11.csv" = validation,
  "mv05ai-abort-rules-2026-08-11.csv" = aborts)
for (name in names(outputs)) write_once(outputs[[name]], name)
message("MV5-AI authorized exactly 150 later label-closed nested-256 groups; calculation=0 labels=0 outcomes=0")
