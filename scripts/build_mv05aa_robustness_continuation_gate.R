#!/usr/bin/env Rscript

options(warn = 2)
for (package in c("digest", "pkgload")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(package, " is required for MV5-AA.", call. = FALSE)
  }
}
pkgload::load_all(".", quiet = TRUE, export_all = TRUE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("usage: build_mv05aa_robustness_continuation_gate.R AUDIT_DIR",
       call. = FALSE)
}
audit_dir <- args[[1L]]
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
file_sha <- function(path) digest::digest(
  file = path, algo = "sha256", serialize = FALSE)
read_public <- function(path) utils::read.csv(
  path, stringsAsFactors = FALSE, check.names = FALSE)
write_once <- function(value, name) {
  path <- file.path(audit_dir, name)
  if (file.exists(path)) stop("Refusing to overwrite: ", path, call. = FALSE)
  write_provenance_csv(value, path)
}

expected_head <- "45c7685"
observed_head <- trimws(system2(
  "git", c("rev-parse", "--short=7", "HEAD"), stdout = TRUE))
if (!identical(observed_head, expected_head)) {
  stop("MV5-AA must apply from prospective gate commit ", expected_head,
       "; observed ", observed_head, ".", call. = FALSE)
}

paths <- c(
  t_code = "R/mv05t_robustness_gate.R",
  t_configuration = "docs/audits/mv05t-configuration-registry-2026-08-10.csv",
  v_spec = "docs/specifications/MV05V_STREAMED_FULL_ROBUSTNESS_PREFREEZE_SPECIFICATION_V1.md",
  v_queue = "docs/audits/mv05v-full-group-queue-2026-08-10.csv",
  v_source_freeze = "docs/audits/mv05v-source-freeze-2026-08-10.csv",
  v_resource = "docs/audits/mv05v-resource-projection-2026-08-10.csv",
  x_resources = "docs/audits/mv05x-pc20-production-resources-2026-08-11.csv",
  z_spec = "docs/specifications/MV05Z_PC20_RETRIEVAL_ROBUSTNESS_EXECUTION_SPECIFICATION_V1.md",
  z_report = "docs/audits/MV05Z_PC20_RETRIEVAL_ROBUSTNESS_EXECUTION_2026-08-11.md",
  z_production = "docs/audits/mv05z-production-summary-2026-08-11.csv",
  z_primary = "docs/audits/mv05z-primary-contrasts-2026-08-11.csv",
  z_macro = "docs/audits/mv05z-macro-estimands-2026-08-11.csv",
  z_intervals = "docs/audits/mv05z-estimand-intervals-2026-08-11.csv",
  z_manifest = "docs/audits/mv05z-artifact-manifest-2026-08-11.csv",
  aa_code = "R/mv05aa_robustness_continuation_gate.R",
  aa_tests = "tests/testthat/test-mv05aa-robustness-continuation-gate.R",
  aa_spec = "docs/specifications/MV05AA_SELECTION_RESISTANT_ROBUSTNESS_CONTINUATION_GATE_SPECIFICATION_V1.md",
  aa_builder = "scripts/build_mv05aa_robustness_continuation_gate.R",
  aa_validator = "scripts/validate_mv05aa_robustness_continuation_gate.R")
if (any(!file.exists(paths))) {
  stop("MV5-AA source set is incomplete: ",
       paste(names(paths)[!file.exists(paths)], collapse = ", "),
       call. = FALSE)
}

source_freeze <- data.frame(
  contract_id = "mv05aa_source_freeze_v1", source_id = names(paths),
  path = unname(paths), sha256 = vapply(paths, file_sha, character(1L)),
  bytes = unname(file.info(paths)$size), accepted_head = expected_head,
  stringsAsFactors = FALSE)

v_queue <- read_public(paths[["v_queue"]])
v_sources <- read_public(paths[["v_source_freeze"]])
v_source_semantic_sha <- .mv05v_digest(stats::setNames(
  v_sources$sha256,
  paste(v_sources$source_class, v_sources$source_id, sep = "\r")))
if (nrow(v_sources) != 176L ||
    !all(grepl("^[0-9a-f]{64}$", v_sources$sha256)) ||
    length(unique(v_queue$source_freeze_sha256)) != 1L ||
    unique(v_queue$source_freeze_sha256) != v_source_semantic_sha) {
  stop("MV5-V frozen source identity is incomplete or drifted.", call. = FALSE)
}
order <- mv05aa_configuration_order_v1(v_queue)

t_configuration <- read_public(paths[["t_configuration"]])
if (!setequal(t_configuration$configuration_id, order$configuration_id)) {
  stop("MV5-T and MV5-V configuration registries disagree.", call. = FALSE)
}

production <- read_public(paths[["z_production"]])
primary <- read_public(paths[["z_primary"]])
macro <- read_public(paths[["z_macro"]])
intervals <- read_public(paths[["z_intervals"]])
evidence <- mv05aa_bind_pc20_evidence_v1(
  production, primary, macro, intervals)
evidence$production_sha256 <- file_sha(paths[["z_production"]])
evidence$primary_sha256 <- file_sha(paths[["z_primary"]])
evidence$macro_sha256 <- file_sha(paths[["z_macro"]])
evidence$interval_sha256 <- file_sha(paths[["z_intervals"]])
evidence$artifact_manifest_sha256 <- file_sha(paths[["z_manifest"]])

x_resources <- read_public(paths[["x_resources"]])
if (nrow(x_resources) != 150L || any(x_resources$disposition != "completed") ||
    any(as.integer(x_resources$exit_status) != 0L) ||
    any(.mv05aa_true(x_resources$labels_opened)) ||
    any(.mv05aa_true(x_resources$outcomes_computed))) {
  stop("MV5-X measured resource precedent is incomplete.", call. = FALSE)
}
resource <- data.frame(
  contract_id = "mv05aa_cosine_resource_envelope_v1",
  precedent_configuration = "cells384_pc20_euclidean_v1",
  precedent_groups = 150L,
  precedent_worker_seconds = sum(as.numeric(x_resources$elapsed_seconds)),
  precedent_max_group_seconds = max(as.numeric(x_resources$elapsed_seconds)),
  precedent_peak_rss_bytes = max(as.numeric(
    x_resources$peak_process_tree_rss_bytes)),
  cosine_max_workers = 1L, cosine_group_cap_seconds = 600,
  cosine_group_rss_cap_bytes = 4294967296,
  cosine_configuration_cap_worker_hours = 8,
  cosine_storage_cap_bytes = 4294967296,
  precedent_is_feasibility_evidence_not_runtime_guarantee = TRUE,
  resource_gate_pass = TRUE, stringsAsFactors = FALSE)

criteria <- mv05aa_continuation_criteria_v1()
criterion_pass <- criteria
criterion_pass$observed_evidence <- c(
  "MV5-V positions cosine second after completed PC20",
  "MV5-Z has 24/24 estimands, 24/24 intervals, and 4/4 primary tests",
  "PC count does not test radial-scale dependence",
  "cosine holds 384 cells and 30 PCs fixed",
  "MV5-V binds 150 coordinate hashes; MV5-T found zero zero-norm cells",
  "MV5-X precedent passed; MV5-V caps retained for one configuration",
  "decision function accepts no estimate, interval, p-value, or subgroup input",
  "authorization fields keep labels/outcomes/clustering/nested configs false")
criterion_pass$passed <- TRUE
decision <- mv05aa_decide_v1(order, evidence, criterion_pass)
queue <- mv05aa_cosine_queue_v1(v_queue, decision)

scope <- data.frame(
  contract_id = "mv05aa_cosine_execution_scope_v1",
  configuration_id = decision$authorized_configuration_id,
  groups = nrow(queue), folds = length(unique(queue$fold_id)),
  seeds = length(unique(queue$seed)),
  representations = length(unique(queue$representation)),
  views = sum(as.integer(queue$view_count)),
  biological_pairs = sum(as.integer(queue$biological_pairs)),
  landscape_request_rows = sum(as.integer(queue$landscape_request_rows)),
  landscape_subchunks = sum(as.integer(queue$landscape_subchunks)),
  energy_request_rows = sum(as.integer(queue$energy_request_rows)),
  assembled_method_rows = sum(as.integer(queue$assembled_method_rows)),
  labels_opened = FALSE, rankings_computed = FALSE,
  outcomes_computed = FALSE, execution_completed = FALSE,
  stringsAsFactors = FALSE)
expected_scope <- c(groups = 150L, folds = 15L, seeds = 5L,
  representations = 2L, views = 13500L, biological_pairs = 70700L,
  landscape_request_rows = 141400L, landscape_subchunks = 720L,
  energy_request_rows = 70700L, assembled_method_rows = 282800L)
if (any(vapply(names(expected_scope), function(name) {
  as.integer(scope[[name]]) != expected_scope[[name]]
}, logical(1L)))) stop("MV5-AA cosine scope count drifted.", call. = FALSE)

validation <- data.frame(
  contract_id = "mv05aa_cosine_validation_plan_v1",
  validation_order = seq_len(12L),
  validation_id = c("source_hashes", "queue_axes", "normalization",
    "point_identity", "ph_and_mst", "exact_landscapes", "matched_energy",
    "method_rows", "atomicity", "clean_repeat", "immutable_resume",
    "independent_reconstruction"),
  required = TRUE, status = "required_before_any_cosine_outcome_prefreeze",
  stringsAsFactors = FALSE)
aborts <- data.frame(
  contract_id = "mv05aa_cosine_abort_rules_v1",
  abort_order = seq_len(10L),
  abort_id = c("source_or_commit_drift", "queue_or_axis_drift",
    "zero_or_nonfinite_norm", "shape_or_point_identity_failure",
    "ph_mst_or_landscape_failure", "energy_or_method_row_failure",
    "atomic_repeat_or_resume_failure", "resource_cap_breach",
    "label_outcome_or_ranking_access", "later_configuration_or_clustering_access"),
  scope = c(rep("current_cosine_configuration", 8L), "whole_sprint",
            "whole_sprint"),
  required_action = "abort_without_substitution_or_automatic_repair",
  stringsAsFactors = FALSE)

write_once(source_freeze, "mv05aa-source-freeze-2026-08-11.csv")
write_once(order, "mv05aa-configuration-order-2026-08-11.csv")
write_once(evidence, "mv05aa-complete-pc20-evidence-binding-2026-08-11.csv")
write_once(criterion_pass, "mv05aa-continuation-criteria-2026-08-11.csv")
write_once(resource, "mv05aa-cosine-resource-envelope-2026-08-11.csv")
write_once(decision, "mv05aa-continuation-decision-2026-08-11.csv")
write_once(scope, "mv05aa-cosine-execution-scope-2026-08-11.csv")
write_once(queue, "mv05aa-cosine-execution-queue-2026-08-11.csv")
write_once(validation, "mv05aa-cosine-validation-plan-2026-08-11.csv")
write_once(aborts, "mv05aa-cosine-abort-rules-2026-08-11.csv")

message("MV5-AA authorized exactly 150 later label-closed cosine groups; ",
        "no calculation or outcome executed.")
