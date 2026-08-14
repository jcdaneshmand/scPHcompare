#!/usr/bin/env Rscript

options(warn = 2)
for (package in c("digest", "pkgload")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(package, " is required for MV5-T gate construction.", call. = FALSE)
  }
}
pkgload::load_all(".", quiet = TRUE, export_all = TRUE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(paste("usage: build_mv05t_robustness_gate.R D1_PRIVATE_ROOT",
             "G_PRIVATE_ROOT AUDIT_DIR EXPECTED_HEAD"), call. = FALSE)
}
d1_root <- normalizePath(args[[1L]], mustWork = TRUE)
g_root <- normalizePath(args[[2L]], mustWork = TRUE)
audit_dir <- args[[3L]]
expected_head <- args[[4L]]
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
head <- trimws(system2("git", c("rev-parse", "HEAD"), stdout = TRUE))
if (!identical(head, expected_head)) {
  stop("MV5-T gate must run at its exact prospective HEAD.", call. = FALSE)
}
file_sha <- function(path) digest::digest(file = path, algo = "sha256",
                                          serialize = FALSE)
sources <- c(
  mv05t_spec = "docs/specifications/MV05T_SELECTION_RESISTANT_ROBUSTNESS_GAP_GATE_SPECIFICATION_V1.md",
  mv05t_code = "R/mv05t_robustness_gate.R",
  mv05t_builder = "scripts/build_mv05t_robustness_gate.R",
  mv05t_tests = "tests/testthat/test-mv05t-robustness-gate.R",
  mv05m_spec = "docs/specifications/MV05M_POST_RETRIEVAL_BENCHMARK_GAP_GATE_SPECIFICATION_V1.md",
  mv05m_audit = "docs/audits/MV05M_POST_RETRIEVAL_BENCHMARK_GAP_GATE_2026-08-10.md",
  mv05s_spec = "docs/specifications/MV05S_PREDICTION_LOCKED_CLUSTERING_OUTCOME_EXECUTION_SPECIFICATION_V1.md",
  mv05s_audit = "docs/audits/MV05S_PREDICTION_LOCKED_CLUSTERING_OUTCOME_EXECUTION_2026-08-10.md",
  d1_manifest = "docs/audits/mv05d1-artifact-manifest-2026-08-07.csv",
  g_manifest = "docs/audits/mv05g-artifact-manifest-2026-08-09.csv",
  d3_resource = "docs/audits/MV05D3_FULL_CELL_PH_PRODUCTION_2026-08-07.md",
  d4_resource = "docs/audits/MV05D4_CELL_LANDSCAPE_DISTANCES_2026-08-07.md",
  h_resource = "docs/audits/MV05H_INTEGRATED_CELL_PH_PRODUCTION_2026-08-09.md",
  i_resource = "docs/audits/MV05I_INTEGRATED_CELL_LANDSCAPE_DISTANCES_2026-08-09.md")
if (any(!file.exists(sources))) stop("MV5-T committed source is missing.",
                                     call. = FALSE)
source_freeze <- data.frame(
  contract_id = "mv05t_source_freeze_v1", source_id = names(sources),
  artifact_locator = unname(sources),
  sha256 = vapply(sources, file_sha, character(1L)),
  bytes = as.numeric(file.info(sources)$size), source_type = "committed_public",
  prospective_head = head, labels_opened = FALSE, outcomes_computed = FALSE,
  admission_executed = FALSE, stringsAsFactors = FALSE)

d1_files <- sort(list.files(d1_root, pattern = "__sct_cell_fold[.]rds$",
                            full.names = TRUE), method = "radix")
g_files <- sort(list.files(file.path(g_root, "results"),
                           pattern = "^mv05g_group__.*[.]rds$", recursive = TRUE,
                           full.names = TRUE), method = "radix")
if (length(d1_files) != 75L || length(g_files) != 75L) {
  stop("MV5-T private coordinate cache file counts drifted.", call. = FALSE)
}
parse_identity <- function(path, source_type) {
  name <- basename(path)
  pattern <- if (source_type == "sct")
    "^(SRA[0-9]+)__([0-9]+)__sct_cell_fold[.]rds$" else
    "^mv05g_group__(SRA[0-9]+)__([0-9]+)[.]rds$"
  match <- regexec(pattern, name)
  values <- regmatches(name, match)[[1L]]
  if (length(values) != 3L) stop("MV5-T private filename drifted.",
                                 call. = FALSE)
  c(fold_study = values[[2L]], seed = values[[3L]])
}
private_rows <- lapply(c(d1_files, g_files), function(path) {
  source_type <- if (path %in% d1_files) "sct" else "integrated"
  identity <- parse_identity(path, source_type)
  data.frame(
    contract_id = "mv05t_private_coordinate_inventory_v1",
    source_type = source_type, fold_study = identity[["fold_study"]],
    seed = as.integer(identity[["seed"]]),
    private_locator = paste0("private_cache:", source_type, "/", basename(path)),
    sha256 = file_sha(path), bytes = unname(file.info(path)$size),
    labels_opened = FALSE, outcomes_computed = FALSE,
    admission_executed = FALSE, stringsAsFactors = FALSE)
})
private_inventory <- do.call(rbind, private_rows)
if (anyDuplicated(private_inventory[c("source_type", "fold_study", "seed")]) ||
    !identical(sort(unique(private_inventory$seed)), 20260805:20260809)) {
  stop("MV5-T private coordinate identity axes drifted.", call. = FALSE)
}
private_freeze <- data.frame(
  contract_id = "mv05t_source_freeze_v1",
  source_id = paste0("private_", private_inventory$source_type, "_",
                     private_inventory$fold_study, "_", private_inventory$seed),
  artifact_locator = private_inventory$private_locator,
  sha256 = private_inventory$sha256, bytes = private_inventory$bytes,
  source_type = "private_hash_only", prospective_head = head,
  labels_opened = FALSE, outcomes_computed = FALSE,
  admission_executed = FALSE, stringsAsFactors = FALSE)
source_freeze <- rbind(source_freeze, private_freeze)
source_freeze$source_freeze_sha256 <- .mv05t_digest(paste(
  source_freeze$artifact_locator, source_freeze$sha256, sep = "\r"))

selected <- c("SRA779509", "SRA703206", "SRA713577")
compatibility <- lapply(selected, function(study) {
  d1_path <- file.path(d1_root, paste0(study, "__20260805__sct_cell_fold.rds"))
  g_path <- file.path(g_root, "results", paste0("mv05g_group__", study,
    "__20260805"), paste0("mv05g_group__", study, "__20260805.rds"))
  d1 <- readRDS(d1_path)
  integrated <- readRDS(g_path)
  sct_views <- d1$payload$cell_views
  integrated_views <- integrated$payload$coordinates
  if (!identical(names(sct_views), names(integrated_views)) ||
      length(sct_views) != 90L) {
    stop("MV5-T selected paired sample axes drifted.", call. = FALSE)
  }
  zero_norm <- 0L
  nested_checks <- 0L
  for (sample_id in names(sct_views)) {
    mv05t_validate_coordinate_pair_v1(sct_views[[sample_id]],
                                       integrated_views[[sample_id]])
    ids192 <- mv05t_nested_point_ids_v1(sample_id, 20260805L,
                                        sct_views[[sample_id]]$point_ids, 192L)
    ids256 <- mv05t_nested_point_ids_v1(sample_id, 20260805L,
                                        sct_views[[sample_id]]$point_ids, 256L)
    if (!all(ids192 %in% ids256)) stop("MV5-T nested point identity failed.",
                                       call. = FALSE)
    nested_checks <- nested_checks + 1L
    zero_norm <- zero_norm + sum(rowSums(sct_views[[sample_id]]$payload^2) <= 0) +
      sum(rowSums(integrated_views[[sample_id]]^2) <= 0)
  }
  data.frame(
    contract_id = "mv05t_coordinate_compatibility_v1", fold_study = study,
    seed = 20260805L, paired_samples = length(sct_views),
    exact_384_by_30_pairs = length(sct_views),
    exact_cell_identity_pairs = length(sct_views),
    nested_192_in_256_checks = nested_checks,
    zero_norm_cells = zero_norm, first20_available = TRUE,
    labels_opened = FALSE, outcomes_computed = FALSE,
    admission_executed = FALSE, stringsAsFactors = FALSE)
})
compatibility <- do.call(rbind, compatibility)
if (any(compatibility$zero_norm_cells != 0L) ||
    any(compatibility$exact_cell_identity_pairs != 90L)) {
  stop("MV5-T coordinate compatibility admission failed.", call. = FALSE)
}

criteria <- mv05t_criteria_v1()
candidates <- mv05t_candidate_registry_v1()
configurations <- mv05t_configuration_registry_v1()
queue <- mv05t_build_admission_queue_v1(
  unique(source_freeze$source_freeze_sha256))
resource <- data.frame(
  contract_id = "mv05t_resource_projection_v1",
  baseline_ph_worker_seconds_per_full_setting = 7720.432,
  baseline_landscape_worker_seconds_per_full_setting = 6267.045,
  baseline_total_worker_seconds_per_full_setting = 13987.477,
  new_full_settings = 4L,
  projected_full_worker_hours = 4 * 13987.477 / 3600,
  conservative_full_storage_bytes = 10180000000,
  full_execution_authorized = FALSE,
  admission_groups = 24L, admission_worker_limit = 1L,
  admission_worker_hour_cap = 2,
  admission_per_group_elapsed_cap_seconds = 600,
  admission_process_rss_cap_bytes = 4294967296,
  admission_storage_cap_bytes = 2147483648,
  labels_opened = FALSE, outcomes_computed = FALSE,
  admission_executed = FALSE, stringsAsFactors = FALSE)
validation <- data.frame(
  contract_id = "mv05t_validation_plan_v1",
  validation_id = c(
    "committed_source_hashes", "private_cache_150_hashes",
    "selected_pair_shapes_and_ids", "nested_cell_identity",
    "first20_coordinate_identity", "cosine_chord_zero_norm",
    "independent_h0_mst", "small_exact_h1_landscape_oracle",
    "matched_energy_oracle", "deterministic_repeat", "immutable_resume",
    "resource_and_public_safety"),
  required = TRUE, execution_state = "prospective_not_run",
  labels_opened = FALSE, outcomes_computed = FALSE,
  admission_executed = FALSE, stringsAsFactors = FALSE)
abort <- data.frame(
  contract_id = "mv05t_abort_rules_v1",
  rule_id = sprintf("MV5T-ABORT-%02d", 1:10),
  trigger = c(
    "committed_or_private_source_hash_drift",
    "coordinate_shape_cell_or_sample_axis_mismatch",
    "nested_subset_or_representation_pair_mismatch",
    "missing_first20_coordinate_or_nonfinite_value",
    "zero_norm_cosine_chord_cell", "label_open_or_outcome_attempt",
    "ph_landscape_or_energy_oracle_failure",
    "partial_stale_or_hash_invalid_artifact_status_pair",
    "time_memory_or_storage_cap_breach",
    "repeat_resume_or_public_safety_failure"),
  disposition = "abort_preserve_completed_immutable_artifacts_review_before_retry",
  automatic_retry = FALSE, labels_opened = FALSE, outcomes_computed = FALSE,
  admission_executed = FALSE, stringsAsFactors = FALSE)

outputs <- list(
  "mv05t-source-freeze-2026-08-10.csv" = source_freeze,
  "mv05t-private-coordinate-inventory-2026-08-10.csv" = private_inventory,
  "mv05t-coordinate-compatibility-2026-08-10.csv" = compatibility,
  "mv05t-selection-criteria-2026-08-10.csv" = criteria,
  "mv05t-candidate-registry-2026-08-10.csv" = candidates,
  "mv05t-configuration-registry-2026-08-10.csv" = configurations,
  "mv05t-admission-queue-2026-08-10.csv" = queue,
  "mv05t-resource-projection-2026-08-10.csv" = resource,
  "mv05t-validation-plan-2026-08-10.csv" = validation,
  "mv05t-abort-rules-2026-08-10.csv" = abort)
for (name in names(outputs)) {
  write_provenance_csv(outputs[[name]], file.path(audit_dir, name))
}
message("MV5-T gate passed: committed_sources=", length(sources),
        " private_files=", nrow(private_inventory),
        " paired_selected_views=", sum(compatibility$paired_samples),
        " candidate_families=", nrow(candidates),
        " admitted_configurations=", nrow(configurations),
        " admission_units=", nrow(queue),
        " outcomes=0 admission=0")
