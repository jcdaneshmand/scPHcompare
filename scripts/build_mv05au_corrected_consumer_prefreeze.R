#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: build_mv05au_corrected_consumer_prefreeze.R OUTPUT_DIR")
out <- normalizePath(args[[1]], mustWork = FALSE); dir.create(out, recursive = TRUE, showWarnings = FALSE)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Unable to resolve the repository root.")
setwd(normalizePath(file.path(
  dirname(gsub("~+~", " ", sub("^--file=", "", script_arg[[1L]]), fixed = TRUE)), ".."
), mustWork = TRUE))
write_csv <- function(x, id) utils::write.csv(x, file.path(out, paste0("mv05au-", id,
  "-2026-08-12.csv")), row.names = FALSE, na = "", quote = TRUE)

consumers <- data.frame(
  contract_id = "mv05au_consumer_inventory_v1",
  consumer = c("verified_corrected_bundle_loader", "average_linkage_dendrogram",
    "pam_partition", "spectral_partition", "ward_d2_dendrogram",
    "descriptive_combined_tree", "legacy_landscape_clustering",
    "betti_euler_curves", "cross_iteration_landscape_curves",
    "cluster_outcome_comparison"),
  input_schema = c(rep("versioned_distance_matrix", 7), "persistence_diagram_or_curve",
    "sampled_landscape_curve", "cluster_assignments_plus_labels"),
  label_free = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE),
  additional_choice = c("none", "none", "k", "k_and_kernel_scale",
    "linkage_justification", "dimension_fusion_policy", "legacy_schema",
    "curve_grid_and_bandwidth", "curve_grid_and_alignment", "labels_and_outcomes"),
  disposition = c("implement_first", "implement_first", rep("closed", 8)),
  reason = c("required verified boundary", "deterministic_no_k_full_hierarchy",
    "requires_k", "requires_k_and_kernel", "linkage_family_not_compared",
    "combined_is_descriptive_only", "must_not_redirect_legacy",
    "matrix_schema_incompatible", "matrix_schema_incompatible",
    "outcome_stage_not_authorized"), stringsAsFactors = FALSE)
write_csv(consumers, "consumer-inventory")

api <- data.frame(contract_id = "mv05au_api_contract_v1",
  api = c("read_corrected_landscape_bundle", "corrected_landscape_average_trees"),
  default_off = TRUE, input = c("verified_sidecar_plus_explicit_iteration_and_view",
    "scph_corrected_landscape_consumer_bundle_v1"),
  output = c("scph_corrected_landscape_consumer_bundle_v1",
    "scph_corrected_landscape_average_trees_v1"),
  primary_dimensions = c("H0;H1", "H0;H1"), combined_status = "descriptive_not_consumed",
  mutates_source = FALSE, writes_legacy = FALSE, stringsAsFactors = FALSE)
write_csv(api, "api-contract")

views <- data.frame(contract_id = "mv05au_view_dimension_policy_v1",
  view_id = rep(c("cell_topology_v1", "gene_topology_v1"), each = 3),
  dimension = rep(c("H0", "H1", "combined"), 2),
  role = rep(c("primary_separate", "primary_separate", "descriptive_only"), 2),
  may_build_tree = rep(c(TRUE, TRUE, FALSE), 2), may_fuse_across_views = FALSE)
write_csv(views, "view-dimension-policy")

stages <- data.frame(contract_id = "mv05au_migration_stages_v1", stage = 0:5,
  action = c("current_artifacts_only", "verified_loader", "separate_average_trees",
    "bounded_label_free_smoke", "partition_prefreeze_if_justified",
    "evaluation_prefreeze_if_justified"),
  authorized_now = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE),
  requires_separate_gate = c(FALSE, TRUE, TRUE, TRUE, TRUE, TRUE))
write_csv(stages, "migration-stages")

validation <- data.frame(contract_id = "mv05au_validation_plan_v1",
  validation_id = sprintf("AU-VAL-%02d", 1:15),
  requirement = c("completion_and_shard_hashes", "matrix_class_and_contract",
    "sidecar_input_key", "sidecar_matrix_key", "explicit_iteration_id",
    "explicit_allowed_view_id", "identical_named_axes", "finite_symmetric_zero_diagonal",
    "H0_H1_separate", "combined_not_consumed", "average_linkage_only",
    "permutation_equivariance", "deterministic_serialization", "source_immutability",
    "default_and_legacy_isolation"), required = TRUE)
write_csv(validation, "validation-plan")

abort <- data.frame(contract_id = "mv05au_abort_rules_v1",
  rule_id = sprintf("AU-ABORT-%02d", 1:14), condition = c("contract_drift",
    "incomplete_or_corrupt_artifact", "sidecar_key_mismatch", "implicit_view",
    "invalid_matrix", "H0_H1_fusion", "combined_primary_use", "implicit_k",
    "label_or_outcome_access", "legacy_redirection", "default_change",
    "source_mutation", "unprofiled_resource", "test_or_check_failure"),
  action = "stop_without_implementation_authorization")
write_csv(abort, "abort-rules")

source_paths <- c("R/corrected_landscape_workflow.R", "R/landscape_public_api.R",
  "R/topological_distance_contract.R", "R/PH_PostProcessing_andAnalysis.R",
  "R/betti_utils.R", "R/cross_iteration_functions.R", "R/mv05n_clustering_gate.R",
  "R/mv05q_clustering_artifacts.R",
  "docs/audits/MV05AT_BROADER_REALISTIC_WORKFLOW_SMOKE_2026-08-12.md",
  "docs/specifications/MV05AU_CORRECTED_MATRIX_CONSUMER_PREFREEZE_SPECIFICATION_V1.md",
  "scripts/build_mv05au_corrected_consumer_prefreeze.R",
  "scripts/validate_mv05au_corrected_consumer_prefreeze.R")
freeze <- data.frame(contract_id = "mv05au_source_freeze_v1", path = source_paths,
  sha256 = vapply(source_paths, function(p) sub(" .*", "", system2("sha256sum", p,
    stdout = TRUE)), character(1)), stringsAsFactors = FALSE)
write_csv(freeze, "source-freeze")

prohibited <- data.frame(contract_id = "mv05au_prohibited_changes_v1",
  counter = c("consumer_code", "consumer_execution", "partitions", "selected_k",
    "combined_primary", "cell_gene_fusion", "labels", "outcomes", "visualizations",
    "legacy_changes", "default_changes", "new_data", "scientific_calculation",
    "optimization", "claims"), value = 0L)
write_csv(prohibited, "prohibited-change-counters")
decision <- data.frame(contract_id = "mv05au_decision_v1", prefreeze_accepted = TRUE,
  decision = "authorize_verified_loader_and_separate_average_trees_implementation_only",
  partitions_authorized = FALSE, evaluation_authorized = FALSE,
  default_change_authorized = FALSE, legacy_change_authorized = FALSE,
  next_sprint = "MV5-AV", stringsAsFactors = FALSE)
write_csv(decision, "continuation-decision")
cat("MV5-AU prefreeze built\n")
