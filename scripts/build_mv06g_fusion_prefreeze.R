#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop("usage: build_mv06g_fusion_prefreeze.R QUEUE COMPLETE_INVENTORY ",
       "COMPLETE_RESUME GROUP_ROOT RUST_LIBRARY OUTPUT_DIR", call. = FALSE)
}
source("R/mv06f_production.R")
source("R/mv06g_fusion_prefreeze.R")
queue <- utils::read.csv(args[[1L]], stringsAsFactors = FALSE,
                         check.names = FALSE)
inventory <- utils::read.csv(args[[2L]], stringsAsFactors = FALSE,
                             check.names = FALSE)
resume <- utils::read.csv(args[[3L]], stringsAsFactors = FALSE,
                          check.names = FALSE)
group_root <- normalizePath(args[[4L]], winslash = "/", mustWork = TRUE)
rust <- normalizePath(args[[5L]], winslash = "/", mustWork = TRUE)
output <- args[[6L]]
dir.create(output, recursive = TRUE, showWarnings = FALSE)
mv06f_validate_queue_v1(queue)
queue <- queue[order(queue$execution_order), , drop = FALSE]
rownames(queue) <- NULL
if (nrow(inventory) != 75L || nrow(resume) != 376L ||
    any(!as.logical(resume$sha256_unchanged)) ||
    any(!as.logical(resume$bytes_unchanged)) ||
    any(!as.logical(resume$mtime_unchanged))) {
  stop("MV6-G requires accepted complete and immutable MV6-F production.",
       call. = FALSE)
}
safe_name <- function(value) gsub("[^A-Za-z0-9_.-]", "_", value)
source_rows <- lapply(seq_len(nrow(queue)), function(index) {
  unit <- queue[index, , drop = FALSE]
  accepted <- inventory[inventory$group_id == unit$group_id, , drop = FALSE]
  directory <- file.path(group_root, safe_name(unit$group_id))
  diagrams <- file.path(directory, "diagrams.rds")
  distances <- file.path(directory, "distances.csv")
  if (nrow(accepted) != 1L || !file.exists(diagrams) ||
      !file.exists(distances) ||
      .mv06f_sha256(diagrams) != accepted$diagrams_sha256 ||
      .mv06f_sha256(distances) != accepted$distances_sha256) {
    stop("An MV6-G source group does not match accepted MV6-F evidence.",
         call. = FALSE)
  }
  data.frame(
    contract_id = "mv06g_source_group_inventory_v1",
    group_id = unit$group_id, fold_id = unit$fold_id, seed = unit$seed,
    execution_order = unit$execution_order,
    diagrams_bytes = file.info(diagrams)$size,
    diagrams_sha256 = .mv06f_sha256(diagrams),
    distances_bytes = file.info(distances)$size,
    distances_sha256 = .mv06f_sha256(distances),
    outcome_label_state = "closed",
    biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
})
source_groups <- do.call(rbind, source_rows)
workload <- mv06g_training_workload_v1(queue)
method_panel <- mv06g_method_panel_v1()
endpoint_plan <- mv06g_endpoint_plan_v1()
contrast_plan <- mv06g_contrast_plan_v1()
inference_plan <- data.frame(
  contract_id = "mv06g_blocked_inference_plan_v1",
  bootstrap_replicates = 2000L, bootstrap_seed = 20260815L,
  bootstrap_unit = "tissue_stratified_heldout_study_block",
  randomization_replicates = 9999L, randomization_seed = 20260816L,
  randomization_unit = "paired_heldout_study_block_sign_flip",
  interval = "percentile_2.5_97.5_type7",
  primary_adjustment = "Holm_across_two_MRR_contrasts",
  minimum_independent_study_blocks = 4L,
  stringsAsFactors = FALSE
)
label_firewall <- data.frame(
  gate_order = 1:10,
  gate = c(
    "complete_mv06f_inventory_hash_locked",
    "complete_mv06f_resume_hash_locked",
    "training_pair_scales_complete_before_labels",
    "all_nine_rankings_complete_before_labels",
    "prediction_manifest_committed_before_label_read",
    "authoritative_metadata_hash_verified_after_lock",
    "no_post_label_rescaling_or_reranking",
    "no_weight_selection_from_outcomes",
    "no_tissue_or_approach_in_prediction_artifacts",
    "independent_evaluation_and_byte_repeat_required"
  ),
  required = TRUE,
  current_state = c(rep("passed", 2L), rep("closed_pending", 8L)),
  stringsAsFactors = FALSE
)
resource_plan <- data.frame(
  contract_id = "mv06g_resource_plan_v1",
  stage = c("stage_1_maximum_group", "stage_2_complete_if_admitted"),
  groups = c(1L, 74L),
  maximum_workers = 1L,
  elapsed_cap_seconds_per_group = 1800,
  peak_process_tree_rss_cap_bytes = 12 * 1024^3,
  private_storage_cap_bytes = 5 * 1024^3,
  complete_worker_cap_seconds = 12 * 3600,
  automatic_retry = FALSE,
  authorized = c(TRUE, FALSE),
  stringsAsFactors = FALSE
)
implementation_paths <- c(
  "R/mv06g_fusion_prefreeze.R",
  "scripts/build_mv06g_fusion_prefreeze.R",
  "scripts/validate_mv06g_fusion_prefreeze.R",
  "scripts/validate_mv06g_fusion_prefreeze_repeat.R",
  "tests/testthat/test-mv06g-fusion-prefreeze.R"
)
implementation_sources <- data.frame(
  path = implementation_paths,
  sha256 = unname(vapply(implementation_paths, .mv06f_sha256, character(1L))),
  stringsAsFactors = FALSE
)
implementation_root <- digest::digest(
  stats::setNames(implementation_sources$sha256, implementation_sources$path),
  algo = "sha256", serialize = TRUE
)
stage1 <- workload[queue$stage == "stage_1_maximum", , drop = FALSE]
contract <- data.frame(
  contract_id = "mv06g_complete_fusion_prefreeze_v1",
  base_revision = "bba0b11",
  queue_root_sha256 = mv06f_queue_root_v1(queue),
  complete_inventory_sha256 = .mv06f_sha256(args[[2L]]),
  complete_resume_sha256 = .mv06f_sha256(args[[3L]]),
  rust_library_sha256 = .mv06f_sha256(rust),
  implementation_root_sha256 = implementation_root,
  groups = 75L, training_biological_pairs = sum(workload$training_biological_pairs),
  training_component_rows = sum(workload$training_component_rows),
  query_biological_pairs = sum(workload$query_biological_pairs),
  query_ranking_rows = sum(workload$query_ranking_rows),
  component_scales = 300L, ranking_methods = nrow(method_panel),
  stage1_group_id = stage1$group_id,
  primary_gene_weight = 0.5,
  metadata_sha256 = "e2b6f1e7d7dc05ba8e80627aad4ea1677c3b718a844397acd06596d7e22a18d0",
  global_panel_scope = "transductive_technical_presence_variance_only",
  stage1_authorized = TRUE, stage2_authorized = FALSE,
  outcome_labels_opened = FALSE, biological_outcomes_computed = FALSE,
  fusion_evaluations = 0L, outcome_jobs = 0L,
  disposition = "prefreeze_pass_stage1_training_scale_sentinel_only",
  stringsAsFactors = FALSE
)
write_stable <- function(value, name) utils::write.csv(
  value, file.path(output, name), row.names = FALSE, na = ""
)
write_stable(contract, "mv06g-contract.csv")
write_stable(source_groups, "mv06g-source-groups.csv")
write_stable(workload, "mv06g-training-workload.csv")
write_stable(method_panel, "mv06g-method-panel.csv")
write_stable(endpoint_plan, "mv06g-endpoint-plan.csv")
write_stable(contrast_plan, "mv06g-contrast-plan.csv")
write_stable(inference_plan, "mv06g-inference-plan.csv")
write_stable(label_firewall, "mv06g-label-firewall.csv")
write_stable(resource_plan, "mv06g-resource-plan.csv")
write_stable(implementation_sources, "mv06g-implementation-sources.csv")
message("Prefroze MV6-G fusion with stage-one training-scale sentinel only.")
