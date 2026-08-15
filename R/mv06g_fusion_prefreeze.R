# Internal MV6-G complete-corpus fusion prefreeze contracts.

mv06g_method_panel_v1 <- function() {
  data.frame(
    method_order = seq_len(9L),
    method_id = c(
      "cell_H0", "cell_H1", "gene_H0", "gene_H1",
      "cell_composite", "fusion_gene_weight_025",
      "fusion_gene_weight_050", "fusion_gene_weight_075",
      "gene_composite"
    ),
    method_role = c(
      rep("component_descriptive", 4L), "component_comparator",
      "weight_sensitivity", "fusion_primary", "weight_sensitivity",
      "component_comparator"
    ),
    gene_weight = c(NA, NA, NA, NA, 0, 0.25, 0.5, 0.75, 1),
    formula_id = c(
      "normalized_cell_H0", "normalized_cell_H1",
      "normalized_gene_H0", "normalized_gene_H1",
      "mean_normalized_cell_H0_H1",
      "convex_cell_gene_composites_w025",
      "convex_cell_gene_composites_w050",
      "convex_cell_gene_composites_w075",
      "mean_normalized_gene_H0_H1"
    ),
    selected_from_outcomes = FALSE,
    stringsAsFactors = FALSE
  )
}

mv06g_training_workload_v1 <- function(queue) {
  required <- c(
    "group_id", "fold_id", "held_out_study", "seed", "training_samples",
    "held_out_samples", "biological_pairs", "landscape_component_rows",
    "stage", "execution_order", "outcome_label_state",
    "biological_outcomes_computed", "fusion_jobs", "clustering_jobs",
    "outcome_jobs"
  )
  if (!is.data.frame(queue) || nrow(queue) != 75L ||
      !all(required %in% names(queue)) || anyDuplicated(queue$group_id) ||
      !identical(sort(as.integer(queue$execution_order)), 1:75) ||
      sum(queue$stage == "stage_1_maximum") != 1L ||
      sum(queue$stage == "stage_2") != 74L ||
      any(queue$outcome_label_state != "closed") ||
      any(as.logical(queue$biological_outcomes_computed)) ||
      any(unlist(queue[c(
        "fusion_jobs", "clustering_jobs", "outcome_jobs"
      )], use.names = FALSE) != 0)) {
    stop("MV6-G requires the complete label-closed MV6-F queue.",
         call. = FALSE)
  }
  queue <- queue[order(queue$execution_order), , drop = FALSE]
  rownames(queue) <- NULL
  training_pairs <- queue$training_samples * (queue$training_samples - 1) / 2
  data.frame(
    contract_id = "mv06g_training_scale_workload_v1",
    group_id = queue$group_id,
    fold_id = queue$fold_id,
    held_out_study = queue$held_out_study,
    seed = queue$seed,
    execution_order = queue$execution_order,
    training_samples = queue$training_samples,
    held_out_samples = queue$held_out_samples,
    training_biological_pairs = as.integer(training_pairs),
    training_component_rows = as.integer(4 * training_pairs),
    component_scales = 4L,
    query_biological_pairs = queue$biological_pairs,
    query_component_rows = queue$landscape_component_rows,
    ranking_methods = nrow(mv06g_method_panel_v1()),
    query_ranking_rows = queue$biological_pairs *
      nrow(mv06g_method_panel_v1()),
    outcome_label_state = "closed",
    biological_outcomes_computed = FALSE,
    fusion_jobs = 0L,
    outcome_jobs = 0L,
    stringsAsFactors = FALSE
  )
}

mv06g_endpoint_plan_v1 <- function() {
  data.frame(
    endpoint_order = 1:2,
    endpoint_id = c(
      "cross_study_tissue_mrr_v1",
      "cross_study_tissue_1nn_balanced_accuracy_v1"
    ),
    role = c("primary", "supportive"),
    aggregation = c(
      "seed_within_sample_then_sample_within_tissue_then_equal_tissue_macro",
      "seed_within_sample_then_sample_within_tissue_then_equal_tissue_macro"
    ),
    inferential_unit = "held_out_biological_sample_blocked_by_study",
    stringsAsFactors = FALSE
  )
}

mv06g_contrast_plan_v1 <- function() {
  data.frame(
    contrast_order = 1:2,
    family_id = "mv06g_F1_primary_mrr",
    endpoint_id = "cross_study_tissue_mrr_v1",
    method_id = "fusion_gene_weight_050",
    comparator_id = c("cell_composite", "gene_composite"),
    null_hypothesis = "paired_macro_difference_equals_zero",
    direction_required = "positive",
    randomization = "two_sided_paired_heldout_study_block_sign_flip",
    multiplicity = "Holm_across_two_primary_MRR_contrasts",
    fusion_benefit_requires_both = TRUE,
    stringsAsFactors = FALSE
  )
}
