.mv06b_assert_record_axis <- function(identity) {
  ids <- sort(c(identity$training_ids, identity$query_ids), method = "radix")
  if (length(ids) != 90L || anyDuplicated(ids) ||
      length(identity$training_ids) + length(identity$query_ids) != 90L) {
    stop("MV6-B requires an exact 90-view training/query axis.", call. = FALSE)
  }
  invisible(ids)
}

#' Summarize one accepted SCT fold cache for the MV6-B gate
#'
#' @keywords internal
mv06b_summarize_sct_group_v1 <- function(record) {
  mv05d1_validate_cell_fold_record_v1(record)
  identity <- record$identity
  .mv06b_assert_record_axis(identity)
  missing <- record$payload$missing_feature_counts
  query_missing <- unname(missing[identity$query_ids])
  if (anyNA(query_missing) || any(query_missing < 0L) ||
      any(missing[identity$training_ids] != 0L)) {
    stop("MV6-B found invalid SCT missing-feature provenance.", call. = FALSE)
  }
  query_complete <- sum(query_missing == 0L)
  query_incomplete <- sum(query_missing > 0L)
  data.frame(
    fold_id = identity$fold_id,
    held_out_study = identity$held_out_study,
    seed = identity$seed,
    training_views = length(identity$training_ids),
    query_views = length(identity$query_ids),
    exact_panel_variance_resolved_views = length(identity$training_ids),
    exact_panel_variance_unresolved_views = query_complete,
    incomplete_panel_views = query_incomplete,
    missing_feature_instances = sum(query_missing),
    maximum_missing_features = if (length(query_missing)) max(query_missing) else 0L,
    full_panel_group = query_incomplete == 0L,
    expression_payload_status = "reconstructable_from_accepted_d0_d1_sources",
    stringsAsFactors = FALSE
  )
}

#' Summarize one accepted integrated-coordinate cache for the MV6-B gate
#'
#' @keywords internal
mv06b_summarize_integrated_group_v1 <- function(record) {
  mv05f_validate_group_record_v1(record)
  identity <- record$identity
  .mv06b_assert_record_axis(identity)
  active <- record$payload$active_features
  if (!identical(sort(names(active), method = "radix"),
                 sort(identity$query_ids, method = "radix"))) {
    stop("MV6-B found an invalid integrated active-feature axis.", call. = FALSE)
  }
  active_count <- vapply(active, length, integer(1L))
  if (any(active_count < 1L) || any(active_count > 500L)) {
    stop("MV6-B found invalid integrated active-feature counts.", call. = FALSE)
  }
  dropped <- 500L - active_count
  data.frame(
    fold_id = identity$fold_id,
    held_out_study = identity$held_out_study,
    seed = identity$seed,
    training_views = length(identity$training_ids),
    query_views = length(identity$query_ids),
    full_panel_views = length(identity$training_ids) + sum(dropped == 0L),
    incomplete_panel_views = sum(dropped > 0L),
    missing_feature_instances = sum(dropped),
    maximum_missing_features = max(dropped),
    full_panel_group = all(dropped == 0L),
    expression_undefined_views = 90L,
    expression_payload_status =
      "undefined_from_accepted_integrated_cell_coordinates",
    stringsAsFactors = FALSE
  )
}

#' Finalize the label-closed MV6-B structural inventory
#'
#' @keywords internal
mv06b_finalize_inventory_v1 <- function(sct_groups, integrated_groups) {
  required_keys <- c("fold_id", "held_out_study", "seed")
  if (!is.data.frame(sct_groups) || !is.data.frame(integrated_groups) ||
      !all(required_keys %in% names(sct_groups)) ||
      !all(required_keys %in% names(integrated_groups)) ||
      nrow(sct_groups) != 75L || nrow(integrated_groups) != 75L) {
    stop("MV6-B requires 75 SCT and 75 integrated group summaries.",
         call. = FALSE)
  }
  key <- function(value) paste(value$fold_id, value$seed, sep = "\r")
  if (anyDuplicated(key(sct_groups)) || anyDuplicated(key(integrated_groups)) ||
      !setequal(key(sct_groups), key(integrated_groups))) {
    stop("MV6-B SCT and integrated group axes differ.", call. = FALSE)
  }
  sct_groups <- sct_groups[order(
    sct_groups$held_out_study, sct_groups$seed, method = "radix"
  ), , drop = FALSE]
  integrated_groups <- integrated_groups[match(
    key(sct_groups), key(integrated_groups)
  ), , drop = FALSE]

  group <- data.frame(
    contract_id = "mv06b_group_structural_inventory_v1",
    fold_id = sct_groups$fold_id,
    held_out_study = sct_groups$held_out_study,
    seed = sct_groups$seed,
    training_views = sct_groups$training_views,
    query_views = sct_groups$query_views,
    sct_incomplete_panel_views = sct_groups$incomplete_panel_views,
    sct_missing_feature_instances = sct_groups$missing_feature_instances,
    sct_maximum_missing_features = sct_groups$maximum_missing_features,
    sct_full_panel_group = sct_groups$full_panel_group,
    integrated_incomplete_panel_views = integrated_groups$incomplete_panel_views,
    integrated_missing_feature_instances = integrated_groups$missing_feature_instances,
    integrated_maximum_missing_features = integrated_groups$maximum_missing_features,
    integrated_full_panel_group = integrated_groups$full_panel_group,
    integrated_expression_payload_defined = FALSE,
    stringsAsFactors = FALSE
  )

  summary <- rbind(
    data.frame(
      contract_id = "mv06b_representation_structural_summary_v1",
      representation = "sct_fold",
      groups = 75L, view_instances = 6750L,
      training_views = sum(sct_groups$training_views),
      query_views = sum(sct_groups$query_views),
      exact_panel_views = sum(sct_groups$exact_panel_variance_resolved_views) +
        sum(sct_groups$exact_panel_variance_unresolved_views),
      incomplete_panel_views = sum(sct_groups$incomplete_panel_views),
      missing_feature_instances = sum(sct_groups$missing_feature_instances),
      maximum_missing_features = max(sct_groups$maximum_missing_features),
      affected_groups = sum(!sct_groups$full_panel_group),
      variance_resolved_views =
        sum(sct_groups$exact_panel_variance_resolved_views),
      variance_unresolved_views =
        sum(sct_groups$exact_panel_variance_unresolved_views),
      expression_undefined_views = 0L,
      strict_axis_complete = FALSE,
      disposition = "incomplete_query_panels_and_unresolved_query_variance",
      stringsAsFactors = FALSE
    ),
    data.frame(
      contract_id = "mv06b_representation_structural_summary_v1",
      representation = "inductive_integration",
      groups = 75L, view_instances = 6750L,
      training_views = sum(integrated_groups$training_views),
      query_views = sum(integrated_groups$query_views),
      exact_panel_views = sum(integrated_groups$full_panel_views),
      incomplete_panel_views = sum(integrated_groups$incomplete_panel_views),
      missing_feature_instances = sum(integrated_groups$missing_feature_instances),
      maximum_missing_features = max(integrated_groups$maximum_missing_features),
      affected_groups = sum(!integrated_groups$full_panel_group),
      variance_resolved_views = 0L,
      variance_unresolved_views = 0L,
      expression_undefined_views = 6750L,
      strict_axis_complete = FALSE,
      disposition = "accepted_artifacts_define_cell_coordinates_not_gene_expression",
      stringsAsFactors = FALSE
    )
  )

  workload <- data.frame(
    contract_id = "mv06b_lower_bound_workload_inventory_v1",
    representation = c("sct_fold", "inductive_integration"),
    view_instances = c(6750L, 6750L),
    h0_h1_diagram_components = c(13500L, 13500L),
    directed_query_training_pairs = c(35350L, 35350L),
    dimension_pair_landscape_distances = c(70700L, 70700L),
    training_fitted_component_scales = c(300L, 300L),
    five_weight_fusion_pair_rows = c(176750L, 176750L),
    runtime_projection_status = c(
      "not_estimable_before_gene_contract_revision_and_bounded_profile",
      "not_estimable_before_gene_representation_definition_and_bounded_profile"
    ),
    stringsAsFactors = FALSE
  )

  decision <- data.frame(
    contract_id = "mv06b_scaleup_gate_decision_v1",
    decision = "stop_contract_revision_required",
    complete_strict_representations = 0L,
    sct_incomplete_panel_views = summary$incomplete_panel_views[[1L]],
    sct_variance_unresolved_views = summary$variance_unresolved_views[[1L]],
    integrated_expression_undefined_views =
      summary$expression_undefined_views[[2L]],
    ph_jobs_executed = 0L, landscape_jobs_executed = 0L,
    fusion_jobs_executed = 0L, biological_outcomes_computed = FALSE,
    outcome_label_state = "closed",
    next_gate = "project_owner_selects_revised_gene_estimand_or_separate_view_disposition",
    stringsAsFactors = FALSE
  )
  list(group = group, summary = summary, workload = workload,
       decision = decision)
}
