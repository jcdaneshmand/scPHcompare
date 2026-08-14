mv05d3_measured_primary_projection_v1 <- function(previous, group_metrics) {
  required <- c(
    "scenario", "normalization_worker_hours",
    "measured_cell_coordinate_worker_hours", "landscape_worker_hours",
    "known_components_lower_bound_worker_hours",
    "planning_cap_with_10_percent_reserve_hours", "outcome_label_state",
    "biological_outcomes_computed"
  )
  if (!is.data.frame(previous) || !all(required %in% names(previous)) ||
      any(previous$outcome_label_state != "closed") ||
      any(as.logical(previous$biological_outcomes_computed))) {
    stop("MV5-D3 requires the label-closed MV5-D1 projection.", call. = FALSE)
  }
  mv05d3_validate_group_metrics_v1(
    group_metrics, expected_groups = nrow(group_metrics)
  )
  primary <- previous[
    previous$scenario == "resource_safe_sct_cell_primary", , drop = FALSE
  ]
  if (nrow(primary) != 1L) {
    stop("MV5-D3 primary projection row is missing or duplicated.",
         call. = FALSE)
  }
  ph_hours <- sum(group_metrics$elapsed_seconds) / 3600
  result <- data.frame(
    contract_id = "mv05d3_measured_primary_projection_v1",
    scenario = "resource_safe_sct_cell_primary_measured_ph",
    normalization_worker_hours = primary$normalization_worker_hours,
    measured_cell_coordinate_worker_hours =
      primary$measured_cell_coordinate_worker_hours,
    measured_cell_ph_worker_hours = ph_hours,
    projected_landscape_worker_hours = primary$landscape_worker_hours,
    measured_plus_projected_total_worker_hours =
      primary$known_components_lower_bound_worker_hours + ph_hours,
    measured_cell_ph_storage_bytes = sum(group_metrics$private_result_bytes),
    planning_cap_with_10_percent_reserve_hours =
      primary$planning_cap_with_10_percent_reserve_hours,
    stringsAsFactors = FALSE
  )
  result$planning_cap_margin_hours <-
    result$planning_cap_with_10_percent_reserve_hours -
    result$measured_plus_projected_total_worker_hours
  result$cap_passes <- result$planning_cap_margin_hours >= 0
  result$full_cell_ph_complete <- TRUE
  result$landscapes_complete <- FALSE
  result$disposition <- ifelse(
    result$cap_passes,
    "cell_ph_complete_landscape_stage_requires_separate_authorization",
    "cell_ph_complete_but_remaining_projection_exceeds_planning_cap"
  )
  result$outcome_label_state <- "closed"
  result$biological_outcomes_computed <- FALSE
  result
}
