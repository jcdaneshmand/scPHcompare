mv05d2_combine_primary_projection_v1 <- function(previous, ph_projection) {
  previous_required <- c(
    "scenario", "normalization_worker_hours",
    "measured_cell_coordinate_worker_hours", "landscape_worker_hours",
    "known_components_lower_bound_worker_hours", "unmeasured_components",
    "planning_cap_with_10_percent_reserve_hours", "outcome_label_state",
    "biological_outcomes_computed"
  )
  ph_required <- c(
    "scenario", "projected_worker_hours", "projected_storage_bytes",
    "outcome_label_state", "biological_outcomes_computed"
  )
  if (!is.data.frame(previous) ||
      !all(previous_required %in% names(previous)) ||
      !is.data.frame(ph_projection) ||
      !all(ph_required %in% names(ph_projection)) ||
      any(previous$outcome_label_state != "closed") ||
      any(ph_projection$outcome_label_state != "closed") ||
      any(as.logical(previous$biological_outcomes_computed)) ||
      any(as.logical(ph_projection$biological_outcomes_computed))) {
    stop("MV5-D2 projections must be complete and label closed.", call. = FALSE)
  }
  primary <- previous[
    previous$scenario == "resource_safe_sct_cell_primary", , drop = FALSE
  ]
  numeric_fields <- c(
    "normalization_worker_hours", "measured_cell_coordinate_worker_hours",
    "landscape_worker_hours", "known_components_lower_bound_worker_hours",
    "planning_cap_with_10_percent_reserve_hours"
  )
  if (nrow(primary) != 1L ||
      !identical(primary$unmeasured_components, "cell_ph") ||
      any(!is.finite(as.numeric(primary[numeric_fields])))) {
    stop("MV5-D1 primary projection is not ready for measured cell PH.",
         call. = FALSE)
  }
  result <- data.frame(
    contract_id = "mv05d2_post_ph_primary_projection_v1",
    scenario = paste0("resource_safe_sct_cell_primary__",
                      ph_projection$scenario),
    ph_assumption = ph_projection$scenario,
    normalization_worker_hours = primary$normalization_worker_hours,
    measured_cell_coordinate_worker_hours =
      primary$measured_cell_coordinate_worker_hours,
    projected_cell_ph_worker_hours = ph_projection$projected_worker_hours,
    landscape_worker_hours = primary$landscape_worker_hours,
    projected_total_worker_hours =
      primary$known_components_lower_bound_worker_hours +
      ph_projection$projected_worker_hours,
    projected_cell_ph_storage_bytes = ph_projection$projected_storage_bytes,
    planning_cap_with_10_percent_reserve_hours =
      primary$planning_cap_with_10_percent_reserve_hours,
    stringsAsFactors = FALSE
  )
  result$planning_cap_margin_hours <-
    result$planning_cap_with_10_percent_reserve_hours -
    result$projected_total_worker_hours
  result$feasibility_complete <- TRUE
  result$cap_passes <- result$planning_cap_margin_hours >= 0
  result$disposition <- ifelse(
    result$cap_passes,
    "bounded_projection_passes_owner_review_required_before_full_ph",
    "bounded_projection_exceeds_planning_cap_do_not_launch_full_ph"
  )
  result$full_production_jobs_launched <- 0L
  result$outcome_label_state <- "closed"
  result$biological_outcomes_computed <- FALSE
  result
}
