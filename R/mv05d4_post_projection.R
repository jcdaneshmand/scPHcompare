mv05d4_measured_primary_projection_v1 <- function(previous, completion) {
  required_previous <- c(
    "normalization_worker_hours", "measured_cell_coordinate_worker_hours",
    "measured_cell_ph_worker_hours", "projected_landscape_worker_hours",
    "planning_cap_with_10_percent_reserve_hours", "outcome_label_state",
    "biological_outcomes_computed"
  )
  required_completion <- c(
    "group_worker_seconds", "full_cell_landscape_distances_complete",
    "outcome_label_state", "biological_outcomes_computed"
  )
  if (!is.data.frame(previous) || nrow(previous) != 1L ||
      !all(required_previous %in% names(previous)) ||
      !is.data.frame(completion) || nrow(completion) != 1L ||
      !all(required_completion %in% names(completion)) ||
      previous$outcome_label_state != "closed" ||
      as.logical(previous$biological_outcomes_computed) ||
      completion$outcome_label_state != "closed" ||
      as.logical(completion$biological_outcomes_computed) ||
      !as.logical(completion$full_cell_landscape_distances_complete)) {
    stop("MV5-D4 projection inputs violate completion or label gates.",
         call. = FALSE)
  }
  landscape_hours <- completion$group_worker_seconds / 3600
  total <- previous$normalization_worker_hours +
    previous$measured_cell_coordinate_worker_hours +
    previous$measured_cell_ph_worker_hours + landscape_hours
  cap <- previous$planning_cap_with_10_percent_reserve_hours
  data.frame(
    contract_id = "mv05d4_measured_primary_projection_v1",
    scenario = "resource_safe_sct_cell_primary_measured_landscapes",
    normalization_worker_hours = previous$normalization_worker_hours,
    measured_cell_coordinate_worker_hours =
      previous$measured_cell_coordinate_worker_hours,
    measured_cell_ph_worker_hours = previous$measured_cell_ph_worker_hours,
    previous_projected_landscape_worker_hours =
      previous$projected_landscape_worker_hours,
    measured_landscape_worker_hours = landscape_hours,
    measured_total_worker_hours = total,
    planning_cap_with_10_percent_reserve_hours = cap,
    planning_cap_margin_hours = cap - total,
    cap_passes = total <= cap,
    full_cell_landscape_distances_complete = TRUE,
    disposition = ifelse(
      total <= cap,
      "cell_primary_precomputation_complete_next_gate_separate",
      "cell_primary_complete_but_planning_cap_exceeded"
    ),
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
}
