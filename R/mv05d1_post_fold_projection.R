# Internal MV5-D1 post-coordinate resource projection.

mv05d1_post_fold_projection_v2 <- function(previous, actual_coordinate_hours) {
  required <- c(
    "scenario", "normalization_worker_hours",
    "cached_sct_fold_worker_hours", "landscape_worker_hours",
    "integrated_reference_mapping_worker_hours", "nominal_cap_hours",
    "planning_cap_with_10_percent_reserve_hours", "outcome_label_state",
    "biological_outcomes_computed"
  )
  actual_coordinate_hours <- as.numeric(actual_coordinate_hours)
  if (!is.data.frame(previous) || !all(required %in% names(previous)) ||
      length(actual_coordinate_hours) != 1L ||
      !is.finite(actual_coordinate_hours) || actual_coordinate_hours <= 0 ||
      any(previous$outcome_label_state != "closed") ||
      any(as.logical(previous$biological_outcomes_computed))) {
    stop("Previous projection or measured coordinate time is invalid.",
         call. = FALSE)
  }
  resource_safe <- startsWith(previous$scenario, "resource_safe_")
  known_lower_bound <- previous$projected_lower_bound_worker_hours
  for (index in which(resource_safe)) {
    known_lower_bound[[index]] <- sum(c(
      previous$normalization_worker_hours[[index]], actual_coordinate_hours,
      previous$landscape_worker_hours[[index]],
      previous$integrated_reference_mapping_worker_hours[[index]]
    ), na.rm = TRUE)
  }
  unmeasured <- rep("none_within_legacy_projection", nrow(previous))
  unmeasured[previous$scenario == "resource_safe_all_planned_views_lower_bound"] <-
    "cell_ph;gene_ph;integrated_reference_mapping"
  unmeasured[previous$scenario == "resource_safe_sct_cell_gene"] <-
    "cell_ph;gene_ph"
  unmeasured[previous$scenario == "resource_safe_sct_cell_primary"] <-
    "cell_ph"
  disposition <- previous$disposition
  disposition[previous$scenario == "resource_safe_all_planned_views_lower_bound"] <-
    "prohibited_scope_and_unmeasured_components"
  disposition[previous$scenario == "resource_safe_sct_cell_gene"] <-
    "deferred_gene_scope_and_unmeasured_ph"
  disposition[previous$scenario == "resource_safe_sct_cell_primary"] <-
    "deferred_pending_measured_production_cell_ph"
  data.frame(
    contract_id = "mv05d1_post_fold_resource_projection_v2",
    scenario = previous$scenario,
    normalization_worker_hours = previous$normalization_worker_hours,
    legacy_compound_fold_projection_hours =
      previous$cached_sct_fold_worker_hours,
    measured_cell_coordinate_worker_hours = ifelse(
      resource_safe, actual_coordinate_hours, NA_real_
    ),
    cell_ph_worker_hours = NA_real_,
    landscape_worker_hours = previous$landscape_worker_hours,
    integrated_reference_mapping_worker_hours =
      previous$integrated_reference_mapping_worker_hours,
    known_components_lower_bound_worker_hours = known_lower_bound,
    unmeasured_components = unmeasured,
    feasibility_complete = !resource_safe,
    nominal_cap_hours = previous$nominal_cap_hours,
    planning_cap_with_10_percent_reserve_hours =
      previous$planning_cap_with_10_percent_reserve_hours,
    cap_passes = ifelse(resource_safe, NA, previous$cap_passes),
    disposition = disposition,
    outcome_label_state = "closed",
    biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
}
