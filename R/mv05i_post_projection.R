mv05i_measured_primary_projection_v1 <- function(
    previous, completion, retrieval_reference, reserve = 1.25) {
  required_previous <- c(
    "measured_coordinate_worker_hours", "measured_ph_worker_hours",
    "measured_coordinate_ph_storage_bytes", "worker_hour_cap",
    "storage_cap_bytes", "measured_peak_rss_bytes", "rss_cap_bytes",
    "measured_max_group_seconds", "group_cap_seconds",
    "outcome_label_state", "biological_outcomes_computed"
  )
  required_completion <- c(
    "group_worker_seconds", "peak_process_tree_rss_bytes",
    "private_total_landscape_stage_bytes",
    "full_cell_landscape_distances_complete",
    "outcome_label_state", "biological_outcomes_computed"
  )
  required_retrieval <- c(
    "elapsed_seconds", "private_result_bytes", "outcome_label_state",
    "biological_outcomes_computed"
  )
  closed <- function(value) {
    all(value$outcome_label_state == "closed") &&
      !any(as.logical(value$biological_outcomes_computed))
  }
  if (!is.data.frame(previous) || nrow(previous) != 1L ||
      !all(required_previous %in% names(previous)) ||
      !is.data.frame(completion) || nrow(completion) != 1L ||
      !all(required_completion %in% names(completion)) ||
      !is.data.frame(retrieval_reference) ||
      !all(required_retrieval %in% names(retrieval_reference)) ||
      !closed(previous) || !closed(completion) || !closed(retrieval_reference) ||
      !as.logical(completion$full_cell_landscape_distances_complete) ||
      !is.numeric(reserve) || length(reserve) != 1L ||
      !is.finite(reserve) || reserve < 1) {
    stop("MV5-I projection inputs violate completion or label gates.",
         call. = FALSE)
  }
  landscape_hours <- completion$group_worker_seconds / 3600
  retrieval_hours <- sum(retrieval_reference$elapsed_seconds) *
    reserve / 3600
  total_hours <- previous$measured_coordinate_worker_hours +
    previous$measured_ph_worker_hours + landscape_hours + retrieval_hours
  measured_storage <- as.numeric(
    previous$measured_coordinate_ph_storage_bytes
  ) + as.numeric(completion$private_total_landscape_stage_bytes)
  retrieval_storage <- sum(as.numeric(
    retrieval_reference$private_result_bytes
  )) * reserve
  total_storage <- measured_storage + retrieval_storage
  peak_rss <- max(
    previous$measured_peak_rss_bytes,
    completion$peak_process_tree_rss_bytes
  )
  maximum_group_seconds <- max(
    previous$measured_max_group_seconds,
    completion$maximum_group_seconds
  )
  passes <- total_hours <= previous$worker_hour_cap &&
    total_storage <= previous$storage_cap_bytes &&
    peak_rss <= previous$rss_cap_bytes &&
    maximum_group_seconds <= previous$group_cap_seconds
  data.frame(
    contract_id = "mv05i_post_distance_projection_summary_v1",
    scenario = "integrated_cell_exact_landscapes_then_label_closed_retrieval",
    measured_coordinate_worker_hours =
      previous$measured_coordinate_worker_hours,
    measured_ph_worker_hours = previous$measured_ph_worker_hours,
    measured_landscape_worker_hours = landscape_hours,
    projected_retrieval_worker_hours = retrieval_hours,
    projected_total_worker_hours = total_hours,
    worker_hour_cap = previous$worker_hour_cap,
    measured_coordinate_ph_landscape_storage_bytes = measured_storage,
    projected_retrieval_storage_bytes = retrieval_storage,
    projected_total_storage_bytes = total_storage,
    storage_cap_bytes = previous$storage_cap_bytes,
    measured_peak_rss_bytes = peak_rss,
    rss_cap_bytes = previous$rss_cap_bytes,
    measured_max_group_seconds = maximum_group_seconds,
    group_cap_seconds = previous$group_cap_seconds,
    full_integrated_landscape_distances_complete = TRUE,
    decision = if (passes)
      "authorize_separate_label_closed_integrated_retrieval_input_sprint"
    else "do_not_authorize_or_narrow_scope",
    next_stage_authorized = passes,
    current_sprint_retrieval_jobs = 0L,
    current_sprint_clustering_jobs = 0L,
    current_sprint_gene_view_jobs = 0L,
    current_sprint_fusion_jobs = 0L,
    current_sprint_new_data_jobs = 0L,
    outcome_label_state = "closed",
    biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
}
