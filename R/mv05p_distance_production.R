.mv05p_required_group_columns <- c(
  "production_group_id", "execution_order", "representation",
  "landscape_request_rows", "energy_pair_rows",
  "shared_pseudobulk_pair_rows", "max_parallel_workers",
  "per_unit_elapsed_cap_seconds", "per_process_rss_cap_bytes",
  "stage_worker_hour_cap", "stage_private_storage_cap_bytes",
  "outcome_label_state", "biological_outcomes_computed",
  "clustering_jobs_executed", "production_executed",
  "source_freeze_sha256"
)

mv05p_validate_group_queue_v1 <- function(groups) {
  missing <- setdiff(.mv05p_required_group_columns, names(groups))
  if (length(missing) || nrow(groups) != 150L ||
      anyDuplicated(groups$production_group_id) ||
      !identical(sort(as.integer(groups$execution_order)), seq_len(150L)) ||
      !setequal(groups$representation,
                c("sct_whole", "inductive_integrated")) ||
      any(table(groups$representation) != 75L) ||
      any(groups$max_parallel_workers != 2L) ||
      any(groups$per_unit_elapsed_cap_seconds != 900) ||
      any(groups$per_process_rss_cap_bytes != 4294967296) ||
      any(groups$stage_worker_hour_cap != 21.6) ||
      any(groups$stage_private_storage_cap_bytes != 10737418240) ||
      any(groups$outcome_label_state != "closed") ||
      any(as.logical(groups$biological_outcomes_computed)) ||
      any(as.integer(groups$clustering_jobs_executed) != 0L) ||
      any(as.logical(groups$production_executed)) ||
      length(unique(groups$source_freeze_sha256)) != 1L) {
    stop("MV5-P group queue differs from the committed MV5-O freeze.",
         call. = FALSE)
  }
  invisible(groups)
}

mv05p_group_projection_v1 <- function(groups, combined, landscape,
                                      total_private_storage_bytes = 1277893355) {
  mv05p_validate_group_queue_v1(groups)
  required_combined <- c("component", "projected_worker_hours")
  required_landscape <- c("representation", "projected_worker_hours",
                          "projected_output_bytes")
  if (!all(required_combined %in% names(combined)) ||
      !all(required_landscape %in% names(landscape))) {
    stop("MV5-P projection inputs are incomplete.", call. = FALSE)
  }
  component <- function(name) {
    value <- combined$projected_worker_hours[combined$component == name]
    if (length(value) != 1L || !is.finite(value)) {
      stop("MV5-P projection component is missing: ", name, call. = FALSE)
    }
    value
  }
  landscape <- landscape[landscape$representation %in%
                           c("sct_whole", "inductive_integrated"), ]
  if (nrow(landscape) != 2L) {
    stop("MV5-P representation projections are incomplete.", call. = FALSE)
  }
  result <- groups[order(groups$execution_order), c(
    "production_group_id", "execution_order", "representation",
    "landscape_request_rows", "energy_pair_rows",
    "shared_pseudobulk_pair_rows"
  )]
  result$core_worker_hours <- 0
  result$raw_projected_bytes <- 0
  for (representation in c("sct_whole", "inductive_integrated")) {
    selected <- result$representation == representation
    landscape_hours <- landscape$projected_worker_hours[
      landscape$representation == representation]
    landscape_bytes <- as.numeric(landscape$projected_output_bytes[
      landscape$representation == representation])
    energy_hours <- component(paste0(representation, "_energy"))
    result$core_worker_hours[selected] <-
      landscape_hours * result$landscape_request_rows[selected] /
        sum(result$landscape_request_rows[selected]) +
      energy_hours * result$energy_pair_rows[selected] /
        sum(result$energy_pair_rows[selected])
    result$raw_projected_bytes[selected] <-
      landscape_bytes * as.numeric(result$landscape_request_rows[selected]) /
        sum(result$landscape_request_rows[selected]) +
      512 * as.numeric(result$energy_pair_rows[selected])
  }
  sct <- result$representation == "sct_whole"
  pseudobulk_hours <- component(
    "shared_pseudobulk_reused_across_representations")
  result$core_worker_hours[sct] <- result$core_worker_hours[sct] +
    pseudobulk_hours * result$shared_pseudobulk_pair_rows[sct] /
      sum(result$shared_pseudobulk_pair_rows[sct])
  result$raw_projected_bytes[sct] <- result$raw_projected_bytes[sct] +
    512 * as.numeric(result$shared_pseudobulk_pair_rows[sct])
  total_hours <- component("all_required_distance_matrices_with_reserve")
  result$projected_worker_seconds <- result$core_worker_hours * 3600 *
    total_hours / sum(result$core_worker_hours)
  result$projected_private_bytes <- result$raw_projected_bytes *
    total_private_storage_bytes / sum(result$raw_projected_bytes)
  if (abs(sum(result$projected_worker_seconds) / 3600 - total_hours) > 1e-10 ||
      abs(sum(result$projected_private_bytes) - total_private_storage_bytes) > 1) {
    stop("MV5-P group projection does not reconcile to the freeze.",
         call. = FALSE)
  }
  result
}

mv05p_launch_budget_v1 <- function(observed_worker_seconds,
                                   projected_incomplete_seconds,
                                   observed_private_bytes,
                                   projected_incomplete_bytes,
                                   worker_hour_cap = 21.6,
                                   storage_cap_bytes = 10737418240) {
  values <- c(observed_worker_seconds, projected_incomplete_seconds,
              observed_private_bytes, projected_incomplete_bytes)
  if (any(!is.finite(values)) || any(values < 0)) {
    stop("MV5-P launch accounting contains invalid values.", call. = FALSE)
  }
  projected_total_worker_hours <-
    (observed_worker_seconds + projected_incomplete_seconds) / 3600
  projected_total_private_bytes <-
    observed_private_bytes + projected_incomplete_bytes
  data.frame(
    projected_total_worker_hours = projected_total_worker_hours,
    worker_hour_cap = worker_hour_cap,
    worker_cap_passed = projected_total_worker_hours <= worker_hour_cap,
    projected_total_private_bytes = projected_total_private_bytes,
    storage_cap_bytes = storage_cap_bytes,
    storage_cap_passed = projected_total_private_bytes <= storage_cap_bytes,
    launch_authorized = projected_total_worker_hours <= worker_hour_cap &&
      projected_total_private_bytes <= storage_cap_bytes,
    outcome_label_state = "closed",
    biological_outcomes_computed = FALSE,
    clustering_jobs_executed = 0L,
    stringsAsFactors = FALSE
  )
}
