# Internal MV6-F stage-two rebind and execution contracts.

.mv06f_stage2_success <- c("completed", "reused_validated")

mv06f_stage2_source_root_v1 <- function(sources) {
  required <- c("path", "sha256", "outcome_label_state",
                "biological_outcomes_computed")
  if (!is.data.frame(sources) || !nrow(sources) ||
      !all(required %in% names(sources)) || anyDuplicated(sources$path) ||
      !all(file.exists(sources$path)) ||
      !identical(tolower(unname(vapply(
        sources$path, .mv06f_sha256, character(1L)
      ))), tolower(unname(sources$sha256))) ||
      any(sources$outcome_label_state != "closed") ||
      any(as.logical(sources$biological_outcomes_computed))) {
    stop("MV6-F stage-two source inventory is stale.", call. = FALSE)
  }
  digest::digest(stats::setNames(sources$sha256, sources$path),
                 algo = "sha256", serialize = TRUE)
}

mv06f_validate_stage2_rebind_contract_v1 <- function(
    queue, contract, sources, resource_plan, rust_library) {
  mv06f_validate_queue_v1(queue)
  required_contract <- c(
    "contract_id", "base_revision", "queue_root_sha256",
    "parent_implementation_root_sha256", "implementation_root_sha256",
    "rust_library_sha256", "groups", "ph_jobs",
    "diagram_dimension_records", "biological_pairs",
    "landscape_component_rows", "stage1_reexecution_required",
    "stage2_authorized", "outcome_label_state",
    "biological_outcomes_computed", "fusion_jobs", "clustering_jobs",
    "outcome_jobs"
  )
  expected_guards <- c(
    group_elapsed_seconds = 1800,
    group_process_tree_rss_bytes = 8 * 1024^3,
    concurrent_process_tree_rss_bytes = 12 * 1024^3,
    production_worker_seconds = 14.4 * 3600,
    private_root_bytes = 10 * 1024^3,
    maximum_workers = 1
  )
  guards <- if (is.data.frame(resource_plan) &&
                 all(c("guard", "value") %in% names(resource_plan))) {
    stats::setNames(as.numeric(resource_plan$value), resource_plan$guard)
  } else numeric()
  source_root <- mv06f_stage2_source_root_v1(sources)
  if (!is.data.frame(contract) || nrow(contract) != 1L ||
      !all(required_contract %in% names(contract)) ||
      contract$contract_id != "mv06f_stage2_rebind_prefreeze_v1" ||
      contract$queue_root_sha256 != mv06f_queue_root_v1(queue) ||
      contract$implementation_root_sha256 != source_root ||
      contract$rust_library_sha256 != .mv06f_sha256(rust_library) ||
      contract$groups != 75L || contract$ph_jobs != 13500L ||
      contract$diagram_dimension_records != 27000L ||
      contract$biological_pairs != 35350L ||
      contract$landscape_component_rows != 141400L ||
      !as.logical(contract$stage1_reexecution_required) ||
      as.logical(contract$stage2_authorized) ||
      contract$outcome_label_state != "closed" ||
      as.logical(contract$biological_outcomes_computed) ||
      any(unlist(contract[c(
        "fusion_jobs", "clustering_jobs", "outcome_jobs"
      )], use.names = FALSE) != 0) ||
      !identical(names(expected_guards), names(guards)) ||
      !identical(as.numeric(expected_guards), as.numeric(guards))) {
    stop("MV6-F stage-two rebind contract is stale or unauthorized.",
         call. = FALSE)
  }
  invisible(TRUE)
}

mv06f_validate_stage2_admission_v1 <- function(admission, contract) {
  required <- c(
    "contract_id", "queue_root_sha256", "implementation_root_sha256",
    "rust_library_sha256", "stage1_scientific_equivalence_passed",
    "stage1_clean_repeat_passed", "stage1_r_oracles_passed",
    "stage1_persim_oracles_passed", "stage1_resume_passed",
    "stage2_authorized", "outcome_label_state",
    "biological_outcomes_computed", "fusion_jobs", "clustering_jobs",
    "outcome_jobs"
  )
  gates <- c(
    "stage1_scientific_equivalence_passed", "stage1_clean_repeat_passed",
    "stage1_r_oracles_passed", "stage1_persim_oracles_passed",
    "stage1_resume_passed", "stage2_authorized"
  )
  if (!is.data.frame(admission) || nrow(admission) != 1L ||
      !all(required %in% names(admission)) ||
      admission$contract_id != "mv06f_stage2_admission_v1" ||
      admission$queue_root_sha256 != contract$queue_root_sha256 ||
      admission$implementation_root_sha256 !=
        contract$implementation_root_sha256 ||
      admission$rust_library_sha256 != contract$rust_library_sha256 ||
      !all(as.logical(admission[gates])) ||
      admission$outcome_label_state != "closed" ||
      as.logical(admission$biological_outcomes_computed) ||
      any(unlist(admission[c(
        "fusion_jobs", "clustering_jobs", "outcome_jobs"
      )], use.names = FALSE) != 0)) {
    stop("MV6-F stage-two admission evidence is incomplete.", call. = FALSE)
  }
  invisible(TRUE)
}

mv06f_stage2_guard_v1 <- function(prior_worker_seconds, active_seconds,
                                  group_elapsed_seconds, group_peak_rss_bytes,
                                  active_rss_bytes, private_root_bytes,
                                  resource_plan) {
  values <- c(
    prior_worker_seconds, active_seconds, group_elapsed_seconds,
    group_peak_rss_bytes, active_rss_bytes, private_root_bytes
  )
  if (length(values) != 6L || any(!is.finite(values)) || any(values < 0)) {
    stop("MV6-F stage-two guard inputs are invalid.", call. = FALSE)
  }
  limits <- stats::setNames(as.numeric(resource_plan$value),
                            resource_plan$guard)
  disposition <- NA_character_
  if (group_elapsed_seconds > limits[["group_elapsed_seconds"]])
    disposition <- "group_elapsed_cap_exceeded"
  if (group_peak_rss_bytes > limits[["group_process_tree_rss_bytes"]])
    disposition <- "group_rss_cap_exceeded"
  if (active_rss_bytes > limits[["concurrent_process_tree_rss_bytes"]])
    disposition <- "concurrent_rss_cap_exceeded"
  if (prior_worker_seconds + active_seconds >
      limits[["production_worker_seconds"]])
    disposition <- "production_worker_cap_exceeded"
  if (private_root_bytes > limits[["private_root_bytes"]])
    disposition <- "private_storage_cap_exceeded"
  list(launch_authorized = is.na(disposition), disposition = disposition)
}

mv06f_validate_stage2_metrics_v1 <- function(metrics, queue, contract,
                                              complete = FALSE) {
  stage2 <- queue[queue$stage == "stage_2", , drop = FALSE]
  required <- c(
    "contract_id", "group_id", "execution_order", "disposition",
    "exit_status", "elapsed_seconds", "charged_worker_seconds",
    "peak_process_tree_rss_bytes", "cumulative_private_bytes",
    "queue_root_sha256", "implementation_root_sha256",
    "rust_library_sha256", "group_directory_complete",
    "outcome_label_state", "biological_outcomes_computed", "fusion_jobs",
    "clustering_jobs", "outcome_jobs"
  )
  if (!is.data.frame(metrics) || !all(required %in% names(metrics)) ||
      nrow(metrics) > 74L || anyDuplicated(metrics$group_id) ||
      any(!metrics$group_id %in% stage2$group_id) ||
      any(!metrics$disposition %in% .mv06f_stage2_success) ||
      any(metrics$exit_status != 0L) ||
      any(!as.logical(metrics$group_directory_complete)) ||
      any(metrics$queue_root_sha256 != contract$queue_root_sha256) ||
      any(metrics$implementation_root_sha256 !=
            contract$implementation_root_sha256) ||
      any(metrics$rust_library_sha256 != contract$rust_library_sha256) ||
      any(metrics$outcome_label_state != "closed") ||
      any(as.logical(metrics$biological_outcomes_computed)) ||
      any(unlist(metrics[c(
        "fusion_jobs", "clustering_jobs", "outcome_jobs"
      )], use.names = FALSE) != 0)) {
    stop("MV6-F stage-two metrics are stale or incomplete.", call. = FALSE)
  }
  expected_order <- stage2$execution_order[match(metrics$group_id,
                                                  stage2$group_id)]
  if (!identical(as.integer(metrics$execution_order),
                 as.integer(expected_order)) ||
      (isTRUE(complete) && nrow(metrics) != 74L)) {
    stop("MV6-F stage-two metric order or completion is invalid.",
         call. = FALSE)
  }
  invisible(metrics)
}
