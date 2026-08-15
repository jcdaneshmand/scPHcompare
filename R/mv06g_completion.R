# Internal MV6-G serial completion and immutable-resume contracts.

.mv06g_completion_success <- c("completed", "reused_validated")

mv06g_completion_source_paths_v1 <- function() {
  c(
    "R/mv06g_completion.R",
    "scripts/build_mv06g_completion_prefreeze.R",
    "scripts/validate_mv06g_completion_prefreeze.R",
    "scripts/validate_mv06g_completion_prefreeze_repeat.R",
    "scripts/run_mv06g_completion_monitor.R",
    "scripts/run_mv06g_completion.R",
    "scripts/validate_mv06g_complete.R",
    "scripts/check_mv06g_complete_resume.R",
    "tests/testthat/test-mv06g-completion.R"
  )
}

mv06g_safe_group_name_v1 <- function(value) {
  gsub("[^A-Za-z0-9_.-]", "_", value)
}

mv06g_completion_root_v1 <- function(sources) {
  if (!is.data.frame(sources) || nrow(sources) != 9L ||
      !identical(sources$path, mv06g_completion_source_paths_v1()) ||
      anyDuplicated(sources$path) || !all(file.exists(sources$path)) ||
      !identical(unname(vapply(sources$path, .mv06f_sha256, character(1L))),
                 unname(sources$sha256))) {
    stop("MV6-G completion source inventory is stale.", call. = FALSE)
  }
  digest::digest(stats::setNames(sources$sha256, sources$path),
                 algo = "sha256", serialize = TRUE)
}

mv06g_validate_completion_policy_v1 <- function(
    policy, queue, parent, rebind_policy, rebind_equivalence, sources,
    rust_library) {
  required <- c(
    "contract_id", "parent_rebind_policy_sha256",
    "rebind_equivalence_sha256", "parent_contract_sha256",
    "queue_root_sha256", "scientific_implementation_root_sha256",
    "execution_implementation_root_sha256", "rust_library_sha256",
    "remaining_groups", "elapsed_cap_seconds_per_group",
    "rss_cap_bytes_per_group", "private_storage_cap_bytes",
    "total_worker_cap_seconds", "maximum_workers", "automatic_retry",
    "remaining_production_authorized", "outcome_label_state",
    "biological_outcomes_computed", "fusion_evaluations", "outcome_jobs",
    "disposition"
  )
  root <- mv06g_completion_root_v1(sources)
  if (!is.data.frame(policy) || nrow(policy) != 1L ||
      !all(required %in% names(policy)) ||
      policy$contract_id != "mv06g_serial_completion_prefreeze_v1" ||
      !is.data.frame(queue) || nrow(queue) != 74L ||
      !identical(as.integer(queue$execution_order), 2:75) ||
      anyDuplicated(queue$group_id) ||
      policy$parent_rebind_policy_sha256 !=
        unique(rebind_policy$.file_sha256) ||
      policy$rebind_equivalence_sha256 !=
        unique(rebind_equivalence$.file_sha256) ||
      policy$parent_contract_sha256 != .mv06f_sha256(parent$.file_path) ||
      policy$queue_root_sha256 != parent$queue_root_sha256 ||
      policy$scientific_implementation_root_sha256 !=
        rebind_policy$production_implementation_root_sha256 ||
      policy$execution_implementation_root_sha256 != root ||
      policy$rust_library_sha256 != .mv06f_sha256(rust_library) ||
      policy$remaining_groups != 74L ||
      policy$elapsed_cap_seconds_per_group != 1800L ||
      policy$rss_cap_bytes_per_group != 12 * 1024^3 ||
      policy$private_storage_cap_bytes != 5 * 1024^3 ||
      policy$total_worker_cap_seconds != 12 * 3600 ||
      policy$maximum_workers != 1L || as.logical(policy$automatic_retry) ||
      !as.logical(policy$remaining_production_authorized) ||
      policy$outcome_label_state != "closed" ||
      as.logical(policy$biological_outcomes_computed) ||
      policy$fusion_evaluations != 0L || policy$outcome_jobs != 0L ||
      policy$disposition != "prefreeze_pass_serial_label_closed_only" ||
      nrow(rebind_equivalence) != 3L ||
      any(!as.logical(rebind_equivalence$sha256_identical)) ||
      any(!as.logical(rebind_equivalence$bytes_identical)) ||
      any(!as.logical(rebind_equivalence$resource_pass))) {
    stop("MV6-G serial completion policy is stale.", call. = FALSE)
  }
  invisible(TRUE)
}

mv06g_completion_guard_v1 <- function(prior_worker_seconds, active_seconds,
                                       peak_rss_bytes, private_bytes, policy) {
  values <- c(prior_worker_seconds, active_seconds, peak_rss_bytes,
              private_bytes)
  if (any(!is.finite(values)) || any(values < 0)) {
    stop("MV6-G completion guard inputs are invalid.", call. = FALSE)
  }
  disposition <- NA_character_
  if (active_seconds > policy$elapsed_cap_seconds_per_group)
    disposition <- "elapsed_cap"
  if (peak_rss_bytes > policy$rss_cap_bytes_per_group)
    disposition <- "rss_cap"
  if (private_bytes > policy$private_storage_cap_bytes)
    disposition <- "storage_cap"
  if (prior_worker_seconds + active_seconds > policy$total_worker_cap_seconds)
    disposition <- "worker_cap"
  list(pass = is.na(disposition), disposition = disposition)
}

mv06g_validate_completion_metric_v1 <- function(metric, queue_row, policy) {
  required <- c(
    "contract_id", "group_id", "execution_order", "disposition",
    "exit_status", "elapsed_seconds", "charged_worker_seconds",
    "peak_process_tree_rss_bytes", "cumulative_worker_seconds",
    "cumulative_private_bytes", "scientific_group_complete",
    "scientific_implementation_root_sha256",
    "execution_implementation_root_sha256", "rust_library_sha256",
    "outcome_label_state", "biological_outcomes_computed",
    "fusion_evaluations", "outcome_jobs"
  )
  if (!is.data.frame(metric) || nrow(metric) != 1L ||
      !all(required %in% names(metric)) ||
      metric$contract_id != "mv06g_serial_completion_metric_v1" ||
      metric$group_id != queue_row$group_id ||
      metric$execution_order != queue_row$execution_order ||
      !metric$disposition %in% .mv06g_completion_success ||
      metric$exit_status != 0L || !as.logical(metric$scientific_group_complete) ||
      metric$peak_process_tree_rss_bytes > policy$rss_cap_bytes_per_group ||
      metric$cumulative_worker_seconds > policy$total_worker_cap_seconds ||
      metric$cumulative_private_bytes > policy$private_storage_cap_bytes ||
      metric$scientific_implementation_root_sha256 !=
        policy$scientific_implementation_root_sha256 ||
      metric$execution_implementation_root_sha256 !=
        policy$execution_implementation_root_sha256 ||
      metric$rust_library_sha256 != policy$rust_library_sha256 ||
      metric$outcome_label_state != "closed" ||
      as.logical(metric$biological_outcomes_computed) ||
      metric$fusion_evaluations != 0L || metric$outcome_jobs != 0L) {
    stop("MV6-G completion metric is stale or failed.", call. = FALSE)
  }
  invisible(TRUE)
}
