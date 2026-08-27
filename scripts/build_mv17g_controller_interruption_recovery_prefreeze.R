#!/usr/bin/env Rscript
options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 11L) {
  stop(
    paste(
      "usage: build_mv17g_controller_interruption_recovery_prefreeze.R",
      "<original-public-prefreeze> <original-private-prefreeze>",
      "<parallel-public-prefreeze> <parallel-private-prefreeze>",
      "<primary-private-root> <attempt1-time> <attempt1-stdout> <attempt1-stderr>",
      "<private-output> <public-output> <implementation-head>"
    ),
    call. = FALSE
  )
}
original_public <- normalizePath(args[[1]], mustWork = TRUE)
original_private <- normalizePath(args[[2]], mustWork = TRUE)
parallel_public <- normalizePath(args[[3]], mustWork = TRUE)
parallel_private <- normalizePath(args[[4]], mustWork = TRUE)
primary_private <- normalizePath(args[[5]], mustWork = TRUE)
attempt1_time <- normalizePath(args[[6]], mustWork = TRUE)
attempt1_stdout <- normalizePath(args[[7]], mustWork = TRUE)
attempt1_stderr <- normalizePath(args[[8]], mustWork = TRUE)
private <- args[[9]]
public <- args[[10]]
implementation_head <- tolower(args[[11]])
if (dir.exists(private) || dir.exists(public) || !grepl("^[0-9a-f]{40}$", implementation_head)) {
  stop("invalid MV17-G controller-recovery target/head", call. = FALSE)
}
dir.create(private, recursive = TRUE)
dir.create(public, recursive = TRUE)

source("R/mv08z_landscape_production.R")
source("R/mv17_null_models.R")
source("R/mv17_calibration.R")
source("R/mv17_full_calibration.R")
source("R/mv17g_parallel_recovery.R")
r <- .mv08z_read_csv
w <- .mv08z_atomic_csv
h <- .mv08z_sha256_file

verify_manifest <- function(root, name) {
  manifest <- r(file.path(root, name))
  paths <- file.path(root, manifest$artifact)
  if (!all(file.exists(paths)) ||
      !identical(unname(as.numeric(file.info(paths)$size)), unname(as.numeric(manifest$bytes))) ||
      !identical(unname(vapply(paths, h, character(1L))), unname(tolower(manifest$sha256)))) {
    stop("MV17-G controller-recovery manifest drift", call. = FALSE)
  }
  manifest
}
original_manifest <- verify_manifest(original_public, "mv17g-artifact-manifest.csv")
parallel_manifest <- verify_manifest(parallel_public, "mv17g-parallel-artifact-manifest.csv")

verify_private <- function(public_root, private_root, binding_name) {
  binding <- r(file.path(public_root, binding_name))
  paths <- file.path(private_root, binding$artifact)
  if (!all(file.exists(paths)) ||
      !identical(unname(as.numeric(file.info(paths)$size)), unname(as.numeric(binding$bytes))) ||
      !identical(unname(vapply(paths, h, character(1L))), unname(tolower(binding$sha256)))) {
    stop("MV17-G controller-recovery private prefreeze drift", call. = FALSE)
  }
  binding
}
original_private_binding <- verify_private(original_public, original_private, "mv17g-private-binding.csv")
parallel_private_binding <- verify_private(parallel_public, parallel_private, "mv17g-parallel-private-binding.csv")

queue <- r(file.path(original_private, "mv17g-primary-grouped-queue.csv"))
contract <- r(file.path(parallel_public, "mv17g-parallel-contract.csv"))
matrix_catalog <- r(file.path(parallel_private, "mv17g-parallel-matrix-catalog.csv"))
matrix_paths <- file.path(primary_private, "matrices", sprintf("%s__%03d.rds", matrix_catalog$view, matrix_catalog$unit_order))
if (nrow(queue) != 1188L || sum(queue$scientific_runs) != 91740L || nrow(matrix_catalog) != 264L ||
    !all(file.exists(matrix_paths)) ||
    !identical(unname(as.numeric(file.info(matrix_paths)$size)), unname(as.numeric(matrix_catalog$bytes))) ||
    !identical(unname(vapply(matrix_paths, h, character(1L))), unname(tolower(matrix_catalog$sha256)))) {
  stop("MV17-G controller-recovery queue/matrix drift", call. = FALSE)
}

scan <- mv17g_checkpoint_scan_v1(queue, primary_private)
stray_partials <- list.files(primary_private, pattern = "[.]partial$", recursive = TRUE, full.names = TRUE)
if (length(stray_partials) || any(scan$state == "partial")) {
  stop("MV17-G controller-recovery checkpoint contains partial evidence", call. = FALSE)
}
artifact_prefix <- mv17g_complete_prefix_v1(scan, require_incomplete = TRUE)
progress <- r(file.path(primary_private, "mv17g-private-progress.csv"))
admitted_prefix <- nrow(progress)
wave_orders <- seq.int(admitted_prefix + 1L, artifact_prefix)
if (!identical(as.integer(progress$job_order), seq_len(admitted_prefix)) ||
    artifact_prefix != 794L || admitted_prefix != 786L ||
    !identical(wave_orders, 787:794)) {
  stop("MV17-G controller-recovery interruption boundary drift", call. = FALSE)
}

validate_result <- function(result, q, matrix_path) {
  expected_seeds <- if (q$null_family == "observed") 0L else seq.int(q$seed_first, length.out = q$replicate_count)
  ok <- identical(result$contract_id, "mv17g_group_result_v1") &&
    identical(result$null_family, q$null_family) &&
    result$seed_first == min(expected_seeds) && result$seed_last == max(expected_seeds) &&
    result$replicate_count == length(expected_seeds) && result$matrix_sha256 == h(matrix_path) &&
    nrow(result$metrics) == length(expected_seeds) * 8L &&
    setequal(unique(result$metrics$seed), expected_seeds) &&
    all(result$metrics$h0_mst_maximum_absolute_error <= 1e-8) && result$finite
  if (!isTRUE(ok)) stop("MV17-G controller-recovery payload drift at order ", q$job_order, call. = FALSE)
  TRUE
}

inspect_one <- function(i, receipt_required) {
  q <- queue[i, , drop = FALSE]
  paths <- mv17g_job_artifacts_v1(q, primary_private)
  if (!all(file.exists(paths)) || file.info(paths[["stdout"]])$size != 0 || file.info(paths[["stderr"]])$size != 0) {
    stop("MV17-G controller-recovery artifact/stream drift at order ", q$job_order, call. = FALSE)
  }
  matrix_path <- file.path(primary_private, "matrices", sprintf("%s__%03d.rds", q$view, q$unit_order))
  validate_result(readRDS(paths[["result"]]), q, matrix_path)
  resource <- try(mv17c_parse_gnu_time_v1(paths[["time"]]), silent = TRUE)
  receipt_valid <- !inherits(resource, "try-error")
  if (receipt_required && !receipt_valid) {
    stop("MV17-G admitted prefix has an invalid resource receipt at order ", q$job_order, call. = FALSE)
  }
  if (receipt_valid && (resource$exit_status != 0L ||
      resource$wall_seconds > contract$child_timeout_seconds ||
      resource$maximum_RSS_bytes > contract$child_RSS_cap_bytes)) {
    stop("MV17-G controller-recovery resource drift at order ", q$job_order, call. = FALSE)
  }
  rows <- do.call(rbind, lapply(names(paths), function(role) {
    data.frame(
      job_order = q$job_order,
      role = role,
      artifact = normalizePath(paths[[role]]),
      bytes = as.numeric(file.info(paths[[role]])$size),
      sha256 = h(paths[[role]]),
      time_receipt_valid = receipt_valid,
      stringsAsFactors = FALSE
    )
  }))
  list(rows = rows, receipt_valid = receipt_valid)
}

admitted <- lapply(seq_len(admitted_prefix), inspect_one, receipt_required = TRUE)
wave <- lapply(wave_orders, inspect_one, receipt_required = FALSE)
admitted_inventory <- do.call(rbind, lapply(admitted, `[[`, "rows"))
wave_inventory <- do.call(rbind, lapply(wave, `[[`, "rows"))
wave_receipt_valid <- vapply(wave, `[[`, logical(1L), "receipt_valid")
valid_wave_orders <- wave_orders[wave_receipt_valid]
incomplete_wave_orders <- wave_orders[!wave_receipt_valid]
if (!identical(valid_wave_orders, c(789L, 791L, 794L)) ||
    !identical(incomplete_wave_orders, c(787L, 788L, 790L, 792L, 793L))) {
  stop("MV17-G interrupted-wave receipt classification drift", call. = FALSE)
}

ps_args <- system2("ps", c("-eo", "args="), stdout = TRUE)
controller_processes <- sum(grepl("run_mv17g_calibration_parallel_recovery[.]R", ps_args))
attempt1_time_bytes <- as.numeric(file.info(attempt1_time)$size)
attempt1_stdout_bytes <- as.numeric(file.info(attempt1_stdout)$size)
attempt1_stderr_bytes <- as.numeric(file.info(attempt1_stderr)$size)
if (controller_processes != 0L || attempt1_time_bytes != 0 || attempt1_stderr_bytes != 0 || attempt1_stdout_bytes < 1) {
  stop("MV17-G controller interruption evidence drift", call. = FALSE)
}
result_paths <- vapply(seq_len(artifact_prefix), function(i) {
  mv17g_job_artifacts_v1(queue[i, , drop = FALSE], primary_private)[["result"]]
}, character(1L))
attempt1_start <- file.info(attempt1_stdout)$mtime
attempt1_last_result <- max(file.info(result_paths)$mtime)
attempt1_observed_seconds <- as.numeric(difftime(attempt1_last_result, attempt1_start, units = "secs"))
if (!is.finite(attempt1_observed_seconds) || attempt1_observed_seconds <= 0 ||
    attempt1_observed_seconds > contract$aggregate_timeout_seconds) {
  stop("MV17-G controller interruption timestamp evidence invalid", call. = FALSE)
}

state <- data.frame(
  contract_id = "mv17g_controller_interruption_state_v2",
  artifact_present_prefix_children = artifact_prefix,
  resource_admitted_prefix_children = admitted_prefix,
  durable_progress_rows = admitted_prefix,
  interrupted_wave_first_order = min(wave_orders),
  interrupted_wave_last_order = max(wave_orders),
  interrupted_wave_children = length(wave_orders),
  receipt_complete_wave_children = length(valid_wave_orders),
  receipt_incomplete_wave_children = length(incomplete_wave_orders),
  runner_restart_job_order = min(wave_orders),
  resume_after_replayed_wave_job_order = max(wave_orders) + 1L,
  pending_children_including_replay = nrow(queue) - admitted_prefix,
  pending_scientific_runs_including_replay = sum(queue$scientific_runs[(admitted_prefix + 1L):nrow(queue)]),
  controller_processes = controller_processes,
  attempt1_GNU_time_bytes = attempt1_time_bytes,
  attempt1_stdout_bytes = attempt1_stdout_bytes,
  attempt1_stderr_bytes = attempt1_stderr_bytes,
  attempt1_timestamp_observed_seconds = attempt1_observed_seconds,
  attempt1_outer_RSS_conservative_upper_bound_bytes = contract$concurrent_child_RSS_cap_bytes,
  owner_authorization_recorded = TRUE,
  manual_replay_exception_children = length(wave_orders),
  automatic_retries = 0L,
  stringsAsFactors = FALSE
)
w(admitted_inventory, file.path(private, "mv17g-controller-admitted-prefix-artifacts.csv"))
w(wave_inventory, file.path(private, "mv17g-controller-interrupted-wave-artifacts.csv"))
w(state, file.path(private, "mv17g-controller-recovery-state.csv"))
private_names <- c(
  "mv17g-controller-admitted-prefix-artifacts.csv",
  "mv17g-controller-interrupted-wave-artifacts.csv",
  "mv17g-controller-recovery-state.csv"
)
private_paths <- file.path(private, private_names)
private_binding <- data.frame(
  contract_id = "mv17g_controller_recovery_private_binding_v2",
  role = c("resource_admitted_prefix_artifacts", "interrupted_wave_artifacts", "interruption_state"),
  artifact = private_names,
  rows = c(nrow(admitted_inventory), nrow(wave_inventory), nrow(state)),
  bytes = as.numeric(file.info(private_paths)$size),
  sha256 = vapply(private_paths, h, character(1L)),
  tracking_state = "private_not_tracked",
  stringsAsFactors = FALSE
)

implementation_files <- c(
  "R/mv17g_parallel_recovery.R",
  "scripts/run_mv17g_calibration_parallel_recovery.R",
  "scripts/run_mv17g_calibration_group_worker.R",
  "scripts/build_mv17g_controller_interruption_recovery_prefreeze.R",
  "scripts/quarantine_mv17g_interrupted_wave.R",
  "scripts/build_mv17g_parallel_calibration_closure.R",
  "tests/testthat/test-mv17g-controller-interruption-recovery.R"
)
implementation <- data.frame(
  contract_id = "mv17g_controller_recovery_implementation_binding_v2",
  file = implementation_files,
  bytes = as.numeric(file.info(implementation_files)$size),
  sha256 = vapply(implementation_files, h, character(1L)),
  stringsAsFactors = FALSE
)
if (anyNA(implementation$bytes)) stop("MV17-G controller-recovery implementation missing", call. = FALSE)

source_paths <- c(
  file.path(original_public, "mv17g-artifact-manifest.csv"),
  file.path(parallel_public, "mv17g-parallel-artifact-manifest.csv"),
  file.path(original_private, original_private_binding$artifact),
  file.path(parallel_private, parallel_private_binding$artifact),
  attempt1_time, attempt1_stdout, attempt1_stderr,
  file.path(primary_private, "mv17g-private-progress.csv")
)
source_binding <- data.frame(
  contract_id = "mv17g_controller_recovery_source_binding_v2",
  role = c(
    "original_public_prefreeze_manifest", "parallel_public_prefreeze_manifest",
    paste0("original_private_", original_private_binding$role),
    paste0("parallel_private_", parallel_private_binding$role),
    "attempt1_GNU_time_incomplete", "attempt1_stdout", "attempt1_stderr", "durable_progress_before_recovery"
  ),
  bytes = as.numeric(file.info(source_paths)$size),
  sha256 = vapply(source_paths, h, character(1L)),
  stringsAsFactors = FALSE
)

recovery_contract <- data.frame(
  contract_id = "mv17g_controller_interruption_recovery_v2",
  implementation_head = implementation_head,
  artifact_present_prefix_children = artifact_prefix,
  resource_admitted_prefix_children = admitted_prefix,
  interrupted_wave_first_order = min(wave_orders),
  interrupted_wave_last_order = max(wave_orders),
  manual_replay_exception_children = length(wave_orders),
  runner_restart_job_order = min(wave_orders),
  resume_after_replayed_wave_job_order = max(wave_orders) + 1L,
  workers = contract$workers,
  threads_per_child = contract$threads_per_child,
  automatic_retries = contract$retries,
  child_timeout_seconds = contract$child_timeout_seconds,
  child_RSS_cap_bytes = contract$child_RSS_cap_bytes,
  concurrent_child_RSS_cap_bytes = contract$concurrent_child_RSS_cap_bytes,
  interrupted_wave_quarantine_required = TRUE,
  preserve_all_interrupted_artifacts = TRUE,
  replay_payload_equality_required = TRUE,
  owner_authorization_recorded = TRUE,
  launch_authorized_after_prefreeze_commit = TRUE,
  real_localization_authorized = FALSE,
  downstream_surfaces = "closed",
  stringsAsFactors = FALSE
)
validation <- data.frame(
  contract_id = "mv17g_controller_interruption_recovery_validation_v2",
  check_id = c(
    "original_prefreeze_bound", "parallel_prefreeze_bound", "private_prefreezes_bound",
    "implementation_head_bound", "queue_1188_children", "queue_91740_runs",
    "artifact_prefix_794", "durable_progress_786", "resource_admitted_prefix_786",
    "interrupted_wave_8", "wave_valid_receipts_3", "wave_incomplete_receipts_5",
    "valid_receipt_orders_exact", "incomplete_receipt_orders_exact",
    "runner_restart_787", "resume_after_wave_795", "artifact_prefix_consecutive",
    "no_partial_evidence", "admitted_artifacts_complete", "admitted_payloads_valid",
    "admitted_resources_pass", "wave_artifacts_complete", "wave_payloads_valid",
    "all_child_streams_empty", "all_264_matrices_bound", "controller_absent",
    "attempt1_GNU_time_incomplete", "attempt1_stderr_empty", "attempt1_stdout_preserved",
    "timestamp_observation_bounded", "eight_workers_one_thread", "zero_automatic_retries",
    "manual_replay_exception_8", "owner_authorization_recorded", "quarantine_required",
    "payload_equality_required", "execution_only_after_commit", "real_localization_closed",
    "downstream_firewall", "aggregate_only_public"
  ),
  passed = c(
    nrow(original_manifest) >= 1L, nrow(parallel_manifest) >= 1L,
    nrow(original_private_binding) == 4L && nrow(parallel_private_binding) == 3L,
    grepl("^[0-9a-f]{40}$", implementation_head), nrow(queue) == 1188L,
    sum(queue$scientific_runs) == 91740L, artifact_prefix == 794L, admitted_prefix == 786L,
    admitted_prefix == nrow(progress), length(wave_orders) == 8L,
    length(valid_wave_orders) == 3L, length(incomplete_wave_orders) == 5L,
    identical(valid_wave_orders, c(789L, 791L, 794L)),
    identical(incomplete_wave_orders, c(787L, 788L, 790L, 792L, 793L)),
    min(wave_orders) == 787L, max(wave_orders) + 1L == 795L,
    artifact_prefix == mv17g_complete_prefix_v1(scan), length(stray_partials) == 0L,
    nrow(admitted_inventory) == admitted_prefix * 4L, TRUE, TRUE,
    nrow(wave_inventory) == length(wave_orders) * 4L, TRUE,
    all(file.info(c(admitted_inventory$artifact, wave_inventory$artifact)[
      c(admitted_inventory$role, wave_inventory$role) %in% c("stdout", "stderr")])$size == 0),
    nrow(matrix_catalog) == 264L, controller_processes == 0L, attempt1_time_bytes == 0,
    attempt1_stderr_bytes == 0, attempt1_stdout_bytes > 0,
    attempt1_observed_seconds > 0 && attempt1_observed_seconds <= contract$aggregate_timeout_seconds,
    contract$workers == 8L && contract$threads_per_child == 1L, contract$retries == 0L,
    recovery_contract$manual_replay_exception_children == 8L,
    recovery_contract$owner_authorization_recorded,
    recovery_contract$interrupted_wave_quarantine_required && recovery_contract$preserve_all_interrupted_artifacts,
    recovery_contract$replay_payload_equality_required,
    recovery_contract$launch_authorized_after_prefreeze_commit,
    !recovery_contract$real_localization_authorized,
    recovery_contract$downstream_surfaces == "closed",
    !any(c("unit_id", "unit_order", "source_path", "artifact") %in% names(state))
  ),
  stringsAsFactors = FALSE
)
if (!all(validation$passed)) stop("MV17-G controller-interruption recovery prefreeze failed", call. = FALSE)

items <- list(
  "mv17g-controller-recovery-contract.csv" = recovery_contract,
  "mv17g-controller-recovery-state-summary.csv" = state,
  "mv17g-controller-recovery-private-binding.csv" = private_binding,
  "mv17g-controller-recovery-source-binding.csv" = source_binding,
  "mv17g-controller-recovery-implementation-binding.csv" = implementation,
  "mv17g-controller-recovery-validation.csv" = validation,
  "mv17g-controller-recovery-decision.csv" = data.frame(
    contract_id = "mv17g_controller_interruption_recovery_decision_v2",
    disposition = "preserve_and_replay_interrupted_wave_787_794_then_continue",
    resource_admitted_prefix_children = admitted_prefix,
    manual_replay_exception_children = length(wave_orders),
    automatic_retries = 0L,
    owner_authorization_recorded = TRUE,
    scientific_estimand_changed = FALSE,
    attempt1_outer_telemetry_complete = FALSE,
    real_localization_authorized = FALSE,
    downstream_surfaces = "closed",
    stringsAsFactors = FALSE
  )
)
for (name in names(items)) w(items[[name]], file.path(public, name))
writeLines(
  c(
    "# MV17-G controller-interruption recovery prefreeze",
    "",
    "The first parallel controller ended without a GNU-time footer during the wave spanning orders 787--794.",
    "Orders 1--786 are the resource-admitted consecutive prefix. All eight interrupted-wave payloads are finite, but five lack complete child GNU-time receipts and therefore are not admitted.",
    "The owner authorized preservation and deterministic replay of the entire interrupted wave. After this prefreeze commits, all 32 wave artifacts must move intact to private quarantine before the unchanged runner restarts at order 787.",
    "Replayed payloads must equal the preserved payloads; automatic retries remain zero and the scientific estimand is unchanged.",
    "The missing controller-level attempt-1 telemetry remains explicit; child receipts and a conservative concurrent-RSS upper bound are retained.",
    "Real localization, labels, outcomes, clustering, fusion, biology, and manuscript claims remain closed."
  ),
  file.path(public, "MV17G_CONTROLLER_INTERRUPTION_RECOVERY_PREFREEZE_2026-08-27.md")
)
files <- sort(list.files(public))
w(
  data.frame(
    contract_id = "mv17g_controller_interruption_recovery_manifest_v2",
    artifact = files,
    bytes = as.numeric(file.info(file.path(public, files))$size),
    sha256 = vapply(file.path(public, files), h, character(1L)),
    stringsAsFactors = FALSE
  ),
  file.path(public, "mv17g-controller-recovery-artifact-manifest.csv")
)
message(
  "Built MV17-G controller-interruption recovery prefreeze; checks=", nrow(validation),
  "; admitted_prefix=", admitted_prefix, "; replay=", length(wave_orders)
)
