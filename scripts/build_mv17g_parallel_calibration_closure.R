#!/usr/bin/env Rscript
options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 17L) {
  stop(
    paste(
      "usage: build_mv17g_parallel_calibration_closure.R",
      "<original-public-prefreeze> <original-private-prefreeze>",
      "<parallel-public-prefreeze> <parallel-private-prefreeze>",
      "<controller-recovery-public-prefreeze> <controller-recovery-private-prefreeze>",
      "<interrupted-wave-quarantine>",
      "<primary-private> <primary-public> <repeat-private> <repeat-public>",
      "<serial-primary-time> <parallel-primary-attempt1-time>",
      "<parallel-primary-recovery-time> <parallel-repeat-time>",
      "<output> <execution-head>"
    ),
    call. = FALSE
  )
}
original_public <- normalizePath(args[[1]], mustWork = TRUE)
original_private <- normalizePath(args[[2]], mustWork = TRUE)
parallel_public <- normalizePath(args[[3]], mustWork = TRUE)
parallel_private <- normalizePath(args[[4]], mustWork = TRUE)
recovery_public <- normalizePath(args[[5]], mustWork = TRUE)
recovery_private <- normalizePath(args[[6]], mustWork = TRUE)
quarantine <- normalizePath(args[[7]], mustWork = TRUE)
primary_private <- normalizePath(args[[8]], mustWork = TRUE)
primary_public <- normalizePath(args[[9]], mustWork = TRUE)
repeat_private <- normalizePath(args[[10]], mustWork = TRUE)
repeat_public <- normalizePath(args[[11]], mustWork = TRUE)
serial_time <- normalizePath(args[[12]], mustWork = TRUE)
primary_attempt1_time <- normalizePath(args[[13]], mustWork = TRUE)
primary_recovery_time <- normalizePath(args[[14]], mustWork = TRUE)
repeat_time <- normalizePath(args[[15]], mustWork = TRUE)
output <- args[[16]]
execution_head <- tolower(args[[17]])
if (dir.exists(output) || !grepl("^[0-9a-f]{40}$", execution_head)) {
  stop("invalid MV17-G parallel closure target/head", call. = FALSE)
}
dir.create(output, recursive = TRUE)

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
    stop("MV17-G parallel closure manifest drift", call. = FALSE)
  }
  manifest
}
original_manifest <- verify_manifest(original_public, "mv17g-artifact-manifest.csv")
parallel_manifest <- verify_manifest(parallel_public, "mv17g-parallel-artifact-manifest.csv")
recovery_manifest <- verify_manifest(recovery_public, "mv17g-controller-recovery-artifact-manifest.csv")
primary_manifest <- verify_manifest(primary_public, "mv17g-artifact-manifest.csv")
repeat_manifest <- verify_manifest(repeat_public, "mv17g-artifact-manifest.csv")

original_private_binding <- r(file.path(original_public, "mv17g-private-binding.csv"))
original_private_paths <- file.path(original_private, original_private_binding$artifact)
parallel_private_binding <- r(file.path(parallel_public, "mv17g-parallel-private-binding.csv"))
parallel_private_paths <- file.path(parallel_private, parallel_private_binding$artifact)
recovery_private_binding <- r(file.path(recovery_public, "mv17g-controller-recovery-private-binding.csv"))
recovery_private_paths <- file.path(recovery_private, recovery_private_binding$artifact)
for (pair in list(
  list(paths = original_private_paths, binding = original_private_binding),
  list(paths = parallel_private_paths, binding = parallel_private_binding),
  list(paths = recovery_private_paths, binding = recovery_private_binding)
)) {
  if (!all(file.exists(pair$paths)) ||
      !identical(unname(as.numeric(file.info(pair$paths)$size)), unname(as.numeric(pair$binding$bytes))) ||
      !identical(unname(vapply(pair$paths, h, character(1L))), unname(tolower(pair$binding$sha256)))) {
    stop("MV17-G parallel closure private prefreeze drift", call. = FALSE)
  }
}

contract <- r(file.path(parallel_public, "mv17g-parallel-contract.csv"))
recovery_contract <- r(file.path(recovery_public, "mv17g-controller-recovery-contract.csv"))
recovery_state <- r(file.path(recovery_public, "mv17g-controller-recovery-state-summary.csv"))
recovery_decision <- r(file.path(recovery_public, "mv17g-controller-recovery-decision.csv"))
recovery_wave <- r(file.path(recovery_private, "mv17g-controller-interrupted-wave-artifacts.csv"))
quarantine_manifest <- r(file.path(quarantine, "mv17g-interrupted-wave-quarantine.csv"))
quarantine_paths <- quarantine_manifest$quarantine_artifact
if (nrow(quarantine_manifest) != 32L || !all(file.exists(quarantine_paths)) ||
    !identical(unname(as.numeric(file.info(quarantine_paths)$size)), unname(as.numeric(quarantine_manifest$bytes))) ||
    !identical(unname(vapply(quarantine_paths, h, character(1L))), unname(tolower(quarantine_manifest$sha256)))) {
  stop("MV17-G interrupted-wave quarantine drift", call. = FALSE)
}
primary_queue <- r(file.path(original_private, "mv17g-primary-grouped-queue.csv"))
repeat_queue <- r(file.path(original_private, "mv17g-repeat-grouped-queue.csv"))
prefix_status <- r(file.path(parallel_private, "mv17g-parallel-prefix-status.csv"))
prefix_inventory <- r(file.path(parallel_private, "mv17g-parallel-prefix-artifacts.csv"))
serial_prefix <- as.integer(prefix_status$serial_prefix_children)

primary_metrics <- r(file.path(primary_private, "mv17g-private-scientific-metrics.csv"))
repeat_metrics <- r(file.path(repeat_private, "mv17g-private-scientific-metrics.csv"))
primary_ledger <- r(file.path(primary_private, "mv17g-private-resource-ledger.csv"))
repeat_ledger <- r(file.path(repeat_private, "mv17g-private-resource-ledger.csv"))
primary_empirical <- r(file.path(primary_private, "mv17g-private-empirical-calibration.csv"))
aggregate_empirical <- r(file.path(primary_public, "mv17g-aggregate-empirical-calibration.csv"))
primary_resource <- r(file.path(primary_public, "mv17g-aggregate-resource.csv"))
repeat_resource <- r(file.path(repeat_public, "mv17g-aggregate-resource.csv"))
primary_status <- r(file.path(primary_public, "mv17g-status.csv"))
repeat_status <- r(file.path(repeat_public, "mv17g-status.csv"))

repeat_units <- unique(repeat_queue[c("view", "unit_order")])
take <- paste(primary_metrics$view, primary_metrics$unit_order) %in% paste(repeat_units$view, repeat_units$unit_order)
metric_keys <- c(
  "view", "unit_order", "null_family", "replicate", "seed",
  "homology_dimension", "summary_id", "value", "h0_mst_maximum_absolute_error"
)
primary_repeat <- primary_metrics[take, metric_keys, drop = FALSE]
repeat_exact_table <- repeat_metrics[metric_keys]
sort_keys <- setdiff(metric_keys, c("value", "h0_mst_maximum_absolute_error"))
primary_repeat <- primary_repeat[do.call(order, c(primary_repeat[sort_keys], list(method = "radix"))), , drop = FALSE]
repeat_exact_table <- repeat_exact_table[do.call(order, c(repeat_exact_table[sort_keys], list(method = "radix"))), , drop = FALSE]
rownames(primary_repeat) <- rownames(repeat_exact_table) <- NULL
repeat_exact <- identical(primary_repeat, repeat_exact_table)

independent_empirical <- mv17g_empirical_table_v1(primary_metrics)
empirical_exact <- isTRUE(all.equal(primary_empirical, independent_empirical, tolerance = 1e-12, check.attributes = FALSE))
independent_aggregate <- mv17g_aggregate_empirical_v1(independent_empirical)
aggregate_exact <- isTRUE(all.equal(aggregate_empirical, independent_aggregate, tolerance = 1e-12, check.attributes = FALSE))
resource_rebuild <- function(ledger) {
  keys <- unique(ledger[c("view", "null_family")])
  do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
    k <- keys[i, , drop = FALSE]
    take <- ledger$view == k$view & ledger$null_family == k$null_family
    data.frame(
      k,
      grouped_children = sum(take),
      scientific_runs = sum(ledger$replicate_count[take]),
      aggregate_child_seconds = sum(ledger$wall_seconds[take]),
      maximum_child_RSS_bytes = max(ledger$maximum_RSS_bytes[take]),
      aggregate_output_bytes = sum(ledger$output_bytes[take]),
      adopted_children = sum(ledger$adopted[take]),
      stringsAsFactors = FALSE
    )
  }))
}
resource_exact <-
  isTRUE(all.equal(primary_resource, resource_rebuild(primary_ledger), tolerance = 1e-12, check.attributes = FALSE)) &&
  isTRUE(all.equal(repeat_resource, resource_rebuild(repeat_ledger), tolerance = 1e-12, check.attributes = FALSE))

serial_outer <- mv17c_parse_gnu_time_v1(serial_time)
if (file.info(primary_attempt1_time)$size != 0) {
  stop("MV17-G interrupted attempt unexpectedly has a GNU-time footer", call. = FALSE)
}
primary_recovery_outer <- mv17c_parse_gnu_time_v1(primary_recovery_time)
repeat_outer <- mv17c_parse_gnu_time_v1(repeat_time)
serial_text <- paste(readLines(serial_time, warn = FALSE), collapse = "\n")
controlled_serial_stop <- grepl("signal 9", serial_text, ignore.case = TRUE) || serial_outer$exit_status %in% c(9L, 137L)

prefix_hashes_exact <- all(vapply(seq_len(serial_prefix), function(i) {
  q <- primary_queue[i, , drop = FALSE]
  paths <- mv17g_job_artifacts_v1(q, primary_private)
  expected <- prefix_inventory[prefix_inventory$job_order == q$job_order, , drop = FALSE]
  expected <- expected[match(names(paths), expected$role), , drop = FALSE]
  identical(unname(as.numeric(file.info(paths)$size)), unname(as.numeric(expected$bytes))) &&
    identical(unname(vapply(paths, h, character(1L))), unname(tolower(expected$sha256)))
}, logical(1L)))
prefix_origin_exact <-
  all(primary_ledger$adopted[seq_len(serial_prefix)]) &&
  all(primary_ledger$execution_origin[seq_len(serial_prefix)] == "serial_prefix_v1") &&
  !any(primary_ledger$adopted[-seq_len(serial_prefix)]) &&
  all(primary_ledger$execution_origin[-seq_len(serial_prefix)] == "parallel_recovery_v1") &&
  !any(repeat_ledger$adopted) && all(repeat_ledger$execution_origin == "parallel_recovery_v1")
replay_result_rows <- recovery_wave[recovery_wave$role == "result", , drop = FALSE]
replay_result_rows <- replay_result_rows[order(replay_result_rows$job_order), , drop = FALSE]
replayed_result_paths <- vapply(replay_result_rows$job_order, function(i) {
  mv17g_job_artifacts_v1(primary_queue[i, , drop = FALSE], primary_private)[["result"]]
}, character(1L))
replay_payloads_exact <-
  identical(as.integer(replay_result_rows$job_order), 787:794) &&
  all(file.exists(replayed_result_paths)) &&
  identical(unname(vapply(replayed_result_paths, h, character(1L))), unname(tolower(replay_result_rows$sha256)))

primary_logs <- list.files(file.path(primary_private, "logs"), pattern = "[.](stdout|stderr)[.]txt$", full.names = TRUE)
repeat_logs <- list.files(file.path(repeat_private, "logs"), pattern = "[.](stdout|stderr)[.]txt$", full.names = TRUE)
all_logs <- c(primary_logs, repeat_logs)
serial_child_seconds <- sum(primary_ledger$wall_seconds[primary_ledger$adopted])
parallel_primary_child_seconds <- sum(primary_ledger$wall_seconds[!primary_ledger$adopted])
parallel_repeat_child_seconds <- sum(repeat_ledger$wall_seconds)
primary_observed_outer_seconds <- serial_outer$wall_seconds +
  recovery_state$attempt1_timestamp_observed_seconds + primary_recovery_outer$wall_seconds

source_paths <- c(
  file.path(original_public, "mv17g-artifact-manifest.csv"),
  file.path(parallel_public, "mv17g-parallel-artifact-manifest.csv"),
  file.path(recovery_public, "mv17g-controller-recovery-artifact-manifest.csv"),
  file.path(quarantine, "mv17g-interrupted-wave-quarantine.csv"),
  original_private_paths,
  parallel_private_paths,
  recovery_private_paths,
  file.path(primary_public, "mv17g-artifact-manifest.csv"),
  file.path(repeat_public, "mv17g-artifact-manifest.csv"),
  file.path(primary_private, "mv17g-private-scientific-metrics.csv"),
  file.path(repeat_private, "mv17g-private-scientific-metrics.csv"),
  file.path(primary_private, "mv17g-private-empirical-calibration.csv"),
  serial_time, primary_attempt1_time, primary_recovery_time, repeat_time
)
source_binding <- data.frame(
  contract_id = "mv17g_parallel_closure_source_binding_v1",
  role = c(
    "original_public_prefreeze_manifest", "parallel_public_prefreeze_manifest",
    "controller_recovery_public_prefreeze_manifest", "interrupted_wave_quarantine_manifest",
    paste0("original_private_", original_private_binding$role),
    paste0("parallel_private_", parallel_private_binding$role),
    paste0("controller_recovery_private_", recovery_private_binding$role),
    "primary_public_manifest", "repeat_public_manifest", "primary_metrics", "repeat_metrics",
    "primary_empirical", "serial_primary_GNU_time", "parallel_primary_attempt1_GNU_time_incomplete",
    "parallel_primary_recovery_GNU_time", "parallel_repeat_GNU_time"
  ),
  bytes = as.numeric(file.info(source_paths)$size),
  sha256 = vapply(source_paths, h, character(1L)),
  stringsAsFactors = FALSE
)

validation <- data.frame(
  contract_id = "mv17g_parallel_closure_validation_v1",
  check_id = c(
    "original_prefreeze_bound", "parallel_prefreeze_bound", "controller_recovery_prefreeze_bound",
    "private_prefreezes_bound", "controller_recovery_private_bound",
    "controller_recovery_artifact_prefix_794", "controller_recovery_admitted_prefix_786",
    "controller_recovery_replay_8", "interrupted_wave_quarantine_bound", "replayed_payloads_exact",
    "controller_recovery_attempt1_telemetry_explicit",
    "execution_head_recorded", "primary_1188_children", "repeat_27_children",
    "primary_91740_runs", "repeat_2085_runs", "primary_733920_metrics", "repeat_16680_metrics",
    "four_H0_H1_summaries", "primary_empirical_7392", "empirical_independent",
    "aggregate_56_rows", "aggregate_independent", "repeat_exact", "resource_aggregate_independent",
    "serial_prefix_hashes_exact", "serial_prefix_origin_exact", "controlled_serial_stop",
    "child_resources_pass", "parallel_outer_resources_pass", "primary_observed_aggregate_pass",
    "serial_outer_covers_prefix",
    "parallel_capacity_covers_primary", "parallel_capacity_covers_repeat", "parallel_outer_bounds_max_child",
    "empty_child_streams", "eight_workers_one_thread_zero_retry", "storage_caps",
    "source_receipts_bound", "real_localization_closed", "downstream_firewall", "aggregate_only_public"
  ),
  passed = c(
    nrow(original_manifest) >= 1L, nrow(parallel_manifest) >= 1L,
    nrow(recovery_manifest) >= 1L,
    nrow(original_private_binding) == 4L && nrow(parallel_private_binding) == 3L,
    nrow(recovery_private_binding) == 3L,
    recovery_contract$artifact_present_prefix_children == 794L &&
      recovery_state$artifact_present_prefix_children == 794L,
    recovery_contract$resource_admitted_prefix_children == 786L &&
      recovery_state$resource_admitted_prefix_children == 786L,
    recovery_contract$manual_replay_exception_children == 8L &&
      recovery_state$manual_replay_exception_children == 8L,
    nrow(quarantine_manifest) == 32L && all(file.exists(quarantine_paths)), replay_payloads_exact,
    file.info(primary_attempt1_time)$size == 0 &&
      !recovery_decision$attempt1_outer_telemetry_complete &&
      recovery_state$attempt1_outer_RSS_conservative_upper_bound_bytes == contract$concurrent_child_RSS_cap_bytes,
    grepl("^[0-9a-f]{40}$", execution_head), nrow(primary_queue) == 1188L, nrow(repeat_queue) == 27L,
    sum(primary_queue$scientific_runs) == 91740L, sum(repeat_queue$scientific_runs) == 2085L,
    nrow(primary_metrics) == 733920L, nrow(repeat_metrics) == 16680L,
    setequal(primary_metrics$homology_dimension, c("H0", "H1")) && setequal(primary_metrics$summary_id, mv17c_summary_registry_v1()$summary_id),
    nrow(primary_empirical) == 7392L, empirical_exact, nrow(aggregate_empirical) == 56L, aggregate_exact,
    repeat_exact, resource_exact, prefix_hashes_exact, prefix_origin_exact, controlled_serial_stop,
    max(c(primary_ledger$maximum_RSS_bytes, repeat_ledger$maximum_RSS_bytes)) <= contract$child_RSS_cap_bytes &&
      max(c(primary_ledger$wall_seconds, repeat_ledger$wall_seconds)) <= contract$child_timeout_seconds,
    primary_recovery_outer$exit_status == 0L && repeat_outer$exit_status == 0L &&
      repeat_outer$wall_seconds <= contract$aggregate_timeout_seconds,
    primary_observed_outer_seconds <= contract$aggregate_timeout_seconds,
    serial_outer$wall_seconds >= serial_child_seconds,
    (recovery_state$attempt1_timestamp_observed_seconds + primary_recovery_outer$wall_seconds) *
      contract$workers >= parallel_primary_child_seconds,
    repeat_outer$wall_seconds * contract$workers >= parallel_repeat_child_seconds,
    max(recovery_state$attempt1_timestamp_observed_seconds, primary_recovery_outer$wall_seconds) >=
      max(primary_ledger$wall_seconds[!primary_ledger$adopted]) &&
      repeat_outer$wall_seconds >= max(repeat_ledger$wall_seconds),
    length(all_logs) == 2L * (1188L + 27L) && all(file.info(all_logs)$size == 0),
    primary_status$workers == 8L && repeat_status$workers == 8L &&
      primary_status$threads_per_child == 1L && repeat_status$threads_per_child == 1L &&
      primary_status$retries == 0L && repeat_status$retries == 0L,
    sum(as.numeric(file.info(list.files(primary_private, recursive = TRUE, full.names = TRUE))$size)) <= contract$private_cap_bytes &&
      sum(as.numeric(file.info(list.files(repeat_private, recursive = TRUE, full.names = TRUE))$size)) <= contract$private_cap_bytes &&
      sum(as.numeric(file.info(list.files(primary_public, full.names = TRUE))$size)) <= contract$public_cap_bytes &&
      sum(as.numeric(file.info(list.files(repeat_public, full.names = TRUE))$size)) <= contract$public_cap_bytes,
    nrow(source_binding) == length(source_paths), !contract$real_localization_authorized,
    !primary_status$labels_opened && !primary_status$outcomes_opened &&
      primary_status$clustering_jobs == 0L && primary_status$fusion_jobs == 0L &&
      primary_status$biology_jobs == 0L && primary_status$manuscript_claim_jobs == 0L,
    !any(c("unit_id", "unit_order", "observed_value") %in% names(aggregate_empirical))
  )
)
if (!all(validation$passed)) stop("MV17-G parallel closure failed", call. = FALSE)

resource_summary <- data.frame(
  contract_id = "mv17g_parallel_closure_resource_v1",
  mode = c("primary", "repeat"),
  grouped_children = c(nrow(primary_queue), nrow(repeat_queue)),
  scientific_runs = c(sum(primary_queue$scientific_runs), sum(repeat_queue$scientific_runs)),
  serial_prefix_children = c(serial_prefix, 0L),
  workers = 8L,
  aggregate_child_seconds = c(sum(primary_ledger$wall_seconds), sum(repeat_ledger$wall_seconds)),
  outer_wall_seconds = c(primary_observed_outer_seconds, repeat_outer$wall_seconds),
  interrupted_attempt_timestamp_observed_seconds = c(recovery_state$attempt1_timestamp_observed_seconds, 0),
  interrupted_attempt_GNU_time_complete = c(FALSE, TRUE),
  maximum_child_RSS_bytes = c(max(primary_ledger$maximum_RSS_bytes), max(repeat_ledger$maximum_RSS_bytes)),
  parallel_outer_maximum_RSS_bytes = c(
    max(recovery_state$attempt1_outer_RSS_conservative_upper_bound_bytes, primary_recovery_outer$maximum_RSS_bytes),
    repeat_outer$maximum_RSS_bytes
  ),
  private_bytes = c(primary_status$private_bytes, repeat_status$private_bytes),
  public_bytes = c(
    sum(as.numeric(file.info(list.files(primary_public, full.names = TRUE))$size)),
    sum(as.numeric(file.info(list.files(repeat_public, full.names = TRUE))$size))
  ),
  stringsAsFactors = FALSE
)
items <- list(
  "mv17g-aggregate-empirical-calibration.csv" = aggregate_empirical,
  "mv17g-resource-summary.csv" = resource_summary,
  "mv17g-source-binding.csv" = source_binding,
  "mv17g-validation.csv" = validation,
  "mv17g-decision.csv" = data.frame(
    contract_id = "mv17g_decision_v1",
    full_H0_H1_calibration_closed = TRUE,
    eight_worker_recovery_closed = TRUE,
    interrupted_controller_artifact_prefix_children = recovery_state$artifact_present_prefix_children,
    interrupted_controller_admitted_prefix_children = recovery_state$resource_admitted_prefix_children,
    interrupted_wave_replay_children = recovery_state$manual_replay_exception_children,
    interrupted_wave_replay_payloads_exact = replay_payloads_exact,
    interrupted_controller_telemetry_complete = FALSE,
    serial_prefix_children = serial_prefix,
    real_localization_prefreeze_eligible = TRUE,
    real_localization_authorized = FALSE,
    labels_outcomes_clustering_fusion_claims = "closed",
    execution_head = execution_head,
    stringsAsFactors = FALSE
  )
)
for (name in names(items)) w(items[[name]], file.path(output, name))
writeLines(
  c(
    "# MV17-G full H0/H1 calibration closure",
    "",
    "All 132 cell and 132 gene units are calibrated separately with 99 fixed replicates per compatible null family.",
    paste0("An immutable serial prefix of ", serial_prefix, " children was adopted exactly; all remaining children and the repeat used eight single-threaded workers."),
    "The parallel primary required two controller segments. The first ended with 786 resource-admitted children and an eight-child interrupted wave; all 32 wave artifacts were preserved, the wave was replayed, and all eight replayed payloads matched exactly.",
    "The missing first-controller GNU-time footer, timestamp-observed duration, and conservative 64-GiB RSS upper bound remain explicit.",
    "Public evidence is aggregate-only. Real localization and every downstream scientific surface remain closed."
  ),
  file.path(output, "MV17G_FULL_CALIBRATION_CLOSURE_2026-08-26.md")
)
files <- sort(list.files(output))
w(
  data.frame(
    contract_id = "mv17g_closure_manifest_v1",
    artifact = files,
    bytes = as.numeric(file.info(file.path(output, files))$size),
    sha256 = vapply(file.path(output, files), h, character(1L)),
    stringsAsFactors = FALSE
  ),
  file.path(output, "mv17g-artifact-manifest.csv")
)
message("Built MV17-G parallel full-calibration closure; checks=", nrow(validation))
