#!/usr/bin/env Rscript
options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8L) {
  stop(
    paste(
      "usage: run_mv17g_calibration_parallel_recovery.R",
      "<original-public-prefreeze> <original-private-prefreeze>",
      "<parallel-public-prefreeze> <parallel-private-prefreeze>",
      "<primary-or-repeat> <matrix-root> <private-output> <public-output>"
    ),
    call. = FALSE
  )
}
original_public <- normalizePath(args[[1]], mustWork = TRUE)
original_private <- normalizePath(args[[2]], mustWork = TRUE)
parallel_public <- normalizePath(args[[3]], mustWork = TRUE)
parallel_private <- normalizePath(args[[4]], mustWork = TRUE)
mode <- match.arg(args[[5]], c("primary", "repeat"))
matrix_root <- normalizePath(args[[6]], mustWork = TRUE)
private_arg <- args[[7]]
if (!dir.exists(private_arg)) dir.create(private_arg, recursive = TRUE)
private <- normalizePath(private_arg, mustWork = TRUE)
public <- args[[8]]
public_partial <- paste0(public, ".partial")
if (dir.exists(public) || dir.exists(public_partial)) stop("MV17-G parallel public output exists", call. = FALSE)

source("R/mv08z_landscape_production.R")
source("R/mv17_null_models.R")
source("R/mv17_calibration.R")
source("R/mv17_localization.R")
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
    stop("MV17-G parallel prefreeze manifest drift", call. = FALSE)
  }
  manifest
}
verify_manifest(original_public, "mv17g-artifact-manifest.csv")
verify_manifest(parallel_public, "mv17g-parallel-artifact-manifest.csv")

original_private_binding <- r(file.path(original_public, "mv17g-private-binding.csv"))
original_private_paths <- file.path(original_private, original_private_binding$artifact)
if (!all(file.exists(original_private_paths)) ||
    !identical(unname(as.numeric(file.info(original_private_paths)$size)), unname(as.numeric(original_private_binding$bytes))) ||
    !identical(unname(vapply(original_private_paths, h, character(1L))), unname(tolower(original_private_binding$sha256)))) {
  stop("MV17-G original private prefreeze drift", call. = FALSE)
}
parallel_private_binding <- r(file.path(parallel_public, "mv17g-parallel-private-binding.csv"))
parallel_private_paths <- file.path(parallel_private, parallel_private_binding$artifact)
if (!all(file.exists(parallel_private_paths)) ||
    !identical(unname(as.numeric(file.info(parallel_private_paths)$size)), unname(as.numeric(parallel_private_binding$bytes))) ||
    !identical(unname(vapply(parallel_private_paths, h, character(1L))), unname(tolower(parallel_private_binding$sha256)))) {
  stop("MV17-G parallel private prefreeze drift", call. = FALSE)
}

decision <- r(file.path(parallel_public, "mv17g-parallel-decision.csv"))
contract <- r(file.path(parallel_public, "mv17g-parallel-contract.csv"))
if (!decision$eight_worker_primary_and_repeat_authorized_after_commit ||
    contract$workers != 8L || contract$threads_per_child != 1L || contract$retries != 0L) {
  stop("MV17-G parallel execution not authorized", call. = FALSE)
}

queue_name <- if (mode == "primary") "mv17g-primary-grouped-queue.csv" else "mv17g-repeat-grouped-queue.csv"
queue <- r(file.path(original_private, queue_name))
expected_jobs <- if (mode == "primary") 1188L else 27L
expected_runs <- if (mode == "primary") 91740L else 2085L
if (nrow(queue) != expected_jobs || sum(queue$scientific_runs) != expected_runs) {
  stop("MV17-G parallel queue cardinality drift", call. = FALSE)
}

matrix_catalog <- r(file.path(parallel_private, "mv17g-parallel-matrix-catalog.csv"))
matrix_paths <- file.path(matrix_root, "matrices", sprintf("%s__%03d.rds", matrix_catalog$view, matrix_catalog$unit_order))
if (nrow(matrix_catalog) != 264L || !all(file.exists(matrix_paths)) ||
    !identical(unname(as.numeric(file.info(matrix_paths)$size)), unname(as.numeric(matrix_catalog$bytes))) ||
    !identical(unname(vapply(matrix_paths, h, character(1L))), unname(tolower(matrix_catalog$sha256)))) {
  stop("MV17-G parallel matrix binding drift", call. = FALSE)
}

dir.create(file.path(private, "jobs"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(private, "logs"), showWarnings = FALSE, recursive = TRUE)
worker <- normalizePath("scripts/run_mv17g_calibration_group_worker.R", mustWork = TRUE)
prefix_status <- r(file.path(parallel_private, "mv17g-parallel-prefix-status.csv"))
serial_prefix <- if (mode == "primary") as.integer(prefix_status$serial_prefix_children) else 0L

validate_result <- function(result, q, matrix_path) {
  expected_seeds <- if (q$null_family == "observed") 0L else seq.int(q$seed_first, length.out = q$replicate_count)
  ok <- identical(result$contract_id, "mv17g_group_result_v1") &&
    identical(result$null_family, q$null_family) &&
    result$seed_first == min(expected_seeds) && result$seed_last == max(expected_seeds) &&
    result$replicate_count == length(expected_seeds) && result$matrix_sha256 == h(matrix_path) &&
    nrow(result$metrics) == length(expected_seeds) * 8L &&
    setequal(unique(result$metrics$seed), expected_seeds) &&
    all(result$metrics$h0_mst_maximum_absolute_error <= 1e-8) && result$finite
  if (!isTRUE(ok)) stop("MV17-G parallel child payload drift at order ", q$job_order, call. = FALSE)
  result
}

collect_one <- function(i) {
  q <- queue[i, , drop = FALSE]
  paths <- mv17g_job_artifacts_v1(q, private)
  if (!all(file.exists(paths)) || file.info(paths[["stdout"]])$size != 0 || file.info(paths[["stderr"]])$size != 0) {
    stop("MV17-G parallel child artifact/stream drift at order ", q$job_order, call. = FALSE)
  }
  resource <- mv17c_parse_gnu_time_v1(paths[["time"]])
  if (resource$exit_status != 0L || resource$wall_seconds > contract$child_timeout_seconds ||
      resource$maximum_RSS_bytes > contract$child_RSS_cap_bytes) {
    stop("MV17-G parallel child resource cap at order ", q$job_order, call. = FALSE)
  }
  matrix_path <- file.path(matrix_root, "matrices", sprintf("%s__%03d.rds", q$view, q$unit_order))
  result <- validate_result(readRDS(paths[["result"]]), q, matrix_path)
  metrics <- result$metrics
  metrics$job_order <- q$job_order
  metrics$view <- q$view
  metrics$unit_order <- q$unit_order
  metrics$null_family <- q$null_family
  ledger <- data.frame(
    job_order = q$job_order,
    view = q$view,
    unit_order = q$unit_order,
    null_family = q$null_family,
    replicate_count = q$replicate_count,
    seed_first = q$seed_first,
    seed_last = q$seed_last,
    wall_seconds = resource$wall_seconds,
    maximum_RSS_bytes = resource$maximum_RSS_bytes,
    output_bytes = as.numeric(file.info(paths[["result"]])$size),
    output_sha256 = h(paths[["result"]]),
    adopted = mode == "primary" && q$job_order <= serial_prefix,
    execution_origin = if (mode == "primary" && q$job_order <= serial_prefix) "serial_prefix_v1" else "parallel_recovery_v1",
    stringsAsFactors = FALSE
  )
  list(
    metrics = metrics[c(
      "job_order", "view", "unit_order", "null_family", "replicate", "seed",
      "homology_dimension", "summary_id", "value", "h0_mst_maximum_absolute_error"
    )],
    ledger = ledger
  )
}

scan <- mv17g_checkpoint_scan_v1(queue, private)
stray_partials <- list.files(private, pattern = "[.]partial$", recursive = TRUE, full.names = TRUE)
if (any(scan$state == "partial") || length(stray_partials)) stop("MV17-G parallel checkpoint contains partial evidence", call. = FALSE)
complete_prefix <- mv17g_complete_prefix_v1(scan, require_incomplete = FALSE)
if (mode == "primary" && complete_prefix < serial_prefix) stop("MV17-G serial prefix missing", call. = FALSE)
if (mode == "primary") {
  prefix_inventory <- r(file.path(parallel_private, "mv17g-parallel-prefix-artifacts.csv"))
  for (i in seq_len(serial_prefix)) {
    q <- queue[i, , drop = FALSE]
    paths <- mv17g_job_artifacts_v1(q, private)
    rows <- prefix_inventory$job_order == q$job_order
    expected <- prefix_inventory[rows, , drop = FALSE]
    expected <- expected[match(names(paths), expected$role), , drop = FALSE]
    if (!identical(unname(as.numeric(file.info(paths)$size)), unname(as.numeric(expected$bytes))) ||
        !identical(unname(vapply(paths, h, character(1L))), unname(tolower(expected$sha256)))) {
      stop("MV17-G serial prefix changed after prefreeze", call. = FALSE)
    }
  }
}

pending <- scan$job_order[scan$state == "absent"]
batches <- mv17g_parallel_batches_v1(pending, contract$workers)
for (batch in batches) {
  rows <- match(batch, queue$job_order)
  wave <- mv17g_run_parallel_wave_v1(queue[rows, , drop = FALSE], private, matrix_root, worker, contract)
  if (any(wave$exit_status != 0L) || any(wave$artifacts != 4L)) {
    stop("MV17-G parallel wave failed; artifacts preserved and retries forbidden", call. = FALSE)
  }
  scan_now <- mv17g_checkpoint_scan_v1(queue, private)
  if (any(scan_now$state == "partial")) stop("MV17-G parallel wave left partial evidence", call. = FALSE)
  prefix_now <- mv17g_complete_prefix_v1(scan_now, require_incomplete = FALSE)
  collected <- lapply(seq_len(prefix_now), collect_one)
  progress <- do.call(rbind, lapply(collected, `[[`, "ledger"))
  w(progress, file.path(private, "mv17g-private-progress.csv"))
}

final_scan <- mv17g_checkpoint_scan_v1(queue, private)
final_partials <- list.files(private, pattern = "[.]partial$", recursive = TRUE, full.names = TRUE)
if (any(final_scan$state != "complete") || length(final_partials) ||
    mv17g_complete_prefix_v1(final_scan, require_incomplete = FALSE) != nrow(queue)) {
  stop("MV17-G parallel execution incomplete", call. = FALSE)
}
collected <- lapply(seq_len(nrow(queue)), collect_one)
metrics <- do.call(rbind, lapply(collected, `[[`, "metrics"))
ledger <- do.call(rbind, lapply(collected, `[[`, "ledger"))
metrics <- metrics[order(metrics$job_order, metrics$replicate, metrics$homology_dimension, metrics$summary_id, method = "radix"), , drop = FALSE]
ledger <- ledger[order(ledger$job_order, method = "radix"), , drop = FALSE]
rownames(metrics) <- rownames(ledger) <- NULL
w(metrics, file.path(private, "mv17g-private-scientific-metrics.csv"))
w(ledger, file.path(private, "mv17g-private-resource-ledger.csv"))
if (mode == "primary") {
  empirical <- mv17g_empirical_table_v1(metrics)
  w(empirical, file.path(private, "mv17g-private-empirical-calibration.csv"))
  aggregate_empirical <- mv17g_aggregate_empirical_v1(empirical)
} else {
  aggregate_empirical <- NULL
}

resource_keys <- unique(ledger[c("view", "null_family")])
resource_public <- do.call(rbind, lapply(seq_len(nrow(resource_keys)), function(i) {
  k <- resource_keys[i, , drop = FALSE]
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
private_files <- sort(list.files(private, recursive = TRUE, full.names = TRUE))
result_binding <- data.frame(
  contract_id = "mv17g_result_private_binding_v1",
  mode = mode,
  files = length(private_files),
  bytes = sum(as.numeric(file.info(private_files)$size)),
  artifact_set_sha256 = digest::digest(
    paste(sort(vapply(private_files, h, character(1L))), collapse = "\n"),
    algo = "sha256", serialize = FALSE
  ),
  tracking_state = "private_not_tracked",
  stringsAsFactors = FALSE
)
status <- data.frame(
  contract_id = "mv17g_status_v1",
  mode = mode,
  grouped_children = nrow(queue),
  scientific_runs = sum(queue$scientific_runs),
  metric_rows = nrow(metrics),
  workers = contract$workers,
  threads_per_child = contract$threads_per_child,
  retries = 0L,
  serial_prefix_children = serial_prefix,
  parallel_children = nrow(queue) - serial_prefix,
  aggregate_child_seconds = sum(ledger$wall_seconds),
  maximum_child_RSS_bytes = max(ledger$maximum_RSS_bytes),
  private_bytes = result_binding$bytes,
  labels_opened = FALSE,
  outcomes_opened = FALSE,
  clustering_jobs = 0L,
  fusion_jobs = 0L,
  biology_jobs = 0L,
  manuscript_claim_jobs = 0L,
  stringsAsFactors = FALSE
)
if (status$aggregate_child_seconds > contract$aggregate_timeout_seconds || status$private_bytes > contract$private_cap_bytes) {
  stop("MV17-G parallel aggregate cap", call. = FALSE)
}

dir.create(public_partial, recursive = TRUE)
if (mode == "primary") w(aggregate_empirical, file.path(public_partial, "mv17g-aggregate-empirical-calibration.csv"))
w(resource_public, file.path(public_partial, "mv17g-aggregate-resource.csv"))
w(result_binding, file.path(public_partial, "mv17g-private-result-binding.csv"))
w(status, file.path(public_partial, "mv17g-status.csv"))
public_files <- sort(list.files(public_partial))
w(
  data.frame(
    contract_id = "mv17g_execution_manifest_v1",
    artifact = public_files,
    bytes = as.numeric(file.info(file.path(public_partial, public_files))$size),
    sha256 = vapply(file.path(public_partial, public_files), h, character(1L)),
    stringsAsFactors = FALSE
  ),
  file.path(public_partial, "mv17g-artifact-manifest.csv")
)
if (sum(as.numeric(file.info(list.files(public_partial, full.names = TRUE))$size)) > contract$public_cap_bytes) {
  stop("MV17-G parallel public cap", call. = FALSE)
}
if (!file.rename(public_partial, public)) stop("MV17-G parallel public promotion failed", call. = FALSE)
