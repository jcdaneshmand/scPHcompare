#!/usr/bin/env Rscript
options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop(
    paste(
      "usage: run_mv17g_gene_geometry_recovery.R",
      "<recovery-public-prefreeze> <recovery-private-prefreeze>",
      "<gene-primary-or-repeat> <matrix-root> <private-output>",
      "<public-output>"
    ),
    call. = FALSE
  )
}

recovery_public <- normalizePath(args[[1L]], mustWork = TRUE)
recovery_private <- normalizePath(args[[2L]], mustWork = TRUE)
mode <- match.arg(args[[3L]], c("gene_primary", "repeat"))
matrix_root <- normalizePath(args[[4L]], mustWork = TRUE)
private_argument <- args[[5L]]
if (!dir.exists(private_argument)) {
  dir.create(private_argument, recursive = TRUE)
}
private <- normalizePath(private_argument, mustWork = TRUE)
public <- args[[6L]]
public_partial <- paste0(public, ".partial")
if (dir.exists(public) || dir.exists(public_partial)) {
  stop("MV17-G geometry-recovery public output exists", call. = FALSE)
}

source("R/mv08z_landscape_production.R")
source("R/mv17_null_models.R")
source("R/mv17_calibration.R")
source("R/mv17_localization.R")
source("R/mv17_full_calibration.R")
source("R/mv17_full_calibration_geometry_v2.R")
source("R/mv17g_parallel_recovery.R")
source("R/mv17g_gene_geometry_recovery.R")
read_csv <- .mv08z_read_csv
write_csv <- .mv08z_atomic_csv
sha256 <- .mv08z_sha256_file

manifest <- read_csv(file.path(
  recovery_public, "mv17g-geometry-recovery-manifest.csv"
))
manifest_paths <- file.path(recovery_public, manifest$artifact)
if (!all(file.exists(manifest_paths)) ||
    !identical(
      unname(as.numeric(file.info(manifest_paths)$size)),
      unname(as.numeric(manifest$bytes))
    ) ||
    !identical(
      unname(vapply(manifest_paths, sha256, character(1L))),
      unname(tolower(manifest$sha256))
    )) {
  stop("MV17-G geometry-recovery public prefreeze drift", call. = FALSE)
}

private_binding <- read_csv(file.path(
  recovery_public, "mv17g-geometry-recovery-private-binding.csv"
))
private_paths <- file.path(recovery_private, private_binding$artifact)
if (!all(file.exists(private_paths)) ||
    !identical(
      unname(as.numeric(file.info(private_paths)$size)),
      unname(as.numeric(private_binding$bytes))
    ) ||
    !identical(
      unname(vapply(private_paths, sha256, character(1L))),
      unname(tolower(private_binding$sha256))
    )) {
  stop("MV17-G geometry-recovery private prefreeze drift", call. = FALSE)
}

contract <- read_csv(file.path(
  recovery_public, "mv17g-geometry-recovery-contract.csv"
))
required_contract <- c(
  "execution_authorized_after_commit", "workers", "threads_per_child",
  "retries", "child_timeout_seconds", "child_RSS_cap_bytes",
  "aggregate_timeout_seconds", "private_cap_bytes", "public_cap_bytes",
  "cell_prefix_children", "gene_primary_children",
  "gene_primary_scientific_runs", "repeat_children",
  "repeat_scientific_runs", "cell_geometry", "gene_geometry",
  "labels_opened", "outcomes_opened"
)
if (nrow(contract) != 1L ||
    !all(required_contract %in% names(contract)) ||
    !isTRUE(contract$execution_authorized_after_commit) ||
    contract$workers != 8L || contract$threads_per_child != 1L ||
    contract$retries != 0L || contract$cell_prefix_children != 528L ||
    contract$gene_primary_children != 660L ||
    contract$gene_primary_scientific_runs != 52404L ||
    contract$repeat_children != 27L ||
    contract$repeat_scientific_runs != 2085L ||
    contract$cell_geometry != "euclidean_shared_pca_v1" ||
    contract$gene_geometry != "euclidean_correlation_chord_v1" ||
    contract$labels_opened || contract$outcomes_opened) {
  stop("MV17-G geometry-recovery execution is not authorized", call. = FALSE)
}

queue_file <- if (mode == "gene_primary") {
  "mv17g-corrected-gene-primary-queue.csv"
} else {
  "mv17g-corrected-repeat-queue.csv"
}
queue <- read_csv(file.path(recovery_private, queue_file))
required_queue <- c(
  "job_order", "view", "unit_order", "null_family", "seed_first",
  "seed_last", "replicate_count", "scientific_runs"
)
if (!all(required_queue %in% names(queue)) || anyDuplicated(queue$job_order)) {
  stop("MV17-G geometry-recovery queue schema drift", call. = FALSE)
}
if (mode == "gene_primary") {
  queue_ok <- nrow(queue) == 660L &&
    identical(as.integer(queue$job_order), 529:1188) &&
    all(queue$view == "gene") &&
    sum(queue$scientific_runs) == 52404L
} else {
  queue_ok <- nrow(queue) == 27L &&
    identical(as.integer(queue$job_order), seq_len(27L)) &&
    setequal(queue$view, c("cell", "gene")) &&
    sum(queue$scientific_runs) == 2085L
}
if (!isTRUE(queue_ok)) {
  stop("MV17-G geometry-recovery queue cardinality drift", call. = FALSE)
}

matrix_catalog <- read_csv(file.path(
  recovery_private, "mv17g-geometry-recovery-matrix-catalog.csv"
))
matrix_paths <- file.path(
  matrix_root, "matrices",
  sprintf("%s__%03d.rds", matrix_catalog$view, matrix_catalog$unit_order)
)
if (nrow(matrix_catalog) != 264L ||
    !setequal(matrix_catalog$view, c("cell", "gene")) ||
    !all(file.exists(matrix_paths)) ||
    !identical(
      unname(as.numeric(file.info(matrix_paths)$size)),
      unname(as.numeric(matrix_catalog$bytes))
    ) ||
    !identical(
      unname(vapply(matrix_paths, sha256, character(1L))),
      unname(tolower(matrix_catalog$sha256))
    )) {
  stop("MV17-G geometry-recovery matrix binding drift", call. = FALSE)
}

dir.create(file.path(private, "jobs"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(private, "logs"), showWarnings = FALSE, recursive = TRUE)
worker <- normalizePath(
  "scripts/run_mv17g_calibration_group_worker_v2.R", mustWork = TRUE
)

validate_result <- function(result, q, matrix_path) {
  expected_seeds <- if (q$null_family == "observed") {
    0L
  } else {
    seq.int(q$seed_first, length.out = q$replicate_count)
  }
  expected_geometry <- mv17g_geometry_id_v2(q$view)
  ok <- identical(result$contract_id, "mv17g_group_result_v2") &&
    identical(result$view, q$view) &&
    identical(result$geometry, expected_geometry) &&
    identical(result$null_family, q$null_family) &&
    result$seed_first == min(expected_seeds) &&
    result$seed_last == max(expected_seeds) &&
    result$replicate_count == length(expected_seeds) &&
    result$matrix_sha256 == sha256(matrix_path) &&
    nrow(result$metrics) == length(expected_seeds) * 8L &&
    setequal(unique(result$metrics$seed), expected_seeds) &&
    identical(unique(result$metrics$geometry), expected_geometry) &&
    all(result$metrics$h0_mst_maximum_absolute_error <= 1e-8) &&
    isTRUE(result$finite) && !result$labels_opened && !result$outcomes_opened
  if (!isTRUE(ok)) {
    stop(
      "MV17-G geometry-recovery child payload drift at order ",
      q$job_order, call. = FALSE
    )
  }
  result
}

collect_one <- function(i) {
  q <- queue[i, , drop = FALSE]
  paths <- mv17g_job_artifacts_v1(q, private)
  if (!all(file.exists(paths)) ||
      file.info(paths[["stdout"]])$size != 0 ||
      file.info(paths[["stderr"]])$size != 0) {
    stop(
      "MV17-G geometry-recovery artifact/stream drift at order ",
      q$job_order, call. = FALSE
    )
  }
  resource <- mv17c_parse_gnu_time_v1(paths[["time"]])
  if (resource$exit_status != 0L ||
      resource$wall_seconds > contract$child_timeout_seconds ||
      resource$maximum_RSS_bytes > contract$child_RSS_cap_bytes) {
    stop(
      "MV17-G geometry-recovery resource cap at order ",
      q$job_order, call. = FALSE
    )
  }
  matrix_path <- file.path(
    matrix_root, "matrices", sprintf("%s__%03d.rds", q$view, q$unit_order)
  )
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
    output_sha256 = sha256(paths[["result"]]),
    execution_origin = "gene_geometry_recovery_v2",
    stringsAsFactors = FALSE
  )
  list(
    metrics = metrics[c(
      "job_order", "view", "unit_order", "null_family", "replicate",
      "seed", "homology_dimension", "summary_id", "value",
      "h0_mst_maximum_absolute_error", "geometry"
    )],
    ledger = ledger
  )
}

scan <- mv17g_checkpoint_scan_v1(queue, private)
partials <- list.files(
  private, pattern = "[.]partial$", recursive = TRUE, full.names = TRUE
)
if (any(scan$state == "partial") || length(partials)) {
  stop("MV17-G geometry-recovery checkpoint has partial evidence", call. = FALSE)
}
prefix <- mv17g_complete_queue_prefix_v2(scan, queue)
pending <- queue$job_order[seq_len(nrow(queue)) > prefix]
batches <- mv17g_parallel_batches_v1(pending, contract$workers)
for (batch in batches) {
  rows <- match(batch, queue$job_order)
  wave <- mv17g_run_parallel_wave_v2(
    queue[rows, , drop = FALSE], private, matrix_root, worker, contract
  )
  if (any(wave$exit_status != 0L) || any(wave$artifacts != 4L)) {
    stop(
      "MV17-G geometry-recovery wave failed; evidence preserved and retries forbidden",
      call. = FALSE
    )
  }
  scan_now <- mv17g_checkpoint_scan_v1(queue, private)
  prefix_now <- mv17g_complete_queue_prefix_v2(scan_now, queue)
  collected <- lapply(seq_len(prefix_now), collect_one)
  progress <- do.call(rbind, lapply(collected, `[[`, "ledger"))
  write_csv(progress, file.path(private, "mv17g-v2-private-progress.csv"))
}

final_scan <- mv17g_checkpoint_scan_v1(queue, private)
final_partials <- list.files(
  private, pattern = "[.]partial$", recursive = TRUE, full.names = TRUE
)
if (any(final_scan$state != "complete") || length(final_partials) ||
    mv17g_complete_queue_prefix_v2(final_scan, queue) != nrow(queue)) {
  stop("MV17-G geometry-recovery execution incomplete", call. = FALSE)
}

collected <- lapply(seq_len(nrow(queue)), collect_one)
metrics <- do.call(rbind, lapply(collected, `[[`, "metrics"))
ledger <- do.call(rbind, lapply(collected, `[[`, "ledger"))
metrics <- metrics[order(
  metrics$job_order, metrics$replicate, metrics$homology_dimension,
  metrics$summary_id, method = "radix"
), , drop = FALSE]
ledger <- ledger[order(ledger$job_order, method = "radix"), , drop = FALSE]
rownames(metrics) <- rownames(ledger) <- NULL
write_csv(metrics, file.path(private, "mv17g-v2-private-scientific-metrics.csv"))
write_csv(ledger, file.path(private, "mv17g-v2-private-resource-ledger.csv"))
if (mode == "gene_primary") {
  empirical <- mv17g_empirical_table_v1(metrics)
  write_csv(
    empirical, file.path(private, "mv17g-v2-private-empirical-calibration.csv")
  )
  aggregate_empirical <- mv17g_aggregate_empirical_v1(empirical)
} else {
  aggregate_empirical <- NULL
}

resource_keys <- unique(ledger[c("view", "null_family")])
resource_public <- do.call(rbind, lapply(seq_len(nrow(resource_keys)), function(i) {
  key <- resource_keys[i, , drop = FALSE]
  take <- ledger$view == key$view & ledger$null_family == key$null_family
  data.frame(
    key,
    grouped_children = sum(take),
    scientific_runs = sum(ledger$replicate_count[take]),
    aggregate_child_seconds = sum(ledger$wall_seconds[take]),
    maximum_child_RSS_bytes = max(ledger$maximum_RSS_bytes[take]),
    aggregate_output_bytes = sum(ledger$output_bytes[take]),
    stringsAsFactors = FALSE
  )
}))
private_files <- sort(list.files(private, recursive = TRUE, full.names = TRUE))
result_binding <- data.frame(
  contract_id = "mv17g_geometry_recovery_private_binding_v2",
  mode = mode,
  files = length(private_files),
  bytes = sum(as.numeric(file.info(private_files)$size)),
  artifact_set_sha256 = digest::digest(
    paste(sort(vapply(private_files, sha256, character(1L))), collapse = "\n"),
    algo = "sha256", serialize = FALSE
  ),
  tracking_state = "private_not_tracked",
  stringsAsFactors = FALSE
)
status <- data.frame(
  contract_id = "mv17g_geometry_recovery_status_v2",
  mode = mode,
  grouped_children = nrow(queue),
  scientific_runs = sum(queue$scientific_runs),
  metric_rows = nrow(metrics),
  workers = contract$workers,
  threads_per_child = contract$threads_per_child,
  retries = 0L,
  aggregate_child_seconds = sum(ledger$wall_seconds),
  maximum_child_RSS_bytes = max(ledger$maximum_RSS_bytes),
  private_bytes = result_binding$bytes,
  cell_geometry = contract$cell_geometry,
  gene_geometry = contract$gene_geometry,
  labels_opened = FALSE,
  outcomes_opened = FALSE,
  clustering_jobs = 0L,
  fusion_jobs = 0L,
  biology_jobs = 0L,
  manuscript_claim_jobs = 0L,
  stringsAsFactors = FALSE
)
if (status$aggregate_child_seconds > contract$aggregate_timeout_seconds ||
    status$private_bytes > contract$private_cap_bytes) {
  stop("MV17-G geometry-recovery aggregate cap", call. = FALSE)
}

dir.create(public_partial, recursive = TRUE)
if (mode == "gene_primary") {
  write_csv(
    aggregate_empirical,
    file.path(public_partial, "mv17g-v2-aggregate-gene-calibration.csv")
  )
}
write_csv(
  resource_public, file.path(public_partial, "mv17g-v2-aggregate-resource.csv")
)
write_csv(
  result_binding, file.path(public_partial, "mv17g-v2-private-result-binding.csv")
)
write_csv(status, file.path(public_partial, "mv17g-v2-status.csv"))
public_files <- sort(list.files(public_partial))
write_csv(
  data.frame(
    contract_id = "mv17g_geometry_recovery_execution_manifest_v2",
    artifact = public_files,
    bytes = as.numeric(file.info(file.path(public_partial, public_files))$size),
    sha256 = vapply(
      file.path(public_partial, public_files), sha256, character(1L)
    ),
    stringsAsFactors = FALSE
  ),
  file.path(public_partial, "mv17g-v2-artifact-manifest.csv")
)
if (sum(as.numeric(file.info(list.files(
  public_partial, full.names = TRUE
))$size)) > contract$public_cap_bytes) {
  stop("MV17-G geometry-recovery public cap", call. = FALSE)
}
if (!file.rename(public_partial, public)) {
  stop("MV17-G geometry-recovery public promotion failed", call. = FALSE)
}
