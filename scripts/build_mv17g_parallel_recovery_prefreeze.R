#!/usr/bin/env Rscript
options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8L) {
  stop(
    paste(
      "usage: build_mv17g_parallel_recovery_prefreeze.R",
      "<original-public-prefreeze> <original-private-prefreeze>",
      "<serial-private-root> <serial-time> <serial-stdout> <serial-stderr>",
      "<private-output> <public-output>"
    ),
    call. = FALSE
  )
}
original_public <- normalizePath(args[[1]], mustWork = TRUE)
original_private <- normalizePath(args[[2]], mustWork = TRUE)
serial_private <- normalizePath(args[[3]], mustWork = TRUE)
serial_time <- normalizePath(args[[4]], mustWork = TRUE)
serial_stdout <- normalizePath(args[[5]], mustWork = TRUE)
serial_stderr <- normalizePath(args[[6]], mustWork = TRUE)
private <- args[[7]]
public <- args[[8]]
if (dir.exists(private) || dir.exists(public)) stop("MV17-G parallel prefreeze output exists", call. = FALSE)
dir.create(private, recursive = TRUE)
dir.create(public, recursive = TRUE)

source("R/mv08z_landscape_production.R")
source("R/mv17_calibration.R")
source("R/mv17g_parallel_recovery.R")
r <- .mv08z_read_csv
w <- .mv08z_atomic_csv
h <- .mv08z_sha256_file
head <- tolower(trimws(Sys.getenv("MV17G_PARALLEL_PREFREEZE_HEAD", unset = "")))
if (!grepl("^[0-9a-f]{40}$", head)) stop("MV17G_PARALLEL_PREFREEZE_HEAD required", call. = FALSE)

verify_manifest <- function(root, name) {
  m <- r(file.path(root, name))
  paths <- file.path(root, m$artifact)
  if (!all(file.exists(paths)) ||
      !identical(unname(as.numeric(file.info(paths)$size)), unname(as.numeric(m$bytes))) ||
      !identical(unname(vapply(paths, h, character(1L))), unname(tolower(m$sha256)))) {
    stop("MV17-G original prefreeze drift", call. = FALSE)
  }
  m
}
original_manifest <- verify_manifest(original_public, "mv17g-artifact-manifest.csv")
original_private_binding <- r(file.path(original_public, "mv17g-private-binding.csv"))
original_private_paths <- file.path(original_private, original_private_binding$artifact)
if (!all(file.exists(original_private_paths)) ||
    !identical(unname(as.numeric(file.info(original_private_paths)$size)), unname(as.numeric(original_private_binding$bytes))) ||
    !identical(unname(vapply(original_private_paths, h, character(1L))), unname(tolower(original_private_binding$sha256)))) {
  stop("MV17-G original private prefreeze drift", call. = FALSE)
}

queue <- r(file.path(original_private, "mv17g-primary-grouped-queue.csv"))
locator <- r(file.path(original_private, "mv17g-private-locator.csv"))
contract_original <- r(file.path(original_public, "mv17g-contract.csv"))
contract_parallel <- mv17g_parallel_recovery_contract_v1(8L)
scan <- mv17g_checkpoint_scan_v1(queue, serial_private)
stray_partials <- list.files(serial_private, pattern = "[.]partial$", recursive = TRUE, full.names = TRUE)
if (length(stray_partials)) stop("MV17-G serial checkpoint contains stray partial files", call. = FALSE)
prefix <- mv17g_complete_prefix_v1(scan, require_incomplete = TRUE)
if (prefix < 1L) stop("MV17-G parallel recovery requires a non-empty serial prefix", call. = FALSE)

matrix_catalog_path <- file.path(serial_private, "mv17g-private-matrix-catalog.csv")
matrix_catalog <- r(matrix_catalog_path)
if (nrow(matrix_catalog) != 264L || nrow(locator) != 264L) stop("MV17-G matrix catalog drift", call. = FALSE)
matrix_paths <- file.path(serial_private, "matrices", sprintf("%s__%03d.rds", matrix_catalog$view, matrix_catalog$unit_order))
if (!all(file.exists(matrix_paths)) ||
    !identical(unname(as.numeric(file.info(matrix_paths)$size)), unname(as.numeric(matrix_catalog$bytes))) ||
    !identical(unname(vapply(matrix_paths, h, character(1L))), unname(tolower(matrix_catalog$sha256)))) {
  stop("MV17-G transition matrix drift", call. = FALSE)
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
  if (!isTRUE(ok)) stop("MV17-G serial prefix payload drift at order ", q$job_order, call. = FALSE)
}

prefix_rows <- lapply(seq_len(prefix), function(i) {
  q <- queue[i, , drop = FALSE]
  paths <- mv17g_job_artifacts_v1(q, serial_private)
  if (!all(file.exists(paths)) || file.info(paths[["stdout"]])$size != 0 || file.info(paths[["stderr"]])$size != 0) {
    stop("MV17-G serial prefix stream/artifact drift", call. = FALSE)
  }
  matrix_path <- file.path(serial_private, "matrices", sprintf("%s__%03d.rds", q$view, q$unit_order))
  validate_result(readRDS(paths[["result"]]), q, matrix_path)
  resource <- mv17c_parse_gnu_time_v1(paths[["time"]])
  if (resource$exit_status != 0L || resource$wall_seconds > contract_original$child_timeout_seconds ||
      resource$maximum_RSS_bytes > contract_original$child_RSS_cap_bytes) {
    stop("MV17-G serial prefix resource drift", call. = FALSE)
  }
  do.call(rbind, lapply(names(paths), function(role) {
    data.frame(
      job_order = q$job_order,
      role = role,
      artifact = normalizePath(paths[[role]]),
      bytes = as.numeric(file.info(paths[[role]])$size),
      sha256 = h(paths[[role]]),
      stringsAsFactors = FALSE
    )
  }))
})
prefix_inventory <- do.call(rbind, prefix_rows)

serial_outer <- mv17c_parse_gnu_time_v1(serial_time)
serial_time_text <- paste(readLines(serial_time, warn = FALSE), collapse = "\n")
controlled_stop <- grepl("signal 15", serial_time_text, ignore.case = TRUE) || serial_outer$exit_status %in% c(15L, 143L)
if (!controlled_stop || file.info(serial_stdout)$size != 0 || file.info(serial_stderr)$size != 0) {
  stop("MV17-G serial controller did not close at the authorized signal-15 boundary", call. = FALSE)
}

meminfo <- readLines("/proc/meminfo", warn = FALSE)
mem_total_kib <- as.numeric(sub("^MemTotal:[[:space:]]+([0-9]+)[[:space:]]+kB$", "\\1", grep("^MemTotal:", meminfo, value = TRUE)))
load_fields <- strsplit(readLines("/proc/loadavg", n = 1L, warn = FALSE), "[[:space:]]+")[[1]]
capacity <- data.frame(
  contract_id = "mv17g_parallel_capacity_v1",
  logical_cpus = parallel::detectCores(logical = TRUE),
  memory_total_bytes = mem_total_kib * 1024,
  load_1_minute_at_prefreeze = as.numeric(load_fields[[1]]),
  workers = contract_parallel$workers,
  threads_per_child = contract_parallel$threads_per_child,
  concurrent_child_RSS_cap_bytes = contract_parallel$concurrent_child_RSS_cap_bytes,
  sufficient = parallel::detectCores(logical = TRUE) >= 16L && mem_total_kib * 1024 >= contract_parallel$concurrent_child_RSS_cap_bytes * 2,
  stringsAsFactors = FALSE
)
if (!capacity$sufficient) stop("MV17-G eight-worker capacity gate failed", call. = FALSE)

prefix_status <- data.frame(
  contract_id = "mv17g_parallel_prefix_status_v1",
  serial_prefix_children = prefix,
  serial_prefix_scientific_runs = sum(queue$scientific_runs[seq_len(prefix)]),
  pending_children = nrow(queue) - prefix,
  pending_scientific_runs = sum(queue$scientific_runs[-seq_len(prefix)]),
  prefix_consecutive = TRUE,
  partial_children = sum(scan$state == "partial"),
  serial_exit_status = serial_outer$exit_status,
  controlled_signal_15_stop = controlled_stop,
  serial_outer_wall_seconds = serial_outer$wall_seconds,
  serial_outer_maximum_RSS_bytes = serial_outer$maximum_RSS_bytes,
  stringsAsFactors = FALSE
)

private_items <- list(
  "mv17g-parallel-prefix-artifacts.csv" = prefix_inventory,
  "mv17g-parallel-matrix-catalog.csv" = matrix_catalog,
  "mv17g-parallel-prefix-status.csv" = prefix_status
)
for (name in names(private_items)) w(private_items[[name]], file.path(private, name))
private_paths <- file.path(private, names(private_items))
private_binding <- data.frame(
  contract_id = "mv17g_parallel_private_binding_v1",
  role = c("serial_prefix_artifacts", "matrix_catalog", "serial_prefix_status"),
  artifact = names(private_items),
  rows = vapply(private_items, nrow, integer(1L)),
  bytes = as.numeric(file.info(private_paths)$size),
  sha256 = vapply(private_paths, h, character(1L)),
  tracking_state = "private_not_tracked",
  stringsAsFactors = FALSE
)

implementation_files <- c(
  "R/mv17g_parallel_recovery.R",
  "scripts/build_mv17g_parallel_recovery_prefreeze.R",
  "scripts/run_mv17g_calibration_parallel_recovery.R",
  "scripts/build_mv17g_parallel_calibration_closure.R",
  "scripts/run_mv17g_calibration_group_worker.R",
  "tests/testthat/test-mv17g-parallel-recovery.R",
  "tests/testthat/test-mv17g-calibration-closure.R",
  "docs/plans/MV17G_EIGHT_WORKER_RECOVERY_SPRINT_PLAN.md"
)
implementation <- data.frame(
  contract_id = "mv17g_parallel_implementation_binding_v1",
  file = implementation_files,
  bytes = as.numeric(file.info(implementation_files)$size),
  sha256 = vapply(implementation_files, h, character(1L)),
  stringsAsFactors = FALSE
)
if (anyNA(implementation$bytes)) stop("MV17-G parallel implementation file missing", call. = FALSE)

source_paths <- c(
  file.path(original_public, "mv17g-artifact-manifest.csv"),
  original_private_paths,
  serial_time,
  serial_stdout,
  serial_stderr,
  matrix_catalog_path
)
source_binding <- data.frame(
  contract_id = "mv17g_parallel_source_binding_v1",
  role = c(
    "original_public_prefreeze_manifest",
    paste0("original_private_", original_private_binding$role),
    "serial_outer_time", "serial_outer_stdout", "serial_outer_stderr", "serial_matrix_catalog"
  ),
  bytes = as.numeric(file.info(source_paths)$size),
  sha256 = vapply(source_paths, h, character(1L)),
  stringsAsFactors = FALSE
)

contract <- cbind(
  contract_parallel,
  original_prefreeze_head = contract_original$execution_head,
  recovery_implementation_head = head,
  serial_prefix_children = prefix,
  primary_children = nrow(queue),
  primary_scientific_runs = sum(queue$scientific_runs),
  parallel_primary_and_repeat_authorized_after_commit = TRUE,
  real_localization_authorized = FALSE,
  labels_authorized = FALSE,
  outcomes_authorized = FALSE,
  clustering_authorized = FALSE,
  fusion_authorized = FALSE,
  biology_authorized = FALSE,
  manuscript_claims_authorized = FALSE
)
validation <- data.frame(
  contract_id = "mv17g_parallel_prefreeze_validation_v1",
  check_id = c(
    "original_prefreeze_bound", "original_private_prefreeze_bound", "implementation_head_bound",
    "serial_controller_stopped", "serial_streams_empty", "prefix_nonempty", "prefix_incomplete",
    "prefix_consecutive", "prefix_no_partials", "prefix_artifacts_complete", "prefix_payloads_valid",
    "prefix_child_resources_pass", "all_264_matrices_bound", "private_transition_bound",
    "eight_workers", "one_thread_per_child", "zero_retries", "capacity_sufficient",
    "concurrent_RSS_bounded", "eight_implementation_files", "execution_only_after_commit",
    "real_localization_closed", "downstream_firewall", "aggregate_only_public"
  ),
  passed = c(
    nrow(original_manifest) >= 1L, nrow(original_private_binding) == 4L, grepl("^[0-9a-f]{40}$", head),
    controlled_stop, file.info(serial_stdout)$size == 0 && file.info(serial_stderr)$size == 0,
    prefix >= 1L, prefix < nrow(queue), prefix == mv17g_complete_prefix_v1(scan),
    !any(scan$state == "partial") && length(stray_partials) == 0L, nrow(prefix_inventory) == prefix * 4L, TRUE,
    max(vapply(seq_len(prefix), function(i) mv17c_parse_gnu_time_v1(mv17g_job_artifacts_v1(queue[i, ], serial_private)[["time"]])$wall_seconds, numeric(1L))) <= contract_original$child_timeout_seconds,
    nrow(matrix_catalog) == 264L, nrow(private_binding) == 3L,
    contract$workers == 8L, contract$threads_per_child == 1L, contract$retries == 0L,
    capacity$sufficient, contract$concurrent_child_RSS_cap_bytes == 8 * contract$child_RSS_cap_bytes,
    nrow(implementation) == 8L, contract$parallel_primary_and_repeat_authorized_after_commit,
    !contract$real_localization_authorized,
    !any(contract[c("labels_authorized", "outcomes_authorized", "clustering_authorized", "fusion_authorized", "biology_authorized", "manuscript_claims_authorized")]),
    !any(c("unit_id", "identity_token", "source_path", "accepted_diagram_sha256", "artifact") %in% names(prefix_status))
  )
)
if (!all(validation$passed)) stop("MV17-G parallel recovery prefreeze failed", call. = FALSE)

items <- list(
  "mv17g-parallel-contract.csv" = contract,
  "mv17g-parallel-capacity.csv" = capacity,
  "mv17g-parallel-prefix-summary.csv" = prefix_status,
  "mv17g-parallel-private-binding.csv" = private_binding,
  "mv17g-parallel-source-binding.csv" = source_binding,
  "mv17g-parallel-implementation-binding.csv" = implementation,
  "mv17g-parallel-validation.csv" = validation,
  "mv17g-parallel-decision.csv" = data.frame(
    contract_id = "mv17g_parallel_decision_v1",
    eight_worker_primary_and_repeat_authorized_after_commit = TRUE,
    serial_prefix_recomputation_forbidden = TRUE,
    retries = 0L,
    real_localization_authorized = FALSE,
    downstream_surfaces = "closed",
    stringsAsFactors = FALSE
  )
)
for (name in names(items)) w(items[[name]], file.path(public, name))
writeLines(
  c(
    "# MV17-G eight-worker recovery prefreeze",
    "",
    paste0("The completed serial prefix contains ", prefix, " atomic children and is privately hash-bound."),
    "After commit, the exact remaining primary queue and the frozen repeat may run with eight single-threaded children.",
    "Completed serial children may be adopted but never recomputed. Retries and all downstream scientific surfaces remain closed."
  ),
  file.path(public, "MV17G_EIGHT_WORKER_RECOVERY_PREFREEZE_2026-08-26.md")
)
files <- sort(list.files(public))
w(
  data.frame(
    contract_id = "mv17g_parallel_prefreeze_manifest_v1",
    artifact = files,
    bytes = as.numeric(file.info(file.path(public, files))$size),
    sha256 = vapply(file.path(public, files), h, character(1L)),
    stringsAsFactors = FALSE
  ),
  file.path(public, "mv17g-parallel-artifact-manifest.csv")
)
message("Built MV17-G eight-worker recovery prefreeze; checks=", nrow(validation), "; prefix=", prefix)
