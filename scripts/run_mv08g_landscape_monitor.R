#!/usr/bin/env Rscript

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1", RAYON_NUM_THREADS = "1")
options(warn = 2)
for (package in c("digest", "processx", "ps")) {
  if (!requireNamespace(package, quietly = TRUE)) stop(package, " required.")
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8L) {
  stop("usage: run_mv08g_landscape_monitor.R PREFREEZE PH475_ROOT PH500_ROOT RUST_LIBRARY PRIVATE_ROOT PUBLIC_DIR EXPECTED_HEAD PRIMARY_PREFREEZE")
}
prefreeze <- args[[1L]]; ph475 <- args[[2L]]; ph500 <- args[[3L]]
rust <- normalizePath(args[[4L]], winslash = "/", mustWork = TRUE)
private_root <- args[[5L]]; public_dir <- args[[6L]]
expected_head <- tolower(trimws(args[[7L]])); primary <- args[[8L]]
head <- tolower(trimws(system2("git", c("rev-parse", "HEAD"), stdout = TRUE)))
if (head != expected_head) stop("MV8-G landscape monitor exact HEAD mismatch.")
source("R/mv08g_panel_sensitivity.R")
sha <- function(path) digest::digest(file = path, algo = "sha256",
                                     serialize = FALSE)
contract <- read.csv(file.path(prefreeze, "mv08g-landscape-contract.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)
within_queue <- read.csv(file.path(prefreeze, "mv08g-landscape-queue.csv"),
                         stringsAsFactors = FALSE, check.names = FALSE)
shift_queue <- read.csv(file.path(prefreeze, "mv08g-matched-shift-queue.csv"),
                        stringsAsFactors = FALSE, check.names = FALSE)
repeat_within <- read.csv(file.path(prefreeze, "mv08g-landscape-repeat-queue.csv"),
                          stringsAsFactors = FALSE, check.names = FALSE)
repeat_shift <- read.csv(file.path(prefreeze,
  "mv08g-matched-shift-repeat-queue.csv"), stringsAsFactors = FALSE,
  check.names = FALSE)
decision <- read.csv(file.path(prefreeze, "mv08g-landscape-decision.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)
if (nrow(contract) != 1L || nrow(within_queue) != 20L || nrow(shift_queue) != 20L ||
    nrow(repeat_within) != 4L || nrow(repeat_shift) != 4L ||
    contract$rust_library_sha256 != sha(rust) ||
    decision$decision != "authorize_20_within_20_matched_and_eight_repeats") {
  stop("MV8-G landscape execution gate is stale.")
}
for (subdir in c("landscape", "matched-shift", "repeat/landscape",
                 "repeat/matched-shift", "logs", "repeat/logs")) {
  dir.create(file.path(private_root, subdir), recursive = TRUE,
             showWarnings = FALSE)
}
dir.create(public_dir, recursive = TRUE, showWarnings = FALSE)
tree_rss <- function(pid) {
  root <- tryCatch(ps::ps_handle(pid), error = function(error) NULL)
  if (is.null(root)) return(0)
  handles <- c(list(root), tryCatch(ps::ps_children(root, recursive = TRUE),
                                   error = function(error) list()))
  sum(vapply(handles, function(handle) tryCatch(
    as.numeric(ps::ps_memory_info(handle)[["rss"]]),
    error = function(error) 0), numeric(1L)))
}
read_empty <- function(path) if (file.exists(path)) read.csv(
  path, stringsAsFactors = FALSE, check.names = FALSE) else data.frame()
write_atomic <- function(value, path) {
  partial <- tempfile(pattern = basename(path), tmpdir = dirname(path))
  write.csv(value, partial, row.names = FALSE, na = "")
  if (!file.rename(partial, path)) {
    unlink(partial); stop("Failed to atomically publish MV8-G landscape ledger.")
  }
}
storage_bytes <- function() {
  files <- list.files(private_root, recursive = TRUE, full.names = TRUE)
  files <- files[!grepl("/logs/|\\\\logs\\\\", files)]
  if (length(files)) sum(as.numeric(file.info(files)$size), na.rm = TRUE) else 0
}
validate_group <- function(path, row) {
  status_path <- file.path(path, "status.csv")
  distances_path <- file.path(path, "distances.csv")
  metrics_path <- file.path(path, "metrics.csv")
  if (!all(file.exists(c(status_path, distances_path, metrics_path)))) return(FALSE)
  status <- read.csv(status_path, stringsAsFactors = FALSE, check.names = FALSE)
  expected_rows <- as.integer(row$component_rows)
  nrow(status) == 1L && status$completion_state == "complete" &&
    status$group_id == row$group_id && status$component_rows == expected_rows &&
    status$distances_sha256 == sha(distances_path) &&
    status$metrics_sha256 == sha(metrics_path)
}
ledger_path <- file.path(private_root, "group-resource-metrics.csv")
repeat_ledger_path <- file.path(private_root, "repeat-group-resource-metrics.csv")
ledger <- read_empty(ledger_path); repeat_ledger <- read_empty(repeat_ledger_path)
run_group <- function(row, scope, repeat_mode = FALSE) {
  scope <- match.arg(scope, c("within475", "matched500_475"))
  base <- if (scope == "within475") "landscape" else "matched-shift"
  root <- file.path(private_root, if (repeat_mode) paste0("repeat/", base) else base)
  output <- file.path(root, row$group_id)
  job_id <- paste(if (repeat_mode) "repeat" else "production", scope,
                  row$group_id, sep = "__")
  metrics <- if (repeat_mode) repeat_ledger else ledger
  hit <- if (nrow(metrics)) metrics[metrics$job_id == job_id, , drop = FALSE]
    else data.frame()
  if (nrow(hit)) {
    if (nrow(hit) != 1L || hit$disposition != "completed" ||
        !validate_group(output, row) ||
        hit$distances_sha256 != sha(file.path(output, "distances.csv"))) {
      stop("MV8-G landscape resume state is stale: ", job_id)
    }
    return(hit)
  }
  if (dir.exists(output)) stop("Unowned MV8-G landscape group exists: ", job_id)
  elapsed_used <- sum(ledger$elapsed_seconds, na.rm = TRUE) +
    sum(repeat_ledger$elapsed_seconds, na.rm = TRUE)
  if (elapsed_used > contract$aggregate_elapsed_cap_seconds ||
      storage_bytes() > contract$aggregate_storage_cap_bytes) {
    stop("MV8-G landscape aggregate resource cap reached before ", job_id)
  }
  script <- if (scope == "within475") "scripts/run_mv08g_landscape_group.R" else
    "scripts/run_mv08g_matched_shift_group.R"
  script_args <- if (scope == "within475") c(
    prefreeze, ph475, rust, row$group_id, root,
    contract$implementation_root_sha256) else c(
      prefreeze, ph500, ph475, rust, row$group_id, root,
      contract$implementation_root_sha256)
  log_dir <- file.path(private_root, if (repeat_mode) "repeat/logs" else "logs")
  stem <- gsub("[^A-Za-z0-9_.-]", "_", job_id)
  started <- Sys.time()
  process <- processx::process$new(Sys.which("Rscript"),
    c("--vanilla", script, script_args),
    stdout = file.path(log_dir, paste0(stem, "__stdout.txt")),
    stderr = file.path(log_dir, paste0(stem, "__stderr.txt")), cleanup_tree = TRUE)
  peak <- 0; cap_failure <- ""
  while (process$is_alive()) {
    Sys.sleep(0.25); peak <- max(peak, tree_rss(process$get_pid()))
    elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    if (elapsed > contract$elapsed_cap_seconds_per_group) {
      cap_failure <- "elapsed_cap_exceeded"; process$kill_tree()
    } else if (peak > contract$rss_cap_bytes_per_group) {
      cap_failure <- "rss_cap_exceeded"; process$kill_tree()
    }
  }
  process$wait(timeout = 5000)
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  status <- process$get_exit_status()
  valid <- identical(status, 0L) && validate_group(output, row)
  disposition <- if (nzchar(cap_failure)) cap_failure else if (valid)
    "completed" else "failed"
  metric <- data.frame(
    contract_id = "mv08g_landscape_resource_metric_v1", job_id = job_id,
    scope = scope, group_id = row$group_id, seed = row$seed,
    view_id = row$view_id, homology_dimension = row$homology_dimension,
    repeat_mode = repeat_mode, disposition = disposition,
    exit_status = if (is.null(status)) NA_integer_ else status,
    elapsed_seconds = elapsed, peak_process_tree_rss_bytes = peak,
    elapsed_cap_seconds = contract$elapsed_cap_seconds_per_group,
    rss_cap_bytes = contract$rss_cap_bytes_per_group,
    component_rows = row$component_rows,
    distances_bytes = if (valid) as.numeric(file.info(
      file.path(output, "distances.csv"))$size) else NA,
    distances_sha256 = if (valid) sha(file.path(output, "distances.csv")) else NA,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE)
  metrics <- if (nrow(metrics)) rbind(metrics, metric) else metric
  if (repeat_mode) {
    repeat_ledger <<- metrics; write_atomic(metrics, repeat_ledger_path)
  } else {
    ledger <<- metrics; write_atomic(metrics, ledger_path)
  }
  if (disposition != "completed") stop("MV8-G landscape stopped under zero-retry: ",
                                        job_id, " (", disposition, ")")
  metric
}
within_metrics <- do.call(rbind, lapply(seq_len(nrow(within_queue)), function(i) {
  run_group(within_queue[i, , drop = FALSE], "within475", FALSE)
}))
shift_metrics <- do.call(rbind, lapply(seq_len(nrow(shift_queue)), function(i) {
  run_group(shift_queue[i, , drop = FALSE], "matched500_475", FALSE)
}))
repeat_rows <- list()
for (index in seq_len(nrow(repeat_within))) {
  prod <- within_metrics[within_metrics$group_id == repeat_within$group_id[[index]], ]
  rep <- run_group(repeat_within[index, , drop = FALSE], "within475", TRUE)
  repeat_rows[[length(repeat_rows) + 1L]] <- data.frame(
    contract_id = "mv08g_landscape_repeat_validation_v1", scope = "within475",
    group_id = prod$group_id, seed = prod$seed, view_id = prod$view_id,
    homology_dimension = prod$homology_dimension,
    production_distances_sha256 = prod$distances_sha256,
    repeat_distances_sha256 = rep$distances_sha256,
    byte_identical = prod$distances_sha256 == rep$distances_sha256 &&
      prod$distances_bytes == rep$distances_bytes,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE)
  prod <- shift_metrics[shift_metrics$group_id == repeat_shift$group_id[[index]], ]
  rep <- run_group(repeat_shift[index, , drop = FALSE], "matched500_475", TRUE)
  repeat_rows[[length(repeat_rows) + 1L]] <- data.frame(
    contract_id = "mv08g_landscape_repeat_validation_v1", scope = "matched500_475",
    group_id = prod$group_id, seed = prod$seed, view_id = prod$view_id,
    homology_dimension = prod$homology_dimension,
    production_distances_sha256 = prod$distances_sha256,
    repeat_distances_sha256 = rep$distances_sha256,
    byte_identical = prod$distances_sha256 == rep$distances_sha256 &&
      prod$distances_bytes == rep$distances_bytes,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE)
}
repeats <- do.call(rbind, repeat_rows)
if (!all(repeats$byte_identical)) stop("MV8-G landscape repeat mismatch.")
aggregate_elapsed <- sum(ledger$elapsed_seconds) + sum(repeat_ledger$elapsed_seconds)
if (aggregate_elapsed > contract$aggregate_elapsed_cap_seconds ||
    storage_bytes() > contract$aggregate_storage_cap_bytes) {
  stop("MV8-G landscape aggregate resource cap exceeded at completion.")
}
decision_out <- data.frame(
  contract_id = "mv08g_landscape_execution_decision_v1",
  decision = "landscapes_complete_await_independent_R_Persim_validation",
  within475_groups = nrow(within_metrics), matched_shift_groups = nrow(shift_metrics),
  repeat_groups = nrow(repeats), aggregate_elapsed_seconds = aggregate_elapsed,
  aggregate_elapsed_cap_seconds = contract$aggregate_elapsed_cap_seconds,
  private_storage_bytes = storage_bytes(),
  aggregate_storage_cap_bytes = contract$aggregate_storage_cap_bytes,
  maximum_process_tree_rss_bytes = max(c(ledger$peak_process_tree_rss_bytes,
    repeat_ledger$peak_process_tree_rss_bytes)), comparison_jobs_authorized = 0L,
  hca_fastq_download_authorized = FALSE, raw_reprocessing_authorized = FALSE,
  label_access_authorized = FALSE, biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE)
write_atomic(within_metrics, file.path(public_dir, "mv08g-within475-inventory.csv"))
write_atomic(shift_metrics, file.path(public_dir, "mv08g-matched-shift-inventory.csv"))
write_atomic(repeats, file.path(public_dir, "mv08g-landscape-repeat-validation.csv"))
write_atomic(decision_out, file.path(public_dir, "mv08g-landscape-decision.csv"))
message("MV8-G landscapes complete: 40 groups and eight exact distance repeats")
