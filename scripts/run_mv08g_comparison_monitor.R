#!/usr/bin/env Rscript

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1")
options(warn = 2)
for (package in c("digest", "processx", "ps")) {
  if (!requireNamespace(package, quietly = TRUE)) stop(package, " required.")
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 9L) {
  stop("usage: run_mv08g_comparison_monitor.R PREFREEZE PRIMARY_PREFREEZE MV07H_LANDSCAPE_ROOT MV08G_LANDSCAPE_ROOT MATCHED_SHIFT_ROOT PRIVATE_ROOT PUBLIC_DIR EXPECTED_HEAD RESULT_SUBDIR")
}
prefreeze <- args[[1L]]; primary <- args[[2L]]; root500 <- args[[3L]]
root475 <- args[[4L]]; shift_root <- args[[5L]]; private_root <- args[[6L]]
public_dir <- args[[7L]]; expected_head <- tolower(trimws(args[[8L]]))
result_subdir <- args[[9L]]
head <- tolower(trimws(system2("git", c("rev-parse", "HEAD"), stdout = TRUE)))
if (head != expected_head) stop("MV8-G comparison monitor exact HEAD mismatch.")
contract <- read.csv(file.path(prefreeze, "mv08g-comparison-contract.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)
decision <- read.csv(file.path(prefreeze, "mv08g-comparison-decision.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)
if (nrow(contract) != 1L || decision$decision !=
      "authorize_one_label_closed_comparison_job" ||
    decision$comparison_jobs_authorized != 1L || contract$workers != 1L ||
    contract$retries != 0L || grepl("[/\\\\]", result_subdir)) {
  stop("MV8-G comparison execution gate is stale or result subdir is unsafe.")
}
dir.create(private_root, recursive = TRUE, showWarnings = FALSE)
dir.create(public_dir, recursive = TRUE, showWarnings = FALSE)
result <- file.path(private_root, result_subdir)
log_dir <- file.path(private_root, "logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
ledger_path <- file.path(private_root, "comparison-resource-metric.csv")
sha <- function(path) digest::digest(file = path, algo = "sha256",
                                     serialize = FALSE)
truth <- function(value) tolower(as.character(value)) %in% c("true", "t", "1")
tree_rss <- function(pid) {
  root <- tryCatch(ps::ps_handle(pid), error = function(error) NULL)
  if (is.null(root)) return(0)
  handles <- c(list(root), tryCatch(ps::ps_children(root, recursive = TRUE),
                                   error = function(error) list()))
  sum(vapply(handles, function(handle) tryCatch(
    as.numeric(ps::ps_memory_info(handle)[["rss"]]),
    error = function(error) 0), numeric(1L)))
}
write_atomic <- function(value, path) {
  partial <- tempfile(pattern = basename(path), tmpdir = dirname(path))
  write.csv(value, partial, row.names = FALSE, na = "")
  if (!file.rename(partial, path)) {
    unlink(partial); stop("Failed to atomically publish MV8-G comparison metric.")
  }
}
validate_result <- function(path) {
  required <- c("mv08g-seed-distance-comparison.csv",
    "mv08g-top10-neighbor-overlap.csv", "mv08g-normalized-matched-shifts.csv",
    "mv08g-candidate-pam-partitions.csv", "mv08g-pam-stability-summary.csv",
    "mv08g-panel-selected-k.csv", "mv08g-panel-selected-k-agreement.csv",
    "mv08g-pam-k-panel-agreement.csv", "mv08g-fixed500k-partitions.csv",
    "mv08g-fixed500k-panel-agreement.csv", "mv08g-component-summary.csv",
    "mv08g-decision.csv")
  manifest_path <- file.path(path, "mv08g-artifact-manifest.csv")
  if (!all(file.exists(file.path(path, c(required,
    "mv08g-artifact-manifest.csv"))))) return(FALSE)
  manifest <- read.csv(manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!setequal(manifest$file, required) ||
      any(vapply(file.path(path, manifest$file), sha, character(1L)) !=
          manifest$sha256) || any(truth(manifest$contains_expression)) ||
      any(truth(manifest$contains_cell_barcode)) ||
      any(truth(manifest$contains_absolute_private_path)) ||
      any(truth(manifest$contains_biological_label))) return(FALSE)
  result_decision <- read.csv(file.path(path, "mv08g-decision.csv"),
                              stringsAsFactors = FALSE, check.names = FALSE)
  nrow(result_decision) == 1L && !truth(result_decision$hca_fastq_download_authorized) &&
    !truth(result_decision$raw_reprocessing_authorized) &&
    !truth(result_decision$label_access_authorized)
}
if (file.exists(ledger_path)) {
  metric <- read.csv(ledger_path, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(metric) != 1L || metric$disposition != "completed" ||
      !validate_result(result) || metric$result_manifest_sha256 !=
        sha(file.path(result, "mv08g-artifact-manifest.csv"))) {
    stop("MV8-G comparison resume state is stale.")
  }
} else {
  if (dir.exists(result) && length(list.files(result, all.files = TRUE,
                                               no.. = TRUE))) {
    stop("Unowned or partial MV8-G comparison output exists; refusing retry.")
  }
  started <- Sys.time()
  process <- processx::process$new(Sys.which("Rscript"), c(
    "--vanilla", "scripts/run_mv08g_comparison.R", primary, root500, root475,
    shift_root, result), stdout = file.path(log_dir, "comparison__stdout.txt"),
    stderr = file.path(log_dir, "comparison__stderr.txt"), cleanup_tree = TRUE)
  peak <- 0; cap_failure <- ""
  while (process$is_alive()) {
    Sys.sleep(0.25); peak <- max(peak, tree_rss(process$get_pid()))
    elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    if (elapsed > contract$elapsed_cap_seconds) {
      cap_failure <- "elapsed_cap_exceeded"; process$kill_tree()
    } else if (peak > contract$rss_cap_bytes) {
      cap_failure <- "rss_cap_exceeded"; process$kill_tree()
    }
  }
  process$wait(timeout = 5000)
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  status <- process$get_exit_status()
  valid <- identical(status, 0L) && validate_result(result)
  bytes <- if (dir.exists(result)) sum(as.numeric(file.info(list.files(
    result, full.names = TRUE))$size), na.rm = TRUE) else 0
  disposition <- if (nzchar(cap_failure)) cap_failure else if (
    bytes > contract$storage_cap_bytes) "storage_cap_exceeded" else if (valid)
      "completed" else "failed"
  metric <- data.frame(
    contract_id = "mv08g_comparison_resource_metric_v1", job_id = "comparison",
    disposition = disposition,
    exit_status = if (is.null(status)) NA_integer_ else status,
    elapsed_seconds = elapsed, peak_process_tree_rss_bytes = peak,
    elapsed_cap_seconds = contract$elapsed_cap_seconds,
    rss_cap_bytes = contract$rss_cap_bytes, result_bytes = bytes,
    storage_cap_bytes = contract$storage_cap_bytes,
    result_manifest_sha256 = if (valid) sha(file.path(result,
      "mv08g-artifact-manifest.csv")) else NA_character_,
    workers = 1L, retries = 0L, outcome_label_state = "closed",
    biological_outcomes_computed = FALSE, stringsAsFactors = FALSE)
  write_atomic(metric, ledger_path)
  if (disposition != "completed") stop("MV8-G comparison stopped under zero-retry: ",
                                        disposition)
}
public_decision <- data.frame(
  contract_id = "mv08g_comparison_execution_decision_v1",
  decision = "comparison_complete_await_independent_reconstruction",
  comparison_jobs = 1L, elapsed_seconds = metric$elapsed_seconds,
  peak_process_tree_rss_bytes = metric$peak_process_tree_rss_bytes,
  result_bytes = metric$result_bytes, hca_fastq_download_authorized = FALSE,
  raw_reprocessing_authorized = FALSE, label_access_authorized = FALSE,
  biological_outcomes_computed = FALSE, stringsAsFactors = FALSE)
write_atomic(metric, file.path(public_dir, "mv08g-comparison-resource-metric.csv"))
write_atomic(public_decision, file.path(public_dir, "mv08g-comparison-decision.csv"))
message("MV8-G comparison complete; awaiting independent reconstruction")
