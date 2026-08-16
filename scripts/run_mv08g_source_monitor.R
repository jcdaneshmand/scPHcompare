#!/usr/bin/env Rscript

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1")
options(warn = 2)
for (package in c("digest", "processx", "ps")) {
  if (!requireNamespace(package, quietly = TRUE)) stop(package, " required.")
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) {
  stop("usage: run_mv08g_source_monitor.R PREFREEZE PRIMARY_CACHE ADDED_CACHE PRIVATE_ROOT PUBLIC_DIR EXPECTED_HEAD REPEAT_SEED")
}
prefreeze <- args[[1L]]; primary_cache <- args[[2L]]; added_cache <- args[[3L]]
private_root <- args[[4L]]; public_dir <- args[[5L]]
expected_head <- tolower(trimws(args[[6L]])); repeat_seed <- as.integer(args[[7L]])
head <- tolower(trimws(system2("git", c("rev-parse", "HEAD"), stdout = TRUE)))
if (head != expected_head) stop("MV8-G source monitor exact HEAD mismatch.")
source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/mv07g_sentinel.R")
source("R/mv08g_panel_sensitivity.R")
queue <- read.csv(file.path(prefreeze, "mv08g-source-queue.csv"),
                  stringsAsFactors = FALSE, check.names = FALSE)
decision <- read.csv(file.path(prefreeze, "mv08g-decision.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)
if (nrow(queue) != 5L || any(queue$workers != 1L) || any(queue$retries != 0L) ||
    !identical(sort(as.integer(queue$seed)), .mv08g_seeds) ||
    decision$decision != "authorize_five_common475_source_bundles_and_one_repeat" ||
    decision$source_jobs_authorized != 5L ||
    decision$source_repeat_jobs_authorized != 1L ||
    !(repeat_seed %in% .mv08g_seeds)) {
  stop("MV8-G source execution differs from the prospective gate.")
}
for (subdir in c("source", "logs", "repeat/source", "repeat/logs")) {
  dir.create(file.path(private_root, subdir), recursive = TRUE,
             showWarnings = FALSE)
}
dir.create(public_dir, recursive = TRUE, showWarnings = FALSE)
sha <- function(path) digest::digest(file = path, algo = "sha256",
                                     serialize = FALSE)
tree_rss <- function(pid) {
  root <- tryCatch(ps::ps_handle(pid), error = function(error) NULL)
  if (is.null(root)) return(0)
  handles <- c(list(root), tryCatch(ps::ps_children(root, recursive = TRUE),
                                   error = function(error) list()))
  sum(vapply(handles, function(handle) tryCatch(
    as.numeric(ps::ps_memory_info(handle)[["rss"]]),
    error = function(error) 0), numeric(1L)))
}
read_ledger <- function(path) if (file.exists(path)) read.csv(
  path, stringsAsFactors = FALSE, check.names = FALSE) else data.frame()
write_atomic_csv <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  partial <- tempfile(pattern = basename(path), tmpdir = dirname(path))
  write.csv(value, partial, row.names = FALSE, na = "")
  if (!file.rename(partial, path)) {
    unlink(partial); stop("Failed to atomically publish MV8-G source ledger.")
  }
}
validate_output <- function(path, seed) {
  record <- readRDS(path)
  mv08g_validate_source_record_v1(record)
  record$identity$seed == seed && record$identity$panel_size == 475L
}
run_unit <- function(row, repeat_mode = FALSE) {
  prefix <- if (repeat_mode) "repeat" else "production"
  output <- file.path(private_root, if (repeat_mode) "repeat/source" else "source",
    paste0("mv08g__", row$seed, "__source.rds"))
  ledger_path <- file.path(private_root,
    if (repeat_mode) "repeat-source-metrics.csv" else "source-metrics.csv")
  ledger <- read_ledger(ledger_path)
  job_id <- paste0(prefix, "__", row$job_id)
  hit <- if (nrow(ledger)) ledger[ledger$job_id == job_id, , drop = FALSE] else
    data.frame()
  if (nrow(hit)) {
    if (nrow(hit) != 1L || hit$disposition != "completed" ||
        !file.exists(output) || hit$output_sha256 != sha(output) ||
        hit$output_bytes != as.numeric(file.info(output)$size) ||
        !validate_output(output, row$seed)) {
      stop("MV8-G source resume state is incomplete or stale: ", job_id)
    }
    return(hit)
  }
  if (file.exists(output)) {
    stop("Unowned MV8-G source output exists; refusing adoption: ", basename(output))
  }
  log_dir <- file.path(private_root, if (repeat_mode) "repeat/logs" else "logs")
  stem <- gsub("[^A-Za-z0-9_.-]", "_", job_id)
  stdout <- file.path(log_dir, paste0(stem, "__stdout.txt"))
  stderr <- file.path(log_dir, paste0(stem, "__stderr.txt"))
  started <- Sys.time()
  process <- processx::process$new(
    command = Sys.which("Rscript"),
    args = c("--vanilla", "scripts/run_mv08g_source_entry.R", prefreeze,
      primary_cache, added_cache, output, as.character(row$seed)),
    stdout = stdout, stderr = stderr, cleanup_tree = TRUE)
  peak <- 0; cap_failure <- ""
  while (process$is_alive()) {
    Sys.sleep(0.25)
    peak <- max(peak, tree_rss(process$get_pid()))
    elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    if (elapsed > row$elapsed_cap_seconds) {
      cap_failure <- "elapsed_cap_exceeded"; process$kill_tree()
    } else if (peak > row$rss_cap_bytes) {
      cap_failure <- "rss_cap_exceeded"; process$kill_tree()
    }
  }
  process$wait(timeout = 5000)
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  status <- process$get_exit_status()
  valid <- identical(status, 0L) && file.exists(output) &&
    tryCatch(validate_output(output, row$seed), error = function(error) FALSE)
  disposition <- if (nzchar(cap_failure)) cap_failure else if (valid)
    "completed" else "failed"
  metric <- data.frame(
    contract_id = "mv08g_source_resource_metric_v1", job_id = job_id,
    seed = row$seed, repeat_mode = repeat_mode, disposition = disposition,
    exit_status = if (is.null(status)) NA_integer_ else status,
    elapsed_seconds = elapsed, peak_process_tree_rss_bytes = peak,
    output_file = basename(output),
    output_bytes = if (file.exists(output)) as.numeric(file.info(output)$size) else NA,
    output_sha256 = if (file.exists(output)) sha(output) else NA_character_,
    elapsed_cap_seconds = row$elapsed_cap_seconds,
    rss_cap_bytes = row$rss_cap_bytes, workers = 1L, retries = 0L,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE)
  ledger <- if (nrow(ledger)) rbind(ledger, metric) else metric
  write_atomic_csv(ledger, ledger_path)
  if (disposition != "completed") {
    stop("MV8-G source job stopped under the zero-retry policy: ", job_id,
         " (", disposition, ")")
  }
  metric
}
production <- do.call(rbind, lapply(seq_len(nrow(queue)), function(index) {
  run_unit(queue[index, , drop = FALSE], FALSE)
}))
repeat_row <- queue[queue$seed == repeat_seed, , drop = FALSE]
repeat_metric <- run_unit(repeat_row, TRUE)
production_row <- production[production$seed == repeat_seed, , drop = FALSE]
repeat_validation <- data.frame(
  contract_id = "mv08g_source_repeat_validation_v1", seed = repeat_seed,
  production_sha256 = production_row$output_sha256,
  repeat_sha256 = repeat_metric$output_sha256,
  byte_identical = production_row$output_sha256 == repeat_metric$output_sha256 &&
    production_row$output_bytes == repeat_metric$output_bytes,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE)
if (!repeat_validation$byte_identical) stop("MV8-G source repeat is not exact.")
public_metrics <- production
public_metrics$output_file <- basename(public_metrics$output_file)
public_decision <- data.frame(
  contract_id = "mv08g_source_execution_decision_v1",
  decision = "source_complete_await_independent_validation",
  source_jobs = nrow(production), source_repeat_jobs = 1L,
  typed_views = 1240L, aggregate_elapsed_seconds = sum(production$elapsed_seconds) +
    repeat_metric$elapsed_seconds,
  maximum_process_tree_rss_bytes = max(c(production$peak_process_tree_rss_bytes,
    repeat_metric$peak_process_tree_rss_bytes)), ph_jobs_authorized = 0L,
  hca_fastq_download_authorized = FALSE, raw_reprocessing_authorized = FALSE,
  label_access_authorized = FALSE, biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE)
write_atomic_csv(public_metrics, file.path(public_dir, "mv08g-source-metrics.csv"))
write_atomic_csv(repeat_validation,
  file.path(public_dir, "mv08g-source-repeat-validation.csv"))
write_atomic_csv(public_decision, file.path(public_dir, "mv08g-source-decision.csv"))
message("MV8-G source execution complete: five bundles plus exact seed ",
        repeat_seed, " repeat; awaiting independent validation")
