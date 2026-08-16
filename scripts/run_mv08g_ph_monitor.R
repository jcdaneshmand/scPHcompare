#!/usr/bin/env Rscript

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1")
options(warn = 2)
for (package in c("digest", "processx", "ps")) {
  if (!requireNamespace(package, quietly = TRUE)) stop(package, " required.")
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop("usage: run_mv08g_ph_monitor.R PREFREEZE SOURCE_ROOT PRIVATE_ROOT PUBLIC_DIR EXPECTED_HEAD PRIMARY_PREFREEZE")
}
prefreeze <- args[[1L]]; source_root <- args[[2L]]; private_root <- args[[3L]]
public_dir <- args[[4L]]; expected_head <- tolower(trimws(args[[5L]]))
primary_prefreeze <- args[[6L]]
head <- tolower(trimws(system2("git", c("rev-parse", "HEAD"), stdout = TRUE)))
if (head != expected_head) stop("MV8-G PH monitor exact HEAD mismatch.")
source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/mv07g_sentinel.R")
source("R/mv08g_panel_sensitivity.R")
queue <- read.csv(file.path(prefreeze, "mv08g-ph-queue.csv"),
                  stringsAsFactors = FALSE, check.names = FALSE)
repeat_queue <- read.csv(file.path(prefreeze, "mv08g-ph-repeat-queue.csv"),
                         stringsAsFactors = FALSE, check.names = FALSE)
fallback <- read.csv(file.path(prefreeze, "mv08g-ph-fallback-policy.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)
decision <- read.csv(file.path(prefreeze, "mv08g-ph-decision.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)
primary_contract <- read.csv(file.path(primary_prefreeze, "mv08g-contract.csv"),
                             stringsAsFactors = FALSE, check.names = FALSE)
ph_contract <- read.csv(file.path(prefreeze, "mv08g-ph-contract.csv"),
                        stringsAsFactors = FALSE, check.names = FALSE)
if (nrow(queue) != 1240L || nrow(repeat_queue) != 4L || nrow(fallback) != 1L ||
    decision$decision != "authorize_1240_PH_jobs_and_four_selected_engine_repeats" ||
    decision$ph_jobs_authorized != 1240L ||
    primary_contract$panel_genes != 475L || nrow(ph_contract) != 1L ||
    ph_contract$aggregate_elapsed_cap_seconds != 172800 ||
    ph_contract$aggregate_storage_cap_bytes != 4 * 1024^3 ||
    any(queue$workers != 1L) ||
    any(queue$retries != 0L)) stop("MV8-G PH execution gate is stale.")
for (subdir in c("ph", "logs", "repeat/ph", "repeat/logs")) {
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
read_csv_or_empty <- function(path) if (file.exists(path)) read.csv(
  path, stringsAsFactors = FALSE, check.names = FALSE) else data.frame()
write_atomic <- function(value, path) {
  partial <- tempfile(pattern = basename(path), tmpdir = dirname(path))
  write.csv(value, partial, row.names = FALSE, na = "")
  if (!file.rename(partial, path)) {
    unlink(partial); stop("Failed to atomically publish MV8-G PH ledger.")
  }
}
private_storage_bytes <- function() {
  files <- c(list.files(file.path(private_root, "ph"), full.names = TRUE),
             list.files(file.path(private_root, "repeat/ph"), full.names = TRUE))
  if (length(files)) sum(as.numeric(file.info(files)$size), na.rm = TRUE) else 0
}
.source_cache <- new.env(parent = emptyenv())
source_for_seed <- function(seed) {
  key <- as.character(seed)
  if (!identical(.source_cache$key, key)) {
    path <- file.path(source_root, paste0("mv08g__", seed, "__source.rds"))
    record <- readRDS(path); mv08g_validate_source_record_v1(record)
    .source_cache$key <- key; .source_cache$record <- record
  }
  .source_cache$record
}
validate_output <- function(path, row) {
  source_record <- source_for_seed(row$seed)
  record <- readRDS(path)
  view <- source_record$views[[row$sample_id]][[row$view_id]]
  mv08g_validate_ph_record_v1(record, view)
  record$identity$seed == row$seed &&
    record$identity$sample_id == row$sample_id &&
    record$identity$view_id == row$view_id
}
attempts_path <- file.path(private_root, "ph-engine-attempts.csv")
selected_path <- file.path(private_root, "ph-selected-metrics.csv")
repeat_attempts_path <- file.path(private_root, "repeat-ph-engine-attempts.csv")
attempts <- read_csv_or_empty(attempts_path)
selected <- read_csv_or_empty(selected_path)
repeat_attempts <- read_csv_or_empty(repeat_attempts_path)

run_attempt <- function(row, engine, repeat_mode = FALSE) {
  engine <- match.arg(engine, c("ripserr", "TDA_ripsDiag_GUDHI"))
  output <- file.path(private_root, if (repeat_mode) "repeat/ph" else "ph",
    paste0("mv08g__", row$seed, "__", row$sample_id, "__", row$view_id,
           "__ph.rds"))
  attempt_id <- paste(if (repeat_mode) "repeat" else "production",
                      row$job_id, engine, sep = "__")
  ledger <- if (repeat_mode) repeat_attempts else attempts
  hit <- if (nrow(ledger)) ledger[ledger$attempt_id == attempt_id, , drop = FALSE]
    else data.frame()
  if (nrow(hit)) {
    completed <- hit$disposition == "completed" && file.exists(output) &&
      hit$output_sha256 == sha(output) &&
      tryCatch(validate_output(output, row), error = function(error) FALSE)
    eligible_rss <- !repeat_mode && engine == "ripserr" &&
      row$view_id == "gene_topology_v1" &&
      hit$disposition == "rss_cap_exceeded" && !file.exists(output)
    if (completed || eligible_rss) return(hit)
    stop("MV8-G PH prior attempt is not resumable: ", attempt_id)
  }
  if (file.exists(output)) {
    stop("Unowned MV8-G PH output exists; refusing adoption: ", basename(output))
  }
  is_fallback <- engine == "TDA_ripsDiag_GUDHI"
  if (is_fallback && row$view_id != "gene_topology_v1") {
    stop("MV8-G exact fallback is gene-view only.")
  }
  script <- if (is_fallback) "scripts/run_mv08g_ph_fallback_entry.R" else
    "scripts/run_mv08g_ph_entry.R"
  source_path <- file.path(source_root, paste0("mv08g__", row$seed,
                                               "__source.rds"))
  script_args <- if (is_fallback) c(source_path, row$sample_id, output) else
    c(source_path, row$sample_id, row$view_id, output)
  cap_elapsed <- if (is_fallback) fallback$fallback_elapsed_cap_seconds else
    row$elapsed_cap_seconds
  cap_rss <- if (is_fallback) fallback$fallback_rss_cap_bytes else row$rss_cap_bytes
  log_dir <- file.path(private_root, if (repeat_mode) "repeat/logs" else "logs")
  stem <- gsub("[^A-Za-z0-9_.-]", "_", attempt_id)
  started <- Sys.time()
  process <- processx::process$new(
    Sys.which("Rscript"), c("--vanilla", script, script_args),
    stdout = file.path(log_dir, paste0(stem, "__stdout.txt")),
    stderr = file.path(log_dir, paste0(stem, "__stderr.txt")),
    cleanup_tree = TRUE)
  peak <- 0; cap_failure <- ""
  while (process$is_alive()) {
    Sys.sleep(0.25); peak <- max(peak, tree_rss(process$get_pid()))
    elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    if (elapsed > cap_elapsed) {
      cap_failure <- "elapsed_cap_exceeded"; process$kill_tree()
    } else if (peak > cap_rss) {
      cap_failure <- "rss_cap_exceeded"; process$kill_tree()
    }
  }
  process$wait(timeout = 5000)
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  status <- process$get_exit_status()
  valid <- identical(status, 0L) && file.exists(output) &&
    tryCatch(validate_output(output, row), error = function(error) FALSE)
  disposition <- if (nzchar(cap_failure)) cap_failure else if (valid)
    "completed" else "failed"
  metric <- data.frame(
    contract_id = "mv08g_ph_engine_attempt_v1", attempt_id = attempt_id,
    job_id = row$job_id, seed = row$seed, sample_id = row$sample_id,
    view_id = row$view_id, repeat_mode = repeat_mode, ph_engine = engine,
    disposition = disposition,
    exit_status = if (is.null(status)) NA_integer_ else status,
    elapsed_seconds = elapsed, peak_process_tree_rss_bytes = peak,
    elapsed_cap_seconds = cap_elapsed, rss_cap_bytes = cap_rss,
    output_file = basename(output),
    output_bytes = if (file.exists(output)) as.numeric(file.info(output)$size) else NA,
    output_sha256 = if (file.exists(output)) sha(output) else NA_character_,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE)
  ledger <- if (nrow(ledger)) rbind(ledger, metric) else metric
  if (repeat_mode) {
    repeat_attempts <<- ledger; write_atomic(ledger, repeat_attempts_path)
  } else {
    attempts <<- ledger; write_atomic(ledger, attempts_path)
  }
  metric
}
run_production <- function(row) {
  if (sum(attempts$elapsed_seconds, na.rm = TRUE) +
      sum(repeat_attempts$elapsed_seconds, na.rm = TRUE) >
      ph_contract$aggregate_elapsed_cap_seconds ||
      private_storage_bytes() > ph_contract$aggregate_storage_cap_bytes) {
    stop("MV8-G PH aggregate resource cap reached before job ", row$job_id)
  }
  hit <- if (nrow(selected)) selected[selected$job_id == row$job_id, , drop = FALSE]
    else data.frame()
  output <- file.path(private_root, "ph", paste0("mv08g__", row$seed, "__",
    row$sample_id, "__", row$view_id, "__ph.rds"))
  if (nrow(hit)) {
    if (nrow(hit) != 1L || hit$disposition != "completed" || !file.exists(output) ||
        hit$output_sha256 != sha(output) || !validate_output(output, row)) {
      stop("MV8-G selected PH resume state is stale: ", row$job_id)
    }
    return(hit)
  }
  primary <- run_attempt(row, "ripserr", FALSE)
  chosen <- if (primary$disposition == "completed") primary else {
    if (row$view_id != "gene_topology_v1" ||
        primary$disposition != "rss_cap_exceeded") {
      stop("MV8-G primary PH stopped under zero-retry policy: ", row$job_id,
           " (", primary$disposition, ")")
    }
    value <- run_attempt(row, "TDA_ripsDiag_GUDHI", FALSE)
    if (value$disposition != "completed") stop("MV8-G exact fallback failed: ",
                                                row$job_id)
    value
  }
  chosen$contract_id <- "mv08g_selected_ph_metric_v1"
  selected <<- if (nrow(selected)) rbind(selected, chosen) else chosen
  write_atomic(selected, selected_path)
  chosen
}
for (index in seq_len(nrow(queue))) run_production(queue[index, , drop = FALSE])
selected <- selected[match(queue$job_id, selected$job_id), , drop = FALSE]
if (anyNA(selected$job_id) || nrow(selected) != 1240L) {
  stop("MV8-G selected PH ledger is incomplete.")
}
repeat_rows <- list()
for (index in seq_len(nrow(repeat_queue))) {
  if (sum(attempts$elapsed_seconds, na.rm = TRUE) +
      sum(repeat_attempts$elapsed_seconds, na.rm = TRUE) >
      ph_contract$aggregate_elapsed_cap_seconds ||
      private_storage_bytes() > ph_contract$aggregate_storage_cap_bytes) {
    stop("MV8-G PH aggregate resource cap reached before repeat ", index)
  }
  spec <- repeat_queue[index, , drop = FALSE]
  row <- queue[queue$seed == spec$seed & queue$sample_id == spec$sample_id &
    queue$view_id == spec$view_id, , drop = FALSE]
  production <- selected[selected$job_id == row$job_id, , drop = FALSE]
  repeated <- run_attempt(row, production$ph_engine, TRUE)
  exact <- repeated$disposition == "completed" &&
    repeated$output_sha256 == production$output_sha256 &&
    repeated$output_bytes == production$output_bytes
  if (!exact) stop("MV8-G PH repeat is not byte-exact: ", row$job_id)
  repeat_rows[[index]] <- data.frame(
    contract_id = "mv08g_ph_repeat_validation_v1", repeat_order = spec$repeat_order,
    job_id = row$job_id, seed = row$seed, sample_id = row$sample_id,
    view_id = row$view_id, selected_engine = production$ph_engine,
    production_sha256 = production$output_sha256,
    repeat_sha256 = repeated$output_sha256, byte_identical = exact,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE)
}
repeats <- do.call(rbind, repeat_rows)
aggregate_elapsed <- sum(attempts$elapsed_seconds, na.rm = TRUE) +
  sum(repeat_attempts$elapsed_seconds, na.rm = TRUE)
if (aggregate_elapsed > ph_contract$aggregate_elapsed_cap_seconds ||
    private_storage_bytes() > ph_contract$aggregate_storage_cap_bytes) {
  stop("MV8-G PH aggregate resource cap exceeded at completion.")
}
public_attempts <- attempts
public_selected <- selected
decision_out <- data.frame(
  contract_id = "mv08g_ph_execution_decision_v1",
  decision = "full_PH_complete_await_independent_validation",
  ph_jobs = nrow(selected), ph_repeat_jobs = nrow(repeats),
  ripserr_selected_records = sum(selected$ph_engine == "ripserr"),
  gudhi_fallback_records = sum(selected$ph_engine == "TDA_ripsDiag_GUDHI"),
  rss_triggered_fallbacks = sum(attempts$disposition == "rss_cap_exceeded"),
  aggregate_elapsed_seconds = aggregate_elapsed,
  maximum_process_tree_rss_bytes = max(c(attempts$peak_process_tree_rss_bytes,
    repeat_attempts$peak_process_tree_rss_bytes)), landscape_jobs_authorized = 0L,
  private_storage_bytes = private_storage_bytes(),
  aggregate_storage_cap_bytes = ph_contract$aggregate_storage_cap_bytes,
  aggregate_elapsed_cap_seconds = ph_contract$aggregate_elapsed_cap_seconds,
  hca_fastq_download_authorized = FALSE, raw_reprocessing_authorized = FALSE,
  label_access_authorized = FALSE, biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE)
write_atomic(public_selected, file.path(public_dir, "mv08g-ph-metrics.csv"))
write_atomic(public_attempts, file.path(public_dir, "mv08g-ph-engine-attempts.csv"))
write_atomic(repeats, file.path(public_dir, "mv08g-ph-repeat-validation.csv"))
write_atomic(decision_out, file.path(public_dir, "mv08g-ph-decision.csv"))
message("MV8-G PH complete: 1,240 records plus four exact repeats; awaiting validation")
