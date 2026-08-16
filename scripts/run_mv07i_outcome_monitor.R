#!/usr/bin/env Rscript

options(warn = 2)
for (package in c("processx", "ps", "digest")) {
  if (!requireNamespace(package, quietly = TRUE)) stop(package, " required.")
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop("usage: run_mv07i_outcome_monitor.R PREFREEZE MV7D MV7E SELECTED ",
       "PRIVATE_ROOT METRICS", call. = FALSE)
}
read_csv <- function(path) utils::read.csv(
  path, stringsAsFactors = FALSE, check.names = FALSE)
prefreeze <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
mv07d <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
mv07e <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
selected <- normalizePath(args[[4L]], winslash = "/", mustWork = TRUE)
private_root <- args[[5L]]
metrics_path <- args[[6L]]
resource <- read_csv(file.path(prefreeze, "mv07i-outcome-resource-resume.csv"))
decision <- read_csv(file.path(prefreeze, "mv07i-outcome-decision.csv"))
if (resource$maximum_workers != 1L || decision$decision !=
      "authorize_MV7I_descriptive_outcome_execution_only") {
  stop("MV7-I outcome monitor admission is stale.", call. = FALSE)
}
named_files <- c(
  metadata_join = "metadata-join.csv", contingency = "contingency-long.csv",
  seed_metrics = "seed-metrics.csv", unit_summaries = "unit-summaries.csv",
  structural_status = "structural-status.csv", provenance = "provenance.csv")
if (file.exists(metrics_path)) {
  metric <- read_csv(metrics_path)
  if (nrow(metric) != 1L || metric$disposition != "completed" ||
      metric$exit_status != 0L) {
    stop("MV7-I outcome monitor preserves a prior failure and refuses retry.")
  }
  artifact_dir <- file.path(private_root, "artifacts")
  status_path <- file.path(artifact_dir, "status.csv")
  paths <- stats::setNames(file.path(artifact_dir, named_files),
                           names(named_files))
  if (!file.exists(status_path) || !all(file.exists(paths))) {
    stop("MV7-I outcome successful checkpoint has missing artifacts.")
  }
  status <- read_csv(status_path)
  expected <- unname(unlist(status[paste0(
    names(named_files), "_sha256")], use.names = FALSE))
  observed <- unname(vapply(paths, function(path) digest::digest(
    file = path, algo = "sha256", serialize = FALSE), character(1L)))
  if (status$completion_state != "complete" ||
      !identical(tolower(observed), tolower(expected))) {
    stop("MV7-I outcome successful checkpoint has hash-drifted artifacts.")
  }
  message("Reused immutable MV7-I outcome checkpoint and artifacts.")
  quit(save = "no", status = 0L)
}
dir.create(private_root, recursive = TRUE, showWarnings = FALSE)
artifact_dir <- file.path(private_root, "artifacts")
log_dir <- file.path(private_root, "logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
tree_rss <- function(pid) {
  root <- tryCatch(ps::ps_handle(pid), error = function(error) NULL)
  if (is.null(root)) return(0)
  handles <- c(list(root), tryCatch(ps::ps_children(root, recursive = TRUE),
                                   error = function(error) list()))
  sum(vapply(handles, function(handle) tryCatch(
    as.numeric(ps::ps_memory_info(handle)[["rss"]]),
    error = function(error) 0), numeric(1L)))
}
started <- Sys.time()
process <- processx::process$new(
  command = Sys.which("Rscript"), args = c(
    "--vanilla", "scripts/run_mv07i_outcome_entry.R", prefreeze, mv07d, mv07e,
    selected, artifact_dir),
  stdout = file.path(log_dir, "stdout.txt"),
  stderr = file.path(log_dir, "stderr.txt"), cleanup_tree = TRUE)
peak <- 0; cap_failure <- NA_character_
while (process$is_alive()) {
  Sys.sleep(0.25)
  peak <- max(peak, tree_rss(process$get_pid()))
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  if (elapsed > resource$elapsed_cap_seconds) {
    cap_failure <- "elapsed_cap_exceeded"; process$kill_tree()
  } else if (peak > resource$process_tree_rss_cap_bytes) {
    cap_failure <- "rss_cap_exceeded"; process$kill_tree()
  }
}
process$wait(timeout = 5000)
elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
exit_status <- process$get_exit_status()
status_path <- file.path(artifact_dir, "status.csv")
complete <- is.na(cap_failure) && identical(exit_status, 0L) &&
  file.exists(status_path) && read_csv(status_path)$completion_state == "complete"
disposition <- if (!is.na(cap_failure)) cap_failure else if (complete)
  "completed" else "failed"
artifact_paths <- if (dir.exists(artifact_dir)) list.files(
  artifact_dir, recursive = TRUE, full.names = TRUE) else character()
metric <- data.frame(
  contract_id = "mv07i_outcome_resource_v1", disposition = disposition,
  exit_status = if (is.null(exit_status)) NA_integer_ else exit_status,
  elapsed_seconds = elapsed, peak_process_tree_rss_bytes = peak,
  elapsed_cap_seconds = resource$elapsed_cap_seconds,
  rss_cap_bytes = resource$process_tree_rss_cap_bytes,
  artifact_bytes = sum(file.info(artifact_paths)$size, na.rm = TRUE),
  p_values_computed = FALSE, method_selection_executed = FALSE,
  outcomes_computed = complete, stringsAsFactors = FALSE)
temporary <- paste0(metrics_path, ".partial.", Sys.getpid())
write.csv(metric, temporary, row.names = FALSE, na = "")
if (file.exists(metrics_path) && !file.remove(metrics_path)) {
  stop("MV7-I outcome resource checkpoint replacement failed.")
}
if (!file.rename(temporary, metrics_path)) {
  stop("MV7-I outcome resource checkpoint publication failed.")
}
if (!complete) stop("MV7-I descriptive outcome production failed: ", disposition)
message("MV7-I descriptive outcome production completed under resource gates.")
