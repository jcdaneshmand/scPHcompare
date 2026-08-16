#!/usr/bin/env Rscript

options(warn = 2)
for (package in c("processx", "ps", "digest")) {
  if (!requireNamespace(package, quietly = TRUE)) stop(package, " required.")
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop("usage: run_mv07i_label_closed_monitor.R MV7I_PREFREEZE ",
       "MV7H_PREFREEZE MV7H_COMPLETE_VALIDATION LANDSCAPE_ROOT PRIVATE_ROOT ",
       "METRICS", call. = FALSE)
}
read_csv <- function(path) utils::read.csv(
  path, stringsAsFactors = FALSE, check.names = FALSE)
mv07i <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
mv07h <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
mv07h_validation <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
landscape_root <- normalizePath(args[[4L]], winslash = "/", mustWork = TRUE)
private_root <- args[[5L]]
metrics_path <- args[[6L]]
resource <- read_csv(file.path(mv07i, "mv07i-resource-resume-contract.csv"))
decision <- read_csv(file.path(mv07i, "mv07i-decision.csv"))
if (resource$maximum_workers != 1L ||
    decision$decision !=
      "authorize_label_closed_matrix_and_clustering_production_only" ||
    as.logical(decision$labels_authorized) || as.logical(decision$outcomes_authorized)) {
  stop("MV7-I monitor admission is stale.", call. = FALSE)
}
if (file.exists(metrics_path)) {
  metric <- read_csv(metrics_path)
  if (nrow(metric) != 1L ||
      !metric$disposition %in% c("completed", "reused_validated") ||
      metric$exit_status != 0L) {
    stop("MV7-I monitor preserves a prior failure and refuses retry.")
  }
  status_path <- file.path(private_root, "artifacts", "status.csv")
  if (!file.exists(status_path) ||
      read_csv(status_path)$completion_state != "complete") {
    stop("MV7-I successful resource checkpoint has stale artifacts.")
  }
  status <- read_csv(status_path)
  named_files <- c(
    matrix_bundle = "matrix-bundle.rds",
    pair_summary = "pair-seed-summary.csv",
    h1_summary = "h1-contribution-summary.csv",
    candidate_partitions = "candidate-pam-partitions.csv",
    stability = "stability-summary.csv",
    selected_partitions = "selected-partitions.csv",
    provenance = "provenance.csv")
  paths <- file.path(private_root, "artifacts", named_files)
  expected_hashes <- unname(unlist(status[paste0(
    names(named_files), "_sha256")], use.names = FALSE))
  if (!all(file.exists(paths))) {
    stop("MV7-I successful resource checkpoint has missing artifacts.")
  }
  observed_hashes <- vapply(paths, function(path) {
    digest::digest(file = path, algo = "sha256", serialize = FALSE)
  }, character(1L))
  if (!identical(tolower(unname(observed_hashes)),
                 tolower(expected_hashes))) {
    stop("MV7-I successful resource checkpoint has hash-drifted artifacts.")
  }
  message("Reused immutable MV7-I resource checkpoint and artifacts.")
  quit(save = "no", status = 0L)
}
dir.create(private_root, recursive = TRUE, showWarnings = FALSE)
artifact_dir <- file.path(private_root, "artifacts")
log_dir <- file.path(private_root, "logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
preexisting <- dir.exists(artifact_dir)
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
    "--vanilla", "scripts/run_mv07i_label_closed_entry.R", mv07i, mv07h,
    mv07h_validation, landscape_root, artifact_dir),
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
status_code <- process$get_exit_status()
status_path <- file.path(artifact_dir, "status.csv")
complete <- is.na(cap_failure) && identical(status_code, 0L) &&
  file.exists(status_path) &&
  read_csv(status_path)$completion_state == "complete"
disposition <- if (!is.na(cap_failure)) cap_failure else if (complete) {
  if (preexisting) "reused_validated" else "completed"
} else "failed"
metric <- data.frame(
  contract_id = "mv07i_label_closed_resource_v1", disposition = disposition,
  exit_status = if (is.null(status_code)) NA_integer_ else status_code,
  elapsed_seconds = elapsed, peak_process_tree_rss_bytes = peak,
  elapsed_cap_seconds = resource$elapsed_cap_seconds,
  rss_cap_bytes = resource$process_tree_rss_cap_bytes,
  artifact_bytes = if (dir.exists(artifact_dir)) sum(file.info(list.files(
    artifact_dir, recursive = TRUE, full.names = TRUE))$size, na.rm = TRUE) else 0,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  labels_loaded = FALSE, outcomes_computed = FALSE, stringsAsFactors = FALSE)
dir.create(dirname(metrics_path), recursive = TRUE, showWarnings = FALSE)
temporary <- paste0(metrics_path, ".partial.", Sys.getpid())
write.csv(metric, temporary, row.names = FALSE, na = "")
if (file.exists(metrics_path) && !file.remove(metrics_path)) {
  stop("MV7-I resource checkpoint replacement failed.")
}
if (!file.rename(temporary, metrics_path)) {
  stop("MV7-I resource checkpoint publication failed.")
}
if (!complete) stop("MV7-I label-closed production failed: ", disposition)
message("MV7-I label-closed production completed under resource gates.")
