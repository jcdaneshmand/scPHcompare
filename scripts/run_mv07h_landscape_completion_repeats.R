#!/usr/bin/env Rscript

options(warn = 2)
for (package in c("processx", "ps", "digest")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(package, " is required for MV7-H landscape repeats.",
         call. = FALSE)
  }
}
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1", RAYON_NUM_THREADS = "1")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) {
  stop("usage: run_mv07h_landscape_completion_repeats.R COMPLETION_PREFREEZE ",
       "PREFREEZE PH_ROOT RUST_LIBRARY PRODUCTION_ROOT REPEAT_ROOT METRICS",
       call. = FALSE)
}
source("R/mv07h_full_topology.R")
read_csv <- function(path) utils::read.csv(
  path, stringsAsFactors = FALSE, check.names = FALSE)
completion_prefreeze <- normalizePath(args[[1L]], winslash = "/",
                                      mustWork = TRUE)
prefreeze <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
ph_root <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
rust <- normalizePath(args[[4L]], winslash = "/", mustWork = TRUE)
production_root <- normalizePath(args[[5L]], winslash = "/", mustWork = TRUE)
repeat_root <- args[[6L]]
metrics_path <- args[[7L]]
selection <- read_csv(file.path(
  completion_prefreeze, "mv07h-landscape-completion-repeat-selection.csv"))
completion_contract <- read_csv(file.path(
  completion_prefreeze, "mv07h-landscape-completion-contract.csv"))
inventory <- read_csv(file.path(
  completion_prefreeze, "mv07h-landscape-completion-input-inventory.csv"))
contract <- read_csv(file.path(prefreeze, "mv07h-contract.csv"))
if (nrow(selection) != 3L || any(selection$stage != "stage_2") ||
    any(selection$workers != 1L) || any(selection$retries != 0L) ||
    any(as.logical(selection$distance_values_used_for_selection)) ||
    completion_contract$independent_repeat_groups != 3L ||
    completion_contract$parent_implementation_root_sha256 !=
      contract$implementation_root_sha256 ||
    completion_contract$rust_library_sha256 != .mv07h_sha256(rust) ||
    any(selection$outcome_label_state != "closed") ||
    any(as.logical(selection$biological_outcomes_computed))) {
  stop("MV7-H landscape repeat admission is stale.", call. = FALSE)
}
safe <- function(value) gsub(":", "_", value, fixed = TRUE)
group_root <- file.path(repeat_root, "landscape")
log_root <- file.path(repeat_root, "logs")
dir.create(group_root, recursive = TRUE, showWarnings = FALSE)
dir.create(log_root, recursive = TRUE, showWarnings = FALSE)
if (length(list.files(group_root, pattern = "__partial__", full.names = TRUE))) {
  stop("MV7-H repeat root contains partial state.", call. = FALSE)
}
tree_rss <- function(pid) {
  root <- tryCatch(ps::ps_handle(pid), error = function(error) NULL)
  if (is.null(root)) return(0)
  handles <- c(list(root), tryCatch(ps::ps_children(root, recursive = TRUE),
                                   error = function(error) list()))
  sum(vapply(handles, function(handle) tryCatch(
    as.numeric(ps::ps_memory_info(handle)[["rss"]]),
    error = function(error) 0
  ), numeric(1L)))
}
write_checkpoint <- function(value) {
  temporary <- paste0(metrics_path, ".partial.", Sys.getpid())
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  dir.create(dirname(metrics_path), recursive = TRUE, showWarnings = FALSE)
  write.csv(value, temporary, row.names = FALSE, na = "")
  if (file.exists(metrics_path) && !file.remove(metrics_path)) {
    stop("MV7-H repeat checkpoint replacement failed.", call. = FALSE)
  }
  if (!file.rename(temporary, metrics_path)) {
    stop("MV7-H repeat checkpoint publication failed.", call. = FALSE)
  }
}
success <- c("completed", "reused_validated")
metrics <- if (file.exists(metrics_path)) read_csv(metrics_path) else NULL
if (!is.null(metrics) && (nrow(metrics) > 3L ||
    any(!metrics$disposition %in% success) ||
    any(metrics$exit_status != 0L) || anyDuplicated(metrics$group_id))) {
  stop("MV7-H repeats preserve prior failure and refuse retry.",
       call. = FALSE)
}
completed <- if (is.null(metrics)) character() else metrics$group_id
for (index in seq_len(nrow(selection))) {
  unit <- selection[index, , drop = FALSE]
  if (unit$group_id %in% completed) next
  stem <- safe(unit$group_id)
  final_dir <- file.path(group_root, stem)
  preexisting <- dir.exists(final_dir)
  started <- Sys.time()
  process <- processx::process$new(
    command = Sys.which("Rscript"), args = c(
      "--vanilla", "scripts/run_mv07h_landscape_group.R", prefreeze,
      ph_root, rust, unit$group_id, group_root,
      contract$implementation_root_sha256
    ), stdout = file.path(log_root, paste0(stem, "__stdout.txt")),
    stderr = file.path(log_root, paste0(stem, "__stderr.txt")),
    cleanup_tree = TRUE
  )
  peak <- 0
  cap_failure <- NA_character_
  while (process$is_alive()) {
    Sys.sleep(0.25)
    peak <- max(peak, tree_rss(process$get_pid()))
    elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    if (elapsed > unit$elapsed_cap_seconds) {
      cap_failure <- "group_elapsed_cap_exceeded"; process$kill_tree()
    } else if (peak > unit$rss_cap_bytes) {
      cap_failure <- "group_rss_cap_exceeded"; process$kill_tree()
    }
  }
  process$wait(timeout = 5000)
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  status_code <- process$get_exit_status()
  status_path <- file.path(final_dir, "status.csv")
  repeat_status <- if (file.exists(status_path)) read_csv(status_path) else NULL
  original <- inventory[inventory$group_id == unit$group_id, , drop = FALSE]
  complete <- is.na(cap_failure) && identical(status_code, 0L) &&
    !is.null(repeat_status) && nrow(repeat_status) == 1L &&
    nrow(original) == 1L &&
    repeat_status$distances_sha256 == original$distances_sha256 &&
    repeat_status$distances_sha256 == .mv07h_sha256(file.path(
      final_dir, "distances.csv"))
  disposition <- if (!is.na(cap_failure)) cap_failure else if (complete) {
    if (preexisting) "reused_validated" else "completed"
  } else "failed"
  metric <- data.frame(
    contract_id = "mv07h_landscape_completion_repeat_resource_v1",
    group_order = unit$group_order, group_id = unit$group_id,
    seed = unit$seed, view_id = unit$view_id,
    homology_dimension = unit$homology_dimension,
    disposition = disposition,
    exit_status = if (is.null(status_code)) NA_integer_ else status_code,
    elapsed_seconds = elapsed, peak_process_tree_rss_bytes = peak,
    elapsed_cap_seconds = unit$elapsed_cap_seconds,
    rss_cap_bytes = unit$rss_cap_bytes,
    distances_sha256 = if (complete)
      repeat_status$distances_sha256 else NA_character_,
    byte_identical_to_production = complete,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    dimension_combination_jobs = 0L, clustering_jobs = 0L,
    label_jobs = 0L, outcome_jobs = 0L, stringsAsFactors = FALSE
  )
  metrics <- if (is.null(metrics)) metric else rbind(metrics, metric)
  write_checkpoint(metrics)
  if (!complete) {
    stop("MV7-H completion repeat failed closed at ", unit$group_id,
         call. = FALSE)
  }
  completed <- c(completed, unit$group_id)
  message("Completed MV7-H completion repeat ", index, "/3")
}
if (nrow(metrics) != 3L || any(!metrics$byte_identical_to_production)) {
  stop("MV7-H completion repeats are incomplete.", call. = FALSE)
}
message("MV7-H landscape completion repeats: 3/3 byte-identical")
