#!/usr/bin/env Rscript

options(warn = 2)
for (package in c("processx", "ps", "digest")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(package, " is required for monitored MV7-H landscapes.",
         call. = FALSE)
  }
}
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1", RAYON_NUM_THREADS = "1")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) {
  stop("usage: run_mv07h_remaining_landscapes.R PREFREEZE PH_ROOT ",
       "RUST_LIBRARY STRESS_VALIDATION STRESS_RESOURCE PRIVATE_ROOT ",
       "METRICS_PATH", call. = FALSE)
}
prefreeze <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
ph_root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
rust <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
stress_validation <- normalizePath(args[[4L]], winslash = "/",
                                   mustWork = TRUE)
stress_resource_path <- normalizePath(args[[5L]], winslash = "/",
                                      mustWork = TRUE)
private_root <- args[[6L]]
metrics_path <- args[[7L]]
source("R/mv07h_full_topology.R")

read_csv <- function(path) utils::read.csv(
  path, stringsAsFactors = FALSE, check.names = FALSE)
queue <- read_csv(file.path(prefreeze, "mv07h-landscape-queue.csv"))
contract <- read_csv(file.path(prefreeze, "mv07h-contract.csv"))
projection <- read_csv(file.path(prefreeze, "mv07h-resource-projection.csv"))
decision <- read_csv(file.path(
  stress_validation, "mv07h-stress-validation-decision.csv"))
stress_resource <- read_csv(stress_resource_path)
stage2 <- queue[queue$stage == "stage_2", , drop = FALSE]
stage2 <- stage2[order(stage2$group_order, method = "radix"), , drop = FALSE]
landscape_cap <- projection$projected_worker_seconds[
  projection$component == "all_landscape_groups"]
landscape_definition <- paste(
  "finite_positive;essential_H0_excluded;all_active_levels;",
  "exact_squared_L2;no_grid;no_level_cap;streamed_groups", sep = "")
if (nrow(queue) != 20L || nrow(stage2) != 19L ||
    !identical(as.integer(stage2$group_order), 2:20) ||
    anyDuplicated(queue$group_id) || any(queue$samples != 124L) ||
    any(queue$component_rows != 7626L) || any(queue$workers != 1L) ||
    any(queue$retries != 0L) ||
    any(queue$outcome_label_state != "closed") ||
    any(as.logical(queue$biological_outcomes_computed)) ||
    any(unlist(queue[c("clustering_jobs", "label_jobs", "outcome_jobs")],
               use.names = FALSE) != 0L) ||
    nrow(contract) != 1L || contract$landscape_groups != 20L ||
    contract$landscape_component_rows != 152520L ||
    contract$landscape_definition != landscape_definition ||
    contract$rust_library_sha256 != .mv07h_sha256(rust) ||
    length(landscape_cap) != 1L || !is.finite(landscape_cap) ||
    nrow(decision) != 1L ||
    decision$decision !=
      "authorize_remaining_19_MV7H_landscape_groups_serially" ||
    decision$remaining_groups_authorized != 19L ||
    decision$outcome_label_state != "closed" ||
    as.logical(decision$biological_outcomes_computed) ||
    any(unlist(decision[c("clustering_jobs", "label_jobs", "outcome_jobs")],
               use.names = FALSE) != 0L) ||
    nrow(stress_resource) != 1L ||
    stress_resource$group_id != queue$group_id[[1L]] ||
    stress_resource$disposition != "completed" ||
    as.logical(stress_resource$repeat_mode) ||
    stress_resource$outcome_label_state != "closed" ||
    as.logical(stress_resource$biological_outcomes_computed)) {
  stop("MV7-H remaining-landscape admission evidence is stale.",
       call. = FALSE)
}

safe_name <- function(value) gsub(":", "_", value, fixed = TRUE)
group_root <- file.path(private_root, "landscape")
log_root <- file.path(private_root, "stage2-logs")
dir.create(group_root, recursive = TRUE, showWarnings = FALSE)
dir.create(log_root, recursive = TRUE, showWarnings = FALSE)
partials <- list.files(group_root, pattern = "__partial__", full.names = TRUE)
if (length(partials)) {
  stop("MV7-H stage two found partial state; quarantine is required.",
       call. = FALSE)
}
observed <- basename(list.dirs(group_root, recursive = FALSE,
                               full.names = TRUE))
if (length(setdiff(observed, safe_name(queue$group_id)))) {
  stop("MV7-H stage two found an unexpected group directory.",
       call. = FALSE)
}

validate_group <- function(unit) {
  directory <- file.path(group_root, safe_name(unit$group_id))
  status_path <- file.path(directory, "status.csv")
  distances_path <- file.path(directory, "distances.csv")
  metrics_group_path <- file.path(directory, "metrics.csv")
  if (!all(file.exists(c(status_path, distances_path, metrics_group_path)))) {
    stop("MV7-H landscape group artifacts are incomplete.", call. = FALSE)
  }
  status <- read_csv(status_path)
  distances <- read_csv(distances_path)
  metrics <- read_csv(metrics_group_path)
  if (nrow(status) != 1L || nrow(metrics) != 1L ||
      nrow(distances) != unit$component_rows ||
      status$group_id != unit$group_id ||
      status$completion_state != "complete" ||
      status$implementation_root_sha256 !=
        contract$implementation_root_sha256 ||
      status$rust_library_sha256 != contract$rust_library_sha256 ||
      status$distances_sha256 != .mv07h_sha256(distances_path) ||
      status$metrics_sha256 != .mv07h_sha256(metrics_group_path) ||
      status$component_rows != unit$component_rows ||
      anyDuplicated(distances$pair_id) ||
      any(distances$group_id != unit$group_id) ||
      any(distances$seed != unit$seed) ||
      any(distances$view_id != unit$view_id) ||
      any(distances$homology_dimension != unit$homology_dimension) ||
      any(!is.finite(distances$squared_distance)) ||
      any(distances$squared_distance < 0) ||
      any(!as.logical(distances$exact)) ||
      any(!as.logical(distances$all_active_levels)) ||
      any(as.logical(distances$level_cap_applied)) ||
      any(distances$outcome_label_state != "closed") ||
      any(as.logical(distances$biological_outcomes_computed)) ||
      any(unlist(distances[c("clustering_jobs", "label_jobs", "outcome_jobs")],
                 use.names = FALSE) != 0L)) {
    stop("MV7-H landscape group is stale or scientifically invalid: ",
         unit$group_id, call. = FALSE)
  }
  list(status = status, metrics = metrics, directory = directory)
}

invisible(validate_group(queue[queue$stage == "stage_1_stress", ,
                               drop = FALSE]))
process_tree_rss <- function(pid) {
  root <- tryCatch(ps::ps_handle(pid), error = function(error) NULL)
  if (is.null(root)) return(0)
  handles <- c(list(root), tryCatch(ps::ps_children(root, recursive = TRUE),
                                   error = function(error) list()))
  sum(vapply(handles, function(handle) tryCatch(
    as.numeric(ps::ps_memory_info(handle)[["rss"]]),
    error = function(error) 0
  ), numeric(1L)))
}
private_bytes <- function() {
  paths <- list.files(private_root, recursive = TRUE, full.names = TRUE,
                      all.files = TRUE, no.. = TRUE)
  if (!length(paths)) return(0)
  info <- file.info(paths)
  sum(info$size[!info$isdir], na.rm = TRUE)
}
write_checkpoint <- function(value) {
  value <- value[order(value$group_order), , drop = FALSE]
  dir.create(dirname(metrics_path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(metrics_path, ".partial.", Sys.getpid())
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  write.csv(value, temporary, row.names = FALSE, na = "")
  if (file.exists(metrics_path) && !file.remove(metrics_path)) {
    stop("MV7-H stage-two checkpoint replacement failed.", call. = FALSE)
  }
  if (!file.rename(temporary, metrics_path)) {
    stop("MV7-H stage-two checkpoint publication failed.", call. = FALSE)
  }
}

success <- c("completed", "reused_validated")
metrics <- if (file.exists(metrics_path)) read_csv(metrics_path) else NULL
if (!is.null(metrics)) {
  required <- c("group_id", "group_order", "disposition", "exit_status",
                "charged_worker_seconds", "distances_sha256")
  if (!all(required %in% names(metrics)) || nrow(metrics) > 19L ||
      anyDuplicated(metrics$group_id) ||
      any(!metrics$group_id %in% stage2$group_id) ||
      any(!metrics$disposition %in% success) ||
      any(metrics$exit_status != 0L)) {
    stop("MV7-H stage two preserves a prior failure and refuses retry.",
         call. = FALSE)
  }
  for (index in seq_len(nrow(metrics))) {
    unit <- stage2[stage2$group_id == metrics$group_id[[index]], , drop = FALSE]
    verified <- validate_group(unit)
    if (metrics$distances_sha256[[index]] !=
        verified$status$distances_sha256) {
      stop("MV7-H stage-two checkpoint is stale.", call. = FALSE)
    }
  }
}
completed_ids <- if (is.null(metrics)) character() else metrics$group_id
prior_worker_seconds <- stress_resource$elapsed_seconds +
  if (is.null(metrics)) 0 else sum(metrics$charged_worker_seconds)

for (index in seq_len(nrow(stage2))) {
  unit <- stage2[index, , drop = FALSE]
  if (unit$group_id %in% completed_ids) next
  if (prior_worker_seconds >= landscape_cap) {
    stop("MV7-H total landscape worker-time cap is exhausted.",
         call. = FALSE)
  }
  stem <- safe_name(unit$group_id)
  final_dir <- file.path(group_root, stem)
  preexisting <- dir.exists(final_dir)
  stdout_path <- file.path(log_root, paste0(stem, "__stdout.txt"))
  stderr_path <- file.path(log_root, paste0(stem, "__stderr.txt"))
  runner_args <- c(
    "--vanilla", "scripts/run_mv07h_landscape_group.R", prefreeze, ph_root,
    rust, unit$group_id, group_root, contract$implementation_root_sha256
  )
  started <- Sys.time()
  process <- processx::process$new(
    command = Sys.which("Rscript"), args = runner_args,
    stdout = stdout_path, stderr = stderr_path, cleanup_tree = TRUE
  )
  peak <- 0
  cap_failure <- NA_character_
  while (process$is_alive()) {
    Sys.sleep(0.25)
    peak <- max(peak, process_tree_rss(process$get_pid()))
    elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    if (elapsed > unit$elapsed_cap_seconds) {
      cap_failure <- "group_elapsed_cap_exceeded"
      process$kill_tree()
    } else if (peak > unit$rss_cap_bytes) {
      cap_failure <- "group_rss_cap_exceeded"
      process$kill_tree()
    } else if (prior_worker_seconds + elapsed > landscape_cap) {
      cap_failure <- "total_landscape_worker_cap_exceeded"
      process$kill_tree()
    }
  }
  process$wait(timeout = 5000)
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  exit_status <- process$get_exit_status()
  verified <- if (is.na(cap_failure) && identical(exit_status, 0L) &&
                  dir.exists(final_dir)) validate_group(unit) else NULL
  complete <- !is.null(verified)
  disposition <- if (!is.na(cap_failure)) cap_failure else if (complete) {
    if (preexisting) "reused_validated" else "completed"
  } else "failed"
  charged <- if (complete) max(elapsed,
    as.numeric(verified$metrics$elapsed_seconds)) else elapsed
  metric <- data.frame(
    contract_id = "mv07h_stage2_landscape_resource_v1",
    group_id = unit$group_id, group_order = unit$group_order,
    seed = unit$seed, view_id = unit$view_id,
    homology_dimension = unit$homology_dimension,
    component_rows = unit$component_rows, disposition = disposition,
    exit_status = if (is.null(exit_status)) NA_integer_ else exit_status,
    elapsed_seconds = elapsed, charged_worker_seconds = charged,
    peak_process_tree_rss_bytes = peak,
    cumulative_landscape_worker_seconds = prior_worker_seconds + charged,
    landscape_worker_cap_seconds = landscape_cap,
    cumulative_private_bytes = private_bytes(),
    distances_sha256 = if (complete)
      verified$status$distances_sha256 else NA_character_,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    clustering_jobs = 0L, label_jobs = 0L, outcome_jobs = 0L,
    stringsAsFactors = FALSE
  )
  metrics <- if (is.null(metrics)) metric else rbind(metrics, metric)
  write_checkpoint(metrics)
  if (!complete) {
    stop("MV7-H stage two stopped at ", unit$group_id, ": ", disposition,
         call. = FALSE)
  }
  prior_worker_seconds <- prior_worker_seconds + charged
  completed_ids <- c(completed_ids, unit$group_id)
  message("Completed MV7-H remaining landscape group ", index, "/19: ",
          unit$group_id)
}
if (nrow(metrics) != 19L || !identical(metrics$group_id, stage2$group_id) ||
    any(!metrics$disposition %in% success) ||
    any(metrics$outcome_label_state != "closed") ||
    any(as.logical(metrics$biological_outcomes_computed))) {
  stop("MV7-H stage-two completion metrics are invalid.", call. = FALSE)
}
message("Completed all 19 authorized MV7-H landscape groups serially.")
