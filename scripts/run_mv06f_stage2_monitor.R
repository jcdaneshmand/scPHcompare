#!/usr/bin/env Rscript

options(warn = 2)
for (package in c("processx", "ps", "digest")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(package, " is required for monitored MV6-F stage two.",
         call. = FALSE)
  }
}
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1", RAYON_NUM_THREADS = "1")
args <- getOption("mv06f.args", commandArgs(trailingOnly = TRUE))
if (length(args) != 13L) {
  stop("usage: run_mv06f_stage2_monitor.R QUEUE CONTRACT SOURCES ",
       "RESOURCE_PLAN ADMISSION CANDIDATE FOLDS RESOURCES PANEL CACHE_DIR ",
       "RUST_LIBRARY PRIVATE_ROOT PUBLIC_METRICS", call. = FALSE)
}
paths <- vapply(args[1:11], normalizePath, character(1L), winslash = "/",
                mustWork = TRUE)
private_root <- args[[12L]]
metrics_path <- args[[13L]]
source("R/mv06f_production.R")
source("R/mv06f_stage2_execution.R")
read_csv <- function(index) utils::read.csv(
  paths[[index]], stringsAsFactors = FALSE, check.names = FALSE
)
queue <- read_csv(1L); contract <- read_csv(2L); sources <- read_csv(3L)
resource_plan <- read_csv(4L); admission <- read_csv(5L)
mv06f_validate_stage2_rebind_contract_v1(
  queue, contract, sources, resource_plan, paths[[11L]]
)
mv06f_validate_stage2_admission_v1(admission, contract)
stage2 <- queue[queue$stage == "stage_2", , drop = FALSE]
stage2 <- stage2[order(stage2$execution_order), , drop = FALSE]

dir.create(private_root, recursive = TRUE, showWarnings = FALSE)
group_root <- file.path(private_root, "groups")
log_root <- file.path(private_root, "logs")
dir.create(group_root, recursive = TRUE, showWarnings = FALSE)
dir.create(log_root, recursive = TRUE, showWarnings = FALSE)
safe_name <- function(value) gsub("[^A-Za-z0-9_.-]", "_", value)
partials <- list.files(group_root, pattern = "\\.partial\\.", full.names = TRUE)
if (length(partials)) {
  stop("MV6-F stage two found partial state; quarantine is required.",
       call. = FALSE)
}
expected_dirs <- safe_name(queue$group_id)
observed_dirs <- basename(list.dirs(group_root, recursive = FALSE,
                                    full.names = TRUE))
if (length(setdiff(observed_dirs, expected_dirs))) {
  stop("MV6-F stage two found an unexpected group directory.",
       call. = FALSE)
}
stage1 <- queue[queue$stage == "stage_1_maximum", , drop = FALSE]
stage1_dir <- file.path(group_root, safe_name(stage1$group_id))
mv06f_validate_group_directory_v1(
  stage1_dir, stage1, contract$queue_root_sha256,
  contract$implementation_root_sha256, contract$rust_library_sha256
)

private_bytes <- function() {
  files <- list.files(private_root, recursive = TRUE, full.names = TRUE,
                      all.files = TRUE, no.. = TRUE)
  if (!length(files)) return(0)
  info <- file.info(files)
  sum(info$size[!info$isdir], na.rm = TRUE)
}
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
write_checkpoint <- function(value) {
  value <- value[order(value$execution_order), , drop = FALSE]
  dir.create(dirname(metrics_path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(metrics_path, ".partial.", Sys.getpid())
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  utils::write.csv(value, temporary, row.names = FALSE, na = "")
  if (file.exists(metrics_path) && !file.remove(metrics_path)) {
    stop("MV6-F stage-two checkpoint replacement failed.", call. = FALSE)
  }
  if (!file.rename(temporary, metrics_path)) {
    stop("MV6-F stage-two checkpoint publication failed.", call. = FALSE)
  }
}

metrics <- NULL
if (file.exists(metrics_path)) {
  metrics <- utils::read.csv(metrics_path, stringsAsFactors = FALSE,
                             check.names = FALSE)
  if (any(!metrics$disposition %in% .mv06f_stage2_success)) {
    stop("MV6-F stage two preserves a prior failure and refuses retry.",
         call. = FALSE)
  }
  mv06f_validate_stage2_metrics_v1(metrics, queue, contract)
  for (index in seq_len(nrow(metrics))) {
    row <- stage2[stage2$group_id == metrics$group_id[[index]], , drop = FALSE]
    mv06f_validate_group_directory_v1(
      file.path(group_root, safe_name(row$group_id)), row,
      contract$queue_root_sha256, contract$implementation_root_sha256,
      contract$rust_library_sha256
    )
  }
}
completed_ids <- if (is.null(metrics)) character() else metrics$group_id
prior_worker_seconds <- if (is.null(metrics)) 0 else
  sum(metrics$charged_worker_seconds)

for (index in seq_len(nrow(stage2))) {
  unit <- stage2[index, , drop = FALSE]
  if (unit$group_id %in% completed_ids) next
  preflight_guard <- mv06f_stage2_guard_v1(
    prior_worker_seconds, 0, 0, 0, 0, private_bytes(), resource_plan
  )
  if (!preflight_guard$launch_authorized) {
    stop("MV6-F stage-two admission closed: ",
         preflight_guard$disposition, call. = FALSE)
  }
  stem <- safe_name(unit$group_id)
  final_dir <- file.path(group_root, stem)
  preexisting <- dir.exists(final_dir)
  stdout_path <- file.path(log_root, paste0(stem, "__stdout.txt"))
  stderr_path <- file.path(log_root, paste0(stem, "__stderr.txt"))
  runner_args <- c(
    "--vanilla", "scripts/run_mv06f_group.R", paths[[1L]], paths[[2L]],
    paths[[3L]], paths[[6L]], paths[[7L]], paths[[8L]], paths[[9L]],
    paths[[10L]], paths[[11L]], unit$group_id, group_root
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
    rss <- process_tree_rss(process$get_pid())
    peak <- max(peak, rss)
    elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    guard <- mv06f_stage2_guard_v1(
      prior_worker_seconds, elapsed, elapsed, peak, rss, private_bytes(),
      resource_plan
    )
    if (!guard$launch_authorized) {
      cap_failure <- guard$disposition
      process$kill_tree()
    }
  }
  process$wait(timeout = 5000)
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  exit_status <- process$get_exit_status()
  complete <- is.na(cap_failure) && identical(exit_status, 0L) &&
    dir.exists(final_dir)
  group_status <- NULL
  group_metrics <- NULL
  if (complete) {
    group_status <- mv06f_validate_group_directory_v1(
      final_dir, unit, contract$queue_root_sha256,
      contract$implementation_root_sha256, contract$rust_library_sha256
    )
    group_metrics <- utils::read.csv(file.path(final_dir, "metrics.csv"),
                                     stringsAsFactors = FALSE)
  }
  disposition <- if (!is.na(cap_failure)) cap_failure else if (complete) {
    if (preexisting) "reused_validated" else "completed"
  } else "failed"
  charged <- if (complete && preexisting) {
    max(elapsed, as.numeric(group_metrics$elapsed_seconds))
  } else elapsed
  metric <- data.frame(
    contract_id = "mv06f_stage2_resource_metric_v1",
    group_id = unit$group_id, fold_id = unit$fold_id, seed = unit$seed,
    execution_order = unit$execution_order,
    biological_pairs = unit$biological_pairs,
    landscape_component_rows = unit$landscape_component_rows,
    disposition = disposition,
    exit_status = if (is.null(exit_status)) NA_integer_ else exit_status,
    elapsed_seconds = elapsed, charged_worker_seconds = charged,
    peak_process_tree_rss_bytes = peak,
    cumulative_private_bytes = private_bytes(),
    queue_root_sha256 = contract$queue_root_sha256,
    implementation_root_sha256 = contract$implementation_root_sha256,
    rust_library_sha256 = contract$rust_library_sha256,
    group_directory_complete = complete,
    diagrams_sha256 = if (complete) group_status$diagrams_sha256 else NA,
    diagram_manifest_sha256 = if (complete)
      group_status$diagram_manifest_sha256 else NA,
    distances_sha256 = if (complete) group_status$distances_sha256 else NA,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    fusion_jobs = 0L, clustering_jobs = 0L, outcome_jobs = 0L,
    stringsAsFactors = FALSE
  )
  metrics <- if (is.null(metrics)) metric else rbind(metrics, metric)
  write_checkpoint(metrics)
  if (!complete) {
    stop("MV6-F stage two stopped at ", unit$group_id, ": ", disposition,
         call. = FALSE)
  }
  prior_worker_seconds <- prior_worker_seconds + charged
  completed_ids <- c(completed_ids, unit$group_id)
  message("Completed MV6-F stage-two group ", index, "/74: ",
          unit$group_id, ".")
}
mv06f_validate_stage2_metrics_v1(metrics, queue, contract, complete = TRUE)
message("Completed all 74 monitored MV6-F stage-two groups.")
