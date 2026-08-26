#!/usr/bin/env Rscript

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1", RAYON_NUM_THREADS = "1")
options(warn = 2)
for (package in c("digest", "processx", "ps")) {
  if (!requireNamespace(package, quietly = TRUE)) stop(package, " required")
}
args <- commandArgs(trailingOnly = TRUE)
if (!(length(args) %in% c(7L, 8L))) stop(paste(
  "usage: run_mv14_cell_landscape_production.R <prefreeze>",
  "<private-bindings> <mv13-private-groups> <rust-library>",
  "<private-output> <public-output> <execution-head> [--resume]"
), call. = FALSE)
prefreeze <- normalizePath(args[[1L]], mustWork = TRUE)
binding_path <- normalizePath(args[[2L]], mustWork = TRUE)
group_root <- normalizePath(args[[3L]], mustWork = TRUE)
rust_library <- normalizePath(args[[4L]], mustWork = TRUE)
private_root <- normalizePath(args[[5L]], mustWork = FALSE)
public_root <- normalizePath(args[[6L]], mustWork = FALSE)
execution_head <- tolower(trimws(args[[7L]]))
resume <- length(args) == 8L && identical(args[[8L]], "--resume")
if (length(args) == 8L && !resume) stop("unknown MV14 runner flag")
environment_head <- tolower(trimws(Sys.getenv("MV14_GIT_HEAD", unset = "")))
if (!grepl("^[0-9a-f]{40}$", execution_head) ||
    execution_head != environment_head) {
  stop("MV14 exact committed HEAD binding failed.", call. = FALSE)
}
if ((dir.exists(private_root) || dir.exists(public_root)) && !resume) {
  stop("MV14 output roots exist; explicit --resume required.", call. = FALSE)
}
if (resume && (!dir.exists(private_root) || !dir.exists(public_root))) {
  stop("MV14 resume requires both roots.", call. = FALSE)
}

source("R/mv14_cell_landscape.R")
.mv14_verify_manifest(prefreeze, "mv14-artifact-manifest.csv")
contract <- .mv14_read_csv(file.path(prefreeze, "mv14-contract.csv"))
groups <- .mv14_read_csv(file.path(prefreeze, "mv14-group-queue.csv"))
queue <- .mv14_read_csv(file.path(prefreeze, "mv14-production-queue.csv"))
inputs <- .mv14_read_csv(file.path(prefreeze, "mv14-input-bindings.csv"))
implementation <- .mv14_read_csv(file.path(
  prefreeze, "mv14-implementation-bindings.csv"
))
decision <- .mv14_read_csv(file.path(prefreeze, "mv14-decision.csv"))
if (nrow(contract) != 1L || contract$execution_head != execution_head ||
    nrow(groups) != 14L || nrow(queue) != 314L ||
    sum(queue$pair_count) != 76372L ||
    !identical(as.integer(queue$production_order), 1:314) ||
    nrow(decision) != 1L ||
    !.mv14_truth(decision$production_authorized_after_commit) ||
    .mv14_sha256_file(binding_path) !=
      inputs$sha256[inputs$role == "private_axis_bindings"] ||
    .mv14_sha256_file(rust_library) != contract$rust_library_sha256 ||
    !all(file.exists(implementation$file)) ||
    !all(vapply(implementation$file, .mv14_sha256_file, character(1L)) ==
           implementation$sha256)) {
  stop("MV14 runner implementation/input drift.", call. = FALSE)
}
for (index in seq_len(nrow(groups))) {
  path <- file.path(group_root, groups$artifact_file[[index]])
  if (!file.exists(path) ||
      as.numeric(file.info(path)$size) != groups$artifact_bytes[[index]] ||
      .mv14_sha256_file(path) != groups$artifact_sha256[[index]]) {
    stop("MV14 source-group rehash failed at group ", index, call. = FALSE)
  }
}
free_disk <- as.numeric(Sys.getenv("MV14_FREE_DISK_BYTES", unset = ""))
available_memory <- as.numeric(Sys.getenv("MV14_AVAILABLE_MEMORY_BYTES", unset = ""))
if (!is.finite(free_disk) || !is.finite(available_memory) ||
    free_disk < contract$minimum_free_disk_bytes ||
    available_memory < contract$minimum_available_memory_bytes) {
  stop("MV14 launch headroom gate failed.", call. = FALSE)
}

dir.create(file.path(private_root, "production"), recursive = TRUE,
           showWarnings = FALSE)
dir.create(file.path(private_root, "logs"), recursive = TRUE,
           showWarnings = FALSE)
dir.create(public_root, recursive = TRUE, showWarnings = FALSE)

tree_rss <- function(pid) {
  root <- tryCatch(ps::ps_handle(pid), error = function(error) NULL)
  if (is.null(root)) return(0)
  handles <- c(list(root), tryCatch(ps::ps_children(root, recursive = TRUE),
                                    error = function(error) list()))
  sum(vapply(handles, function(handle) tryCatch(
    as.numeric(ps::ps_memory_info(handle)[["rss"]]),
    error = function(error) 0), numeric(1L)))
}
classify_stderr <- function(path) {
  lines <- if (file.exists(path)) readLines(path, warn = FALSE) else character()
  text <- trimws(paste(lines, collapse = "\n"))
  if (!nzchar(text)) return("empty")
  if (length(lines) == 1L && grepl(
    "^Completed MV14 group_[0-9]+/chunk_[0-9]+; pairs=[0-9]+$",
    text, perl = TRUE
  )) "expected_completion" else "unexpected"
}
chunk_paths <- function(row) {
  root <- file.path(private_root, "production",
                    .mv14_safe_group(row$group_order),
                    .mv14_safe_chunk(row$chunk_order))
  c(distance = file.path(root, "distances.csv"),
    status = file.path(root, "status.csv"))
}
ledger_path <- file.path(public_root, "mv14-resource-ledger.csv")
completion_path <- file.path(public_root, "mv14-chunk-completions.csv")
progress_path <- file.path(public_root, "mv14-progress.csv")
ledger <- if (file.exists(ledger_path)) .mv14_read_csv(ledger_path) else data.frame()
completed <- if (file.exists(completion_path)) {
  .mv14_read_csv(completion_path)
} else data.frame()
if (resume) {
  if (nrow(completed) > nrow(queue) || nrow(ledger) != nrow(completed) ||
      !identical(as.integer(completed$production_order), seq_len(nrow(completed))) ||
      !identical(as.integer(ledger$production_order), seq_len(nrow(ledger))) ||
      any(ledger$disposition != "completed")) {
    stop("MV14 resume is not an exact completed prefix.", call. = FALSE)
  }
  if (nrow(completed)) for (index in seq_len(nrow(completed))) {
    paths <- chunk_paths(queue[index, , drop = FALSE])
    status <- if (all(file.exists(paths))) {
      .mv14_read_csv(paths[["status"]])
    } else data.frame()
    if (nrow(status) != 1L || status$completion_state != "complete" ||
        status$execution_head != execution_head ||
        status$pair_subset_sha256 != queue$pair_subset_sha256[[index]] ||
        .mv14_sha256_file(paths[["distance"]]) !=
          completed$distances_sha256[[index]] ||
        .mv14_sha256_file(paths[["status"]]) !=
          completed$status_sha256[[index]]) {
      stop("MV14 completed-prefix drift at order ", index, call. = FALSE)
    }
  }
} else if (any(file.exists(c(ledger_path, completion_path, progress_path)))) {
  stop("MV14 public output unexpectedly pre-exists.", call. = FALSE)
}

publish_progress <- function(state) {
  value <- data.frame(
    contract_id = "mv14_progress_v1", execution_head = execution_head,
    state = state, completed_chunks = nrow(completed), total_chunks = nrow(queue),
    completed_pairs = if (nrow(completed)) sum(completed$pair_count) else 0L,
    total_pairs = sum(queue$pair_count),
    aggregate_child_seconds = if (nrow(ledger)) sum(ledger$elapsed_seconds) else 0,
    private_bytes = .mv14_private_bytes(private_root), workers = 1L,
    retries = 0L, comparison_jobs = 0L, clustering_jobs = 0L,
    fusion_jobs = 0L, label_jobs = 0L, outcome_jobs = 0L,
    manuscript_claim_jobs = 0L, outcome_label_state = "closed",
    biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
  )
  .mv14_atomic_csv(value, progress_path)
}

run_child <- function(index, row) {
  paths <- chunk_paths(row)
  if (dir.exists(dirname(paths[[1L]])) || any(file.exists(paths))) {
    stop("MV14 unowned chunk output at order ", index, call. = FALSE)
  }
  stdout <- file.path(private_root, "logs", sprintf("chunk_%04d.stdout", index))
  stderr <- file.path(private_root, "logs", sprintf("chunk_%04d.stderr", index))
  if (any(file.exists(c(stdout, stderr)))) {
    stop("MV14 ambiguous child logs at order ", index, call. = FALSE)
  }
  dir.create(dirname(dirname(paths[[1L]])), recursive = TRUE, showWarnings = FALSE)
  started <- Sys.time()
  process <- processx::process$new(
    Sys.which("Rscript"), c(
      "--vanilla", "scripts/run_mv14_cell_landscape_chunk.R", prefreeze,
      binding_path, group_root, rust_library, row$group_order, row$chunk_order,
      file.path(private_root, "production"), execution_head
    ), stdout = stdout, stderr = stderr, cleanup_tree = TRUE
  )
  peak <- 0; cap_failure <- ""
  while (process$is_alive()) {
    Sys.sleep(0.25)
    peak <- max(peak, tree_rss(process$get_pid()))
    elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    if (elapsed > contract$child_elapsed_cap_seconds) {
      cap_failure <- "elapsed_cap_exceeded"; process$kill_tree()
    } else if (peak > contract$child_rss_cap_bytes) {
      cap_failure <- "rss_cap_exceeded"; process$kill_tree()
    }
  }
  process$wait(timeout = 5000)
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  exit_status <- process$get_exit_status()
  stderr_class <- classify_stderr(stderr)
  output_ok <- all(file.exists(paths)) &&
    !any(file.exists(paste0(paths, ".partial")))
  status <- if (output_ok) .mv14_read_csv(paths[["status"]]) else data.frame()
  valid <- identical(exit_status, 0L) && output_ok &&
    stderr_class %in% c("empty", "expected_completion") && !nzchar(cap_failure) &&
    nrow(status) == 1L && status$completion_state == "complete" &&
    status$execution_head == execution_head && status$pair_count == row$pair_count &&
    status$pair_subset_sha256 == row$pair_subset_sha256 &&
    status$distances_sha256 == .mv14_sha256_file(paths[["distance"]])
  disposition <- if (nzchar(cap_failure)) cap_failure else if (valid) "completed" else "failed"
  metric <- data.frame(
    contract_id = "mv14_resource_ledger_v1", execution_head = execution_head,
    production_order = index, group_order = row$group_order,
    chunk_order = row$chunk_order, pair_count = row$pair_count,
    disposition = disposition,
    exit_status = if (is.null(exit_status)) NA_integer_ else exit_status,
    elapsed_seconds = elapsed, peak_process_tree_rss_bytes = peak,
    elapsed_cap_seconds = contract$child_elapsed_cap_seconds,
    rss_cap_bytes = contract$child_rss_cap_bytes,
    distances_bytes = if (output_ok) file.info(paths[["distance"]])$size else NA,
    distances_sha256 = if (output_ok) .mv14_sha256_file(paths[["distance"]]) else NA,
    status_bytes = if (output_ok) file.info(paths[["status"]])$size else NA,
    status_sha256 = if (output_ok) .mv14_sha256_file(paths[["status"]]) else NA,
    stdout_bytes = file.info(stdout)$size, stderr_bytes = file.info(stderr)$size,
    stderr_class = stderr_class, scientific_engine_version = 2L,
    workers = 1L, retries = 0L, outcome_label_state = "closed",
    biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
  )
  ledger <<- if (nrow(ledger)) rbind(ledger, metric) else metric
  .mv14_atomic_csv(ledger, ledger_path)
  if (!valid) stop("MV14 child failed closed at order ", index, call. = FALSE)
  completion <- data.frame(
    contract_id = "mv14_chunk_completion_v1", execution_head = execution_head,
    production_order = index, group_order = row$group_order,
    chunk_order = row$chunk_order, pair_count = row$pair_count,
    pair_subset_sha256 = row$pair_subset_sha256,
    distances_bytes = metric$distances_bytes,
    distances_sha256 = metric$distances_sha256,
    status_bytes = metric$status_bytes, status_sha256 = metric$status_sha256,
    exact = TRUE, all_active_levels = TRUE, grid_points = 0L,
    level_cap_applied = FALSE, scientific_engine_version = 2L,
    workers = 1L, retries = 0L, outcome_label_state = "closed",
    biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
  )
  completed <<- if (nrow(completed)) rbind(completed, completion) else completion
  .mv14_atomic_csv(completed, completion_path)
}

start <- nrow(completed) + 1L
publish_progress("running")
if (start <= nrow(queue)) for (index in seq.int(start, nrow(queue))) {
  if ((nrow(ledger) && sum(ledger$elapsed_seconds) >=
       contract$aggregate_elapsed_cap_seconds) ||
      .mv14_private_bytes(private_root) >= contract$private_storage_cap_bytes) {
    publish_progress("stopped_before_next_chunk_resource_cap")
    stop("MV14 aggregate resource cap reached; evidence preserved.")
  }
  run_child(index, queue[index, , drop = FALSE])
  publish_progress("running")
}
publish_progress("landscape_production_complete_closure_pending")
cat("MV14 full cell-landscape production completed 314/314 chunks; independent closure required\n")
