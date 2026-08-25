#!/usr/bin/env Rscript

# Run the committed MV8-ZT landscape queue serially. Scientific outputs and
# logs remain private; public receipts contain only aggregate chunk/resource
# evidence. There are no automatic retries or mixed-engine fallbacks.

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1", RAYON_NUM_THREADS = "1")
options(warn = 2)
for (package in c("digest", "processx", "ps")) {
  if (!requireNamespace(package, quietly = TRUE)) stop(package, " required", call. = FALSE)
}
args <- commandArgs(trailingOnly = TRUE)
if (!(length(args) %in% c(9L, 10L))) stop(paste(
  "usage: run_mv08zt_full_landscape_production.R <mv08zt-prefreeze>",
  "<mv08z-prefreeze> <private-bindings> <mv08s-private> <mv08v-private>",
  "<rust-library> <private-output> <public-output> <expected-head> [--resume]"
), call. = FALSE)

zt_root <- normalizePath(args[[1L]], mustWork = TRUE)
z_root <- normalizePath(args[[2L]], mustWork = TRUE)
bindings_path <- normalizePath(args[[3L]], mustWork = TRUE)
s_root <- normalizePath(args[[4L]], mustWork = TRUE)
v_root <- normalizePath(args[[5L]], mustWork = TRUE)
rust_library <- normalizePath(args[[6L]], mustWork = TRUE)
private_root <- normalizePath(args[[7L]], mustWork = FALSE)
public_root <- normalizePath(args[[8L]], mustWork = FALSE)
expected_head <- tolower(trimws(args[[9L]]))
resume <- length(args) == 10L && identical(args[[10L]], "--resume")
if (length(args) == 10L && !resume) stop("unknown MV8-ZT runner flag", call. = FALSE)
environment_head <- tolower(trimws(Sys.getenv("MV08ZT_GIT_HEAD", unset = "")))
if (!grepl("^[0-9a-f]{40}$", expected_head) || expected_head != environment_head) {
  stop("MV8-ZT exact committed HEAD binding failed", call. = FALSE)
}
if ((dir.exists(private_root) || dir.exists(public_root)) && !resume) {
  stop("MV8-ZT output roots exist; explicit --resume required", call. = FALSE)
}
if (resume && (!dir.exists(private_root) || !dir.exists(public_root))) {
  stop("MV8-ZT resume requires both existing roots", call. = FALSE)
}

source("R/mv08z_landscape_production.R")
read_csv <- .mv08z_read_csv
sha_file <- .mv08z_sha256_file
truth <- .mv08z_truth
atomic_csv <- .mv08z_atomic_csv
.mv08z_verify_manifest(zt_root, "mv08zt-artifact-manifest.csv")
.mv08z_verify_manifest(z_root, "mv08z-artifact-manifest.csv")

contract <- read_csv(file.path(zt_root, "mv08zt-contract.csv"))
queue <- read_csv(file.path(zt_root, "mv08zt-production-queue.csv"))
implementation <- read_csv(file.path(zt_root, "mv08zt-implementation-bindings.csv"))
z_inputs <- read_csv(file.path(z_root, "mv08z-input-manifest.csv"))
bindings <- read_csv(bindings_path)
if (nrow(contract) != 1L || nrow(queue) != 628L ||
    !identical(as.integer(queue$production_order), 1:628) ||
    sum(as.integer(queue$pair_count)) != 152744L ||
    !all(file.exists(implementation$file)) ||
    !all(vapply(implementation$file, sha_file, character(1L)) == implementation$sha256) ||
    sha_file(bindings_path) != z_inputs$sha256[z_inputs$role == "private_unit_bindings"] ||
    sha_file(rust_library) != contract$rust_library_sha256 ||
    as.numeric(file.info(rust_library)$size) != contract$rust_library_bytes ||
    contract$scientific_engine_version != 2L ||
    !truth(contract$fresh_production_authorized_after_commit)) {
  stop("MV8-ZT implementation or input binding drift", call. = FALSE)
}

root_for <- function(role) {
  if (role == "mv08s_private_v3") return(s_root)
  if (role == "mv08v_recovery_private_v2") return(v_root)
  stop("MV8-ZT unknown PH source role", call. = FALSE)
}
binding_key <- paste(bindings$source_role, bindings$output_file, sep = "\r")
source_rows <- bindings[!duplicated(binding_key), , drop = FALSE]
source_paths <- vapply(seq_len(nrow(source_rows)), function(index) {
  file.path(root_for(source_rows$source_role[[index]]), source_rows$output_file[[index]])
}, character(1L))
if (!all(file.exists(source_paths)) ||
    !all(as.numeric(file.info(source_paths)$size) == as.numeric(source_rows$output_bytes)) ||
    !all(vapply(source_paths, sha_file, character(1L)) == source_rows$output_sha256)) {
  stop("MV8-ZT PH input rehash failed", call. = FALSE)
}

free_disk <- as.numeric(Sys.getenv("MV08ZT_FREE_DISK_BYTES", unset = ""))
available_memory <- as.numeric(Sys.getenv("MV08ZT_AVAILABLE_MEMORY_BYTES", unset = ""))
if (!is.finite(free_disk) || !is.finite(available_memory) ||
    free_disk < contract$minimum_free_disk_bytes ||
    available_memory < contract$minimum_available_memory_bytes) {
  stop("MV8-ZT launch headroom gate failed", call. = FALSE)
}

dir.create(private_root, recursive = TRUE, showWarnings = FALSE)
dir.create(public_root, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(private_root, "production"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(private_root, "logs"), recursive = TRUE, showWarnings = FALSE)

tree_rss <- function(pid) {
  root <- tryCatch(ps::ps_handle(pid), error = function(error) NULL)
  if (is.null(root)) return(0)
  handles <- c(list(root), tryCatch(ps::ps_children(root, recursive = TRUE),
                                    error = function(error) list()))
  sum(vapply(handles, function(handle) tryCatch(
    as.numeric(ps::ps_memory_info(handle)[["rss"]]),
    error = function(error) 0), numeric(1L)))
}
private_bytes <- function() {
  files <- list.files(private_root, recursive = TRUE, full.names = TRUE,
                      all.files = TRUE, no.. = TRUE)
  files <- files[!file.info(files)$isdir]
  sum(as.numeric(file.info(files)$size))
}
classify_stderr <- function(path) {
  lines <- if (file.exists(path)) readLines(path, warn = FALSE) else character()
  text <- trimws(paste(lines, collapse = "\n"))
  if (!nzchar(text)) return("empty")
  if (length(lines) == 1L && grepl(
    "^Completed MV8-Z group_[0-9]+/chunk_[0-9]+; pairs=[0-9]+$",
    text, perl = TRUE
  )) "expected_completion" else "unexpected"
}

ledger_path <- file.path(public_root, "mv08zt-resource-ledger.csv")
completion_path <- file.path(public_root, "mv08zt-chunk-completions.csv")
progress_path <- file.path(public_root, "mv08zt-progress.csv")
ledger <- if (file.exists(ledger_path)) read_csv(ledger_path) else data.frame()
completed <- if (file.exists(completion_path)) read_csv(completion_path) else data.frame()

chunk_paths <- function(row) {
  root <- file.path(private_root, "production",
                    .mv08z_safe_group(row$group_order),
                    .mv08z_safe_chunk(row$chunk_order))
  c(distance = file.path(root, "distances.csv"), status = file.path(root, "status.csv"))
}
if (resume) {
  if (nrow(completed) > nrow(queue) ||
      !identical(as.integer(completed$production_order), seq_len(nrow(completed))) ||
      !identical(as.integer(completed$global_chunk_order),
                 as.integer(queue$global_chunk_order[seq_len(nrow(completed))])) ||
      nrow(ledger) != nrow(completed) ||
      !identical(as.integer(ledger$production_order), seq_len(nrow(ledger))) ||
      any(ledger$disposition != "completed")) {
    stop("MV8-ZT resume is not an exact completed strict prefix", call. = FALSE)
  }
  if (nrow(completed)) for (index in seq_len(nrow(completed))) {
    paths <- chunk_paths(queue[index, , drop = FALSE])
    status <- read_csv(paths[["status"]])
    if (!all(file.exists(paths)) || nrow(status) != 1L ||
        sha_file(paths[["distance"]]) != completed$distances_sha256[[index]] ||
        sha_file(paths[["status"]]) != completed$status_sha256[[index]] ||
        status$completion_state != "complete" || status$mode != "production" ||
        status$execution_head != expected_head ||
        status$pair_subset_sha256 != queue$pair_subset_sha256[[index]]) {
      stop("MV8-ZT completed prefix artifact drift at order ", index, call. = FALSE)
    }
  }
  expected_logs <- unlist(lapply(seq_len(nrow(completed)), function(index) {
    file.path(private_root, "logs", sprintf("chunk_%04d.%s", index, c("stdout", "stderr")))
  }))
  actual_logs <- list.files(file.path(private_root, "logs"), full.names = TRUE)
  if (!setequal(normalizePath(actual_logs, mustWork = TRUE),
                normalizePath(expected_logs, mustWork = TRUE))) {
    stop("MV8-ZT resume has ambiguous logs beyond completed prefix", call. = FALSE)
  }
} else if (any(file.exists(c(ledger_path, completion_path, progress_path)))) {
  stop("MV8-ZT public output unexpectedly pre-exists", call. = FALSE)
}

publish_progress <- function(state) {
  value <- data.frame(
    contract_id = "mv08zt_progress_v1", execution_head = expected_head,
    state = state, completed_chunks = nrow(completed), total_chunks = nrow(queue),
    completed_pairs = if (nrow(completed)) sum(completed$pair_count) else 0L,
    total_pairs = sum(queue$pair_count),
    aggregate_child_seconds = if (nrow(ledger)) sum(ledger$elapsed_seconds) else 0,
    private_bytes = private_bytes(), workers = 1L, retries = 0L,
    comparison_jobs = 0L, clustering_jobs = 0L, fusion_jobs = 0L,
    label_jobs = 0L, outcome_jobs = 0L, adoption_jobs = 0L,
    manuscript_claim_jobs = 0L, outcome_label_state = "closed",
    biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
  )
  atomic_csv(value, progress_path)
}

run_child <- function(index, row) {
  paths <- chunk_paths(row)
  if (any(file.exists(paths)) || dir.exists(dirname(paths[[1L]]))) {
    stop("MV8-ZT unowned chunk output exists at order ", index, call. = FALSE)
  }
  stdout <- file.path(private_root, "logs", sprintf("chunk_%04d.stdout", index))
  stderr <- file.path(private_root, "logs", sprintf("chunk_%04d.stderr", index))
  if (any(file.exists(c(stdout, stderr)))) {
    stop("MV8-ZT ambiguous prior logs at order ", index, call. = FALSE)
  }
  dir.create(dirname(dirname(paths[[1L]])), recursive = TRUE, showWarnings = FALSE)
  Sys.setenv(MV08ZT_PREFREEZE = zt_root, MV08ZT_GIT_HEAD = expected_head)
  started <- Sys.time()
  process <- processx::process$new(
    Sys.which("Rscript"), c(
      "--vanilla", "scripts/run_mv08z_landscape_chunk.R", z_root,
      bindings_path, s_root, v_root, rust_library, row$group_order,
      row$chunk_order, file.path(private_root, "production"), expected_head,
      "production", zt_root
    ), stdout = stdout, stderr = stderr, cleanup_tree = TRUE
  )
  peak <- 0
  cap_failure <- ""
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
  output_ok <- all(file.exists(paths)) && !any(file.exists(paste0(paths, ".partial")))
  status <- if (output_ok) read_csv(paths[["status"]]) else data.frame()
  valid <- identical(exit_status, 0L) && output_ok &&
    stderr_class %in% c("empty", "expected_completion") && !nzchar(cap_failure) &&
    nrow(status) == 1L && status$completion_state == "complete" &&
    status$execution_head == expected_head && status$mode == "production" &&
    status$scientific_engine_version == 2L &&
    status$pair_count == row$pair_count &&
    status$pair_subset_sha256 == row$pair_subset_sha256 &&
    status$distances_sha256 == sha_file(paths[["distance"]])
  disposition <- if (nzchar(cap_failure)) cap_failure else if (valid) "completed" else "failed"
  metric <- data.frame(
    contract_id = "mv08zt_resource_ledger_v1", execution_head = expected_head,
    production_order = index, global_chunk_order = row$global_chunk_order,
    group_order = row$group_order, chunk_order = row$chunk_order,
    pair_count = row$pair_count, disposition = disposition,
    exit_status = if (is.null(exit_status)) NA_integer_ else exit_status,
    elapsed_seconds = elapsed, peak_process_tree_rss_bytes = peak,
    elapsed_cap_seconds = contract$child_elapsed_cap_seconds,
    rss_cap_bytes = contract$child_rss_cap_bytes,
    distances_bytes = if (output_ok) as.numeric(file.info(paths[["distance"]])$size) else NA_real_,
    distances_sha256 = if (output_ok) sha_file(paths[["distance"]]) else NA_character_,
    status_bytes = if (output_ok) as.numeric(file.info(paths[["status"]])$size) else NA_real_,
    status_sha256 = if (output_ok) sha_file(paths[["status"]]) else NA_character_,
    stdout_bytes = as.numeric(file.info(stdout)$size),
    stderr_bytes = as.numeric(file.info(stderr)$size),
    stderr_class = stderr_class, scientific_engine_version = 2L,
    workers = 1L, retries = 0L,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
  ledger <<- if (nrow(ledger)) rbind(ledger, metric) else metric
  atomic_csv(ledger, ledger_path)
  if (!valid) stop("MV8-ZT chunk failed closed at production order ", index, call. = FALSE)
  completion <- data.frame(
    contract_id = "mv08zt_chunk_completion_v1", execution_head = expected_head,
    production_order = index, global_chunk_order = row$global_chunk_order,
    group_order = row$group_order, chunk_order = row$chunk_order,
    pair_count = row$pair_count, pair_subset_sha256 = row$pair_subset_sha256,
    distances_bytes = metric$distances_bytes,
    distances_sha256 = metric$distances_sha256,
    status_bytes = metric$status_bytes, status_sha256 = metric$status_sha256,
    exact = TRUE, all_active_levels = TRUE, grid_points = 0L,
    level_cap_applied = FALSE, scientific_engine_version = 2L,
    workers = 1L, retries = 0L,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
  completed <<- if (nrow(completed)) rbind(completed, completion) else completion
  atomic_csv(completed, completion_path)
}

start <- nrow(completed) + 1L
publish_progress("running")
if (start <= nrow(queue)) for (index in seq.int(start, nrow(queue))) {
  if ((nrow(ledger) && sum(ledger$elapsed_seconds) >=
       contract$aggregate_elapsed_cap_seconds) ||
      private_bytes() >= contract$private_storage_cap_bytes) {
    publish_progress("stopped_before_next_chunk_resource_cap")
    stop("MV8-ZT aggregate resource cap reached; evidence preserved", call. = FALSE)
  }
  run_child(index, queue[index, , drop = FALSE])
  publish_progress("running")
}
publish_progress("landscape_production_complete_closure_pending")
cat("MV8-ZT full landscape production completed 628/628 chunks; MV8-ZG closure required\n")
