#!/usr/bin/env Rscript

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1", RAYON_NUM_THREADS = "1")
options(warn = 2)
for (package in c("digest", "processx", "ps")) {
  if (!requireNamespace(package, quietly = TRUE)) stop(package, " required", call. = FALSE)
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 12L) stop(paste(
  "usage: run_mv08za_landscape_sentinel.R <mv08z-prefreeze> <mv08za-prefreeze>",
  "<private-bindings> <private-sentinel> <mv08s-private> <mv08v-private>",
  "<rust-library> <private-output> <public-output> <execution-head>",
  "<helper-recovery-prefreeze> <traversal-recovery-prefreeze>"
), call. = FALSE)
z_root <- normalizePath(args[[1L]], mustWork = TRUE)
za_root <- normalizePath(args[[2L]], mustWork = TRUE)
bindings <- normalizePath(args[[3L]], mustWork = TRUE)
sentinel_private <- normalizePath(args[[4L]], mustWork = TRUE)
s_root <- normalizePath(args[[5L]], mustWork = TRUE)
v_root <- normalizePath(args[[6L]], mustWork = TRUE)
rust_library <- normalizePath(args[[7L]], mustWork = TRUE)
private_root <- normalizePath(args[[8L]], mustWork = FALSE)
public_root <- normalizePath(args[[9L]], mustWork = FALSE)
execution_head <- tolower(args[[10L]])
recovery_root <- normalizePath(args[[11L]], mustWork = TRUE)
traversal_root <- normalizePath(args[[12L]], mustWork = TRUE)
Sys.setenv(MV08ZB_RECOVERY_PREFREEZE = recovery_root)
if (!grepl("^[0-9a-f]{40}$", execution_head) ||
    execution_head != tolower(Sys.getenv("MV08ZA_GIT_HEAD", unset = "")) ||
    dir.exists(private_root) || dir.exists(public_root)) {
  stop("MV8-ZA requires exact mechanical head and fresh roots", call. = FALSE)
}
source("R/mv08z_landscape_production.R")
.mv08z_verify_manifest(z_root, "mv08z-artifact-manifest.csv")
.mv08z_verify_manifest(za_root, "mv08za-artifact-manifest.csv")
.mv08z_verify_manifest(recovery_root, "mv08zb-artifact-manifest.csv")
.mv08z_verify_manifest(traversal_root, "mv08zc-artifact-manifest.csv")
z_inputs <- .mv08z_read_csv(file.path(z_root, "mv08z-input-manifest.csv"))
z_resource <- .mv08z_read_csv(file.path(z_root, "mv08z-resource-policy.csv"))
sentinel <- .mv08z_read_csv(file.path(z_root, "mv08z-sentinel-selection.csv"))
implementation <- .mv08z_read_csv(file.path(traversal_root, "mv08zc-implementation-bindings.csv"))
prior_implementation <- .mv08z_read_csv(file.path(recovery_root, "mv08zb-implementation-bindings.csv"))
decision <- .mv08z_read_csv(file.path(za_root, "mv08za-decision.csv"))
if (nrow(sentinel) != 1L || nrow(decision) != 1L ||
    decision$authorized_child_processes != 3L ||
    .mv08z_sha256_file(bindings) != z_inputs$sha256[z_inputs$role == "private_unit_bindings"] ||
    .mv08z_sha256_file(sentinel_private) !=
      z_inputs$sha256[z_inputs$role == "private_sentinel_selection"] ||
    .mv08z_sha256_file(rust_library) !=
      z_inputs$sha256[z_inputs$role == "admitted_private_rust_library"] ||
    !all(file.exists(implementation$file)) ||
    !all(vapply(implementation$file, .mv08z_sha256_file, character(1L)) ==
           implementation$sha256) ||
    .mv08z_sha256_file("scripts/run_mv08z_landscape_oracle.R") !=
      prior_implementation$sha256[prior_implementation$role == "oracle_worker"]) {
  stop("MV8-ZA binding drift", call. = FALSE)
}

dir.create(private_root, recursive = TRUE)
dir.create(public_root, recursive = TRUE)
dir.create(file.path(private_root, "logs"))
tree_rss <- function(pid) {
  root <- tryCatch(ps::ps_handle(pid), error = function(error) NULL)
  if (is.null(root)) return(0)
  handles <- c(list(root), tryCatch(ps::ps_children(root, recursive = TRUE),
                                    error = function(error) list()))
  sum(vapply(handles, function(handle) tryCatch(
    as.numeric(ps::ps_memory_info(handle)[["rss"]]),
    error = function(error) 0), numeric(1L)))
}
classify_stderr <- function(path, expected_pattern) {
  text <- if (file.exists(path)) paste(readLines(path, warn = FALSE), collapse = "\n") else ""
  if (!nzchar(trimws(text))) return("empty")
  if (grepl(expected_pattern, trimws(text), perl = TRUE) &&
      length(readLines(path, warn = FALSE)) == 1L) "expected_completion" else "unexpected"
}
ledger <- list()
run_child <- function(stage, script, script_args, output_path, expected_pattern) {
  policy <- z_resource[z_resource$stage == stage, , drop = FALSE]
  if (nrow(policy) != 1L) stop("MV8-ZA missing resource policy", call. = FALSE)
  stdout <- file.path(private_root, "logs", paste0(stage, ".stdout"))
  stderr <- file.path(private_root, "logs", paste0(stage, ".stderr"))
  started <- Sys.time()
  process <- processx::process$new(
    Sys.which("Rscript"), c("--vanilla", script, script_args),
    stdout = stdout, stderr = stderr, cleanup_tree = TRUE
  )
  peak <- 0
  cap_failure <- ""
  while (process$is_alive()) {
    Sys.sleep(0.25)
    peak <- max(peak, tree_rss(process$get_pid()))
    elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    if (elapsed > policy$elapsed_cap_seconds) {
      cap_failure <- "elapsed_cap_exceeded"; process$kill_tree()
    } else if (peak > policy$rss_cap_bytes) {
      cap_failure <- "rss_cap_exceeded"; process$kill_tree()
    }
  }
  process$wait(timeout = 5000)
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  status <- process$get_exit_status()
  stderr_class <- classify_stderr(stderr, expected_pattern)
  valid <- identical(status, 0L) && file.exists(output_path) &&
    stderr_class %in% c("empty", "expected_completion") && !nzchar(cap_failure)
  row <- data.frame(
    contract_id = "mv08za_resource_ledger_v1", execution_head = execution_head,
    stage = stage, disposition = if (valid) "completed" else
      if (nzchar(cap_failure)) cap_failure else "failed",
    exit_status = if (is.null(status)) NA_integer_ else status,
    elapsed_seconds = elapsed, peak_process_tree_rss_bytes = peak,
    elapsed_cap_seconds = policy$elapsed_cap_seconds,
    rss_cap_bytes = policy$rss_cap_bytes,
    output_bytes = if (file.exists(output_path)) as.numeric(file.info(output_path)$size) else NA_real_,
    output_sha256 = if (file.exists(output_path)) .mv08z_sha256_file(output_path) else NA_character_,
    stdout_bytes = as.numeric(file.info(stdout)$size),
    stderr_bytes = as.numeric(file.info(stderr)$size), stderr_class = stderr_class,
    workers = 1L, retries = 0L, outcome_label_state = "closed",
    biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
  )
  ledger[[length(ledger) + 1L]] <<- row
  .mv08z_atomic_csv(do.call(rbind, ledger),
                    file.path(public_root, "mv08za-resource-ledger.csv"))
  if (!valid) stop("MV8-ZA child failed closed: ", stage, call. = FALSE)
}
group_order <- as.integer(sentinel$group_order)
chunk_order <- as.integer(sentinel$chunk_order)
primary_root <- file.path(private_root, "primary")
repeat_root <- file.path(private_root, "repeat")
primary_distance <- file.path(primary_root, .mv08z_safe_group(group_order),
                              .mv08z_safe_chunk(chunk_order), "distances.csv")
repeat_distance <- file.path(repeat_root, .mv08z_safe_group(group_order),
                             .mv08z_safe_chunk(chunk_order), "distances.csv")
common <- c(z_root, bindings, s_root, v_root, rust_library,
            group_order, chunk_order)
run_child("sentinel_primary_chunk", "scripts/run_mv08z_landscape_chunk.R",
          c(common, primary_root, execution_head, "sentinel_primary"),
          primary_distance, "^Completed MV8-Z group_[0-9]+/chunk_[0-9]+; pairs=250$")
run_child("sentinel_repeat_chunk", "scripts/run_mv08z_landscape_chunk.R",
          c(common, repeat_root, execution_head, "sentinel_repeat"),
          repeat_distance, "^Completed MV8-Z group_[0-9]+/chunk_[0-9]+; pairs=250$")
oracle_output <- file.path(private_root, "oracle.csv")
run_child("sentinel_canonical_R_oracle", "scripts/run_mv08z_landscape_oracle.R",
          c(z_root, sentinel_private, s_root, v_root, rust_library,
            oracle_output, execution_head, "maximum_burden"), oracle_output,
          "^MV8-Z sentinel oracle passed$")
progress <- data.frame(
  contract_id = "mv08za_progress_v1", execution_head = execution_head,
  state = "sentinel_execution_complete_closure_pending",
  completed_child_processes = 3L, Rust_chunks = 2L,
  pairs_per_Rust_chunk = 250L, canonical_R_oracle_pairs = 1L,
  workers = 1L, retries = 0L, production_pairs = 0L,
  comparison_jobs = 0L, clustering_jobs = 0L, fusion_jobs = 0L,
  label_jobs = 0L, outcome_jobs = 0L, outcome_label_state = "closed",
  biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
)
.mv08z_atomic_csv(progress, file.path(public_root, "mv08za-progress.csv"))
cat("MV8-ZA sentinel execution complete; closure pending\n")
