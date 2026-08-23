#!/usr/bin/env Rscript

# Recovery-only overlay for MV8-P jobs 124:129. The stopped MV8-P v1 roots are
# read-only prerequisites. This runner writes only fresh recovery roots and
# never calls PH or any downstream analysis.

args <- commandArgs(trailingOnly = TRUE)
if (!(length(args) %in% 9:10)) stop(paste(
  "usage: run_mv08pr_source_production_recovery.R <mv08p-audit-dir>",
  "<mv08pr-prefreeze-dir> <original-private-dir> <original-public-dir>",
  "<primary-raw-dir> <added-raw-dir> <hca-count-root>",
  "<recovery-private-dir> <recovery-public-dir> [--resume]"), call. = FALSE)

audit_dir <- normalizePath(args[[1L]], mustWork = TRUE)
recovery_audit_dir <- normalizePath(args[[2L]], mustWork = TRUE)
original_private <- normalizePath(args[[3L]], mustWork = TRUE)
original_public <- normalizePath(args[[4L]], mustWork = TRUE)
primary_raw <- normalizePath(args[[5L]], mustWork = TRUE)
added_raw <- normalizePath(args[[6L]], mustWork = TRUE)
hca_root <- normalizePath(args[[7L]], mustWork = TRUE)
recovery_private <- normalizePath(args[[8L]], mustWork = FALSE)
recovery_public <- normalizePath(args[[9L]], mustWork = FALSE)
resume <- length(args) == 10L && identical(args[[10L]], "--resume")
if (length(args) == 10L && !resume) stop("optional argument must be --resume", call. = FALSE)
if (!resume && (dir.exists(recovery_private) || dir.exists(recovery_public))) {
  stop("refusing to overwrite MV8-PR recovery roots", call. = FALSE)
}
if (resume && (!dir.exists(recovery_private) || !dir.exists(recovery_public))) {
  stop("MV8-PR resume requires both recovery roots", call. = FALSE)
}
for (package in c("digest", "processx")) {
  if (!requireNamespace(package, quietly = TRUE)) stop(package, " is required", call. = FALSE)
}

sha_file <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
atomic_csv <- function(x, path) {
  partial <- paste0(path, ".partial")
  utils::write.csv(x, partial, row.names = FALSE, quote = TRUE, na = "")
  if (!file.rename(partial, path)) stop("failed to publish ", basename(path), call. = FALSE)
}
tree_rss <- function(pid) {
  process <- tryCatch(processx::process$new("ps", c("-eo", "pid=,ppid=,rss="),
    stdout = "|", stderr = "|"), error = function(e) NULL)
  if (is.null(process)) return(NA_real_)
  process$wait(timeout = 10000)
  lines <- process$read_all_output_lines()
  fields <- strsplit(trimws(lines), "[[:space:]]+")
  fields <- fields[lengths(fields) == 3L]
  if (!length(fields)) return(NA_real_)
  table <- do.call(rbind, lapply(fields, as.numeric))
  ids <- as.integer(table[, 1L]); parents <- as.integer(table[, 2L]); rss <- table[, 3L] * 1024
  keep <- pid; prior <- integer()
  while (!identical(sort(keep), sort(prior))) {
    prior <- keep; keep <- unique(c(keep, ids[parents %in% keep]))
  }
  sum(rss[ids %in% keep], na.rm = TRUE)
}
classify_stderr <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) return("empty")
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")
  expected <- grepl("could not find glmGamPoi installed", text, fixed = TRUE) &&
    grepl("Falling back to native (slower) implementation", text, fixed = TRUE) &&
    !grepl("Error in |Execution halted", text)
  if (expected) "known_glmGamPoi_native_fallback" else "unexpected_nonempty"
}

contract <- utils::read.csv(file.path(recovery_audit_dir, "mv08pr-contract.csv"),
  check.names = FALSE, stringsAsFactors = FALSE)
evidence <- utils::read.csv(file.path(recovery_audit_dir, "mv08pr-original-evidence.csv"),
  check.names = FALSE, stringsAsFactors = FALSE)
queue_path <- file.path(audit_dir, "mv08p-remaining-source-queue.csv")
queue <- utils::read.csv(queue_path, check.names = FALSE, stringsAsFactors = FALSE)
original_resource_path <- file.path(original_public, "mv08p-source-production-resource.csv")
original_progress_path <- file.path(original_public, "mv08p-source-production-progress.csv")
original_resource <- utils::read.csv(original_resource_path, check.names = FALSE, stringsAsFactors = FALSE)
original_progress <- utils::read.csv(original_progress_path, check.names = FALSE, stringsAsFactors = FALSE)

if (nrow(contract) != 1L || contract$failed_job_order != 124L ||
    contract$recovery_first_job_order != 124L || contract$recovery_last_job_order != 129L ||
    contract$recovery_jobs != 6L || contract$explicit_retry_jobs != 1L ||
    contract$future_globals_max_size_bytes != 2 * 1024^3 || contract$workers != 1L ||
    contract$elapsed_cap_seconds != 1800L || contract$rss_cap_bytes != 12 * 1024^3 ||
    contract$recovery_state != "authorized_after_commit" ||
    contract$topology_execution_state != "closed" || contract$outcome_label_state != "closed" ||
    isTRUE(contract$biological_outcomes_computed)) stop("MV8-PR recovery contract drift", call. = FALSE)
if (nrow(queue) != 129L || !identical(as.integer(queue$job_order), seq_len(129L)) ||
    nrow(original_resource) != 124L ||
    !identical(as.integer(original_resource$job_order), seq_len(124L)) ||
    !all(original_resource$disposition[1:123] == "completed") ||
    original_resource$disposition[[124L]] != "child_failed" ||
    original_resource$unit_id[[124L]] != queue$unit_id[[124L]] ||
    nrow(original_progress) != 1L || original_progress$state != "stopped" ||
    original_progress$completed_jobs != 123L || original_progress$last_job_order != 124L) {
  stop("MV8-PR stopped-run prerequisite drift", call. = FALSE)
}

evidence_paths <- c(
  file.path(original_public, "mv08p-source-production-resource.csv"),
  file.path(original_public, "mv08p-source-production-progress.csv"),
  file.path(original_private, "logs", paste0(queue$unit_id[[124L]], "__primary-stderr.txt")),
  file.path(original_private, "logs", paste0(queue$unit_id[[124L]], "__primary-stdout.txt")))
if (nrow(evidence) != 4L || !all(file.exists(evidence_paths)) ||
    !identical(as.numeric(file.info(evidence_paths)$size), as.numeric(evidence$bytes)) ||
    !identical(tolower(vapply(evidence_paths, sha_file, character(1L))), tolower(evidence$sha256))) {
  stop("MV8-PR original stopped evidence drift", call. = FALSE)
}
for (i in 1:123) {
  cache_path <- file.path(original_private, "cache", queue$output_file[[i]])
  audit_path <- file.path(original_private, "worker-audit", paste0(queue$unit_id[[i]], "__primary.csv"))
  if (!all(file.exists(c(cache_path, audit_path))) ||
      !identical(tolower(sha_file(cache_path)), tolower(original_resource$cache_sha256[[i]])) ||
      !identical(tolower(sha_file(audit_path)), tolower(original_resource$worker_audit_sha256[[i]]))) {
    stop("MV8-PR accepted-prefix private evidence drift at job ", i, call. = FALSE)
  }
}

dir.create(recovery_private, recursive = TRUE, showWarnings = FALSE)
dir.create(recovery_public, recursive = TRUE, showWarnings = FALSE)
resource_path <- file.path(recovery_public, "mv08pr-source-production-resource.csv")
progress_path <- file.path(recovery_public, "mv08pr-source-production-progress.csv")
resource <- NULL
if (resume) {
  if (!file.exists(resource_path)) stop("MV8-PR resume ledger absent", call. = FALSE)
  resource <- utils::read.csv(resource_path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!nrow(resource) || any(resource$disposition != "completed") ||
      anyDuplicated(resource$job_order) ||
      !identical(as.integer(resource$job_order), seq.int(124L, 123L + nrow(resource))) ||
      any(resource$job_order >= 129L)) stop("MV8-PR resume ledger is not a completed strict prefix", call. = FALSE)
  for (i in seq_len(nrow(resource))) {
    order <- 123L + i
    cache_path <- file.path(recovery_private, "cache", queue$output_file[[order]])
    audit_path <- file.path(recovery_private, "worker-audit", paste0(queue$unit_id[[order]], "__primary.csv"))
    if (!all(file.exists(c(cache_path, audit_path))) ||
        !identical(tolower(sha_file(cache_path)), tolower(resource$cache_sha256[[i]])) ||
        !identical(tolower(sha_file(audit_path)), tolower(resource$worker_audit_sha256[[i]]))) {
      stop("MV8-PR resume private evidence drift at job ", order, call. = FALSE)
    }
  }
}
write_progress <- function(resource, state) {
  completed <- if (is.null(resource)) 0L else sum(resource$disposition == "completed")
  last_order <- if (is.null(resource)) 123L else tail(resource$job_order, 1L)
  out <- data.frame(contract_id = "mv08pr_source_production_progress_v1", state = state,
    authorized_jobs = 6L, completed_jobs = completed, remaining_jobs = 6L - completed,
    overall_completed_jobs = 123L + completed, last_job_order = last_order,
    persistence_computed = FALSE, landscapes_computed = FALSE, clustering_computed = FALSE,
    fusion_computed = FALSE, outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE)
  atomic_csv(out, progress_path)
}
write_progress(resource, "running")
worker <- normalizePath("scripts/run_mv08o_residual_source_worker.R", mustWork = TRUE)
start_order <- if (is.null(resource)) 124L else 124L + nrow(resource)

for (i in seq.int(start_order, 129L)) {
  job <- queue[i, , drop = FALSE]
  is_hca <- job$dataset_scope == "external8"
  source_kind <- if (is_hca) "h5" else "raw"
  if (is_hca) {
    candidates <- list.files(hca_root, pattern = "^filtered_feature_bc_matrix\\.h5$",
      recursive = TRUE, full.names = TRUE)
    normalized <- gsub("\\\\", "/", candidates)
    source_path <- candidates[grepl(paste0("/", job$unit_id, "/"), normalized, fixed = TRUE)]
    if (length(source_path) != 1L) stop("expected one exact-reference HCA matrix", call. = FALSE)
  } else {
    source_dir <- if (job$source_tier == "primary90") primary_raw else added_raw
    source_path <- file.path(source_dir, paste0(job$unit_id, "__raw.rds"))
  }
  if (!file.exists(source_path) ||
      !identical(tolower(sha_file(source_path)), tolower(job$source_sha256))) {
    stop("MV8-PR source identity drift at order ", i, call. = FALSE)
  }
  cache_path <- file.path(recovery_private, "cache", job$output_file)
  audit_path <- file.path(recovery_private, "worker-audit", paste0(job$unit_id, "__primary.csv"))
  stdout_path <- file.path(recovery_private, "logs", paste0(job$unit_id, "__primary-stdout.txt"))
  stderr_path <- file.path(recovery_private, "logs", paste0(job$unit_id, "__primary-stderr.txt"))
  if (any(file.exists(c(cache_path, audit_path, stdout_path, stderr_path)))) {
    stop("refusing to overwrite partial MV8-PR job artifacts at order ", i, call. = FALSE)
  }
  dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(audit_path), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(stdout_path), recursive = TRUE, showWarnings = FALSE)
  child_env <- Sys.getenv()
  child_env[["MV08_FUTURE_GLOBALS_MAX_SIZE_BYTES"]] <- as.character(2 * 1024^3)
  started <- Sys.time(); peak <- 0
  process <- processx::process$new(Sys.which("Rscript"), c("--vanilla", worker, audit_dir,
    source_kind, source_path, job$unit_id, cache_path, audit_path, "primary", queue_path),
    stdout = stdout_path, stderr = stderr_path, cleanup_tree = TRUE, env = child_env)
  disposition <- "running"
  while (process$is_alive()) {
    Sys.sleep(0.25)
    elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    sampled_rss <- tree_rss(process$get_pid())
    if (is.finite(sampled_rss)) peak <- max(peak, sampled_rss)
    if (elapsed > job$elapsed_cap_seconds || peak > job$rss_cap_bytes) {
      disposition <- if (elapsed > job$elapsed_cap_seconds) "elapsed_cap_exceeded" else "rss_cap_exceeded"
      process$kill_tree(); break
    }
  }
  process$wait(timeout = 10000)
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  status <- process$get_exit_status()
  if (disposition == "running") disposition <- if (identical(status, 0L)) "completed" else "child_failed"
  audit <- if (disposition == "completed" && file.exists(audit_path) && file.exists(cache_path))
    utils::read.csv(audit_path, check.names = FALSE, stringsAsFactors = FALSE) else NULL
  required <- if (is.null(audit)) logical() else audit$dataset_scope == "internal124" |
    audit$panel_id == "common475" | (audit$panel_id == "exact500" &
      audit$representation_id == "sct_pearson_residual_all_qc_fit_selected384")
  row <- data.frame(contract_id = "mv08pr_source_production_resource_v1",
    recovery_sequence = i - 123L, job_order = i, unit_id = job$unit_id,
    dataset_scope = job$dataset_scope, source_tier = job$source_tier, fit_cells = job$fit_cells,
    attempt_number = ifelse(i == 124L, 2L, 1L), disposition = disposition,
    exit_status = status, elapsed_seconds = elapsed, peak_rss_bytes = peak,
    elapsed_cap_seconds = job$elapsed_cap_seconds, rss_cap_bytes = job$rss_cap_bytes,
    future_globals_max_size_bytes = 2 * 1024^3,
    cache_bytes = if (file.exists(cache_path)) file.info(cache_path)$size else NA_real_,
    cache_sha256 = if (file.exists(cache_path)) sha_file(cache_path) else NA_character_,
    worker_audit_sha256 = if (file.exists(audit_path)) sha_file(audit_path) else NA_character_,
    stderr_bytes = if (file.exists(stderr_path)) file.info(stderr_path)$size else NA_real_,
    stderr_class = if (file.exists(stderr_path)) classify_stderr(stderr_path) else "absent",
    worker_rows = if (is.null(audit)) 0L else nrow(audit),
    all_required_geometry_valid = if (is.null(audit)) FALSE else
      all(audit$values_finite[required]) && all(audit$zero_variance_gene_count[required] == 0L) &&
      all(audit$correlation_chord_valid[required]),
    persistence_computed = FALSE, landscapes_computed = FALSE, clustering_computed = FALSE,
    fusion_computed = FALSE, outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE)
  resource <- if (is.null(resource)) row else rbind(resource, row)
  atomic_csv(resource, resource_path)
  write_progress(resource, if (disposition == "completed") "running" else "stopped")
  message("MV8-PR ", i, "/129 ", job$unit_id, ": ", disposition)
  if (disposition != "completed" || row$stderr_class == "unexpected_nonempty" ||
      elapsed > job$elapsed_cap_seconds || peak > job$rss_cap_bytes ||
      row$worker_rows != ifelse(is_hca, 4L, 10L) || !row$all_required_geometry_valid) {
    stop("MV8-PR recovery stopped; all evidence preserved", call. = FALSE)
  }
}
write_progress(resource, "source_production_complete_closure_pending")
cat("MV8-PR recovery completed jobs 124-129; merged closure pending; topology remains closed\n")
