#!/usr/bin/env Rscript

# Execute only MV8-P's 129 authorized source fits. Jobs run serially in frozen
# ascending-size order, each in a monitored child process. Private matrices and
# logs remain under the supplied private root; public output contains hashes,
# resource evidence, and aggregate progress only. No PH is called here.

args <- commandArgs(trailingOnly = TRUE)
if (!(length(args) %in% 6:7)) stop(paste(
  "usage: run_mv08p_full_source_production.R <mv08p-audit-dir> <primary-raw-dir>",
  "<added-raw-dir> <hca-count-root> <private-output-dir> <public-output-dir> [--resume]"), call. = FALSE)
audit_dir <- normalizePath(args[[1L]], mustWork = TRUE)
primary_raw <- normalizePath(args[[2L]], mustWork = TRUE)
added_raw <- normalizePath(args[[3L]], mustWork = TRUE)
hca_root <- normalizePath(args[[4L]], mustWork = TRUE)
private_dir <- normalizePath(args[[5L]], mustWork = FALSE)
public_dir <- normalizePath(args[[6L]], mustWork = FALSE)
resume <- length(args) == 7L && identical(args[[7L]], "--resume")
if (length(args) == 7L && !resume) stop("optional argument must be --resume", call. = FALSE)
if (!resume && (dir.exists(private_dir) || dir.exists(public_dir))) {
  stop("refusing to overwrite MV8-P production roots", call. = FALSE)
}
if (resume && (!dir.exists(private_dir) || !dir.exists(public_dir))) {
  stop("MV8-P resume requires both existing output roots", call. = FALSE)
}
dir.create(private_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(public_dir, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
options(warn = 2)
for (package in c("processx", "ps", "digest")) {
  if (!requireNamespace(package, quietly = TRUE)) stop(package, " is required", call. = FALSE)
}

sha_file <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
atomic_csv <- function(x, path) {
  partial <- paste0(path, ".partial")
  utils::write.csv(x, partial, row.names = FALSE, quote = TRUE, na = "")
  if (!file.rename(partial, path)) stop("failed to publish ", basename(path), call. = FALSE)
}
tree_rss <- function(pid) {
  root <- tryCatch(ps::ps_handle(pid), error = function(e) NULL)
  if (is.null(root)) return(0)
  handles <- c(list(root), tryCatch(ps::ps_children(root, recursive = TRUE), error = function(e) list()))
  sum(vapply(handles, function(handle) {
    tryCatch(as.numeric(ps::ps_memory_info(handle)[["rss"]]), error = function(e) 0)
  }, numeric(1L)))
}
classify_stderr <- function(path) {
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")
  if (!nzchar(text)) return("empty")
  known <- grepl("vst.flavor.*glmGamPoi", text) &&
    grepl("Falling back to native \\(slower\\) implementation", text) &&
    !grepl("(^|\\n)(Error|Execution halted)", text)
  if (known) "known_glmGamPoi_native_fallback" else "unexpected_nonempty"
}
write_progress <- function(resource, state) {
  completed <- if (is.null(resource)) 0L else sum(resource$disposition == "completed")
  last_order <- if (is.null(resource) || !nrow(resource)) 0L else max(resource$job_order)
  progress <- data.frame(
    contract_id = "mv08p_source_production_progress_v1", state = state,
    authorized_jobs = 129L, completed_jobs = completed, remaining_jobs = 129L - completed,
    last_job_order = last_order, persistence_computed = FALSE, landscapes_computed = FALSE,
    clustering_computed = FALSE, fusion_computed = FALSE, outcome_label_state = "closed",
    biological_outcomes_computed = FALSE, stringsAsFactors = FALSE)
  atomic_csv(progress, file.path(public_dir, "mv08p-source-production-progress.csv"))
}

queue_path <- file.path(audit_dir, "mv08p-remaining-source-queue.csv")
queue <- utils::read.csv(queue_path, check.names = FALSE, stringsAsFactors = FALSE)
if (nrow(queue) != 129L || !identical(as.integer(queue$job_order), seq_len(129L)) ||
    any(diff(queue$fit_cells) < 0L) ||
    !all(queue$authorization_state == "authorized_after_mv08p_commit") ||
    !all(queue$execution_policy == "serial_one_worker_ascending_fit_cells_no_retry") ||
    !all(queue$workers == 1L & queue$retries == 0L & queue$elapsed_cap_seconds == 1800L &
      queue$rss_cap_bytes == 12 * 1024^3) ||
    any(queue$ph_authorized | queue$landscapes_authorized | queue$comparisons_authorized |
      queue$clustering_authorized | queue$fusion_authorized | queue$labels_authorized |
      queue$outcomes_authorized | queue$biological_outcomes_computed)) {
  stop("MV8-P production queue drift", call. = FALSE)
}
worker <- normalizePath("scripts/run_mv08o_residual_source_worker.R", mustWork = TRUE)
resource_path <- file.path(public_dir, "mv08p-source-production-resource.csv")
resource <- NULL
if (resume) {
  if (!file.exists(resource_path)) stop("MV8-P resume resource ledger is absent", call. = FALSE)
  resource <- utils::read.csv(resource_path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!nrow(resource) || any(resource$disposition != "completed") || anyDuplicated(resource$job_order) ||
      !identical(as.integer(resource$job_order), seq_len(nrow(resource))) ||
      any(resource$job_order >= 129L)) stop("MV8-P resume ledger is not a completed strict prefix", call. = FALSE)
  for (i in seq_len(nrow(resource))) {
    cache_path <- file.path(private_dir, "cache", queue$output_file[[i]])
    audit_path <- file.path(private_dir, "worker-audit", paste0(queue$unit_id[[i]], "__primary.csv"))
    if (!all(file.exists(c(cache_path, audit_path))) ||
        !identical(tolower(sha_file(cache_path)), tolower(resource$cache_sha256[[i]])) ||
        !identical(tolower(sha_file(audit_path)), tolower(resource$worker_audit_sha256[[i]]))) {
      stop("MV8-P resume private evidence drift at job ", i, call. = FALSE)
    }
  }
}
start_index <- if (is.null(resource)) 1L else nrow(resource) + 1L
write_progress(resource, "running")

for (i in seq.int(start_index, nrow(queue))) {
  job <- queue[i, , drop = FALSE]
  is_hca <- job$dataset_scope == "external8"
  source_kind <- if (is_hca) "h5" else "raw"
  if (is_hca) {
    candidates <- list.files(hca_root, pattern = "^filtered_feature_bc_matrix\\.h5$",
      recursive = TRUE, full.names = TRUE)
    normalized_candidates <- gsub("\\\\", "/", candidates)
    source_path <- candidates[grepl(paste0("/", job$unit_id, "/"), normalized_candidates, fixed = TRUE)]
    if (length(source_path) != 1L) stop("expected one exact-reference HCA matrix for ", job$unit_id, call. = FALSE)
  } else {
    source_dir <- if (job$source_tier == "primary90") primary_raw else added_raw
    source_path <- file.path(source_dir, paste0(job$unit_id, "__raw.rds"))
  }
  if (!file.exists(source_path)) stop("required source absent for ", job$unit_id, call. = FALSE)
  cache_path <- file.path(private_dir, "cache", job$output_file)
  audit_path <- file.path(private_dir, "worker-audit", paste0(job$unit_id, "__primary.csv"))
  stdout_path <- file.path(private_dir, "logs", paste0(job$unit_id, "__primary-stdout.txt"))
  stderr_path <- file.path(private_dir, "logs", paste0(job$unit_id, "__primary-stderr.txt"))
  if (any(file.exists(c(cache_path, audit_path, stdout_path, stderr_path)))) {
    stop("refusing to overwrite partial MV8-P job artifacts at order ", i, call. = FALSE)
  }
  dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(audit_path), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(stdout_path), recursive = TRUE, showWarnings = FALSE)
  started <- Sys.time(); peak <- 0
  process <- processx::process$new(Sys.which("Rscript"), c("--vanilla", worker, audit_dir,
    source_kind, source_path, job$unit_id, cache_path, audit_path, "primary", queue_path),
    stdout = stdout_path, stderr = stderr_path, cleanup_tree = TRUE)
  disposition <- "running"
  while (process$is_alive()) {
    Sys.sleep(0.25)
    elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    peak <- max(peak, tree_rss(process$get_pid()))
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
  required <- if (is.null(audit)) logical() else
    audit$dataset_scope == "internal124" | audit$panel_id == "common475" |
      (audit$panel_id == "exact500" &
        audit$representation_id == "sct_pearson_residual_all_qc_fit_selected384")
  row <- data.frame(
    contract_id = "mv08p_source_production_resource_v1", job_order = i,
    unit_id = job$unit_id, dataset_scope = job$dataset_scope, source_tier = job$source_tier,
    fit_cells = job$fit_cells, disposition = disposition, exit_status = status,
    elapsed_seconds = elapsed, peak_rss_bytes = peak,
    elapsed_cap_seconds = job$elapsed_cap_seconds, rss_cap_bytes = job$rss_cap_bytes,
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
  message("MV8-P ", i, "/129 ", job$unit_id, ": ", disposition)
  if (disposition != "completed" || row$stderr_class == "unexpected_nonempty" ||
      elapsed > job$elapsed_cap_seconds || peak > job$rss_cap_bytes ||
      row$worker_rows != ifelse(is_hca, 4L, 10L) || !row$all_required_geometry_valid) {
    stop("MV8-P source production stopped; partial evidence preserved", call. = FALSE)
  }
}
write_progress(resource, "source_production_complete_closure_pending")
cat("MV8-P source production completed 129/129; private caches require closure rehash; topology remains closed\n")
