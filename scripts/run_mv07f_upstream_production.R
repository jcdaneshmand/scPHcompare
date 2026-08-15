#!/usr/bin/env Rscript

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1")
options(warn = 2)
for (package in c("processx", "ps", "digest", "Matrix")) {
  if (!requireNamespace(package, quietly = TRUE)) stop(package, " required.")
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop("usage: run_mv07f_upstream_production.R PREFREEZE_DIR SOURCE_ROOT ",
    "RETAINED_SHA256 PRIVATE_ROOT OUTPUT_DIR EXPECTED_PREFREEZE_HEAD",
    call. = FALSE)
}
prefreeze_dir <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
source_root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
retained_sha <- tolower(trimws(args[[3L]])); private_root <- args[[4L]]
output_dir <- args[[5L]]; prefreeze_head <- tolower(trimws(args[[6L]]))
if (!grepl("^[0-9a-f]{64}$", retained_sha) ||
    !grepl("^[0-9a-f]{40}$", prefreeze_head)) stop("Invalid frozen hash.")
source("R/provenance_utils.R")
source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/mv05_resource_safe_execution.R")
source("R/mv07f_upstream_production.R")
source("R/mv07f_validation_utils.R")
sha <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
head <- tolower(trimws(system2("git", c("rev-parse", "HEAD"), stdout = TRUE)))
ancestor <- system2("git", c("merge-base", "--is-ancestor", prefreeze_head, head),
                    stdout = FALSE, stderr = FALSE)
if (!identical(ancestor, 0L)) stop("MV7-F prefreeze is not an execution ancestor.")
source_freeze <- read.csv(file.path(prefreeze_dir, "mv07f-source-freeze.csv"),
  stringsAsFactors = FALSE, check.names = FALSE)
if (!all(source_freeze$accepted_head == prefreeze_head) ||
    any(vapply(seq_len(nrow(source_freeze)), function(i)
      sha(source_freeze$artifact_locator[[i]]) != source_freeze$sha256[[i]],
      logical(1L)))) stop("MV7-F source drift detected.")
queue <- read.csv(file.path(prefreeze_dir, "mv07f-upstream-queue.csv"),
  stringsAsFactors = FALSE, check.names = FALSE)
caps <- read.csv(file.path(prefreeze_dir, "mv07f-resource-caps.csv"),
  stringsAsFactors = FALSE, check.names = FALSE)
if (nrow(queue) != 204L || sum(queue$stage == "raw") != 34L ||
    sum(queue$stage == "sct") != 170L) stop("MV7-F queue drift.")
if (dir.exists(output_dir) && length(list.files(output_dir, all.files=TRUE,
                                                no..=TRUE))) {
  stop("MV7-F output directory must be empty.", call. = FALSE)
}

dirs <- file.path(private_root, c("raw", "sct", "mv07f-raw-audit",
  "mv07f-sct-audit", "mv07f-raw-log", "mv07f-sct-log", "receipts",
  "overlap"))
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
names(dirs) <- c("raw", "sct", "raw_audit", "sct_audit", "raw_log",
                 "sct_log", "receipts", "overlap")
receipt <- data.frame(
  contract_id = "mv07f_execution_receipt_v1", execution_head = head,
  prefreeze_head = prefreeze_head,
  queue_sha256 = sha(file.path(prefreeze_dir, "mv07f-upstream-queue.csv")),
  source_freeze_sha256 = sha(file.path(prefreeze_dir, "mv07f-source-freeze.csv")),
  retained_metadata_sha256 = retained_sha, receipt_before_source_parse = TRUE,
  label_access = FALSE, panel_fit = FALSE, pca = FALSE, ph = FALSE,
  landscape = FALSE, outcomes = FALSE, stringsAsFactors = FALSE)
receipt_path <- file.path(dirs[["receipts"]], paste0(head, ".csv"))
if (file.exists(receipt_path)) {
  old <- read.csv(receipt_path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!mv07f_provenance_record_matches_v1(old, receipt)) {
    stop("Existing execution receipt drifted.")
  }
} else {
  write_provenance_csv(receipt, receipt_path)
}

all_paths <- list.files(source_root, pattern = "[.]rds$", recursive = TRUE,
  full.names = TRUE, ignore.case = TRUE)
all_ids <- tools::file_path_sans_ext(basename(all_paths))
raw_queue <- queue[queue$stage == "raw", , drop = FALSE]
sct_queue <- queue[queue$stage == "sct", , drop = FALSE]
source_paths <- vapply(raw_queue$sample_id, function(id) {
  hit <- all_paths[all_ids == id]
  if (length(hit) != 1L) stop("Expected one individual source for ", id)
  hit
}, character(1L))

tree_rss <- function(pid) {
  root <- tryCatch(ps::ps_handle(pid), error = function(e) NULL)
  if (is.null(root)) return(0)
  handles <- c(list(root), tryCatch(ps::ps_children(root, recursive = TRUE),
                                   error = function(e) list()))
  sum(vapply(handles, function(handle) tryCatch(
    as.numeric(ps::ps_memory_info(handle)[["rss"]]), error = function(e) 0),
    numeric(1L)))
}
run_monitored <- function(script_args, stdout_path, stderr_path,
                          elapsed_cap = 1800, rss_cap = 8 * 1024^3) {
  started <- Sys.time(); peak <- 0
  process <- processx::process$new(Sys.which("Rscript"),
    c("--vanilla", script_args), stdout = stdout_path, stderr = stderr_path,
    cleanup_tree = TRUE)
  disposition <- "running"
  while (process$is_alive()) {
    Sys.sleep(0.25)
    elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    peak <- max(peak, tree_rss(process$get_pid()))
    if (elapsed > elapsed_cap || peak > rss_cap) {
      disposition <- if (elapsed > elapsed_cap) "elapsed_cap_exceeded" else
        "rss_cap_exceeded"
      process$kill_tree(); break
    }
  }
  process$wait(timeout = 5000)
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  status <- process$get_exit_status()
  if (identical(disposition, "running")) {
    disposition <- if (identical(status, 0L)) "completed" else "child_failed"
  }
  list(disposition = disposition, exit_status = status,
       elapsed_seconds = elapsed, peak_rss_bytes = peak)
}

raw_metrics <- vector("list", nrow(raw_queue))
for (i in seq_len(nrow(raw_queue))) {
  row <- raw_queue[i, , drop = FALSE]; id <- row$sample_id
  raw_metrics[[i]] <- run_monitored(c("scripts/run_mv05d0_raw_shard_entry.R",
    source_paths[[i]], dirs[["raw"]],
    file.path(dirs[["raw_audit"]], paste0(id, ".csv")), id,
    as.character(row$expected_post_qc_cells), retained_sha,
    dirs[["overlap"]], source_root),
    file.path(dirs[["raw_log"]], paste0(id, "-stdout.txt")),
    file.path(dirs[["raw_log"]], paste0(id, "-stderr.txt")))
  if (raw_metrics[[i]]$disposition != "completed") {
    stop("MV7-F raw child failed for ", id, ": ",
         raw_metrics[[i]]$disposition, call. = FALSE)
  }
  message("MV7-F raw ", i, "/", nrow(raw_queue), ": ", id)
}

selection_rows <- vector("list", nrow(sct_queue))
for (i in seq_len(nrow(sct_queue))) {
  row <- sct_queue[i, , drop = FALSE]; id <- row$sample_id
  raw <- readRDS(file.path(dirs[["raw"]], paste0(id, "__raw.rds")))
  mv05d0_validate_raw_sample_cache_v2(raw)
  selected <- select_matched_cells(colnames(raw$counts), n = 384L,
                                   seed = row$seed)
  selection_rows[[i]] <- data.frame(
    contract_id = "mv07f_matched_cell_selection_v1", sample_id = id,
    seed = as.integer(row$seed), eligible_cells = ncol(raw$counts),
    selected_cells = 384L,
    selected_cell_sha256 = attr(selected, "selected_cell_sha256"),
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE)
  rm(raw, selected); invisible(gc())
}
selection <- do.call(rbind, selection_rows)
selection_path <- file.path(private_root, "mv07f-selection.csv")
if (file.exists(selection_path)) {
  old <- read.csv(selection_path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!identical(old, selection)) stop("Existing MV7-F selection drifted.")
} else {
  partial <- tempfile("mv07f-selection-", tmpdir = private_root)
  write_provenance_csv(selection, partial)
  if (!file.rename(partial, selection_path)) stop("Atomic selection publish failed.")
}

sct_metrics <- vector("list", nrow(sct_queue))
for (i in seq_len(nrow(sct_queue))) {
  row <- sct_queue[i, , drop = FALSE]; id <- row$sample_id
  key <- paste0(id, "__", row$seed)
  sct_metrics[[i]] <- run_monitored(c("scripts/run_mv05d0_sct_cache_entry.R",
    file.path(dirs[["raw"]], paste0(id, "__raw.rds")), selection_path,
    dirs[["sct"]], file.path(dirs[["sct_audit"]], paste0(key, ".csv")),
    id, as.character(row$seed)),
    file.path(dirs[["sct_log"]], paste0(key, "-stdout.txt")),
    file.path(dirs[["sct_log"]], paste0(key, "-stderr.txt")))
  if (sct_metrics[[i]]$disposition != "completed") {
    stop("MV7-F SCT child failed for ", key, ": ",
         sct_metrics[[i]]$disposition, call. = FALSE)
  }
  message("MV7-F SCT ", i, "/", nrow(sct_queue), ": ", key)
}

raw_rows <- vector("list", nrow(raw_queue))
for (i in seq_len(nrow(raw_queue))) {
  row <- raw_queue[i, , drop = FALSE]; id <- row$sample_id
  audit <- read.csv(file.path(dirs[["raw_audit"]], paste0(id, ".csv")),
                    stringsAsFactors = FALSE, check.names = FALSE)
  raw_rows[[i]] <- data.frame(
    contract_id = "mv07f_raw_production_v1", sample_id = id,
    expected_post_qc_cells = row$expected_post_qc_cells,
    observed_genes = audit$genes, observed_cells = audit$cells,
    individual_source_file = audit$individual_source_file,
    individual_source_sha256 = audit$individual_source_sha256,
    counts_sha256 = audit$counts_sha256,
    private_cache_file = audit$private_raw_cache_file,
    private_cache_bytes = audit$private_raw_cache_size_bytes,
    private_cache_sha256 = audit$private_raw_cache_sha256,
    disposition = audit$disposition,
    elapsed_seconds = raw_metrics[[i]]$elapsed_seconds,
    peak_rss_bytes = raw_metrics[[i]]$peak_rss_bytes,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE)
}
raw_result <- do.call(rbind, raw_rows)
sct_rows <- vector("list", nrow(sct_queue))
for (i in seq_len(nrow(sct_queue))) {
  row <- sct_queue[i, , drop = FALSE]; id <- row$sample_id
  key <- paste0(id, "__", row$seed)
  audit <- read.csv(file.path(dirs[["sct_audit"]], paste0(key, ".csv")),
                    stringsAsFactors = FALSE, check.names = FALSE)
  cache <- readRDS(file.path(dirs[["sct"]], audit$private_cache_file))
  matrix <- mv05d0_sct_matrix_from_cache_v1(cache)
  sct_rows[[i]] <- data.frame(
    contract_id = "mv07f_sct_production_v1", sample_id = id,
    seed = as.integer(row$seed), selected_cells = audit$selected_cells,
    selected_cell_sha256 = audit$selected_cell_sha256,
    normalization_cache_key = audit$normalization_cache_key,
    payload_sha256 = audit$payload_sha256, sct_genes = nrow(matrix),
    sct_cells = ncol(matrix), finite_sct_payload = all(is.finite(matrix@x)),
    private_cache_file = audit$private_cache_file,
    private_cache_bytes = audit$private_cache_size_bytes,
    private_cache_sha256 = audit$private_cache_sha256,
    disposition = audit$disposition,
    elapsed_seconds = sct_metrics[[i]]$elapsed_seconds,
    peak_rss_bytes = sct_metrics[[i]]$peak_rss_bytes,
    seurat_version = audit$seurat_version,
    sctransform_version = audit$sctransform_version,
    matrix_version = audit$matrix_version, r_version = audit$r_version,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE)
  rm(cache, matrix); invisible(gc())
}
sct_result <- do.call(rbind, sct_rows)

cache_paths <- c(file.path(dirs[["raw"]], raw_result$private_cache_file),
                 file.path(dirs[["sct"]], sct_result$private_cache_file))
if (length(unique(cache_paths)) != 204L || any(!file.exists(cache_paths))) {
  stop("MV7-F final cache axis is incomplete.", call. = FALSE)
}
unexpected_raw <- setdiff(list.files(dirs[["raw"]], full.names = TRUE),
                          file.path(dirs[["raw"]], raw_result$private_cache_file))
unexpected_sct <- setdiff(list.files(dirs[["sct"]], full.names = TRUE),
                          file.path(dirs[["sct"]], sct_result$private_cache_file))
if (length(unexpected_raw) || length(unexpected_sct) ||
    any(!sct_result$finite_sct_payload)) stop("MV7-F partial or invalid state.")
worker_seconds <- sum(vapply(raw_metrics, `[[`, numeric(1L), "elapsed_seconds")) +
  sum(vapply(sct_metrics, `[[`, numeric(1L), "elapsed_seconds"))
peak_rss <- max(c(vapply(raw_metrics, `[[`, numeric(1L), "peak_rss_bytes"),
                  vapply(sct_metrics, `[[`, numeric(1L), "peak_rss_bytes")))
storage <- sum(as.numeric(file.info(cache_paths)$size))
if (worker_seconds > caps$elapsed_cap_seconds[caps$scope == "aggregate_worker"] ||
    storage > caps$storage_cap_bytes[caps$scope == "aggregate_storage"] ||
    peak_rss > max(caps$rss_cap_bytes, na.rm = TRUE)) {
  stop("MV7-F aggregate resource cap exceeded.", call. = FALSE)
}
resource <- data.frame(
  contract_id = "mv07f_resource_summary_v1", raw_worker_seconds =
    sum(vapply(raw_metrics, `[[`, numeric(1L), "elapsed_seconds")),
  sct_worker_seconds = sum(vapply(sct_metrics, `[[`, numeric(1L), "elapsed_seconds")),
  total_worker_seconds = worker_seconds, peak_process_tree_rss_bytes = peak_rss,
  unique_cache_storage_bytes = storage, aggregate_elapsed_cap_seconds = 14400,
  rss_cap_bytes = 8 * 1024^3, storage_cap_bytes = 4 * 1024^3,
  resource_gate_passed = TRUE, stringsAsFactors = FALSE)
summary <- data.frame(
  contract_id = "mv07f_production_summary_v1", execution_head = head,
  prefreeze_head = prefreeze_head, raw_jobs = nrow(raw_result),
  sct_jobs = nrow(sct_result), selections = nrow(selection),
  raw_built = sum(raw_result$disposition == "built_atomic"),
  raw_reused = sum(raw_result$disposition == "reuse_validated"),
  sct_built = sum(sct_result$disposition == "built_atomic"),
  sct_reused = sum(sct_result$disposition == "reuse_validated"),
  partial_files = 0L, label_access = FALSE, panel_fit_jobs = 0L,
  pca_jobs = 0L, ph_jobs = 0L, landscape_jobs = 0L, outcome_jobs = 0L,
  stringsAsFactors = FALSE)
manifest <- rbind(
  data.frame(cache_kind = "raw", sample_id = raw_result$sample_id,
    seed = NA_integer_, private_cache_file = raw_result$private_cache_file,
    private_cache_bytes = raw_result$private_cache_bytes,
    private_cache_sha256 = raw_result$private_cache_sha256,
    payload_sha256 = raw_result$counts_sha256, stringsAsFactors = FALSE),
  data.frame(cache_kind = "sct", sample_id = sct_result$sample_id,
    seed = sct_result$seed, private_cache_file = sct_result$private_cache_file,
    private_cache_bytes = sct_result$private_cache_bytes,
    private_cache_sha256 = sct_result$private_cache_sha256,
    payload_sha256 = sct_result$payload_sha256, stringsAsFactors = FALSE))
manifest$contract_id <- "mv07f_cache_manifest_v1"
manifest <- manifest[c("contract_id", "cache_kind", "sample_id", "seed",
  "private_cache_file", "private_cache_bytes", "private_cache_sha256",
  "payload_sha256")]
snapshot <- data.frame(private_cache_file = manifest$private_cache_file,
  cache_kind = manifest$cache_kind, bytes = as.numeric(file.info(cache_paths)$size),
  sha256 = vapply(cache_paths, sha, character(1L)),
  mtime_numeric = as.numeric(file.info(cache_paths)$mtime), stringsAsFactors = FALSE)
snapshot_path <- file.path(private_root, "mv07f-accepted-cache-snapshot.csv")
if (file.exists(snapshot_path)) {
  old <- read.csv(snapshot_path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!identical(old, snapshot)) stop("Existing accepted cache snapshot drifted.")
} else {
  write_provenance_csv(snapshot, snapshot_path)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write_provenance_csv(receipt, file.path(output_dir, "mv07f-execution-provenance.csv"))
write_provenance_csv(raw_result, file.path(output_dir, "mv07f-raw-production.csv"))
write_provenance_csv(sct_result, file.path(output_dir, "mv07f-sct-production.csv"))
write_provenance_csv(selection, file.path(output_dir, "mv07f-selection-summary.csv"))
write_provenance_csv(manifest, file.path(output_dir, "mv07f-cache-manifest.csv"))
write_provenance_csv(resource, file.path(output_dir, "mv07f-resource-summary.csv"))
write_provenance_csv(summary, file.path(output_dir, "mv07f-production-summary.csv"))
message("MV7-F upstream production complete: 34 raw + 170 SCT caches")
