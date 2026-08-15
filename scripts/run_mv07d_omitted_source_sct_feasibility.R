#!/usr/bin/env Rscript

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1")
options(warn = 2)
for (package in c("processx", "ps", "digest")) {
  if (!requireNamespace(package, quietly = TRUE)) stop(package, " required.")
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop("usage: run_mv07d_omitted_source_sct_feasibility.R SENTINEL_CSV ",
       "INDIVIDUAL_SOURCE_ROOT RETAINED_METADATA_SHA256 PRIVATE_ROOT OUTPUT_CSV",
       call. = FALSE)
}
sentinel_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
source_root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
retained_sha <- tolower(trimws(args[[3L]]))
private_root <- args[[4L]]; output_path <- args[[5L]]
if (!grepl("^[0-9a-f]{64}$", retained_sha) || file.exists(output_path)) {
  stop("Retained source SHA is invalid or output already exists.", call. = FALSE)
}
dirs <- file.path(private_root, c("raw", "raw-audit", "raw-log", "overlap",
                                  "sct", "sct-audit", "sct-log"))
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
names(dirs) <- c("raw", "raw_audit", "raw_log", "overlap", "sct",
                 "sct_audit", "sct_log")
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
source("R/provenance_utils.R")
source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/mv05_resource_safe_execution.R")

sentinels <- read.csv(sentinel_path, stringsAsFactors = FALSE, check.names = FALSE)
required <- c("sample_id", "tissue", "post_qc_cells", "selection_boundary",
              "seed", "selected_cells", "accepted_head", "ph_authorized")
if (nrow(sentinels) != 6L || !all(required %in% names(sentinels)) ||
    anyDuplicated(sentinels$sample_id) || any(sentinels$seed != 20260805L) ||
    any(sentinels$selected_cells != 384L) ||
    any(tolower(sentinels$ph_authorized) == "true")) {
  stop("Sentinel prefreeze is invalid.", call. = FALSE)
}
accepted_head <- unique(tolower(sentinels$accepted_head))
if (length(accepted_head) != 1L || !grepl("^[0-9a-f]{40}$", accepted_head)) {
  stop("Sentinel accepted HEAD is invalid.", call. = FALSE)
}
ancestor_status <- system2("git", c("merge-base", "--is-ancestor", accepted_head,
                                     "HEAD"), stdout = FALSE, stderr = FALSE)
if (!identical(ancestor_status, 0L)) {
  stop("Sentinel prefreeze HEAD is not an ancestor of execution HEAD.",
       call. = FALSE)
}

all_paths <- list.files(source_root, pattern = "\\.rds$", recursive = TRUE,
                        full.names = TRUE, ignore.case = TRUE)
base_ids <- tools::file_path_sans_ext(basename(all_paths))
source_paths <- vapply(sentinels$sample_id, function(id) {
  hit <- all_paths[base_ids == id]
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
      process$kill_tree()
      break
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

raw_metrics <- vector("list", nrow(sentinels))
for (i in seq_len(nrow(sentinels))) {
  row <- sentinels[i, , drop = FALSE]; stem <- row$sample_id
  raw_metrics[[i]] <- run_monitored(c(
    "scripts/run_mv05d0_raw_shard_entry.R", source_paths[[i]], dirs[["raw"]],
    file.path(dirs[["raw_audit"]], paste0(stem, ".csv")), stem,
    as.character(row$post_qc_cells), retained_sha, dirs[["overlap"]], source_root),
    file.path(dirs[["raw_log"]], paste0(stem, "-stdout.txt")),
    file.path(dirs[["raw_log"]], paste0(stem, "-stderr.txt")))
  if (raw_metrics[[i]]$disposition != "completed") {
    stop("Raw feasibility failed for ", stem, ": ", raw_metrics[[i]]$disposition,
         call. = FALSE)
  }
  message("MV7-D raw source ", i, "/", nrow(sentinels), ": ", stem)
}

selection_rows <- vector("list", nrow(sentinels))
for (i in seq_len(nrow(sentinels))) {
  id <- sentinels$sample_id[[i]]
  raw_path <- file.path(dirs[["raw"]], paste0(id, "__raw.rds"))
  raw <- readRDS(raw_path); mv05d0_validate_raw_sample_cache_v2(raw)
  selected <- select_matched_cells(colnames(raw$counts), n = 384L,
                                   seed = sentinels$seed[[i]])
  selection_rows[[i]] <- data.frame(
    contract_id = "mv05d0_matched_cell_selection_summary_v1", sample_id = id,
    seed = sentinels$seed[[i]], eligible_cells = ncol(raw$counts),
    selected_cells = length(selected),
    selected_cell_sha256 = attr(selected, "selected_cell_sha256"),
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE)
  rm(raw); invisible(gc())
}
selection <- do.call(rbind, selection_rows)
selection_path <- file.path(private_root, "mv07d-selection.csv")
write_provenance_csv(selection, selection_path)

sct_metrics <- vector("list", nrow(sentinels))
for (i in seq_len(nrow(sentinels))) {
  row <- sentinels[i, , drop = FALSE]; stem <- row$sample_id
  raw_path <- file.path(dirs[["raw"]], paste0(stem, "__raw.rds"))
  sct_metrics[[i]] <- run_monitored(c(
    "scripts/run_mv05d0_sct_cache_entry.R", raw_path, selection_path,
    dirs[["sct"]], file.path(dirs[["sct_audit"]], paste0(stem, ".csv")),
    stem, as.character(row$seed)),
    file.path(dirs[["sct_log"]], paste0(stem, "-stdout.txt")),
    file.path(dirs[["sct_log"]], paste0(stem, "-stderr.txt")))
  if (sct_metrics[[i]]$disposition != "completed") {
    stop("SCT feasibility failed for ", stem, ": ", sct_metrics[[i]]$disposition,
         call. = FALSE)
  }
  message("MV7-D SCT source ", i, "/", nrow(sentinels), ": ", stem)
}

rows <- vector("list", nrow(sentinels))
for (i in seq_len(nrow(sentinels))) {
  row <- sentinels[i, , drop = FALSE]; id <- row$sample_id
  raw_audit <- read.csv(file.path(dirs[["raw_audit"]], paste0(id, ".csv")),
                        stringsAsFactors = FALSE, check.names = FALSE)
  sct_audit <- read.csv(file.path(dirs[["sct_audit"]], paste0(id, ".csv")),
                        stringsAsFactors = FALSE, check.names = FALSE)
  cache <- readRDS(file.path(dirs[["sct"]], sct_audit$private_cache_file))
  matrix <- mv05d0_sct_matrix_from_cache_v1(cache)
  rows[[i]] <- data.frame(
    contract_id = "mv07d_omitted_source_sct_feasibility_v1",
    sample_id = id, tissue = row$tissue,
    selection_boundary = row$selection_boundary,
    expected_post_qc_cells = row$post_qc_cells,
    observed_source_cells = raw_audit$cells,
    observed_source_genes = raw_audit$genes,
    individual_source_file = raw_audit$individual_source_file,
    individual_source_sha256 = raw_audit$individual_source_sha256,
    seed = row$seed, selected_cells = sct_audit$selected_cells,
    selected_cell_sha256 = sct_audit$selected_cell_sha256,
    sct_genes = nrow(matrix), sct_cells = ncol(matrix),
    sct_payload_sha256 = sct_audit$payload_sha256,
    raw_elapsed_seconds = raw_metrics[[i]]$elapsed_seconds,
    raw_peak_rss_bytes = raw_metrics[[i]]$peak_rss_bytes,
    sct_elapsed_seconds = sct_metrics[[i]]$elapsed_seconds,
    sct_peak_rss_bytes = sct_metrics[[i]]$peak_rss_bytes,
    elapsed_cap_seconds = 1800, rss_cap_bytes = 8 * 1024^3,
    finite_sct_payload = all(is.finite(matrix@x)),
    source_sct_feasible = raw_audit$cells == row$post_qc_cells &&
      sct_audit$selected_cells == 384L && ncol(matrix) == 384L &&
      all(is.finite(matrix@x)),
    ph_computed = FALSE, landscape_computed = FALSE,
    biological_outcomes_computed = FALSE, accepted_head = accepted_head,
    stringsAsFactors = FALSE)
  rm(cache, matrix); invisible(gc())
}
result <- do.call(rbind, rows)
if (!all(result$source_sct_feasible) || any(result$ph_computed) ||
    any(result$landscape_computed) || any(result$biological_outcomes_computed) ||
    any(result$raw_peak_rss_bytes > result$rss_cap_bytes) ||
    any(result$sct_peak_rss_bytes > result$rss_cap_bytes) ||
    any(result$raw_elapsed_seconds > result$elapsed_cap_seconds) ||
    any(result$sct_elapsed_seconds > result$elapsed_cap_seconds)) {
  stop("MV7-D feasibility acceptance failed.", call. = FALSE)
}
write_provenance_csv(result, output_path)
message("MV7-D omitted-source/SCT feasibility: 6/6 pass; PH and outcomes remain closed")
