#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) stop(paste(
  "usage: run_mv08s_ph_entry.R <mv08s-prefreeze> <job-id>",
  "<source-rds> <common-panel> <ripserr|gudhi> <output-rds>"
), call. = FALSE)
prefreeze <- normalizePath(args[[1L]], mustWork = TRUE)
job_id <- args[[2L]]
source_path <- normalizePath(args[[3L]], mustWork = TRUE)
panel_path <- normalizePath(args[[4L]], mustWork = TRUE)
engine <- match.arg(args[[5L]], c("ripserr", "gudhi"))
output_path <- normalizePath(args[[6L]], mustWork = FALSE)
if (file.exists(output_path)) stop("refusing to overwrite MV8-S PH output", call. = FALSE)
for (package in c("digest", "ripserr")) {
  if (!requireNamespace(package, quietly = TRUE)) stop(package, " required", call. = FALSE)
}
if (engine == "gudhi" && !requireNamespace("TDA", quietly = TRUE)) {
  stop("TDA required for MV8-S GUDHI fallback", call. = FALSE)
}
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
options(warn = 2)
source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/mv07g_sentinel.R")
source("R/mv07h_full_topology.R")
source("R/mv08s_ph_sentinel.R")

sha_file <- function(path) digest::digest(
  file = path, algo = "sha256", serialize = FALSE
)
queue <- utils::read.csv(
  file.path(prefreeze, "mv08s-ph-sentinel-queue.csv"),
  check.names = FALSE, stringsAsFactors = FALSE
)
row <- queue[queue$job_id == job_id, , drop = FALSE]
if (nrow(row) != 1L ||
    row$authorization_state != "authorized_after_mv08s_commit" ||
    row$outcome_label_state != "closed" ||
    isTRUE(row$biological_outcomes_computed) || row$workers != 1L ||
    row$retries != 0L) {
  stop("MV8-S PH job is not authorized", call. = FALSE)
}
panel <- utils::read.csv(panel_path, check.names = FALSE, stringsAsFactors = FALSE)
source_record <- readRDS(source_path)
if (row$execution_role == "source_produced_gene_ph") {
  binding <- utils::read.csv(
    file.path(prefreeze, "mv08s-source-cache-bindings.csv"),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  binding <- binding[binding$unit_id == row$unit_id, , drop = FALSE]
  if (nrow(binding) != 1L || sha_file(source_path) != binding$cache_sha256) {
    stop("MV8-S source-produced cache file drift", call. = FALSE)
  }
  mv08s_validate_residual_cache_v1(source_record, binding)
  view <- mv08s_residual_gene_view_v1(source_record, row, panel)
  source_cache_key <- source_record$cache_key
} else {
  binding <- utils::read.csv(
    file.path(prefreeze, "mv08s-external-input-bindings.csv"),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  binding <- binding[binding$unit_id == row$unit_id, , drop = FALSE]
  mv08s_validate_baseline_record_v1(source_record, binding)
  view <- source_record$views[[row$view_kind]]
  source_cache_key <- source_record$cache_key
}
validate_topology_view(view)
if (engine == "gudhi" && view$view_id != "gene_topology_v1") {
  stop("MV8-S exact GUDHI resource fallback is gene-view only", call. = FALSE)
}
result <- if (engine == "ripserr") {
  run_topology_view_ph(view, max_dim = 1L, threshold = -1, field = 2L)
} else {
  mv07h_run_topology_view_ph_gudhi_v1(
    view, max_dim = 1L, threshold = -1, field = 2L
  )
}
record <- mv08s_new_ph_record_v1(row, source_cache_key, view, result)
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
partial <- tempfile(pattern = paste0(basename(output_path), "."),
                    tmpdir = dirname(output_path))
on.exit(if (file.exists(partial)) unlink(partial), add = TRUE)
saveRDS(record, partial, compress = FALSE, version = 3)
if (!file.rename(partial, output_path)) {
  stop("failed to atomically publish MV8-S PH record", call. = FALSE)
}
cat("MV8-S PH job=", job_id, "; engine=", engine,
    "; points=", length(view$point_ids), "; checks=pass\n", sep = "")
