#!/usr/bin/env Rscript

# Execute one authorized MV8-V source-produced gene PH record. The parent
# runner supplies the exact source cache and owns monitoring/resource policy.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) stop(paste(
  "usage: run_mv08v_ph_entry.R <mv08u-prefreeze> <job-id>",
  "<source-rds> <common-panel> <ripserr|gudhi> <output-rds>"
), call. = FALSE)
prefreeze <- normalizePath(args[[1L]], mustWork = TRUE)
job_id <- args[[2L]]
source_path <- normalizePath(args[[3L]], mustWork = TRUE)
panel_path <- normalizePath(args[[4L]], mustWork = TRUE)
engine <- match.arg(args[[5L]], c("ripserr", "gudhi"))
output_path <- normalizePath(args[[6L]], mustWork = FALSE)
if (file.exists(output_path)) stop("refusing to overwrite MV8-V PH output", call. = FALSE)
for (package in c("digest", "ripserr")) {
  if (!requireNamespace(package, quietly = TRUE)) stop(package, " required", call. = FALSE)
}
if (engine == "gudhi" && !requireNamespace("TDA", quietly = TRUE)) {
  stop("TDA required for exact GUDHI fallback", call. = FALSE)
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
  file.path(prefreeze, "mv08u-full-ph-queue.csv"),
  check.names = FALSE, stringsAsFactors = FALSE
)
row <- queue[queue$job_id == job_id, , drop = FALSE]
if (nrow(row) != 1L ||
    row$authorization_state != "authorized_after_mv08u_commit" ||
    row$view_kind != "gene_topology_v1" ||
    row$execution_role != "source_produced_gene_ph" ||
    row$outcome_label_state != "closed" ||
    isTRUE(row$biological_outcomes_computed) || row$workers != 1L ||
    row$retries != 0L || row$landscape_authorized ||
    row$comparison_authorized || row$labels_authorized ||
    row$outcomes_authorized) {
  stop("MV8-V PH job is not authorized", call. = FALSE)
}
if (sha_file(source_path) != row$source_cache_sha256) {
  stop("MV8-V source cache SHA-256 drift", call. = FALSE)
}
panel <- utils::read.csv(panel_path, check.names = FALSE,
                         stringsAsFactors = FALSE)
cache <- readRDS(source_path)
mv08s_validate_residual_cache_v1(cache)
view <- mv08s_residual_gene_view_v1(cache, row, panel)
validate_topology_view(view)
if (engine == "gudhi" && row$fallback_trigger != "rss_cap_exceeded_only") {
  stop("MV8-V GUDHI attempt lacks the frozen resource trigger", call. = FALSE)
}
result <- if (engine == "ripserr") {
  run_topology_view_ph(view, max_dim = 1L, threshold = -1, field = 2L)
} else {
  mv07h_run_topology_view_ph_gudhi_v1(
    view, max_dim = 1L, threshold = -1, field = 2L
  )
}
record <- mv08s_new_ph_record_v1(row, cache$cache_key, view, result)
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
partial <- paste0(output_path, ".partial")
if (file.exists(partial)) stop("MV8-V partial output pre-exists", call. = FALSE)
on.exit(if (file.exists(partial)) unlink(partial), add = TRUE)
saveRDS(record, partial, compress = FALSE, version = 3)
if (!file.rename(partial, output_path)) {
  stop("failed to atomically publish MV8-V PH record", call. = FALSE)
}
cat("MV8-V job=", job_id, "; engine=", engine,
    "; points=", length(view$point_ids), "; checks=pass\n", sep = "")
