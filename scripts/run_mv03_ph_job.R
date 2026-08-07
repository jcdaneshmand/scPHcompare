#!/usr/bin/env Rscript

options(warn = 2)

source(file.path("R", "toy_baseline.R"), local = FALSE)
source(file.path("R", "dual_view_topology.R"), local = FALSE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: run_mv03_ph_job.R <typed-view-rds> <result-rds>",
       call. = FALSE)
}

view_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
result_path <- args[[2L]]
dir.create(dirname(result_path), recursive = TRUE, showWarnings = FALSE)
view <- readRDS(view_path)
validate_topology_view(view)
if (!isTRUE(view$scientific_eligible)) {
  stop("MV-03 PH jobs require a scientifically eligible typed view.",
       call. = FALSE)
}

started <- proc.time()
result <- run_topology_view_ph(view, max_dim = 1L, threshold = -1, field = 2L)
timing <- proc.time() - started
result$execution <- list(
  elapsed_seconds = unname(timing[["elapsed"]]),
  user_cpu_seconds = unname(timing[["user.self"]]),
  system_cpu_seconds = unname(timing[["sys.self"]]),
  completed_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
)
temporary_path <- paste0(result_path, ".partial")
saveRDS(result, temporary_path, compress = FALSE)
if (file.exists(result_path)) {
  stop("Refusing to overwrite an existing MV-03 PH result.", call. = FALSE)
}
if (!file.rename(temporary_path, result_path)) {
  stop("Failed to atomically publish MV-03 PH result.", call. = FALSE)
}
