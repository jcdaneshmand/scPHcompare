#!/usr/bin/env Rscript

options(digits = 17)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: landscape_reference_worker.R <exact|adaptive> <intervals>",
       call. = FALSE)
}
method <- match.arg(args[[1]], c("exact", "adaptive"))
intervals <- as.integer(args[[2]])
if (is.na(intervals) || intervals < 1L) stop("intervals must be positive.")

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE)
project_dir <- dirname(dirname(script_path))
source(file.path(project_dir, "R", "landscape_contract.R"), local = .GlobalEnv)
source(file.path(project_dir, "R", "landscape_reference.R"), local = .GlobalEnv)

births <- seq_len(intervals) / (20 * intervals)
deaths <- 3 - births
first <- cbind(0, births, deaths)
second <- cbind(0, births + 0.01, deaths + 0.015)
started <- proc.time()[["elapsed"]]
result <- landscape_reference_distance(
  first, second, method = method,
  exact_max_intervals = intervals + 1L,
  abs_tol = 1e-8, rel_tol = 1e-8
)
elapsed <- proc.time()[["elapsed"]] - started
rss_bytes <- if (.Platform$OS.type == "unix") {
  status <- readLines("/proc/self/status", warn = FALSE)
  value <- sub("^VmHWM:[[:space:]]+([0-9]+)[[:space:]]+kB$", "\\1",
               grep("^VmHWM:", status, value = TRUE))
  as.numeric(value) * 1024
} else NA_real_
cat(
  result$distances[["combined"]], elapsed, rss_bytes,
  sep = "\t"
)
cat("\n")
