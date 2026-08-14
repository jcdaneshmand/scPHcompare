#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% 4:5) {
  stop(paste("usage: run_mv04_diagram_sensitivities.R MANIFEST INTERVALS",
             "OUTPUT ENVIRONMENT [bottleneck|wasserstein_p1|both]"))
}
suppressPackageStartupMessages(library(TDA))
manifest <- read.csv(args[[1L]], stringsAsFactors = FALSE, check.names = FALSE)
intervals <- read.csv(args[[2L]], stringsAsFactors = FALSE, check.names = FALSE)
selected_strata <- Sys.getenv("MV04_STRATA", unset = "")
if (nzchar(selected_strata)) {
  selected_strata <- strsplit(selected_strata, ",", fixed = TRUE)[[1L]]
  if (!all(selected_strata %in% manifest$stratum_id)) {
    stop("MV04_STRATA contains an unknown stratum ID.")
  }
  manifest <- manifest[manifest$stratum_id %in% selected_strata, , drop = FALSE]
}
selected_samples <- Sys.getenv("MV04_SAMPLE_IDS", unset = "")
if (nzchar(selected_samples)) {
  selected_samples <- strsplit(selected_samples, ",", fixed = TRUE)[[1L]]
  manifest <- manifest[manifest$sample_id %in% selected_samples, , drop = FALSE]
  if (any(table(manifest$stratum_id) < 2L)) {
    stop("MV04_SAMPLE_IDS must retain at least two samples per selected stratum.")
  }
}

diagram <- function(id) {
  part <- intervals[intervals$diagram_id == id, , drop = FALSE]
  cbind(part$homology_dimension, part$birth, part$death)
}

started_all <- proc.time()[["elapsed"]]
mode <- if (length(args) == 5L) args[[5L]] else "both"
methods <- switch(mode, bottleneck = "bottleneck",
                  wasserstein_p1 = "wasserstein_p1",
                  both = c("bottleneck", "wasserstein_p1"),
                  stop("Unknown sensitivity mode: ", mode))
rows <- list()
index <- 0L
for (stratum in sort(unique(manifest$stratum_id))) {
  part <- manifest[manifest$stratum_id == stratum, , drop = FALSE]
  part <- part[order(part$sample_id), , drop = FALSE]
  combinations <- utils::combn(seq_len(nrow(part)), 2L)
  diagrams <- lapply(part$diagram_id, diagram)
  for (dimension in 0:1) for (method in methods) {
    for (column in seq_len(ncol(combinations))) {
      left <- combinations[1L, column]
      right <- combinations[2L, column]
      started <- proc.time()[["elapsed"]]
      value <- if (method == "bottleneck") {
        TDA::bottleneck(diagrams[[left]], diagrams[[right]], dimension = dimension)
      } else {
        TDA::wasserstein(diagrams[[left]], diagrams[[right]], p = 1,
                         dimension = dimension)
      }
      index <- index + 1L
      rows[[index]] <- data.frame(
        method_id = if (method == "bottleneck") "bottleneck_v1" else "wasserstein_p1_v1",
        stratum_id = stratum, homology_dimension = paste0("H", dimension),
        first_sample_id = part$sample_id[[left]],
        second_sample_id = part$sample_id[[right]],
        first_diagram_id = part$diagram_id[[left]],
        second_diagram_id = part$diagram_id[[right]], distance = value,
        pair_seconds = proc.time()[["elapsed"]] - started,
        status = "completed", stringsAsFactors = FALSE
      )
      write.csv(do.call(rbind, rows), args[[3L]], row.names = FALSE)
    }
  }
}
result <- do.call(rbind, rows)
stratum_sizes <- table(manifest$stratum_id)
expected_rows <- sum(stratum_sizes * (stratum_sizes - 1L) / 2L) *
  2L * length(methods)
stopifnot(nrow(result) == expected_rows, all(is.finite(result$distance)),
          all(result$distance >= 0))
status_lines <- readLines("/proc/self/status", warn = FALSE)
peak_line <- status_lines[grepl("^VmHWM:", status_lines)]
peak_rss_bytes <- if (length(peak_line)) {
  as.numeric(gsub("[^0-9]", "", peak_line[[1L]])) * 1024
} else NA_real_
environment <- data.frame(
  tda_version = as.character(packageVersion("TDA")),
  r_version = as.character(getRversion()),
  elapsed_seconds = proc.time()[["elapsed"]] - started_all,
  peak_rss_bytes = peak_rss_bytes, method_scope = mode,
  row_count = nrow(result), stringsAsFactors = FALSE
)
write.csv(environment, args[[4L]], row.names = FALSE)
message("Completed ", nrow(result), " sensitivity distances.")
