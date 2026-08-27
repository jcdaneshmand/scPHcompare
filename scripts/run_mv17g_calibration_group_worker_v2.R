#!/usr/bin/env Rscript
options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop(
    paste(
      "usage: run_mv17g_calibration_group_worker_v2.R",
      "<matrix-rds> <cell-or-gene> <null-family> <seed-first>",
      "<replicate-count> <output-rds>"
    ),
    call. = FALSE
  )
}

matrix_path <- normalizePath(args[[1L]], mustWork = TRUE)
view <- match.arg(args[[2L]], c("cell", "gene"))
family <- args[[3L]]
seed_first <- as.integer(args[[4L]])
replicate_count <- as.integer(args[[5L]])
output <- args[[6L]]
if (file.exists(output) || is.na(seed_first) || is.na(replicate_count) ||
    replicate_count < 1L) {
  stop("invalid MV17-G worker-v2 contract", call. = FALSE)
}

source("R/mv17_null_models.R")
source("R/mv17_calibration.R")
source("R/mv17_localization.R")
source("R/mv17_full_calibration.R")
source("R/mv17_full_calibration_geometry_v2.R")

x <- readRDS(matrix_path)
seeds <- if (identical(family, "observed")) {
  0L
} else {
  seq.int(seed_first, length.out = replicate_count)
}
metrics <- mv17g_group_metrics_v2(x, view, family, seeds)
geometry <- mv17g_geometry_id_v2(view)
if (!identical(unique(metrics$geometry), geometry)) {
  stop("MV17-G worker-v2 geometry binding failed", call. = FALSE)
}

out <- list(
  contract_id = "mv17g_group_result_v2",
  view = view,
  geometry = geometry,
  null_family = family,
  seed_first = min(seeds),
  seed_last = max(seeds),
  replicate_count = length(seeds),
  points = nrow(x),
  coordinates = ncol(x),
  matrix_sha256 = digest::digest(
    matrix_path, algo = "sha256", file = TRUE, serialize = FALSE
  ),
  metrics = metrics,
  finite = all(is.finite(metrics$value)),
  labels_opened = FALSE,
  outcomes_opened = FALSE
)
if (!out$finite) {
  stop("non-finite MV17-G worker-v2 result", call. = FALSE)
}
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
partial <- paste0(output, ".partial")
saveRDS(out, partial, version = 3)
if (!file.rename(partial, output)) {
  stop("MV17-G worker-v2 atomic promotion failed", call. = FALSE)
}
