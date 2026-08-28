#!/usr/bin/env Rscript
options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop(
    paste(
      "usage: run_mv17gr_exact_h1_profile_worker.R",
      "<matrix-rds> <null-family> <seed> <engine> <output-rds>"
    ), call. = FALSE
  )
}
matrix_path <- normalizePath(args[[1L]], mustWork = TRUE)
null_family <- args[[2L]]
seed <- as.integer(args[[3L]])
engine <- args[[4L]]
output <- args[[5L]]
if (file.exists(output) || file.exists(paste0(output, ".partial")) ||
    is.na(seed)) {
  stop("invalid MV17-GR worker target/seed", call. = FALSE)
}

source("R/mv17_null_models.R")
source("R/mv17_calibration.R")
source("R/mv17_localization.R")
source("R/mv17_full_calibration_geometry_v2.R")
source("R/mv17gr_exact_h1_resource_profile.R")

x <- readRDS(matrix_path)
coordinates <- mv17gr_materialize_gene_case_v1(x, null_family, seed)
result <- mv17gr_run_exact_h1_v1(coordinates, engine)
result$null_family <- null_family
result$seed <- seed
result$matrix_sha256 <- digest::digest(
  matrix_path, algo = "sha256", file = TRUE, serialize = FALSE
)
result$points <- nrow(coordinates)
result$coordinates <- ncol(coordinates)
result$geometry <- "euclidean_correlation_chord_v1"
result$labels_opened <- FALSE
result$outcomes_opened <- FALSE
if (!result$finite) stop("non-finite MV17-GR result", call. = FALSE)
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
partial <- paste0(output, ".partial")
saveRDS(result, partial, version = 3)
if (!file.rename(partial, output)) {
  stop("MV17-GR worker atomic promotion failed", call. = FALSE)
}
