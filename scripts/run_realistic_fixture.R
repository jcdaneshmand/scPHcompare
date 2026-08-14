#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args) >= 1L) args[[1]] else file.path("tmp", "realistic-fixture")
seed <- if (length(args) >= 2L) as.integer(args[[2]]) else 20260805L
max_dimension <- if (length(args) >= 3L) as.integer(args[[3]]) else 0L

result <- scPHcompare::run_realistic_fixture(
  output_dir,
  seed = seed,
  num_cores = 1L,
  max_dimension = max_dimension
)
print(result$manifest)
cat("Fixture artifacts:", normalizePath(output_dir, winslash = "/"), "\n")
