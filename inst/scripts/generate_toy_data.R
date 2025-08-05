#!/usr/bin/env Rscript

# Script to generate toy sparse expression matrices and metadata
# for scPHcompare package

set.seed(123)

suppressPackageStartupMessages({
  library(Matrix)
})

# Determine package root based on script location
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args)
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", args[file_arg]))
} else {
  # Fallback when running via source()
  normalizePath(sys.frames()[[1]]$ofile)
}

pkg_root <- normalizePath(file.path(dirname(script_path), "..", ".."))

out_dir <- file.path(pkg_root, "inst", "extdata", "toy")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

samples <- data.frame(
  sample = c("sample1", "sample2", "sample3"),
  sra    = c("SRR000001", "SRR000002", "SRR000003"),
  tissue = c("brain", "heart", "liver"),
  approach = c("scRNA-seq", "scRNA-seq", "snRNA-seq"),
  stringsAsFactors = FALSE
)

metadata <- data.frame(
  `File Path` = character(),
  `Sample Name` = character(),
  SRA = character(),
  Tissue = character(),
  Approach = character(),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(samples))) {
  m <- rsparsematrix(nrow = 100, ncol = 50, density = 0.05)
  sample <- samples$sample[i]
  file_rel <- file.path("inst", "extdata", "toy", paste0(sample, ".sparse.RData"))
  file_abs <- file.path(pkg_root, file_rel)
  save(m, file = file_abs)
  metadata[i, ] <- list(file_rel, sample, samples$sra[i], samples$tissue[i], samples$approach[i])
}

metadata_path <- file.path(pkg_root, "inst", "extdata", "toy", "metadata.csv")
write.csv(metadata, metadata_path, row.names = FALSE)

message("Toy data and metadata generated in ", out_dir)
