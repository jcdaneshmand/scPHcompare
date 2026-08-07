#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("usage: assemble_mv04_distance_bundle.R PAIRS MANIFEST OUTPUT_DIR AUDIT_DIR")
}

pairs_path <- args[[1L]]
manifest_path <- args[[2L]]
output_dir <- args[[3L]]
audit_dir <- args[[4L]]

source("R/topological_distance_contract.R")
pairs <- read.csv(pairs_path, check.names = FALSE, stringsAsFactors = FALSE)
manifest <- read.csv(manifest_path, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(nrow(pairs) == 408L, nrow(manifest) == 56L,
          all(pairs$status == "completed"),
          all(tolower(pairs$exact) == "true"),
          all(tolower(pairs$all_active_levels) == "true"),
          all(pairs$absolute_error_estimate == 0))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
bundles <- list()
contributions <- list()
matrix_checks <- list()

for (stratum in sort(unique(manifest$stratum_id))) {
  manifest_part <- manifest[manifest$stratum_id == stratum, , drop = FALSE]
  pair_part <- pairs[pairs$stratum_id == stratum, , drop = FALSE]
  sample_ids <- sort(manifest_part$sample_id, method = "radix")
  h0 <- distance_pairs_to_matrix_v1(
    pair_part[pair_part$homology_dimension == "H0", , drop = FALSE], sample_ids
  )
  h1 <- distance_pairs_to_matrix_v1(
    pair_part[pair_part$homology_dimension == "H1", , drop = FALSE], sample_ids
  )
  diagram_ids <- manifest_part$diagram_id[match(sample_ids, manifest_part$sample_id)]
  bundle <- new_mv04_distance_bundle(
    h0, h1, stratum, "full_l2_exact_critical_pairs_v1", diagram_ids
  )
  bundles[[stratum]] <- bundle
  contributions[[stratum]] <- mv04_distance_contributions(bundle)
  matrix_checks[[stratum]] <- data.frame(
    stratum_id = stratum,
    sample_count = length(sample_ids),
    pair_count = nrow(pair_part) / 2L,
    h0_symmetric = identical(h0, t(h0)),
    h1_symmetric = identical(h1, t(h1)),
    h0_zero_diagonal = all(diag(h0) == 0),
    h1_zero_diagonal = all(diag(h1) == 0),
    h0_scale = bundle$scales$H0$scale,
    h1_scale = bundle$scales$H1$scale,
    cache_key = bundle$cache_key,
    stringsAsFactors = FALSE
  )
  saveRDS(bundle, file.path(output_dir, paste0(stratum, ".rds")), version = 3)
}

write.csv(pairs, file.path(audit_dir, "mv04-landscape-pairs-2026-08-05.csv"),
          row.names = FALSE)
write.csv(do.call(rbind, contributions),
          file.path(audit_dir, "mv04-h0-h1-contributions-2026-08-05.csv"),
          row.names = FALSE)
write.csv(do.call(rbind, matrix_checks),
          file.path(audit_dir, "mv04-matrix-validation-2026-08-05.csv"),
          row.names = FALSE)
saveRDS(bundles, file.path(output_dir, "mv04-eligible-distance-bundle.rds"),
        version = 3)
message("Assembled ", length(bundles), " immutable MV-04 distance strata.")
