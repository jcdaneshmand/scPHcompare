#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(
    "usage: validate_mv05d0_v2_reproducibility.R CACHE_A CACHE_B OUTPUT_CSV",
    call. = FALSE
  )
}
path_a <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
path_b <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
output_path <- args[[3L]]
source("R/provenance_utils.R")
source("R/mv05_resource_safe_execution.R")
a <- readRDS(path_a)
b <- readRDS(path_b)
mv05d0_validate_normalization_cache_record_v2(a)
mv05d0_validate_normalization_cache_record_v2(b)
matrix_a <- mv05d0_sct_matrix_from_cache_v1(a)
matrix_b <- mv05d0_sct_matrix_from_cache_v1(b)
result <- data.frame(
  contract_id = "mv05d0_v2_reproducibility_v1",
  sample_id = a$identity$sample_id, seed = a$identity$seed,
  normalization_cache_key = a$cache_key,
  cache_key_identical = identical(a$cache_key, b$cache_key),
  payload_sha256 = a$payload_sha256,
  payload_hash_identical = identical(a$payload_sha256, b$payload_sha256),
  file_a_sha256 = digest::digest(
    file = path_a, algo = "sha256", serialize = FALSE
  ),
  file_b_sha256 = digest::digest(
    file = path_b, algo = "sha256", serialize = FALSE
  ),
  file_hash_identical = identical(
    digest::digest(file = path_a, algo = "sha256", serialize = FALSE),
    digest::digest(file = path_b, algo = "sha256", serialize = FALSE)
  ),
  file_a_bytes = unname(file.info(path_a)$size),
  file_b_bytes = unname(file.info(path_b)$size),
  matrix_exact_numeric_identity = identical(matrix_a, matrix_b),
  r_version = a$identity$runtime$r_version,
  seurat_version = a$identity$runtime$seurat_version,
  sctransform_version = a$identity$runtime$sctransform_version,
  matrix_version = a$identity$runtime$matrix_version,
  outcome_label_state = "closed",
  biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
if (!all(result[c("cache_key_identical", "payload_hash_identical",
                  "file_hash_identical", "matrix_exact_numeric_identity")])) {
  stop("Independent MV5-D0 v2 caches are not reproducible.", call. = FALSE)
}
write_provenance_csv(result, output_path)
message("Validated byte-exact MV5-D0 v2 cache reproducibility.")
