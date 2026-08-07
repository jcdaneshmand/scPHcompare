#!/usr/bin/env Rscript

options(warn = 2)

if (!requireNamespace("Matrix", quietly = TRUE)) {
  stop("Matrix is required to inspect historical integrated matrices.",
       call. = FALSE)
}
if (!requireNamespace("digest", quietly = TRUE)) {
  stop("digest is required to bind extracted caches to their payloads.",
       call. = FALSE)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) {
  stop(
    "Usage: extract_mv03_integrated_pilot.R <cohort> <representation> ",
    "<integrated-list-rds> <residual-cache-rds> <cache-rds> ",
    "<object-audit-csv> <sample-audit-csv>",
    call. = FALSE
  )
}

cohort <- args[[1L]]
representation <- args[[2L]]
source_path <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
residual_cache_path <- normalizePath(args[[4L]], winslash = "/", mustWork = TRUE)
cache_path <- args[[5L]]
object_audit_path <- args[[6L]]
sample_audit_path <- args[[7L]]

dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(object_audit_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(sample_audit_path), recursive = TRUE, showWarnings = FALSE)

residual_cache <- readRDS(residual_cache_path)
if (!identical(residual_cache$contract_id,
               "mv03_sct_whole_pearson_residual_cache_v1") ||
    !identical(residual_cache$cohort, cohort)) {
  stop("Residual cache has an incompatible contract or cohort.", call. = FALSE)
}

started <- proc.time()[["elapsed"]]
source <- readRDS(source_path)
read_seconds <- proc.time()[["elapsed"]] - started
if (!is.list(source) || is.null(names(source)) || anyDuplicated(names(source))) {
  stop("Integrated source must be a uniquely named list.", call. = FALSE)
}

sample_ids <- names(residual_cache$samples)
if (!all(sample_ids %in% names(source))) {
  stop("Integrated source is missing one or more frozen pilot samples.",
       call. = FALSE)
}

matrix_values <- function(x) {
  if (inherits(x, "sparseMatrix")) {
    return(methods::slot(x, "x"))
  }
  as.numeric(x)
}

pilot <- vector("list", length(sample_ids))
names(pilot) <- sample_ids
sample_rows <- vector("list", length(sample_ids))
for (index in seq_along(sample_ids)) {
  sample_id <- sample_ids[[index]]
  value <- source[[sample_id]]
  if ((!is.matrix(value) && !inherits(value, "Matrix")) ||
      is.null(rownames(value)) || is.null(colnames(value)) ||
      anyDuplicated(rownames(value)) || anyDuplicated(colnames(value))) {
    stop("Integrated sample is not a uniquely named matrix: ", sample_id,
         call. = FALSE)
  }
  canonical_ids <- sub("^.*__", "", colnames(value))
  residual_ids <- colnames(residual_cache$samples[[sample_id]])
  unique_canonical <- !anyDuplicated(canonical_ids)
  paired_set <- unique_canonical && setequal(canonical_ids, residual_ids)
  if (!paired_set) {
    stop("Integrated/residual cell pairing is not bijective for ", sample_id,
         call. = FALSE)
  }
  value <- value[, match(residual_ids, canonical_ids), drop = FALSE]
  colnames(value) <- residual_ids
  finite <- all(is.finite(matrix_values(value)))
  if (!finite) {
    stop("Integrated sample contains nonfinite values: ", sample_id,
         call. = FALSE)
  }
  pilot[[index]] <- value
  sample_rows[[index]] <- data.frame(
    cohort = cohort,
    representation = representation,
    sample_id = sample_id,
    genes = nrow(value),
    cells = ncol(value),
    residual_cells = length(residual_ids),
    unique_canonical_cell_ids = unique_canonical,
    paired_cell_set_exact = paired_set,
    paired_cell_order_after_reindex = identical(colnames(value), residual_ids),
    finite_values = finite,
    extraction_eligible = paired_set && finite &&
      ncol(value) == ncol(residual_cache$samples[[sample_id]]),
    stringsAsFactors = FALSE
  )
}
sample_audit <- do.call(rbind, sample_rows)
if (!all(sample_audit$extraction_eligible)) {
  utils::write.csv(sample_audit, sample_audit_path, row.names = FALSE)
  stop("At least one integrated pilot sample is ineligible.", call. = FALSE)
}

cache <- list(
  contract_id = "mv03_integrated_corrected_values_cache_v1",
  cohort = cohort,
  representation = representation,
  extraction = "integrated_assay_data_corrected_values",
  source_path = source_path,
  source_size_bytes = unname(file.info(source_path)$size),
  paired_residual_cache_sha256 = digest::digest(
    file = residual_cache_path, algo = "sha256", serialize = FALSE
  ),
  samples = pilot
)
temporary_cache <- paste0(cache_path, ".partial")
saveRDS(cache, temporary_cache, compress = FALSE)
if (file.exists(cache_path)) {
  stop("Refusing to overwrite an existing MV-03 integrated cache.",
       call. = FALSE)
}
if (!file.rename(temporary_cache, cache_path)) {
  stop("Failed to atomically publish the MV-03 integrated cache.",
       call. = FALSE)
}
cache_sha256 <- digest::digest(file = cache_path, algo = "sha256",
                               serialize = FALSE)

object_audit <- data.frame(
  cohort = cohort,
  representation = representation,
  source_file = basename(source_path),
  source_size_bytes = unname(file.info(source_path)$size),
  source_modified_utc = format(
    file.info(source_path)$mtime, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"
  ),
  read_seconds = read_seconds,
  list_samples = length(source),
  pilot_samples = length(pilot),
  eligible_pilot_samples = sum(sample_audit$extraction_eligible),
  feature_count_min = min(sample_audit$genes),
  feature_count_max = max(sample_audit$genes),
  cache_file = basename(cache_path),
  cache_size_bytes = unname(file.info(cache_path)$size),
  cache_sha256 = cache_sha256,
  stringsAsFactors = FALSE
)

utils::write.csv(object_audit, object_audit_path, row.names = FALSE)
utils::write.csv(sample_audit, sample_audit_path, row.names = FALSE)
message("Extracted eligible paired integrated values for ", length(pilot),
        " ", cohort, " pilot samples.")
