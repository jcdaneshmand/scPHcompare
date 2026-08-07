#!/usr/bin/env Rscript

options(warn = 2)

if (!requireNamespace("SeuratObject", quietly = TRUE)) {
  stop("SeuratObject is required to inspect the historical Seurat object.",
       call. = FALSE)
}
if (!requireNamespace("digest", quietly = TRUE)) {
  stop("digest is required to bind extracted caches to their payloads.",
       call. = FALSE)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop(
    "Usage: extract_mv03_sct_residuals.R <cohort> <seurat-rds> ",
    "<manifest-csv> <cache-rds> <object-audit-csv> <sample-audit-csv>",
    call. = FALSE
  )
}

cohort <- args[[1L]]
object_path <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
manifest_path <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
cache_path <- args[[4L]]
object_audit_path <- args[[5L]]
sample_audit_path <- args[[6L]]

manifest <- utils::read.csv(
  manifest_path, stringsAsFactors = FALSE, check.names = FALSE
)
required_manifest <- c("cohort", "sample_id", "filtered_cells")
if (!all(required_manifest %in% names(manifest))) {
  stop("Pilot manifest has an incompatible schema.", call. = FALSE)
}
manifest <- manifest[manifest$cohort == cohort, , drop = FALSE]
if (!nrow(manifest) || anyDuplicated(manifest$sample_id)) {
  stop("Cohort manifest is empty or has duplicate sample IDs.", call. = FALSE)
}

dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(object_audit_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(sample_audit_path), recursive = TRUE, showWarnings = FALSE)

started <- proc.time()[["elapsed"]]
object <- readRDS(object_path)
read_seconds <- proc.time()[["elapsed"]] - started

if (!inherits(object, "Seurat")) {
  stop("Historical source is not a Seurat object.", call. = FALSE)
}
if (!"SCT" %in% names(object@assays)) {
  stop("Historical Seurat object has no SCT assay.", call. = FALSE)
}
if (!"orig.ident" %in% colnames(object@meta.data)) {
  stop("Historical Seurat object has no orig.ident metadata.", call. = FALSE)
}

sct <- object[["SCT"]]
sct_slots <- methods::slotNames(sct)
if (!"scale.data" %in% sct_slots) {
  stop("SCT assay has no scale.data slot; Pearson residuals unavailable.",
       call. = FALSE)
}
residuals <- methods::slot(sct, "scale.data")
if (!is.matrix(residuals) && !inherits(residuals, "Matrix")) {
  stop("SCT scale.data is not matrix-like.", call. = FALSE)
}
if (!nrow(residuals) || !ncol(residuals) || is.null(rownames(residuals)) ||
    is.null(colnames(residuals)) || anyDuplicated(rownames(residuals)) ||
    anyDuplicated(colnames(residuals))) {
  stop("SCT scale.data lacks complete unique axes.", call. = FALSE)
}

metadata <- object@meta.data
if (!identical(rownames(metadata), colnames(object))) {
  stop("Seurat metadata rows do not exactly match object cells.", call. = FALSE)
}
cell_index <- match(colnames(residuals), rownames(metadata))
if (anyNA(cell_index)) {
  stop("SCT residual cells are not all represented in object metadata.",
       call. = FALSE)
}
residual_sample <- as.character(metadata$orig.ident[cell_index])

model_count <- 0L
model_classes <- ""
if ("SCTModel.list" %in% sct_slots) {
  models <- methods::slot(sct, "SCTModel.list")
  model_count <- length(models)
  model_classes <- paste(sort(unique(vapply(
    models, function(x) paste(class(x), collapse = "/"), character(1)
  ))), collapse = ";")
}

sample_rows <- vector("list", nrow(manifest))
pilot <- vector("list", nrow(manifest))
names(pilot) <- manifest$sample_id
for (index in seq_len(nrow(manifest))) {
  sample_id <- manifest$sample_id[[index]]
  columns <- which(residual_sample == sample_id)
  extracted <- length(columns) > 0L
  if (extracted) {
    value <- residuals[, columns, drop = FALSE]
    cell_ids <- colnames(value)
    canonical_ids <- sub("^.*__", "", cell_ids)
    if (anyDuplicated(canonical_ids)) {
      stop("Canonical pilot cell IDs are duplicated for ", sample_id,
           call. = FALSE)
    }
    colnames(value) <- canonical_ids
    pilot[[index]] <- value
    finite <- all(is.finite(as.numeric(value)))
  } else {
    pilot[[index]] <- NULL
    finite <- FALSE
  }
  sample_rows[[index]] <- data.frame(
    cohort = cohort,
    sample_id = sample_id,
    manifest_filtered_cells = manifest$filtered_cells[[index]],
    residual_cells = length(columns),
    residual_genes = if (extracted) nrow(pilot[[index]]) else 0L,
    exact_manifest_cell_match = length(columns) == manifest$filtered_cells[[index]],
    unique_canonical_cell_ids = extracted &&
      !anyDuplicated(colnames(pilot[[index]])),
    finite_residual_values = finite,
    extraction_eligible = extracted && finite &&
      length(columns) == manifest$filtered_cells[[index]],
    stringsAsFactors = FALSE
  )
}
sample_audit <- do.call(rbind, sample_rows)
if (!all(sample_audit$extraction_eligible)) {
  utils::write.csv(sample_audit, sample_audit_path, row.names = FALSE)
  stop("At least one pilot sample lacks eligible SCT Pearson residuals.",
       call. = FALSE)
}

cache <- list(
  contract_id = "mv03_sct_whole_pearson_residual_cache_v1",
  cohort = cohort,
  extraction = "SCT_assay_scale.data_Pearson_residuals",
  source_path = object_path,
  source_size_bytes = unname(file.info(object_path)$size),
  genes = rownames(residuals),
  samples = pilot
)
temporary_cache <- paste0(cache_path, ".partial")
saveRDS(cache, temporary_cache, compress = FALSE)
if (file.exists(cache_path)) {
  stop("Refusing to overwrite an existing MV-03 residual cache.", call. = FALSE)
}
if (!file.rename(temporary_cache, cache_path)) {
  stop("Failed to atomically publish the MV-03 residual cache.", call. = FALSE)
}
cache_sha256 <- digest::digest(file = cache_path, algo = "sha256",
                               serialize = FALSE)

object_audit <- data.frame(
  cohort = cohort,
  source_file = basename(object_path),
  source_size_bytes = unname(file.info(object_path)$size),
  source_modified_utc = format(
    file.info(object_path)$mtime, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"
  ),
  read_seconds = read_seconds,
  seurat_class = paste(class(object), collapse = "/"),
  assays = paste(names(object@assays), collapse = ";"),
  sct_assay_class = paste(class(sct), collapse = "/"),
  sct_assay_slots = paste(sct_slots, collapse = ";"),
  sct_model_count = model_count,
  sct_model_classes = model_classes,
  object_cells = ncol(object),
  metadata_samples = length(unique(as.character(metadata$orig.ident))),
  residual_genes = nrow(residuals),
  residual_cells = ncol(residuals),
  residual_cells_match_object = identical(colnames(residuals), colnames(object)),
  extraction_semantics = "SCT_assay_scale.data_Pearson_residuals",
  pilot_samples = nrow(sample_audit),
  eligible_pilot_samples = sum(sample_audit$extraction_eligible),
  cache_file = basename(cache_path),
  cache_size_bytes = unname(file.info(cache_path)$size),
  cache_sha256 = cache_sha256,
  stringsAsFactors = FALSE
)

utils::write.csv(object_audit, object_audit_path, row.names = FALSE)
utils::write.csv(sample_audit, sample_audit_path, row.names = FALSE)
message("Extracted eligible SCT Pearson residuals for ", nrow(sample_audit),
        " ", cohort, " pilot samples.")
