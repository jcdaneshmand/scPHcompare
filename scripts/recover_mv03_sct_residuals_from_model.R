#!/usr/bin/env Rscript

options(warn = 2)

for (package in c("Seurat", "SeuratObject", "digest")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(package, " is required for SCT residual recovery.", call. = FALSE)
  }
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) {
  stop(
    "Usage: recover_mv03_sct_residuals_from_model.R <cohort> ",
    "<seurat-rds> <integrated-cache-rds> <manifest-csv> <cache-rds> ",
    "<object-audit-csv> <sample-audit-csv>",
    call. = FALSE
  )
}

cohort <- args[[1L]]
object_path <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
integrated_path <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
manifest_path <- normalizePath(args[[4L]], winslash = "/", mustWork = TRUE)
cache_path <- args[[5L]]
object_audit_path <- args[[6L]]
sample_audit_path <- args[[7L]]

dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(object_audit_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(sample_audit_path), recursive = TRUE, showWarnings = FALSE)

manifest <- utils::read.csv(
  manifest_path, stringsAsFactors = FALSE, check.names = FALSE
)
manifest <- manifest[manifest$cohort == cohort, , drop = FALSE]
integrated <- readRDS(integrated_path)
if (!identical(integrated$cohort, cohort) || !nrow(manifest)) {
  stop("Integrated cache or manifest has an incompatible cohort.",
       call. = FALSE)
}
sample_ids <- sort(names(integrated$samples), method = "radix")
desired_features <- Reduce(
  intersect, lapply(integrated$samples, rownames)
)
desired_features <- sort(desired_features, method = "radix")
integrated_sha256 <- digest::digest(
  file = integrated_path, algo = "sha256", serialize = FALSE
)
rm(integrated)
invisible(gc())

read_started <- proc.time()[["elapsed"]]
object <- readRDS(object_path)
read_seconds <- proc.time()[["elapsed"]] - read_started
if (!inherits(object, "Seurat") || !"SCT" %in% names(object@assays) ||
    !"RNA" %in% names(object@assays) ||
    !"orig.ident" %in% colnames(object@meta.data)) {
  stop("Historical object lacks the Seurat/SCT/RNA/orig.ident contract.",
       call. = FALSE)
}
pilot_cells <- rownames(object@meta.data)[
  as.character(object@meta.data$orig.ident) %in% sample_ids
]
if (!length(pilot_cells)) {
  stop("No pilot cells were found in the historical object.", call. = FALSE)
}
object <- subset(object, cells = pilot_cells)
invisible(gc())

available_features <- rownames(object[["RNA"]])
missing_rna <- setdiff(desired_features, available_features)
if (length(missing_rna)) {
  stop("Integrated features are absent from the RNA assay: ",
       paste(utils::head(missing_rna, 10L), collapse = ", "), call. = FALSE)
}

recovery_started <- proc.time()[["elapsed"]]
object <- Seurat::GetResidual(
  object = object,
  features = desired_features,
  assay = "SCT",
  umi.assay = "RNA",
  replace.value = FALSE,
  na.rm = TRUE,
  verbose = TRUE
)
recovery_seconds <- proc.time()[["elapsed"]] - recovery_started
residuals <- methods::slot(object[["SCT"]], "scale.data")
if (!all(desired_features %in% rownames(residuals))) {
  stop("GetResidual did not materialize every requested feature.",
       call. = FALSE)
}
residuals <- residuals[desired_features, , drop = FALSE]

metadata <- object@meta.data
cell_index <- match(colnames(residuals), rownames(metadata))
if (anyNA(cell_index)) {
  stop("Recovered residual cells do not map to metadata.", call. = FALSE)
}
residual_sample <- as.character(metadata$orig.ident[cell_index])
pilot <- vector("list", length(sample_ids))
names(pilot) <- sample_ids
sample_rows <- vector("list", length(sample_ids))
for (index in seq_along(sample_ids)) {
  sample_id <- sample_ids[[index]]
  columns <- which(residual_sample == sample_id)
  value <- residuals[, columns, drop = FALSE]
  canonical_ids <- sub("^.*__", "", colnames(value))
  if (anyDuplicated(canonical_ids)) {
    stop("Recovered canonical cell IDs are duplicated for ", sample_id,
         call. = FALSE)
  }
  colnames(value) <- canonical_ids
  finite <- all(is.finite(as.numeric(value)))
  manifest_count <- manifest$filtered_cells[manifest$sample_id == sample_id]
  eligible <- length(manifest_count) == 1L && ncol(value) == manifest_count &&
    finite
  pilot[[index]] <- value
  sample_rows[[index]] <- data.frame(
    cohort = cohort,
    sample_id = sample_id,
    manifest_filtered_cells = manifest_count,
    recovered_cells = ncol(value),
    recovered_genes = nrow(value),
    exact_manifest_cell_match = length(manifest_count) == 1L &&
      ncol(value) == manifest_count,
    unique_canonical_cell_ids = !anyDuplicated(canonical_ids),
    finite_residual_values = finite,
    extraction_eligible = eligible,
    stringsAsFactors = FALSE
  )
}
sample_audit <- do.call(rbind, sample_rows)
if (!all(sample_audit$extraction_eligible)) {
  utils::write.csv(sample_audit, sample_audit_path, row.names = FALSE)
  stop("At least one model-recovered residual sample is ineligible.",
       call. = FALSE)
}

cache <- list(
  contract_id = "mv03_sct_whole_pearson_residual_cache_v1",
  cohort = cohort,
  extraction = "SCT_model_recovered_Pearson_residuals_for_integrated_features",
  source_path = object_path,
  source_size_bytes = unname(file.info(object_path)$size),
  paired_integrated_cache_sha256 = integrated_sha256,
  genes = desired_features,
  samples = pilot
)
temporary_cache <- paste0(cache_path, ".partial")
saveRDS(cache, temporary_cache, compress = FALSE)
if (file.exists(cache_path)) {
  stop("Refusing to overwrite an existing recovered residual cache.",
       call. = FALSE)
}
if (!file.rename(temporary_cache, cache_path)) {
  stop("Failed to atomically publish recovered residual cache.",
       call. = FALSE)
}
cache_sha256 <- digest::digest(file = cache_path, algo = "sha256",
                               serialize = FALSE)

sct <- object[["SCT"]]
models <- methods::slot(sct, "SCTModel.list")
object_audit <- data.frame(
  cohort = cohort,
  source_file = basename(object_path),
  source_size_bytes = unname(file.info(object_path)$size),
  read_seconds = read_seconds,
  recovery_seconds = recovery_seconds,
  recovery_api = "Seurat::GetResidual",
  seurat_version = as.character(utils::packageVersion("Seurat")),
  sctransform_version = as.character(utils::packageVersion("sctransform")),
  sct_assay_class = paste(class(sct), collapse = "/"),
  sct_model_count = length(models),
  subset_pilot_cells = ncol(object),
  requested_features = length(desired_features),
  recovered_features = nrow(residuals),
  extraction_semantics =
    "SCT_model_recovered_Pearson_residuals_for_integrated_features",
  pilot_samples = length(pilot),
  eligible_pilot_samples = sum(sample_audit$extraction_eligible),
  cache_file = basename(cache_path),
  cache_size_bytes = unname(file.info(cache_path)$size),
  cache_sha256 = cache_sha256,
  stringsAsFactors = FALSE
)
utils::write.csv(object_audit, object_audit_path, row.names = FALSE)
utils::write.csv(sample_audit, sample_audit_path, row.names = FALSE)
message("Recovered eligible SCT Pearson residuals for ", length(pilot),
        " ", cohort, " pilot samples and ", length(desired_features),
        " integrated features.")
