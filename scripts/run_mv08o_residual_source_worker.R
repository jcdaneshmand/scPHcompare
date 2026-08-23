#!/usr/bin/env Rscript

# One MV8-O source/view worker.  It fits SCT on every frozen-QC cell in one
# authorized source, materializes the exact500 data/residual layers, then
# evaluates only the immutable selected-384 axes.  It never calls PH.

args <- commandArgs(trailingOnly = TRUE)
if (!(length(args) %in% 7:8)) {
  stop(paste(
    "usage: run_mv08o_residual_source_worker.R <mv08n-audit-dir>",
    "<raw|h5> <source-path> <unit-id> <output-rds> <output-csv> <run-role> [source-queue.csv]"
  ), call. = FALSE)
}

audit_dir <- normalizePath(args[[1L]], mustWork = TRUE)
source_kind <- args[[2L]]
source_path <- normalizePath(args[[3L]], mustWork = TRUE)
unit_id <- args[[4L]]
output_rds <- normalizePath(args[[5L]], mustWork = FALSE)
output_csv <- normalizePath(args[[6L]], mustWork = FALSE)
run_role <- args[[7L]]
source_queue_path <- if (length(args) == 8L) normalizePath(args[[8L]], mustWork = TRUE) else
  file.path(audit_dir, "mv08n-residual-source-queue.csv")
if (!(source_kind %in% c("raw", "h5")) || !nzchar(unit_id) ||
    !(run_role %in% c("primary", "repeat")) || file.exists(output_rds) ||
    file.exists(output_csv)) {
  stop("invalid MV8-O worker arguments or refusing to overwrite output", call. = FALSE)
}
dir.create(dirname(output_rds), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(output_csv), recursive = TRUE, showWarnings = FALSE)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
options(warn = 1)
future_globals_max_size_env <- Sys.getenv("MV08_FUTURE_GLOBALS_MAX_SIZE_BYTES", unset = "")
future_globals_max_size_bytes <- if (nzchar(future_globals_max_size_env))
  suppressWarnings(as.numeric(future_globals_max_size_env)) else NA_real_
if (nzchar(future_globals_max_size_env)) {
  if (!is.finite(future_globals_max_size_bytes) ||
      future_globals_max_size_bytes != 2 * 1024^3) {
    stop("MV08 future globals limit must be the authorized bounded 2 GiB value", call. = FALSE)
  }
  options(future.globals.maxSize = future_globals_max_size_bytes)
}
for (package in c("digest", "Matrix", "Seurat")) {
  if (!requireNamespace(package, quietly = TRUE)) stop(package, " is required", call. = FALSE)
}
source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/mv05_resource_safe_execution.R")

sha_file <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
sha_object <- function(x) digest::digest(x, algo = "sha256", serialize = TRUE)
atomic_csv <- function(x, path) {
  partial <- paste0(path, ".partial")
  utils::write.csv(x, partial, row.names = FALSE, quote = TRUE, na = "")
  if (!file.rename(partial, path)) stop("failed to publish ", basename(path), call. = FALSE)
}
atomic_rds <- function(x, path) {
  partial <- tempfile(pattern = paste0(basename(path), "."), tmpdir = dirname(path))
  on.exit(if (file.exists(partial)) unlink(partial), add = TRUE)
  saveRDS(x, partial, compress = FALSE, version = 3)
  if (!file.rename(partial, path)) stop("failed to publish ", basename(path), call. = FALSE)
}
read_csv <- function(name) utils::read.csv(file.path(audit_dir, name), check.names = FALSE,
                                            stringsAsFactors = FALSE)
geometry <- function(values) {
  sds <- apply(values, 1L, stats::sd)
  zero <- sum(!is.finite(sds) | sds <= sqrt(.Machine$double.eps))
  out <- list(finite = all(is.finite(values)), zero_variance = zero,
              minimum_sd = suppressWarnings(min(sds, na.rm = TRUE)),
              valid = FALSE, distance_sha256 = NA_character_)
  if (!out$finite || zero != 0L) return(out)
  corr <- stats::cor(t(values), method = "pearson")
  corr[corr < -1] <- -1; corr[corr > 1] <- 1
  distance <- sqrt(2 * (1 - corr)); distance[distance < 0] <- 0; diag(distance) <- 0
  out$valid <- all(is.finite(distance)) && isSymmetric(distance) && all(diag(distance) == 0) &&
    all(distance[upper.tri(distance)] >= 0)
  if (out$valid) out$distance_sha256 <- sha_object(distance)
  out
}

queue <- utils::read.csv(source_queue_path, check.names = FALSE, stringsAsFactors = FALSE)
row <- queue[queue$unit_id == unit_id, , drop = FALSE]
if (nrow(row) != 1L || !(row$authorization_state %in% c("source_view_sentinel_authorized", "authorized_after_mv08p_commit")) ||
    row$fit_cells < 384L || row$selected_cells_per_axis != 384L || row$workers != 1L ||
    row$retries != 0L || row$elapsed_cap_seconds != 1800L ||
    row$rss_cap_bytes != 12 * 1024^3 || row$outcome_label_state != "closed" ||
    isTRUE(row$biological_outcomes_computed)) {
  stop("unit is not an authorized MV8-O source sentinel", call. = FALSE)
}
if (!identical(tolower(sha_file(source_path)), tolower(row$source_sha256[[1L]]))) {
  stop("source hash drift", call. = FALSE)
}
exact <- utils::read.csv("docs/audits/mv07h-prefreeze-evidence-v4/mv07h-panel.csv",
                         check.names = FALSE, stringsAsFactors = FALSE)
exact <- exact[order(exact$panel_order), , drop = FALSE]
if (nrow(exact) != 500L || !identical(as.integer(exact$panel_order), seq_len(500L)) ||
    anyDuplicated(exact$feature_id)) stop("exact500 axis drift", call. = FALSE)
exact_sha <- unique(tolower(exact$panel_sha256))
if (length(exact_sha) != 1L) stop("exact500 SHA drift", call. = FALSE)
common <- utils::read.csv("docs/audits/mv08e-reference-reconciliation-evidence/mv08e-common475-panel.csv",
                          check.names = FALSE, stringsAsFactors = FALSE)
common <- common[order(common$panel_order_475), , drop = FALSE]
common_index <- match(common$feature_id, exact$feature_id)
if (nrow(common) != 475L || !identical(as.integer(common$panel_order_475), seq_len(475L)) ||
    anyNA(common_index) || anyDuplicated(common_index)) {
  stop("common475 axis drift", call. = FALSE)
}

if (source_kind == "raw") {
  raw <- readRDS(source_path)
  mv05d0_validate_raw_sample_cache_v2(raw)
  if (!identical(raw$sample_id, unit_id)) stop("raw sample identity mismatch", call. = FALSE)
  counts <- if (inherits(raw$counts, "dgCMatrix")) raw$counts else
    methods::as(raw$counts, "dgCMatrix")
  selection <- utils::read.csv("docs/audits/mv07fp-prefreeze-evidence-v4/mv07fp-cache-manifest.csv",
                               check.names = FALSE, stringsAsFactors = FALSE)
  selected_rows <- selection[selection$sample_id == unit_id, , drop = FALSE]
  selected_rows <- selected_rows[order(selected_rows$seed), , drop = FALSE]
  if (nrow(selected_rows) != 5L || !identical(as.integer(selected_rows$seed), 20260805:20260809) ||
      any(selected_rows$selected_cells != 384L) ||
      any(selected_rows$outcome_label_state != "closed") ||
      any(selected_rows$biological_outcomes_computed)) {
    stop("internal selected-cell axis drift", call. = FALSE)
  }
  selected_by_seed <- lapply(selected_rows$seed, function(seed) {
    select_matched_cells(colnames(counts), n = 384L, seed = seed)
  })
  names(selected_by_seed) <- as.character(selected_rows$seed)
  if (!identical(unname(vapply(selected_by_seed, attr, character(1L), which = "selected_cell_sha256")),
                 as.character(selected_rows$selected_cell_sha256))) {
    stop("internal selected-cell hash drift", call. = FALSE)
  }
  if (ncol(counts) != row$fit_cells) stop("internal frozen-QC count drift", call. = FALSE)
} else {
  raw_ids <- Seurat::Read10X_h5(source_path, use.names = FALSE, unique.features = TRUE)
  raw_names <- Seurat::Read10X_h5(source_path, use.names = TRUE, unique.features = TRUE)
  if (is.list(raw_ids) || is.list(raw_names) || !identical(dim(raw_ids), dim(raw_names))) {
    stop("HCA input must be one Gene Expression matrix", call. = FALSE)
  }
  counts <- as(raw_ids, "dgCMatrix")
  total <- Matrix::colSums(counts)
  nfeature <- Matrix::colSums(counts > 0)
  mito <- grepl("^MT-", rownames(raw_names))
  ribo <- grepl("^(RPS|RPL)", rownames(raw_names))
  p_mito <- 100 * Matrix::colSums(counts[mito, , drop = FALSE]) / total
  p_ribo <- 100 * Matrix::colSums(counts[ribo, , drop = FALSE]) / total
  eligible <- nfeature >= 500 & nfeature <= 9000 & p_mito <= 20 & p_ribo > 5 & total > 0
  if (sum(eligible) != row$fit_cells) stop("HCA frozen-QC count drift", call. = FALSE)
  counts <- counts[, eligible, drop = FALSE]
  selected <- select_matched_cells(colnames(counts), n = 384L, seed = 20260805L)
  external_axis <- utils::read.csv(
    "docs/audits/mv08n-pearson-residual-migration-prefreeze-v1/mv08n-external-source-axis.csv",
    check.names = FALSE, stringsAsFactors = FALSE)
  expected <- external_axis[external_axis$unit_id == unit_id, , drop = FALSE]
  if (nrow(expected) != 1L || expected$selection_seed != 20260805L ||
      expected$selected_cells != 384L || expected$qc_eligible_cells != ncol(counts) ||
      expected$outcome_label_state != "closed" || isTRUE(expected$biological_outcomes_computed)) {
    stop("HCA selected-cell preflight drift", call. = FALSE)
  }
  observed_selection_sha <- attr(selected, "selected_cell_sha256")
  if (expected$selected_axis_state == "frozen_mv08k") {
    if (!nzchar(expected$selected_cell_sha256) ||
        !identical(observed_selection_sha, expected$selected_cell_sha256[[1L]])) {
      stop("HCA frozen selected-cell axis drift", call. = FALSE)
    }
  } else if (expected$selected_axis_state != "freeze_in_source_preflight" ||
             nzchar(expected$selected_cell_sha256)) {
    stop("HCA selected-cell preflight policy drift", call. = FALSE)
  }
  selected_by_seed <- list("20260805" = selected)
}

panel_index <- match(exact$feature_id, rownames(counts))
if (anyNA(panel_index) || anyDuplicated(panel_index)) {
  stable_panel <- sub("\\.[0-9]+$", "", sub("^.*-", "", exact$feature_id))
  stable_counts <- sub("\\.[0-9]+$", "", rownames(counts))
  panel_index <- match(stable_panel, stable_counts)
}
if (anyNA(panel_index) || anyDuplicated(panel_index)) stop("exact500 feature mapping drift", call. = FALSE)
panel_source_ids <- rownames(counts)[panel_index]
if (min(Matrix::rowSums(counts[panel_index, , drop = FALSE] > 0)) <= 0L) {
  stop("exact500 raw feature absence after frozen QC", call. = FALSE)
}

started <- proc.time()[["elapsed"]]
object <- Seurat::CreateSeuratObject(counts = counts, min.cells = 0L, min.features = 0L,
                                     project = unit_id)
sct_warnings <- character()
object <- withCallingHandlers(
  Seurat::SCTransform(object, assay = "RNA", new.assay.name = "SCT",
    variable.features.n = 3000L, return.only.var.genes = FALSE, min_cells = 5L,
    verbose = FALSE),
  warning = function(warning) {
    sct_warnings <<- c(sct_warnings, conditionMessage(warning))
    invokeRestart("muffleWarning")
  }
)
sct_seconds <- proc.time()[["elapsed"]] - started
sct_data <- Seurat::GetAssayData(object, assay = "SCT", layer = "data")
data_index <- match(panel_source_ids, rownames(sct_data))
if (anyNA(data_index) || anyDuplicated(data_index)) stop("SCT data lost exact500 feature", call. = FALSE)
residual_started <- proc.time()[["elapsed"]]
object <- Seurat::GetResidual(object, features = panel_source_ids, assay = "SCT", umi.assay = "RNA",
                              replace.value = FALSE, na.rm = TRUE, verbose = FALSE)
residual_seconds <- proc.time()[["elapsed"]] - residual_started
scale_data <- methods::slot(object[["SCT"]], "scale.data")
residual_index <- match(panel_source_ids, rownames(scale_data))
if (anyNA(residual_index) || anyDuplicated(residual_index)) stop("GetResidual lost exact500 feature", call. = FALSE)

representations <- lapply(selected_by_seed, function(selected) {
  data <- as.matrix(sct_data[data_index, selected, drop = FALSE])
  residual <- as.matrix(scale_data[residual_index, selected, drop = FALSE])
  rownames(data) <- exact$feature_id; rownames(residual) <- exact$feature_id
  list(sct_data_all_qc_fit_selected384 = data,
       sct_pearson_residual_all_qc_fit_selected384 = residual)
})
panel_ids <- if (row$dataset_scope == "external8") c("exact500", "common475") else "exact500"
summary_rows <- list()
for (seed in names(representations)) for (representation_id in names(representations[[seed]])) {
  for (panel_id in panel_ids) {
    value <- representations[[seed]][[representation_id]]
    if (panel_id == "common475") value <- value[common_index, , drop = FALSE]
    geo <- geometry(value)
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(contract_id = "mv08o_residual_source_sentinel_v1", run_role = run_role,
      unit_id = unit_id, dataset_scope = row$dataset_scope, seed = as.integer(seed),
      representation_id = representation_id, panel_id = panel_id, fit_cells = ncol(counts), selected_cells = ncol(value),
      selected_cell_sha256 = attr(selected_by_seed[[seed]], "selected_cell_sha256"),
      panel_genes = nrow(value), panel_sha256 = exact_sha, values_finite = geo$finite,
      zero_variance_gene_count = geo$zero_variance, minimum_gene_sd = geo$minimum_sd,
      correlation_chord_valid = geo$valid, matrix_sha256 = sha_object(value),
      distance_sha256 = geo$distance_sha256, sct_seconds = sct_seconds,
      sct_warning_count = length(sct_warnings),
      sct_warning_sha256 = sha_object(sort(unique(sct_warnings), method = "radix")),
      residual_seconds = residual_seconds, total_elapsed_seconds = proc.time()[["elapsed"]] - started,
      persistence_computed = FALSE, landscapes_computed = FALSE, clustering_computed = FALSE,
      fusion_computed = FALSE, outcome_label_state = "closed", biological_outcomes_computed = FALSE,
      stringsAsFactors = FALSE)
  }
}
summary <- do.call(rbind, summary_rows)
atomic_csv(summary, output_csv)
required_geometry <- if (row$dataset_scope == "internal124") {
  rep(TRUE, nrow(summary))
} else {
  summary$panel_id == "common475" |
    (summary$panel_id == "exact500" &
      summary$representation_id == "sct_pearson_residual_all_qc_fit_selected384")
}
if (nrow(summary) != row$selected_axes * 2L * length(panel_ids) ||
    !all(summary$values_finite[required_geometry]) ||
    !all(summary$zero_variance_gene_count[required_geometry] == 0L) ||
    !all(summary$correlation_chord_valid[required_geometry]) ||
    any(summary$selected_cells != 384L) ||
    any(summary$panel_genes != ifelse(summary$panel_id == "exact500", 500L, 475L))) {
  stop("MV8-O source/view geometry gate failed", call. = FALSE)
}
identity <- list(contract_id = "mv08o_residual_source_cache_v1", unit_id = unit_id,
  dataset_scope = row$dataset_scope, source_kind = source_kind, source_sha256 = sha_file(source_path),
  exact500_panel_sha256 = exact_sha, fit_cells = ncol(counts), selection_sha256 =
    vapply(selected_by_seed, attr, character(1L), which = "selected_cell_sha256"),
  normalization = list(method = "SCTransform", min_cells = 5L, variable_features_n = 3000L,
                       return_only_var_genes = FALSE, workers = 1L,
                       future_globals_max_size_bytes = future_globals_max_size_bytes,
                       warning_count = length(sct_warnings),
                       warning_sha256 = sha_object(sort(unique(sct_warnings), method = "radix"))),
  representation_ids = names(representations[[1L]]), run_role = run_role,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE)
cache <- list(contract_id = "mv08o_residual_source_cache_v1", identity = identity,
  panels = exact[c("panel_order", "feature_id", "gene")], representations = representations,
  payload_sha256 = NULL, cache_key = NULL,
  downstream_execution = list(ph_jobs = 0L, landscape_jobs = 0L, comparison_jobs = 0L,
                              clustering_jobs = 0L, fusion_jobs = 0L, label_jobs = 0L,
                              biological_outcome_jobs = 0L))
cache$payload_sha256 <- sha_object(cache[c("identity", "panels", "representations", "downstream_execution")])
cache$cache_key <- paste0("mv08o_residual_source_cache_v1:", cache$payload_sha256)
atomic_rds(cache, output_rds)
cat("MV8-O worker unit=", unit_id, "; role=", run_role, "; rows=", nrow(summary),
    "; fit_cells=", ncol(counts), "; checks=pass\n", sep = "")
