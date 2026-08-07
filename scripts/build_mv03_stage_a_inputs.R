#!/usr/bin/env Rscript

options(warn = 2)

if (!requireNamespace("digest", quietly = TRUE)) {
  stop("digest is required for MV-03 identity binding.", call. = FALSE)
}

source(file.path("R", "toy_baseline.R"), local = FALSE)
source(file.path("R", "dual_view_topology.R"), local = FALSE)
source(file.path("R", "mv03_pilot.R"), local = FALSE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop(
    "Usage: build_mv03_stage_a_inputs.R <residual-cache> ",
    "<integrated-cache> <manifest-csv> <audit-dir> <view-cache-dir>",
    call. = FALSE
  )
}

residual_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
integrated_path <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
manifest_path <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
audit_dir <- args[[4L]]
view_dir <- args[[5L]]
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(view_dir, recursive = TRUE, showWarnings = FALSE)

residual_cache <- readRDS(residual_path)
integrated_cache <- readRDS(integrated_path)
if (!identical(residual_cache$cohort, "large") ||
    !identical(integrated_cache$cohort, "large")) {
  stop("Stage A requires large-cohort caches.", call. = FALSE)
}

fit_scope_id <- "descriptive_all_pilot_samples"
seed <- 20260805L
panel_fit <- fit_mv03_gene_panel(
  residual_cache$samples, integrated_cache$samples,
  cohort = "large", fit_scope_id = fit_scope_id, panel_size = 500L
)

panel_table <- data.frame(
  order = seq_along(panel_fit$panel),
  gene = panel_fit$panel,
  cohort = "large",
  fit_scope_id = fit_scope_id,
  panel_sha256 = panel_fit$panel_sha256,
  stringsAsFactors = FALSE
)
utils::write.csv(
  panel_table,
  file.path(audit_dir, "mv03-gene-panel-large-seed-independent-2026-08-05.csv"),
  row.names = FALSE
)
utils::write.csv(
  panel_fit$eligibility,
  file.path(audit_dir, "mv03-gene-eligibility-large-2026-08-05.csv"),
  row.names = FALSE
)
utils::write.csv(
  panel_fit$ranks,
  file.path(audit_dir, "mv03-sct-variance-ranks-large-2026-08-05.csv"),
  row.names = FALSE
)

sample_ids <- sort(names(residual_cache$samples), method = "radix")
selected_cells <- lapply(sample_ids, function(sample_id) {
  select_matched_cells(
    colnames(residual_cache$samples[[sample_id]]), n = 384L, seed = seed
  )
})
names(selected_cells) <- sample_ids
subsample_table <- do.call(rbind, lapply(sample_ids, function(sample_id) {
  ids <- selected_cells[[sample_id]]
  data.frame(
    cohort = "large",
    sample_id = sample_id,
    seed = seed,
    selected_order = seq_along(ids),
    cell_id = as.character(ids),
    eligible_cell_count = attr(ids, "eligible_cell_count"),
    selected_cell_sha256 = attr(ids, "selected_cell_sha256"),
    stringsAsFactors = FALSE
  )
}))
utils::write.csv(
  subsample_table,
  file.path(audit_dir, "mv03-matched-cells-large-seed-20260805.csv"),
  row.names = FALSE
)

prepared <- prepare_mv03_sources(
  panel_fit$residual_samples, panel_fit$panel, cohort = "large",
  representation = "sct_whole", fit_scope_id = fit_scope_id, seed = seed,
  selected_cells = selected_cells
)
pca_model <- fit_cell_topology_pca(
  prepared$sources, n_components = 30L, pca_seed = seed
)

standardization_table <- data.frame(
  gene = panel_fit$panel,
  center = unname(prepared$center),
  scale = unname(prepared$scale),
  cohort = "large",
  representation = "sct_whole",
  fit_scope_id = fit_scope_id,
  seed = seed,
  standardization_id = prepared$standardization_id,
  stringsAsFactors = FALSE
)
utils::write.csv(
  standardization_table,
  file.path(audit_dir, "mv03-standardization-large-sct-whole-seed-20260805.csv"),
  row.names = FALSE
)

rotation_table <- data.frame(
  gene = rownames(pca_model$rotation), pca_model$rotation,
  check.names = FALSE, stringsAsFactors = FALSE
)
utils::write.csv(
  rotation_table,
  file.path(audit_dir, "mv03-pca-loadings-large-sct-whole-seed-20260805.csv"),
  row.names = FALSE
)
utils::write.csv(
  data.frame(
    component = seq_along(pca_model$singular_values),
    singular_value = pca_model$singular_values,
    pca_cache_key = pca_model$cache_key,
    algorithm = pca_model$algorithm,
    fit_sample_count = length(pca_model$fit_sample_ids),
    fit_scope_id = fit_scope_id,
    stringsAsFactors = FALSE
  ),
  file.path(audit_dir, "mv03-pca-summary-large-sct-whole-seed-20260805.csv"),
  row.names = FALSE
)

manifest <- utils::read.csv(
  manifest_path, stringsAsFactors = FALSE, check.names = FALSE
)
stage_a <- manifest[
  manifest$cohort == "large" &
    manifest$pilot_role %in% c("stress_cell_minimum", "stress_cell_maximum"),
  , drop = FALSE
]
if (nrow(stage_a) != 2L) {
  stop("Frozen manifest must contain exactly two Stage A stress samples.",
       call. = FALSE)
}

job_rows <- list()
for (sample_id in sort(stage_a$sample_id, method = "radix")) {
  source_object <- prepared$sources[[sample_id]]
  views <- list(
    cell_topology_v1 = construct_cell_topology_view(source_object, pca_model),
    gene_topology_v1 = construct_gene_topology_view(source_object)
  )
  for (view_id in names(views)) {
    path <- file.path(view_dir, paste(sample_id, view_id, sep = "__"))
    path <- paste0(path, ".rds")
    saveRDS(views[[view_id]], path, compress = FALSE)
    job_rows[[length(job_rows) + 1L]] <- data.frame(
      stage = "A",
      cohort = "large",
      representation = "sct_whole",
      sample_id = sample_id,
      seed = seed,
      view_id = view_id,
      point_count = length(views[[view_id]]$point_ids),
      coordinate_count = length(views[[view_id]]$coordinate_ids),
      view_cache_key = views[[view_id]]$cache_key,
      view_rds = gsub("\\\\", "/", path),
      stringsAsFactors = FALSE
    )
  }
}
jobs <- do.call(rbind, job_rows)
utils::write.csv(
  jobs,
  file.path(audit_dir, "mv03-stage-a-job-manifest-2026-08-05.csv"),
  row.names = FALSE
)
message("Built ", nrow(jobs), " scientifically eligible Stage A view jobs.")
