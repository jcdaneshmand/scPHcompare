#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("usage: build_mv05ab_execution_scaffold.R HISTORICAL_DIR AUDIT_DIR",
       call. = FALSE)
}
historical_dir <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
audit_dir <- args[[2L]]
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/topological_distance_contract.R")
source("R/mv05_benchmark_contract.R")
source("R/mv05_benchmark_execution.R")
source("R/mv05_inductive_mapping.R")

large_path <- file.path(historical_dir, "joined_metadata_cellcounts.csv")
large_raw <- utils::read.csv(
  large_path, stringsAsFactors = FALSE, check.names = TRUE
)
large <- data.frame(
  cohort = "large", sample_id = large_raw$orig.ident,
  study = large_raw$SRA, tissue = tolower(large_raw$Tissue.x),
  approach = large_raw$Approach.x,
  filtered_cells = large_raw$Number_of_Cells_After_Filtering,
  stringsAsFactors = FALSE
)
loso <- mv05_build_loso_manifest_v1(large)
loso_summary <- mv05_loso_manifest_summary_v1(loso)

loso_path <- file.path(
  audit_dir, "mv05a-loso-execution-manifest-2026-08-06.csv"
)
loso_summary_path <- file.path(
  audit_dir, "mv05a-loso-execution-summary-2026-08-06.csv"
)
utils::write.csv(loso$table, loso_path, row.names = FALSE)
utils::write.csv(loso_summary, loso_summary_path, row.names = FALSE)

genes <- paste0("g", 1:4)
cells <- paste0("c", 1:5)
first <- matrix(
  c(1, 2, 4, 7, 11, 2, 5, 3, 8, 13, 4, 1, 6, 9, 15, 3, 8, 2, 12, 5),
  nrow = 4, byrow = TRUE, dimnames = list(genes, cells)
)
second <- first + matrix(
  c(0, 1, 0, 1, 0, 1, 0, 2, 0, 1, 0, 2, 1, 0, 1, 2, 0, 1, 1, 0),
  nrow = 4, byrow = TRUE
)
third <- first + matrix(
  c(1, 0, 2, 0, 1, 0, 1, 0, 2, 0, 2, 0, 1, 0, 2, 0, 1, 2, 0, 1),
  nrow = 4, byrow = TRUE
)
matrices <- list(sample_b = second, sample_a = first, sample_c = third)
sources <- lapply(names(matrices), function(id) new_dual_view_source(
  matrices[[id]], sample_id = id, cohort = "fixture",
  representation = "sct_fixture", fit_scope_id = "fixture_fold:training",
  subsample_seed = 20260805L, standardization_id = "fixture_scaling_v1",
  contract_profile = "analytical_fixture", expected_genes = 4L,
  expected_cells = 5L, expected_pcs = 2L
))
pca <- fit_cell_topology_pca(sources, n_components = 2L)
cell_views <- lapply(sources, construct_cell_topology_view, pca_model = pca)
bundles <- list(
  mv05_cell_energy_baseline_v1(cell_views),
  mv05_pseudobulk_baseline_v1(sources),
  mv05_gene_correlation_baseline_v1(sources)
)
baseline_validation <- do.call(rbind, lapply(bundles, function(bundle) {
  mv05_validate_baseline_bundle_v1(bundle)
  values <- bundle$distance_matrix[upper.tri(bundle$distance_matrix)]
  data.frame(
    method_id = bundle$method_id, formula_id = bundle$formula_id,
    sample_count = length(bundle$sample_ids), pair_count = length(values),
    minimum_distance = min(values), maximum_distance = max(values),
    exact_zero_diagonal = all(diag(bundle$distance_matrix) == 0),
    symmetric = identical(bundle$distance_matrix, t(bundle$distance_matrix)),
    fit_scope_id = bundle$fit_scope_id,
    distance_sha256 = bundle$distance_sha256,
    cache_key = bundle$cache_key,
    analytical_fixture_only = TRUE,
    biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
}))
baseline_path <- file.path(
  audit_dir, "mv05a-baseline-analytical-validation-2026-08-06.csv"
)
utils::write.csv(baseline_validation, baseline_path, row.names = FALSE)

shape_seed <- 20260805L
shape_genes <- paste0("shape_gene_", seq_len(500L))
shape_matrices <- .with_preserved_seed(shape_seed, lapply(seq_len(3L), function(i) {
  value <- matrix(stats::rnorm(500L * 384L), nrow = 500L)
  value[seq_len(50L), ] <- value[seq_len(50L), ] + (i - 1L) * 0.15
  rownames(value) <- shape_genes
  colnames(value) <- paste0("shape_sample_", i, "_cell_", seq_len(384L))
  value
}))
names(shape_matrices) <- paste0("shape_sample_", seq_len(3L))
shape_sources <- lapply(names(shape_matrices), function(id) new_dual_view_source(
  shape_matrices[[id]], sample_id = id, cohort = "synthetic_shape_fixture",
  representation = "sct_shape_fixture",
  fit_scope_id = "synthetic_shape_fixture:training",
  subsample_seed = shape_seed, standardization_id = "shape_fixture_scaling_v1"
))
pca_started <- proc.time()[["elapsed"]]
shape_pca <- fit_cell_topology_pca(shape_sources, n_components = 30L)
pca_elapsed <- proc.time()[["elapsed"]] - pca_started
shape_views <- lapply(
  shape_sources, construct_cell_topology_view, pca_model = shape_pca
)
timed_baseline <- function(method_id, code) {
  started <- proc.time()[["elapsed"]]
  bundle <- force(code)
  elapsed <- proc.time()[["elapsed"]] - started
  values <- bundle$distance_matrix[upper.tri(bundle$distance_matrix)]
  data.frame(
    method_id = method_id, genes = 500L, cells_per_sample = 384L,
    pcs = 30L, samples = 3L, pairs = length(values),
    pca_fit_elapsed_seconds = pca_elapsed,
    baseline_elapsed_seconds = elapsed,
    bundle_bytes = as.numeric(utils::object.size(bundle)),
    minimum_distance = min(values), maximum_distance = max(values),
    exact_zero_diagonal = all(diag(bundle$distance_matrix) == 0),
    symmetric = identical(bundle$distance_matrix, t(bundle$distance_matrix)),
    distance_sha256 = bundle$distance_sha256,
    cache_key = bundle$cache_key,
    synthetic_scientific_shape_fixture_only = TRUE,
    biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
}
shape_validation <- do.call(rbind, list(
  timed_baseline(
    "cell_distribution_energy_shared_pca_v1",
    mv05_cell_energy_baseline_v1(shape_views)
  ),
  timed_baseline(
    "pseudobulk_shared_panel_euclidean_v1",
    mv05_pseudobulk_baseline_v1(shape_sources)
  ),
  timed_baseline(
    "gene_correlation_frobenius_v1",
    mv05_gene_correlation_baseline_v1(shape_sources)
  )
))
shape_path <- file.path(
  audit_dir, "mv05a-baseline-scientific-shape-feasibility-2026-08-06.csv"
)
utils::write.csv(shape_validation, shape_path, row.names = FALSE)

stopifnot(
  nrow(loso$table) == 124L * 18L,
  nrow(loso_summary) == 18L,
  all(loso_summary$total_samples == 124L),
  sum(loso_summary$held_out_samples) == 124L,
  !any(c("tissue", "approach") %in% names(loso$table)),
  all(loso$table$outcome_label_state == "closed"),
  nrow(baseline_validation) == 3L,
  all(baseline_validation$exact_zero_diagonal),
  all(baseline_validation$symmetric),
  all(!baseline_validation$biological_outcomes_computed),
  nrow(shape_validation) == 3L,
  all(shape_validation$exact_zero_diagonal),
  all(shape_validation$symmetric),
  all(!shape_validation$biological_outcomes_computed)
)
message("Built MV5-A execution scaffold without computing biological outcomes.")
