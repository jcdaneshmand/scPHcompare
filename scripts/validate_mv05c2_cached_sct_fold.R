#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(
    "usage: validate_mv05c2_cached_sct_fold.R CACHE_DIR MV05C_JOB_DIR OUTPUT_CSV",
    call. = FALSE
  )
}
cache_dir <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
job_dir <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
output_path <- args[[3L]]

source("R/provenance_utils.R")
source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/topological_distance_contract.R")
source("R/mv03_pilot.R")
source("R/mv05_benchmark_execution.R")
source("R/mv05_resource_safe_execution.R")

seed <- 20260805L
cache_paths <- sort(list.files(
  cache_dir, pattern = paste0("__", seed, "__sct\\.rds$"), full.names = TRUE
), method = "radix")
records <- lapply(cache_paths, readRDS)
names(records) <- vapply(
  records, function(record) record$identity$sample_id, character(1L)
)
records <- records[order(names(records), method = "radix")]
if (length(records) != 6L) {
  stop("Cached SCT fold validation requires six sample caches.", call. = FALSE)
}

job_paths <- sort(list.files(
  job_dir, pattern = paste0("__", seed, "\\.rds$"), full.names = TRUE
), method = "radix")
if (length(job_paths) != 2L) {
  stop("Cached SCT fold validation requires both MV5-C LOSO directions.",
       call. = FALSE)
}
rows <- list()
for (job_path in job_paths) {
  reference <- readRDS(job_path)
  started <- proc.time()[["elapsed"]]
  observed <- mv05c2_run_cached_sct_fold_v1(
    records = records,
    training_ids = reference$training_sample_ids,
    fold_id = reference$fold_id,
    fit_scope_id = reference$fit_scope_id,
    seed = reference$seed,
    cohort = "mv05c_one_tissue_feasibility_v1",
    panel_size = 500L, n_components = 30L
  )
  elapsed <- proc.time()[["elapsed"]] - started
  if (!identical(observed$panel, reference$sct_fold$panel)) {
    stop("Cached fold feature panel differs from MV5-C.", call. = FALSE)
  }
  for (sample_id in sort(names(observed$cell_diagrams), method = "radix")) {
    observed_diagram <- observed$cell_diagrams[[sample_id]]$diagram
    expected_diagram <- reference$sct_fold$cell_diagrams[[sample_id]]$diagram
    same_axes <- identical(dim(observed_diagram), dim(expected_diagram)) &&
      identical(dimnames(observed_diagram), dimnames(expected_diagram))
    same_nonfinite <- same_axes && identical(
      is.finite(observed_diagram), is.finite(expected_diagram)
    )
    finite <- is.finite(observed_diagram) & is.finite(expected_diagram)
    maximum_difference <- if (same_nonfinite && any(finite)) {
      max(abs(observed_diagram[finite] - expected_diagram[finite]))
    } else if (same_nonfinite) {
      0
    } else {
      Inf
    }
    rows[[length(rows) + 1L]] <- data.frame(
      contract_id = "mv05c2_cached_sct_fold_equivalence_v1",
      fold_id = reference$fold_id, seed = reference$seed,
      artifact_role = "cell_topology_diagram", sample_id = sample_id,
      method_id = "cell_topology_v1",
      compared_values = length(observed_diagram),
      maximum_absolute_difference = maximum_difference,
      exact_numeric_identity = same_nonfinite && maximum_difference == 0,
      cached_fold_seconds = elapsed,
      outcome_label_state = "closed",
      biological_outcomes_computed = FALSE,
      stringsAsFactors = FALSE
    )
  }
  observed_gene <- observed$gene_diagrams
  expected_gene <- reference$sct_fold$gene_diagrams
  if (xor(is.null(observed_gene), is.null(expected_gene))) {
    stop("Cached fold gene eligibility differs from MV5-C.", call. = FALSE)
  }
  if (!is.null(observed_gene)) {
    for (sample_id in sort(names(observed_gene), method = "radix")) {
      observed_diagram <- observed_gene[[sample_id]]$diagram
      expected_diagram <- expected_gene[[sample_id]]$diagram
      same_axes <- identical(dim(observed_diagram), dim(expected_diagram)) &&
        identical(dimnames(observed_diagram), dimnames(expected_diagram))
      same_nonfinite <- same_axes && identical(
        is.finite(observed_diagram), is.finite(expected_diagram)
      )
      finite <- is.finite(observed_diagram) & is.finite(expected_diagram)
      maximum_difference <- if (same_nonfinite && any(finite)) {
        max(abs(observed_diagram[finite] - expected_diagram[finite]))
      } else if (same_nonfinite) {
        0
      } else {
        Inf
      }
      rows[[length(rows) + 1L]] <- data.frame(
        contract_id = "mv05c2_cached_sct_fold_equivalence_v1",
        fold_id = reference$fold_id, seed = reference$seed,
        artifact_role = "gene_topology_diagram", sample_id = sample_id,
        method_id = "gene_topology_v1",
        compared_values = length(observed_diagram),
        maximum_absolute_difference = maximum_difference,
        exact_numeric_identity = same_nonfinite && maximum_difference == 0,
        cached_fold_seconds = elapsed,
        outcome_label_state = "closed",
        biological_outcomes_computed = FALSE,
        stringsAsFactors = FALSE
      )
    }
  }
  expected_baselines <- reference$sct_fold$baselines
  if (!setequal(names(observed$baselines), names(expected_baselines))) {
    stop("Cached fold baseline availability differs from MV5-C.", call. = FALSE)
  }
  for (method in sort(names(observed$baselines), method = "radix")) {
    first <- observed$baselines[[method]]$distance_matrix
    second <- expected_baselines[[method]]$distance_matrix
    same_axes <- identical(dim(first), dim(second)) &&
      identical(dimnames(first), dimnames(second))
    maximum_difference <- if (same_axes) max(abs(first - second)) else Inf
    rows[[length(rows) + 1L]] <- data.frame(
      contract_id = "mv05c2_cached_sct_fold_equivalence_v1",
      fold_id = reference$fold_id, seed = reference$seed,
      artifact_role = "matched_baseline_matrix", sample_id = "",
      method_id = observed$baselines[[method]]$method_id,
      compared_values = length(first),
      maximum_absolute_difference = maximum_difference,
      exact_numeric_identity = same_axes && maximum_difference == 0,
      cached_fold_seconds = elapsed,
      outcome_label_state = "closed",
      biological_outcomes_computed = FALSE,
      stringsAsFactors = FALSE
    )
  }
}
result <- do.call(rbind, rows)
if (!all(result$exact_numeric_identity)) {
  stop("Cached SCT fold differs numerically from immutable MV5-C.", call. = FALSE)
}
write_provenance_csv(result, output_path)
message("Validated cached SCT fold equivalence in both LOSO directions.")
