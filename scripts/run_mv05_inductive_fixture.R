#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("usage: run_mv05_inductive_fixture.R OUTPUT_CSV", call. = FALSE)
}
output_csv <- args[[1L]]
dir.create(dirname(output_csv), recursive = TRUE, showWarnings = FALSE)

source("R/provenance_utils.R")
source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/topological_distance_contract.R")
source("R/mv05_benchmark_contract.R")
source("R/mv05_benchmark_execution.R")
source("R/mv05_inductive_mapping.R")

seed <- 20260805L
genes <- c(
  paste0("RPL", seq_len(30L)), paste0("RPS", seq_len(30L)),
  paste0("MT-MV05", seq_len(5L)), paste0("GENE", seq_len(455L))
)
counts <- .with_preserved_seed(seed, matrix(
  stats::rpois(length(genes) * 120L, lambda = 3) + 1L,
  nrow = length(genes),
  dimnames = list(genes, paste0("cell", seq_len(120L)))
))
counts[71:120, 81:120] <- counts[71:120, 81:120] + 1L
reference_counts <- Matrix::Matrix(counts[, seq_len(80L)], sparse = TRUE)
query_counts <- Matrix::Matrix(counts[, 81:120], sparse = TRUE)
colnames(reference_counts) <- paste0("reference_cell", seq_len(ncol(reference_counts)))
colnames(query_counts) <- paste0("query_cell", seq_len(ncol(query_counts)))

reference <- Seurat::CreateSeuratObject(reference_counts, project = "mv05_reference")
query <- Seurat::CreateSeuratObject(query_counts, project = "mv05_query")
reference <- Seurat::SCTransform(
  reference, variable.features.n = 200L, return.only.var.genes = FALSE,
  verbose = FALSE, seed.use = seed
)
query <- Seurat::SCTransform(
  query, variable.features.n = 200L, return.only.var.genes = FALSE,
  verbose = FALSE, seed.use = seed
)
features <- Seurat::VariableFeatures(reference)
features <- features[features %in% rownames(Seurat::GetAssayData(
  query, assay = "SCT", layer = "data"
))]
features <- sort(features, method = "radix")
reference <- Seurat::RunPCA(
  reference, assay = "SCT", features = features, npcs = 10L,
  approx = FALSE, seed.use = seed, verbose = FALSE
)

run_once <- function() mv05_try_inductive_mapping_v1(
  reference = reference, query = query, features = features,
  dimensions = 1:10, fold_id = "synthetic_two_study_fixture:query",
  fit_scope_id = "synthetic_two_study_fixture:reference_training",
  reference_sample_ids = "synthetic_reference_study",
  held_out_sample_id = "synthetic_query_study", seed = seed,
  k_anchor = 3L, k_score = 10L, k_weight = 20L, verbose = FALSE
)
first <- run_once()
second <- run_once()
completed <- identical(first$status, "completed") &&
  identical(second$status, "completed")
first_hash <- if (completed) first$result$query_embedding_sha256 else ""
second_hash <- if (completed) second$result$query_embedding_sha256 else ""
manifest <- data.frame(
  contract_id = "mv05_inductive_mapping_v1",
  fixture_id = "synthetic_label_closed_two_study_v1",
  status = if (completed && identical(first_hash, second_hash))
    "completed_deterministic" else if (completed) "nondeterministic" else
      "held_out_mapping_unavailable",
  first_reason = first$reason,
  second_reason = second$reason,
  reference_cells = ncol(reference), query_cells = ncol(query),
  reference_selected_features = length(features), dimensions = 10L,
  first_anchor_count = if (completed) first$result$anchor_count else NA_integer_,
  second_anchor_count = if (completed) second$result$anchor_count else NA_integer_,
  first_embedding_sha256 = first_hash, second_embedding_sha256 = second_hash,
  reference_identity_sha256 = if (completed)
    first$result$reference_identity_sha256 else "",
  first_elapsed_seconds = first$elapsed_seconds,
  second_elapsed_seconds = second$elapsed_seconds,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  seurat_version = as.character(utils::packageVersion("Seurat")),
  seurat_object_version = as.character(utils::packageVersion("SeuratObject")),
  r_version = paste(R.version$major, R.version$minor, sep = "."),
  stringsAsFactors = FALSE
)
utils::write.csv(manifest, output_csv, row.names = FALSE, na = "")
message("MV5-B inductive fixture status: ", manifest$status)
