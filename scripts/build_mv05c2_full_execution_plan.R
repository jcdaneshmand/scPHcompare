#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(
    "usage: build_mv05c2_full_execution_plan.R METADATA MV05A_LOSO_MANIFEST AUDIT_DIR",
    call. = FALSE
  )
}
metadata_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
loso_path <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
audit_dir <- args[[3L]]
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

source("R/provenance_utils.R")
source("R/mv05_resource_safe_execution.R")

frozen_studies <- c(
  "SRA550660", "SRA621638", "SRA628554", "SRA645804", "SRA667709",
  "SRA703206", "SRA713577", "SRA716608", "SRA728025", "SRA740174",
  "SRA749327", "SRA760933", "SRA779509", "SRA784242", "SRA826293"
)
seeds <- 20260805:20260809
metadata <- utils::read.csv(
  metadata_path, stringsAsFactors = FALSE, check.names = TRUE
)
loso <- utils::read.csv(loso_path, stringsAsFactors = FALSE, check.names = FALSE)
base <- unique(loso[c("sample_id", "sample_study")])
candidate <- base[base$sample_study %in% frozen_studies, , drop = FALSE]
candidate <- candidate[order(
  candidate$sample_study, candidate$sample_id, method = "radix"
), , drop = FALSE]
filtered <- stats::setNames(
  as.integer(metadata$Number_of_Cells_After_Filtering), metadata$orig.ident
)
candidate$filtered_cells <- unname(filtered[candidate$sample_id])
if (nrow(candidate) != 90L || length(unique(candidate$sample_study)) != 15L ||
    anyNA(candidate$filtered_cells) || any(candidate$filtered_cells < 384L)) {
  stop("Frozen MV5-D candidate set is not 90 eligible samples/15 studies.",
       call. = FALSE)
}
names(candidate)[names(candidate) == "sample_study"] <- "study"
candidate$contract_id <- "mv05c2_label_closed_candidate_manifest_v1"
candidate$outcome_label_state <- "closed"
candidate$biological_outcomes_computed <- FALSE
candidate <- candidate[c(
  "contract_id", "sample_id", "study", "filtered_cells",
  "outcome_label_state", "biological_outcomes_computed"
)]

fold_rows <- list()
for (held_out in frozen_studies) {
  fold_rows[[length(fold_rows) + 1L]] <- data.frame(
    fold_id = paste0("large_loso_v1:", held_out),
    fit_scope_id = paste0("large_loso_v1:", held_out, ":training"),
    held_out_study = held_out,
    sample_id = candidate$sample_id,
    execution_role = ifelse(
      candidate$study == held_out, "held_out_query", "training_reference"
    ),
    outcome_label_state = "closed", stringsAsFactors = FALSE
  )
}
fold_table <- do.call(rbind, fold_rows)

normalization <- expand.grid(
  sample_id = candidate$sample_id, seed = seeds,
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
)
normalization <- normalization[order(
  normalization$seed, normalization$sample_id, method = "radix"
), , drop = FALSE]
normalization$contract_id <- "mv05c2_normalization_plan_v1"
normalization$selected_cells <- 384L
normalization$reuse_scope <- "sample_seed_across_all_folds"
normalization$planned_private_cache_file <- paste0(
  gsub("[^A-Za-z0-9_.-]", "_", normalization$sample_id), "__",
  normalization$seed, "__sct.rds"
)
normalization$outcome_label_state <- "closed"
normalization$biological_outcomes_computed <- FALSE
normalization <- normalization[c(
  "contract_id", "sample_id", "seed", "selected_cells", "reuse_scope",
  "planned_private_cache_file", "outcome_label_state",
  "biological_outcomes_computed"
)]

fold_seed <- do.call(rbind, lapply(frozen_studies, function(held_out) {
  data.frame(
    contract_id = "mv05c2_fold_seed_execution_plan_v1",
    fold_id = paste0("large_loso_v1:", held_out),
    fit_scope_id = paste0("large_loso_v1:", held_out, ":training"),
    held_out_study = held_out, seed = seeds,
    training_samples = sum(candidate$study != held_out),
    held_out_samples = sum(candidate$study == held_out),
    cached_normalization_inputs = 90L,
    execution_unit = "fold_seed_representation_then_sample_diagram_chunks",
    maximum_heavy_workers = 2L,
    elapsed_cap_seconds = 1800,
    rss_cap_bytes = 8 * 1024^3,
    outcome_label_state = "closed",
    biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
}))

strata <- data.frame(
  representation = c(
    "sct_fold", "sct_fold", "inductive_integrated",
    "inductive_integrated", "sct_fold", "sct_fold"
  ),
  view_id = c(
    "cell_topology_v1", "cell_topology_v1", "cell_topology_v1",
    "cell_topology_v1", "gene_topology_v1", "gene_topology_v1"
  ),
  homology_dimension = rep(c("H0", "H1"), 3L),
  stringsAsFactors = FALSE
)
pairs <- mv05c2_build_query_training_pairs_v1(fold_table, seeds, strata)
chunked <- mv05c2_assign_pair_chunks_v1(pairs, max_pairs = 250L)
chunks <- mv05c2_pair_chunk_summary_v1(chunked)
scope_groups <- split(
  chunked,
  interaction(
    chunked$fold_id, chunked$seed, chunked$representation,
    chunked$view_id, chunked$homology_dimension, drop = TRUE, lex.order = TRUE
  )
)
scope <- do.call(rbind, lapply(scope_groups, function(group) data.frame(
  contract_id = "mv05c2_pair_scope_summary_v1",
  fold_id = group$fold_id[[1L]], fit_scope_id = group$fit_scope_id[[1L]],
  seed = group$seed[[1L]], representation = group$representation[[1L]],
  view_id = group$view_id[[1L]],
  homology_dimension = group$homology_dimension[[1L]],
  held_out_samples = length(unique(group$query_sample_id)),
  training_samples = length(unique(group$training_sample_id)),
  query_training_pairs = nrow(group),
  full_matrix_pairs_90_samples = choose(90L, 2L),
  excluded_full_matrix_pairs = choose(90L, 2L) - nrow(group),
  supports_primary_retrieval = TRUE,
  supports_full_matrix_clustering = FALSE,
  supports_within_study_pair_contrasts = FALSE,
  outcome_label_state = "closed",
  biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)))
rownames(scope) <- NULL

projection <- data.frame(
  contract_id = "mv05c2_architectural_reduction_projection_v1",
  candidate_samples = 90L, folds = 15L, seeds = 5L,
  naive_sample_normalizations = 90L * 15L * 5L,
  cached_sample_seed_normalizations = nrow(normalization),
  normalization_operation_reduction_factor =
    (90L * 15L * 5L) / nrow(normalization),
  naive_landscape_pair_rows = 1802250L,
  query_training_landscape_pair_rows = nrow(pairs),
  landscape_pair_row_reduction_factor = 1802250 / nrow(pairs),
  pair_chunk_size = 250L, pair_chunks = nrow(chunks),
  full_matrix_clustering_disposition = "deferred_not_supported_by_pair_scope",
  within_study_contrast_disposition = "deferred_not_supported_by_pair_scope",
  outcome_label_state = "closed",
  biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
write_provenance_csv(
  candidate,
  file.path(audit_dir, "mv05c2-label-closed-candidate-manifest-2026-08-06.csv")
)
write_provenance_csv(
  normalization,
  file.path(audit_dir, "mv05c2-normalization-plan-2026-08-06.csv")
)
write_provenance_csv(
  fold_seed,
  file.path(audit_dir, "mv05c2-fold-seed-execution-plan-2026-08-06.csv")
)
write_provenance_csv(
  scope,
  file.path(audit_dir, "mv05c2-pair-scope-summary-2026-08-06.csv")
)
write_provenance_csv(
  chunks,
  file.path(audit_dir, "mv05c2-pair-chunk-manifest-2026-08-06.csv")
)
write_provenance_csv(
  projection,
  file.path(audit_dir, "mv05c2-architectural-reduction-projection-2026-08-06.csv")
)
message(
  "Built label-closed 90-sample resource-safe plan with ", nrow(pairs),
  " requested landscape pairs in ", nrow(chunks), " chunks."
)
