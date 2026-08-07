#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(paste(
    "usage: build_mv05_statistical_plan.R HISTORICAL_DIR PILOT_MANIFEST",
    "AUDIT_DIR"
  ), call. = FALSE)
}
if (!requireNamespace("digest", quietly = TRUE)) stop("digest is required.")

historical_dir <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
pilot_path <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
audit_dir <- args[[3L]]
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

source("R/topological_distance_contract.R")
source("R/mv05_benchmark_contract.R")

large_path <- file.path(historical_dir, "joined_metadata_cellcounts.csv")
bone_path <- file.path(historical_dir, "joined_metadata_cellcounts_bonemarrow.csv")
if (!file.exists(large_path) || !file.exists(bone_path)) {
  stop("Required existing-data metadata files are missing.", call. = FALSE)
}

large_raw <- read.csv(large_path, stringsAsFactors = FALSE, check.names = TRUE)
bone_raw <- read.csv(bone_path, stringsAsFactors = FALSE, check.names = TRUE)
large <- data.frame(
  cohort = "large", sample_id = large_raw$orig.ident,
  study = large_raw$SRA, tissue = tolower(large_raw$Tissue.x),
  approach = large_raw$Approach.x,
  filtered_cells = large_raw$Number_of_Cells_After_Filtering,
  stringsAsFactors = FALSE
)
bone <- data.frame(
  cohort = "bone", sample_id = bone_raw$orig.ident,
  study = bone_raw$SRA, tissue = tolower(bone_raw$Tissue),
  approach = bone_raw$Approach,
  filtered_cells = bone_raw$Number_of_Cells_After_Filtering,
  stringsAsFactors = FALSE
)
metadata <- .mv05_validate_metadata(rbind(large, bone))
pilot <- read.csv(pilot_path, stringsAsFactors = FALSE, check.names = FALSE)
pilot_metadata <- .mv05_validate_metadata(pilot)

contract <- new_mv05_benchmark_contract_v1()
feasibility <- mv05_design_feasibility_v1(metadata)
pilot_feasibility <- mv05_design_feasibility_v1(pilot_metadata)
loso <- mv05_loso_fold_summary_v1(metadata)

cohort_groups <- split(metadata, metadata$cohort)
cohort_summary <- do.call(rbind, lapply(cohort_groups, function(values) {
  tissue_part <- feasibility[feasibility$cohort == values$cohort[[1L]], ]
  eligible_tissues <- tissue_part$tissue[tissue_part$cross_study_tissue_eligible]
  data.frame(
    cohort = values$cohort[[1L]], sample_count = nrow(values),
    study_count = length(unique(values$study)),
    tissue_count = length(unique(values$tissue)),
    approach_count = length(unique(values$approach)),
    cross_study_eligible_tissues = length(eligible_tissues),
    cross_study_eligible_samples = sum(values$tissue %in% eligible_tissues),
    confirmatory_role = if (values$cohort[[1L]] == "large")
      "cross-study tissue candidate subset" else "technical approach only",
    stringsAsFactors = FALSE
  )
}))
rownames(cohort_summary) <- NULL

integration_code <- readLines("R/Integration_flexible.R", warn = FALSE)
integration_audit <- data.frame(
  code_path = "R/Integration_flexible.R",
  find_integration_anchors = any(grepl("FindIntegrationAnchors\\(", integration_code)),
  integrate_data = any(grepl("IntegrateData\\(", integration_code)),
  find_transfer_anchors = any(grepl("FindTransferAnchors\\(", integration_code)),
  map_query = any(grepl("MapQuery\\(", integration_code)),
  transfer_data = any(grepl("TransferData\\(", integration_code)),
  disposition = paste(
    "current integration is transductive; confirmatory LOSO integration",
    "requires an audited training-reference/held-out mapping implementation"
  ),
  stringsAsFactors = FALSE
)

legacy_audit <- data.frame(
  code_object = c(
    "apply_all_clustering_methods", "assign_ph_clusters",
    "enhanced_cluster_comparison_with_pvals", "random p-value fraction",
    "perform_hierarchical_clustering_ph", "perform_spectral_clustering"
  ),
  observed_behavior = c(
    "derives k_tissue/k_study/k_approach directly from evaluation labels",
    "copies one sample cluster assignment to every cell",
    "evaluates copied sample assignments over cell rows with ARI/NMI/Jaccard/VI/purity",
    "uses mean(null >= observed), permitting exact zero and unblocked nulls",
    "defaults include ward.D2/complete/mcquitty on arbitrary dissimilarities",
    "uses a Gaussian kernel and kernlab specc without frozen seed/eigengap provenance"
  ),
  mv05_disposition = c(
    "oracle-k historical sensitivity only", "prohibited for primary metrics",
    "prohibited for MV-05 primary metrics", "replace with blocked paired inference",
    "average linkage only", "replace with spectral_gaussian_laplacian_v1 before use"
  ),
  stringsAsFactors = FALSE
)

contract_summary <- data.frame(
  contract_id = contract$contract_id, cache_key = contract$cache_key,
  outer_split = contract$outer_split,
  confirmatory_tissue_rule = contract$confirmatory_tissue_rule,
  subsample_seeds = paste(contract$subsample_seeds, collapse = ";"),
  bootstrap_replicates = contract$bootstrap_replicates,
  monte_carlo_replicates = contract$monte_carlo_replicates,
  pilot_policy = contract$pilot_policy,
  integration_policy = contract$integration_policy,
  large_metadata_sha256 = digest::digest(file = large_path, algo = "sha256"),
  bone_metadata_sha256 = digest::digest(file = bone_path, algo = "sha256"),
  pilot_manifest_sha256 = digest::digest(file = pilot_path, algo = "sha256"),
  stringsAsFactors = FALSE
)

outputs <- list(
  "mv05-contract-summary-2026-08-06.csv" = contract_summary,
  "mv05-endpoint-registry-2026-08-06.csv" = contract$endpoint_registry,
  "mv05-method-registry-2026-08-06.csv" = contract$method_registry,
  "mv05-label-use-registry-2026-08-06.csv" = contract$label_use_registry,
  "mv05-multiplicity-registry-2026-08-06.csv" = contract$multiplicity_registry,
  "mv05-existing-data-cohort-feasibility-2026-08-06.csv" = cohort_summary,
  "mv05-existing-data-tissue-feasibility-2026-08-06.csv" = feasibility,
  "mv05-pilot-tissue-feasibility-2026-08-06.csv" = pilot_feasibility,
  "mv05-loso-fold-feasibility-2026-08-06.csv" = loso,
  "mv05-integration-induction-audit-2026-08-06.csv" = integration_audit,
  "mv05-legacy-evaluation-disposition-2026-08-06.csv" = legacy_audit
)
for (name in names(outputs)) {
  write.csv(outputs[[name]], file.path(audit_dir, name), row.names = FALSE)
}

manifest_inputs <- c(
  "R/mv05_benchmark_contract.R",
  "scripts/build_mv05_statistical_plan.R",
  "tests/testthat/test-mv05-benchmark-contract.R",
  "docs/specifications/MV05_STATISTICAL_BENCHMARK_PLAN_V1.md",
  file.path(audit_dir, names(outputs))
)
manifest <- data.frame(
  artifact_path = gsub("\\\\", "/", manifest_inputs),
  sha256 = vapply(manifest_inputs, digest::digest, character(1L),
                  file = TRUE, algo = "sha256"),
  bytes = as.numeric(file.info(manifest_inputs)$size),
  artifact_role = c(
    "executable_contract", "evidence_builder", "contract_tests",
    "frozen_specification", rep("generated_design_evidence", length(outputs))
  ),
  biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
write.csv(
  manifest,
  file.path(audit_dir, "mv05-artifact-manifest-2026-08-06.csv"),
  row.names = FALSE
)

stopifnot(
  nrow(large) == 124L, nrow(bone) == 25L,
  sum(feasibility$cross_study_tissue_eligible) == 5L,
  sum(metadata$cohort == "large" & metadata$tissue %in%
        feasibility$tissue[feasibility$cross_study_tissue_eligible]) == 90L,
  sum(pilot_feasibility$cross_study_tissue_eligible) == 0L,
  length(unique(large$study)) == 18L,
  nrow(loso) == length(unique(large$study)),
  all(loso$no_study_overlap),
  !integration_audit$find_transfer_anchors,
  !integration_audit$map_query
)
message("Built MV-05 statistical-plan evidence without computing outcomes.")
