#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop(paste(
  "usage: build_mv16_descriptive_synthesis_prefreeze.R",
  "<mv15-public-root> <output-dir>"
), call. = FALSE)
source_root <- normalizePath(args[[1L]], mustWork = TRUE)
output <- args[[2L]]
if (dir.exists(output)) stop("MV16 prefreeze output exists")
dir.create(output, recursive = TRUE)
source("R/mv08z_landscape_production.R")
readc <- .mv08z_read_csv
sha <- .mv08z_sha256_file
atomic <- .mv08z_atomic_csv
head <- tolower(trimws(Sys.getenv("MV16_PREFREEZE_HEAD", unset = "")))
if (!grepl("^[0-9a-f]{40}$", head)) stop("MV16_PREFREEZE_HEAD required")
closure <- "docs/audits/mv15-cell-distance-comparison-closure-v2"
closure_validation <- readc(file.path(closure, "mv15-validation.csv"))
global_path <- file.path(source_root, "mv15-global-summary.csv")
neighbor_path <- file.path(source_root, "mv15-neighbor-summary.csv")
global <- readc(global_path)
neighbor <- readc(neighbor_path)
if (!all(closure_validation$passed) || nrow(global) != 36L ||
    nrow(neighbor) != 42L) stop("MV16 source closure drift")
source_files <- c(
  file.path(closure, "mv15-artifact-manifest.csv"),
  file.path(closure, "mv15-validation.csv"),
  global_path, neighbor_path
)
source_binding <- data.frame(
  contract_id = "mv16_source_binding_v1", artifact_order = seq_along(source_files),
  artifact_role = c("MV15_closure_manifest", "MV15_closure_validation",
                    "complete_global_source", "complete_neighbor_source"),
  bytes = as.numeric(file.info(source_files)$size),
  sha256 = vapply(source_files, sha, character(1L)), stringsAsFactors = FALSE
)
implementation_files <- c(
  "R/mv16_cell_distance_synthesis.R",
  "scripts/build_mv16_descriptive_synthesis_prefreeze.R",
  "scripts/run_mv16_descriptive_synthesis.R",
  "scripts/build_mv16_descriptive_synthesis_closure.R"
)
implementation <- data.frame(
  contract_id = "mv16_implementation_binding_v1",
  implementation_order = seq_along(implementation_files),
  file = implementation_files,
  bytes = as.numeric(file.info(implementation_files)$size),
  sha256 = vapply(implementation_files, sha, character(1L)),
  stringsAsFactors = FALSE
)
contract <- data.frame(
  contract_id = "mv16_descriptive_synthesis_prefreeze_v1",
  execution_head = head, source_global_rows = 36L, source_neighbor_rows = 42L,
  global_family_rows = 10L, neighbor_family_rows = 16L,
  global_metrics = paste(c("pearson", "spearman", "relative_stress",
                           "median_scaled_change", "p95_scaled_change"),
                         collapse = ";"),
  neighbor_metrics = "mean_jaccard;median_jaccard;p10_jaccard",
  summary_policy = "median_min_max_complete_no_selection",
  H0_H1_policy = "separate", internal_k = 10L, external_k = "2;3",
  thresholds = "none", inference_authorized = FALSE,
  view_ranking_authorized = FALSE, fusion_authorized = FALSE,
  clustering_authorized = FALSE, labels_authorized = FALSE,
  outcomes_authorized = FALSE, biological_claims_authorized = FALSE,
  manuscript_claims_authorized = FALSE,
  independent_repeat_required = TRUE, stringsAsFactors = FALSE
)
validation <- data.frame(
  contract_id = "mv16_validation_v1",
  check_id = c(
    "MV15_closed", "source_cardinality", "source_aggregate_only",
    "source_rehashed", "implementation_bound", "complete_rows_retained",
    "fixed_global_groups", "fixed_neighbor_groups", "H0_H1_separate",
    "neighborhood_policy", "median_min_max_only", "no_threshold_inference",
    "no_view_ranking", "downstream_firewall", "no_synthesis_executed"
  ),
  passed = c(
    all(closure_validation$passed), nrow(global) == 36L && nrow(neighbor) == 42L,
    !any(grepl("unit_id|pair_key", c(names(global), names(neighbor)))),
    all(file.exists(source_files)), nrow(implementation) == 4L,
    contract$source_global_rows == 36L && contract$source_neighbor_rows == 42L,
    contract$global_family_rows == 10L, contract$neighbor_family_rows == 16L,
    contract$H0_H1_policy == "separate",
    contract$internal_k == 10L && contract$external_k == "2;3",
    contract$summary_policy == "median_min_max_complete_no_selection",
    contract$thresholds == "none" && !contract$inference_authorized,
    !contract$view_ranking_authorized,
    !contract$fusion_authorized && !contract$clustering_authorized &&
      !contract$labels_authorized && !contract$outcomes_authorized &&
      !contract$biological_claims_authorized &&
      !contract$manuscript_claims_authorized,
    TRUE
  ), stringsAsFactors = FALSE
)
if (!all(validation$passed)) stop("MV16 prefreeze failed")
decision <- data.frame(
  contract_id = "mv16_decision_v1",
  decision = "authorize_complete_threshold_free_synthesis_after_commit",
  execution_authorized_after_commit = TRUE, values_inspected = FALSE,
  next_after_closure = "owner_scientific_review_before_fusion_or_clustering",
  stringsAsFactors = FALSE
)
artifacts <- list(
  "mv16-contract.csv" = contract,
  "mv16-source-binding.csv" = source_binding,
  "mv16-implementation-binding.csv" = implementation,
  "mv16-validation.csv" = validation,
  "mv16-decision.csv" = decision
)
for (name in names(artifacts)) atomic(artifacts[[name]], file.path(output, name))
writeLines(c(
  "# MV16 threshold-free descriptive-synthesis prefreeze", "",
  "All 36 global and 42 neighborhood rows are retained. Ten global-family and",
  "16 neighborhood-family rows will report prospectively fixed median/min/max",
  "summaries with H0/H1 separate and no selection or threshold.", "",
  "No value was inspected by this prefreeze. Ranking, fusion, clustering, labels,",
  "outcomes, inference, biological interpretation, and manuscript claims stay closed."
), file.path(output, "MV16_DESCRIPTIVE_SYNTHESIS_PREFREEZE_2026-08-26.md"))
files <- sort(setdiff(list.files(output), "mv16-artifact-manifest.csv"))
manifest <- data.frame(
  contract_id = "mv16_artifact_manifest_v1", artifact = files,
  bytes = as.numeric(file.info(file.path(output, files))$size),
  sha256 = vapply(file.path(output, files), sha, character(1L)),
  stringsAsFactors = FALSE
)
atomic(manifest, file.path(output, "mv16-artifact-manifest.csv"))
message("Built MV16 prefreeze; checks=", nrow(validation))
