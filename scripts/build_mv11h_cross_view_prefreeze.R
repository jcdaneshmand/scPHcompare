#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop(paste(
  "usage: build_mv11h_cross_view_prefreeze.R <gene-partitions>",
  "<cell-partitions> <output> <execution-head>"
), call. = FALSE)
gene_path <- normalizePath(args[[1L]], mustWork = TRUE)
cell_path <- normalizePath(args[[2L]], mustWork = TRUE)
output <- args[[3L]]
execution_head <- tolower(trimws(args[[4L]]))
if (dir.exists(output) || !grepl("^[0-9a-f]{40}$", execution_head)) {
  stop("invalid MV11-H output or execution head", call. = FALSE)
}
source("R/mv08z_landscape_production.R")
sha <- .mv08z_sha256_file; atomic <- .mv08z_atomic_csv
sources <- data.frame(
  contract_id = "mv11h_private_source_binding_v1",
  source_id = c("historical_gene_partitions", "historical_cell_partitions"),
  private_path_not_published = c(
    "tmp/mv10e-full-clustering-private-v1/mv10e-sample-partitions.csv",
    "tmp/mv11e-cell-benchmark-private-v1/mv11e-sample-partitions.csv"
  ),
  bytes = c(56572844, 15569137),
  sha256 = c(
    "4a3e3719065cdec4e8b78b5ff40c686b5d3d27d77f1918487aaa0ffd0b84c710",
    "ba62457a402143c5ace6191f65fd2e820482421dc0c94f1b3d03d0b74189a02d"
  ),
  independently_closed_by = c("MV10-F", "MV11-F"),
  closure_manifest = c(
    "docs/audits/mv10f-full-clustering-closure-v1/mv10f-artifact-manifest.csv",
    "docs/audits/mv11f-cell-benchmark-closure-v1/mv11f-artifact-manifest.csv"
  ),
  closure_manifest_bytes = c(565, 688),
  closure_manifest_sha256 = c(
    "ba4dcb369bf3bba071afd8d2cf225893221dbf2bae6b872f46e663a9067666bb",
    "99cd32bf08de82fe76435d21020f66b207e35368c47f0fb84c34658c6d96cb31"
  ), stringsAsFactors = FALSE
)
source_paths <- c(gene_path, cell_path)
closure_paths <- sources$closure_manifest
if (!identical(as.numeric(file.info(source_paths)$size), sources$bytes) ||
    !identical(unname(vapply(source_paths, sha, character(1L))),
               sources$sha256) || !all(file.exists(closure_paths)) ||
    !identical(as.numeric(file.info(closure_paths)$size),
               sources$closure_manifest_bytes) ||
    !identical(unname(vapply(closure_paths, sha, character(1L))),
               sources$closure_manifest_sha256)) {
  stop("MV11-H source or closure binding drift", call. = FALSE)
}
implementation_files <- c(
  "R/mv05_benchmark_contract.R", "R/mv05n_clustering_gate.R",
  "R/mv10_clustering_benchmark.R", "R/mv11_cell_benchmark.R",
  "R/mv11g_cross_view_agreement.R",
  "scripts/run_mv11i_cross_view_agreement.R",
  "scripts/build_mv11j_cross_view_closure.R"
)
if (!all(file.exists(implementation_files))) {
  stop("MV11-H implementation stack is incomplete", call. = FALSE)
}
implementation <- data.frame(
  contract_id = "mv11h_implementation_binding_v1",
  file = implementation_files,
  bytes = as.numeric(file.info(implementation_files)$size),
  sha256 = vapply(implementation_files, sha, character(1L)),
  stringsAsFactors = FALSE
)
contract <- data.frame(
  contract_id = "mv11h_cross_view_prefreeze_v1",
  execution_head = execution_head,
  source_gene_rows = 167400L, source_cell_rows = 55800L,
  selected_rows_per_view = 12400L, samples = 124L, seeds = 5L,
  homology_dimensions = "H0;H1_separate", methods = 5L,
  common_k = "2;3", comparison_units = 100L,
  public_seed_rows = 100L, public_summary_rows = 20L,
  metric = "adjusted_rand_index", exact_identity_reported = TRUE,
  elapsed_cap_seconds = 600, public_storage_cap_bytes = 1024^2,
  precommit_real_source_dry_run = TRUE,
  precommit_values_emitted_or_inspected = FALSE,
  execution_authorized_after_commit = TRUE,
  labels_allowed = FALSE, outcomes_allowed = FALSE,
  view_ranking_allowed = FALSE, fusion_allowed = FALSE,
  inference_allowed = FALSE, biological_claims_allowed = FALSE,
  manuscript_claims_allowed = FALSE, stringsAsFactors = FALSE
)
decision <- data.frame(
  contract_id = "mv11h_decision_v1",
  fixed_execution_authorized_after_commit = TRUE,
  labels_authorized = FALSE, outcomes_authorized = FALSE,
  view_ranking_authorized = FALSE, fusion_authorized = FALSE,
  biological_claims_authorized = FALSE,
  manuscript_claims_authorized = FALSE,
  next_action = paste(
    "commit this prefreeze; execute fixed symmetric K=2/K=3 agreement once;",
    "independently repeat both outputs before descriptive review"
  ), stringsAsFactors = FALSE
)
checks <- c(
  two_sources_bound = nrow(sources) == 2L,
  source_hashes_exact = all(vapply(source_paths, sha, character(1L)) ==
                              sources$sha256),
  source_bytes_exact = all(as.numeric(file.info(source_paths)$size) ==
                             sources$bytes),
  two_closures_bound = all(file.exists(closure_paths)),
  closure_hashes_exact = all(vapply(closure_paths, sha, character(1L)) ==
                               sources$closure_manifest_sha256),
  implementation_bound = all(implementation$bytes > 0),
  common_axes_fixed = contract$selected_rows_per_view == 12400L,
  one_hundred_units = contract$comparison_units == 100L,
  one_hundred_seed_rows = contract$public_seed_rows == 100L,
  twenty_summary_rows = contract$public_summary_rows == 20L,
  ARI_only = contract$metric == "adjusted_rand_index",
  exact_identity = contract$exact_identity_reported,
  common_k_only = contract$common_k == "2;3",
  H0_H1_separate = contract$homology_dimensions == "H0;H1_separate",
  structural_dry_run_disclosed = contract$precommit_real_source_dry_run,
  dry_run_values_uninspected =
    !contract$precommit_values_emitted_or_inspected,
  positive_resource_caps = contract$elapsed_cap_seconds > 0 &&
    contract$public_storage_cap_bytes > 0,
  execution_after_commit = contract$execution_authorized_after_commit,
  labels_outcomes_closed = !contract$labels_allowed &&
    !contract$outcomes_allowed,
  ranking_fusion_closed = !contract$view_ranking_allowed &&
    !contract$fusion_allowed,
  inference_claims_closed = !contract$inference_allowed &&
    !contract$biological_claims_allowed &&
    !contract$manuscript_claims_allowed
)
validation <- data.frame(
  contract_id = "mv11h_validation_v1", check_id = names(checks),
  passed = unname(checks), stringsAsFactors = FALSE
)
if (!all(validation$passed)) stop("MV11-H validation failed", call. = FALSE)
dir.create(output, recursive = TRUE)
atomic(contract, file.path(output, "mv11h-contract.csv"))
atomic(sources, file.path(output, "mv11h-source-bindings.csv"))
atomic(implementation, file.path(output, "mv11h-implementation-bindings.csv"))
atomic(decision, file.path(output, "mv11h-decision.csv"))
atomic(validation, file.path(output, "mv11h-validation.csv"))
writeLines(c(
  "# MV11-H common-K cross-view agreement prefreeze", "",
  "This contract freezes symmetric ARI agreement at common K=2 and K=3 for",
  "matched historical cell/gene partitions. It discloses a value-uninspected",
  "100/20-cardinality dry run. It permits no labels, outcomes, view ranking,",
  "fusion, inference, biology, or manuscript claims.", "",
  paste0("Validation: ", sum(validation$passed), "/", nrow(validation),
         " checks pass.")
), file.path(output, "MV11H_CROSS_VIEW_PREFREEZE_2026-08-25.md"))
files <- sort(setdiff(list.files(output), "mv11h-artifact-manifest.csv"))
manifest <- data.frame(
  contract_id = "mv11h_artifact_manifest_v1", artifact = files,
  bytes = as.numeric(file.info(file.path(output, files))$size),
  sha256 = vapply(file.path(output, files), sha, character(1L)),
  stringsAsFactors = FALSE
)
atomic(manifest, file.path(output, "mv11h-artifact-manifest.csv"))
cat("Completed MV11-H prefreeze; checks=", nrow(validation), "\n", sep = "")
