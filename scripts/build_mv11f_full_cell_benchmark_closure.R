#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) stop(paste(
  "usage: build_mv11f_full_cell_benchmark_closure.R <prefreeze>",
  "<matrix-bundle> <full-private> <full-public> <output> <execution-head>"
), call. = FALSE)
prefreeze <- normalizePath(args[[1L]], mustWork = TRUE)
bundle_path <- normalizePath(args[[2L]], mustWork = TRUE)
private <- normalizePath(args[[3L]], mustWork = TRUE)
public <- normalizePath(args[[4L]], mustWork = TRUE)
output <- args[[5L]]
execution_head <- tolower(trimws(args[[6L]]))
if (dir.exists(output)) stop("MV11-F output already exists", call. = FALSE)

source("R/mv05_benchmark_contract.R")
source("R/mv05n_clustering_gate.R")
source("R/mv08z_landscape_production.R")
source("R/mv10_clustering_benchmark.R")
source("R/mv11_cell_benchmark.R")
readc <- .mv08z_read_csv; sha <- .mv08z_sha256_file
atomic <- .mv08z_atomic_csv; truth <- .mv08z_truth
.mv08z_verify_manifest(prefreeze, "mv11b-artifact-manifest.csv")
.mv08z_verify_manifest(public, "mv11e-artifact-manifest.csv")
contract <- readc(file.path(prefreeze, "mv11b-contract.csv"))
queue <- readc(file.path(prefreeze, "mv11b-workload-queue.csv"))
source_binding <- readc(file.path(prefreeze, "mv11b-source-binding.csv"))
receipt <- readc(file.path(public, "mv11e-terminal-receipt.csv"))
ledger <- readc(file.path(public, "mv11e-resource-ledger.csv"))
saved <- list(
  assignments = readc(file.path(private, "mv11e-sample-partitions.csv")),
  quality = readc(file.path(public, "mv11e-partition-quality.csv")),
  stability = readc(file.path(public, "mv11e-seed-stability.csv")),
  primary_k = readc(file.path(public, "mv11e-primary-k-selection.csv")),
  agreement = readc(file.path(public, "mv11e-method-agreement.csv"))
)
bundle <- readRDS(bundle_path)
partition_rows <- list(); quality_rows <- list()
for (i in seq_len(nrow(queue))) {
  binding <- queue[i, , drop = FALSE]
  matrix <- mv11_cell_matrix_v1(bundle, binding$seed,
                                binding$homology_dimension)
  grid <- mv10_partition_grid_v1(matrix)
  metadata <- data.frame(
    execution_head = execution_head,
    catalog_order = binding$catalog_order,
    stack_id = binding$stack_id,
    representation_id = binding$representation_id,
    panel_id = binding$panel_id,
    seed = as.integer(binding$seed),
    homology_dimension = binding$homology_dimension,
    source_distances_sha256 = binding$source_distances_sha256,
    stringsAsFactors = FALSE
  )
  partition_rows[[i]] <- cbind(
    metadata[rep(1L, nrow(grid$partitions)), , drop = FALSE], grid$partitions
  )
  quality_rows[[i]] <- cbind(
    metadata[rep(1L, nrow(grid$quality)), , drop = FALSE], grid$quality
  )
}
repeat_values <- list(
  assignments = do.call(rbind, partition_rows),
  quality = do.call(rbind, quality_rows)
)
rownames(repeat_values$assignments) <- NULL
rownames(repeat_values$quality) <- NULL
repeat_values$stability <- mv10_seed_stability_v1(repeat_values$assignments)
repeat_values$primary_k <- mv11_select_primary_k_v1(repeat_values$assignments)
repeat_values$agreement <- mv10_method_agreement_v1(repeat_values$assignments)
temp <- tempfile("mv11f-repeat-"); dir.create(temp)
file_names <- c(
  assignments = "mv11e-sample-partitions.csv",
  quality = "mv11e-partition-quality.csv",
  stability = "mv11e-seed-stability.csv",
  primary_k = "mv11e-primary-k-selection.csv",
  agreement = "mv11e-method-agreement.csv"
)
for (name in names(file_names)) {
  atomic(repeat_values[[name]], file.path(temp, file_names[[name]]))
}
saved_paths <- c(
  assignments = file.path(private, file_names[["assignments"]]),
  setNames(file.path(public, file_names[names(file_names) != "assignments"]),
           names(file_names)[names(file_names) != "assignments"])
)
repeat_paths <- setNames(file.path(temp, file_names), names(file_names))
exact <- vapply(names(file_names), function(name) {
  sha(saved_paths[[name]]) == sha(repeat_paths[[name]])
}, logical(1L))
checks <- c(
  prefreeze_head_exact = execution_head == contract$execution_head,
  source_hash_exact = sha(bundle_path) == source_binding$sha256,
  source_bytes_exact = as.numeric(file.info(bundle_path)$size) ==
    source_binding$bytes,
  terminal_complete = identical(receipt$completion_state, "complete"),
  ten_matrices = receipt$matrices == 10L,
  four_hundred_fifty_fits = receipt$partition_fits == 450L,
  fifty_five_thousand_eight_hundred_assignments =
    nrow(saved$assignments) == 55800L,
  four_hundred_fifty_quality = nrow(saved$quality) == 450L,
  ninety_stability = nrow(saved$stability) == 90L,
  two_primary_k = nrow(saved$primary_k) == 2L,
  nine_hundred_agreement = nrow(saved$agreement) == 900L,
  ten_ledger_rows = nrow(ledger) == 10L,
  all_children_complete = all(ledger$disposition == "completed"),
  stderr_empty = receipt$stderr_bytes == 0,
  one_worker = receipt$workers == 1L && all(ledger$workers == 1L),
  zero_retries = receipt$retries == 0L && all(ledger$retries == 0L),
  elapsed_caps = all(ledger$elapsed_seconds <= ledger$elapsed_cap_seconds),
  rss_caps = all(ledger$peak_process_tree_rss_bytes <=
                   ledger$process_tree_rss_cap_bytes),
  aggregate_elapsed_cap = receipt$elapsed_seconds <=
    receipt$aggregate_elapsed_cap_seconds,
  private_storage_cap = receipt$private_bytes <=
    receipt$private_storage_cap_bytes,
  assignments_exact_repeat = exact[["assignments"]],
  quality_exact_repeat = exact[["quality"]],
  stability_exact_repeat = exact[["stability"]],
  primary_k_exact_repeat = exact[["primary_k"]],
  agreement_exact_repeat = exact[["agreement"]],
  H0_H1_separate = !truth(receipt$H0_H1_combined),
  cell_gene_separate = !truth(receipt$cell_gene_combined),
  labels_outcomes_closed = !truth(receipt$labels_used) &&
    !truth(receipt$outcomes_used),
  cross_view_closed = !truth(receipt$cross_view_comparison_performed),
  claims_closed = !truth(receipt$biological_claims) &&
    !truth(receipt$manuscript_claims)
)
validation <- data.frame(
  contract_id = "mv11f_validation_v1", check_id = names(checks),
  passed = unname(checks), stringsAsFactors = FALSE
)
if (!all(validation$passed)) {
  stop("MV11-F closure failed: ",
       paste(validation$check_id[!validation$passed], collapse = ", "),
       call. = FALSE)
}
artifact_repeat <- data.frame(
  contract_id = "mv11f_artifact_repeat_v1",
  artifact_id = names(file_names),
  rows = vapply(saved, nrow, integer(1L)),
  saved_sha256 = vapply(saved_paths, sha, character(1L)),
  repeat_sha256 = vapply(repeat_paths, sha, character(1L)),
  exact_repeat = unname(exact), stringsAsFactors = FALSE
)
summary <- data.frame(
  contract_id = "mv11f_cell_benchmark_closure_summary_v1",
  execution_head = execution_head, matrices = 10L, partition_fits = 450L,
  private_assignment_rows = 55800L, quality_rows = 450L,
  stability_rows = 90L, primary_k_rows = 2L, agreement_rows = 900L,
  selected_H0_k = saved$primary_k$selected_k[
    saved$primary_k$homology_dimension == "H0"
  ],
  selected_H1_k = saved$primary_k$selected_k[
    saved$primary_k$homology_dimension == "H1"
  ],
  exact_repeat = TRUE, labels_used = FALSE, outcomes_used = FALSE,
  cross_view_comparison_performed = FALSE,
  result_interpretation_state = "label_closed_cell_benchmark_only",
  stringsAsFactors = FALSE
)
decision <- data.frame(
  contract_id = "mv11f_decision_v1",
  historical_cell_benchmark_closed = TRUE,
  common_k_cross_view_prefreeze_eligible_next = TRUE,
  labels_authorized = FALSE, outcomes_authorized = FALSE,
  cross_view_comparison_authorized_now = FALSE,
  fusion_authorized = FALSE, biological_claims_authorized = FALSE,
  manuscript_claims_authorized = FALSE,
  next_action = paste(
    "commit this closure; prospectively freeze common K=2 and K=3 historical",
    "cell-versus-gene partition agreement without labels or view ranking"
  ), stringsAsFactors = FALSE
)
dir.create(output, recursive = TRUE)
atomic(validation, file.path(output, "mv11f-validation.csv"))
atomic(artifact_repeat, file.path(output, "mv11f-artifact-repeat.csv"))
atomic(summary, file.path(output, "mv11f-closure-summary.csv"))
atomic(decision, file.path(output, "mv11f-decision.csv"))
readme <- c(
  "# MV11-F matched historical cell benchmark closure", "",
  "All ten H0/H1 cell matrices completed the frozen five-method K=2:10 grid.",
  "All five aggregate artifacts independently repeated byte-for-byte. Labels,",
  "outcomes, cross-view comparison, fusion, inference, and claims remained",
  "closed. The selected cell-native PAM K values are recorded as label-closed",
  "method diagnostics, not as biological conclusions.", "",
  paste0("Validation: ", sum(validation$passed), "/", nrow(validation),
         " checks pass.")
)
writeLines(readme,
           file.path(output, "MV11F_CELL_BENCHMARK_CLOSURE_2026-08-25.md"))
files <- sort(setdiff(list.files(output), "mv11f-artifact-manifest.csv"))
manifest <- data.frame(
  contract_id = "mv11f_artifact_manifest_v1", artifact = files,
  bytes = as.numeric(file.info(file.path(output, files))$size),
  sha256 = vapply(file.path(output, files), sha, character(1L)),
  stringsAsFactors = FALSE
)
atomic(manifest, file.path(output, "mv11f-artifact-manifest.csv"))
cat("Completed MV11-F full cell benchmark closure; checks=",
    nrow(validation), "\n", sep = "")
