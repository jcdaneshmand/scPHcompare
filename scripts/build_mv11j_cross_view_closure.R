#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) stop(paste(
  "usage: build_mv11j_cross_view_closure.R <prefreeze> <gene-partitions>",
  "<cell-partitions> <production> <output> <execution-head>"
), call. = FALSE)
prefreeze <- normalizePath(args[[1L]], mustWork = TRUE)
gene_path <- normalizePath(args[[2L]], mustWork = TRUE)
cell_path <- normalizePath(args[[3L]], mustWork = TRUE)
production <- normalizePath(args[[4L]], mustWork = TRUE)
output <- args[[5L]]
execution_head <- tolower(trimws(args[[6L]]))
if (dir.exists(output)) stop("MV11-J output already exists", call. = FALSE)
source("R/mv05_benchmark_contract.R")
source("R/mv05n_clustering_gate.R")
source("R/mv08z_landscape_production.R")
source("R/mv10_clustering_benchmark.R")
source("R/mv11_cell_benchmark.R")
source("R/mv11g_cross_view_agreement.R")
readc <- .mv08z_read_csv; sha <- .mv08z_sha256_file
atomic <- .mv08z_atomic_csv; truth <- .mv08z_truth
.mv08z_verify_manifest(prefreeze, "mv11h-artifact-manifest.csv")
.mv08z_verify_manifest(production, "mv11i-artifact-manifest.csv")
contract <- readc(file.path(prefreeze, "mv11h-contract.csv"))
sources <- readc(file.path(prefreeze, "mv11h-source-bindings.csv"))
receipt <- readc(file.path(production, "mv11i-terminal-receipt.csv"))
saved_seed <- readc(file.path(production, "mv11i-seed-agreement.csv"))
saved_summary <- readc(file.path(production, "mv11i-agreement-summary.csv"))
gene <- readc(gene_path); cell <- readc(cell_path)
repeat_result <- mv11g_cross_view_agreement_v1(gene, cell)
temp <- tempfile("mv11j-repeat-"); dir.create(temp)
atomic(repeat_result$seed_agreement, file.path(temp, "seed.csv"))
atomic(repeat_result$summary, file.path(temp, "summary.csv"))
repeat_seed <- readc(file.path(temp, "seed.csv"))
repeat_summary <- readc(file.path(temp, "summary.csv"))
checks <- c(
  execution_head = execution_head == contract$execution_head,
  source_hashes = all(vapply(c(gene_path, cell_path), sha, character(1L)) ==
                        sources$sha256),
  terminal_complete = identical(receipt$completion_state, "complete"),
  source_rows = receipt$source_gene_rows == 167400L &&
    receipt$source_cell_rows == 55800L,
  selected_rows = receipt$selected_rows_per_view == 12400L,
  one_hundred_units = receipt$comparison_units == 100L,
  one_hundred_seed_rows = nrow(saved_seed) == 100L,
  twenty_summary_rows = nrow(saved_summary) == 20L,
  finite_ARI = all(is.finite(saved_seed$adjusted_rand)),
  bounded_ARI = all(saved_seed$adjusted_rand >= -1 &
                      saved_seed$adjusted_rand <= 1),
  five_seeds_per_summary = all(saved_summary$seeds == 5L),
  seed_exact_repeat = sha(file.path(temp, "seed.csv")) ==
    sha(file.path(production, "mv11i-seed-agreement.csv")),
  summary_exact_repeat = sha(file.path(temp, "summary.csv")) ==
    sha(file.path(production, "mv11i-agreement-summary.csv")),
  seed_value_identity = identical(saved_seed, repeat_seed),
  summary_value_identity = identical(saved_summary, repeat_summary),
  elapsed_cap = receipt$elapsed_seconds <= receipt$elapsed_cap_seconds,
  no_sample_ids_public = !"sample_id" %in% names(saved_seed) &&
    !"sample_id" %in% names(saved_summary),
  labels_outcomes_closed = !truth(receipt$labels_used) &&
    !truth(receipt$outcomes_used),
  ranking_fusion_closed = !truth(receipt$view_ranking_performed) &&
    !truth(receipt$fusion_performed),
  inference_claims_closed = !truth(receipt$inference_performed) &&
    !truth(receipt$biological_claims) && !truth(receipt$manuscript_claims)
)
validation <- data.frame(
  contract_id = "mv11j_validation_v1", check_id = names(checks),
  passed = unname(checks), stringsAsFactors = FALSE
)
if (!all(validation$passed)) {
  stop("MV11-J closure failed: ",
       paste(validation$check_id[!validation$passed], collapse = ", "),
       call. = FALSE)
}
repeat_artifacts <- data.frame(
  contract_id = "mv11j_artifact_repeat_v1",
  artifact_id = c("seed_agreement", "summary"),
  rows = c(nrow(saved_seed), nrow(saved_summary)),
  saved_sha256 = c(
    sha(file.path(production, "mv11i-seed-agreement.csv")),
    sha(file.path(production, "mv11i-agreement-summary.csv"))
  ),
  repeat_sha256 = c(sha(file.path(temp, "seed.csv")),
                    sha(file.path(temp, "summary.csv"))),
  exact_repeat = TRUE, stringsAsFactors = FALSE
)
decision <- data.frame(
  contract_id = "mv11j_decision_v1",
  common_k_cross_view_agreement_closed = TRUE,
  descriptive_review_eligible_next = TRUE,
  labels_authorized = FALSE, outcomes_authorized = FALSE,
  view_ranking_authorized = FALSE, fusion_authorized = FALSE,
  biological_claims_authorized = FALSE,
  manuscript_claims_authorized = FALSE,
  next_action = paste(
    "commit this closure; prospectively freeze a descriptive review of the",
    "20 symmetric agreement summaries without thresholds or view ranking"
  ), stringsAsFactors = FALSE
)
dir.create(output, recursive = TRUE)
atomic(validation, file.path(output, "mv11j-validation.csv"))
atomic(repeat_artifacts, file.path(output, "mv11j-artifact-repeat.csv"))
atomic(decision, file.path(output, "mv11j-decision.csv"))
writeLines(c(
  "# MV11-J common-K cross-view agreement closure", "",
  "All 100 seed-level and 20 summary rows repeated exactly. The result remains",
  "a symmetric label-closed partition-agreement diagnostic; it does not rank",
  "cell versus gene topology or authorize fusion, biology, or claims.", "",
  paste0("Validation: ", sum(validation$passed), "/", nrow(validation),
         " checks pass.")
), file.path(output, "MV11J_CROSS_VIEW_CLOSURE_2026-08-25.md"))
files <- sort(setdiff(list.files(output), "mv11j-artifact-manifest.csv"))
manifest <- data.frame(
  contract_id = "mv11j_artifact_manifest_v1", artifact = files,
  bytes = as.numeric(file.info(file.path(output, files))$size),
  sha256 = vapply(file.path(output, files), sha, character(1L)),
  stringsAsFactors = FALSE
)
atomic(manifest, file.path(output, "mv11j-artifact-manifest.csv"))
cat("Completed MV11-J cross-view closure; checks=", nrow(validation), "\n",
    sep = "")
