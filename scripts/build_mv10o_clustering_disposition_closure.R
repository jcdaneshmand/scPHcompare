#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) stop(paste(
  "usage: build_mv10o_clustering_disposition_closure.R <prefreeze>",
  "<mv10h-synthesis> <mv10n-output> <output-dir> <execution-head>"
), call. = FALSE)
prefreeze <- normalizePath(args[[1L]], mustWork = TRUE)
synthesis <- normalizePath(args[[2L]], mustWork = TRUE)
production <- normalizePath(args[[3L]], mustWork = TRUE)
output <- args[[4L]]
execution_head <- tolower(trimws(args[[5L]]))
if (dir.exists(output) && length(list.files(output, all.files = TRUE,
                                            no.. = TRUE))) {
  stop("refusing to overwrite MV10-O closure")
}
if (!dir.exists(output)) dir.create(output, recursive = TRUE)
source("R/mv08z_landscape_production.R")
source("R/mv10m_clustering_disposition.R")
readc <- .mv08z_read_csv; sha <- .mv08z_sha256_file; atomic <- .mv08z_atomic_csv
.mv08z_verify_manifest(prefreeze, "mv10m-artifact-manifest.csv")
.mv08z_verify_manifest(production, "mv10n-artifact-manifest.csv")
contract <- readc(file.path(prefreeze, "mv10m-contract.csv"))
schemas <- readc(file.path(prefreeze, "mv10m-source-schema.csv"))
receipt <- readc(file.path(production, "mv10n-terminal-receipt.csv"))
if (execution_head != contract$execution_head ||
    receipt$completion_state != "complete") stop("MV10-O prerequisite drift")
values <- lapply(schemas$file, function(name) readc(file.path(synthesis, name)))
names(values) <- schemas$source_id
repeat_result <- do.call(mv10m_build_methodological_disposition_v1, values)
names_map <- c(
  selected_primary_seed_metrics = "mv10n-selected-primary-seed-metrics.csv",
  primary_summary = "mv10n-primary-summary.csv",
  representation_context = "mv10n-representation-context.csv",
  method_sensitivity = "mv10n-method-sensitivity.csv",
  disposition = "mv10n-disposition.csv"
)
rehash <- do.call(rbind, lapply(seq_along(names_map), function(index) {
  id <- names(names_map)[[index]]
  observed <- readc(file.path(production, names_map[[index]]))
  expected <- repeat_result[[id]]
  data.frame(
    contract_id = "mv10o_rehash_v1", output_order = index,
    artifact_id = id, rows = nrow(observed),
    expected_rows = nrow(expected),
    exact_column_names = identical(names(observed), names(expected)),
    numeric_repeat = isTRUE(all.equal(observed, expected, tolerance = 1e-14,
                                      check.attributes = FALSE)),
    stringsAsFactors = FALSE
  )
}))
summary <- repeat_result$primary_summary
disposition <- repeat_result$disposition
validation <- data.frame(
  contract_id = "mv10o_validation_v1",
  check_id = c("prefreeze_manifest", "production_manifest", "execution_head",
               "terminal_complete", "five_outputs", "ten_seed_rows",
               "two_summary_rows", "six_representation_rows",
               "twenty_four_sensitivity_rows", "one_disposition_row",
               "five_independent_rehashes", "all_numeric_repeats",
               "H0_H1_separate", "frozen_PAM", "selected_K_2_3",
               "structural_rule_applied", "stability_descriptive",
               "silhouette_descriptive", "representation_no_ranking",
               "method_no_ranking", "internal_only", "labels_closed",
               "outcomes_closed", "biology_closed", "claims_closed"),
  passed = c(TRUE, TRUE, execution_head == contract$execution_head,
             receipt$completion_state == "complete", receipt$output_tables == 5L,
             receipt$selected_seed_rows == 10L, receipt$primary_summary_rows == 2L,
             receipt$representation_context_rows == 6L,
             receipt$method_sensitivity_rows == 24L,
             receipt$disposition_rows == 1L, nrow(rehash) == 5L,
             all(rehash$numeric_repeat & rehash$exact_column_names),
             disposition$H0_H1_remain_separate,
             disposition$primary_method == "pam_dissimilarity_v1",
             disposition$selected_H0_k == 2L && disposition$selected_H1_k == 3L,
             all(summary$structurally_nondegenerate) ==
               disposition$all_selected_partitions_structurally_nondegenerate,
             grepl("descriptive", disposition$stability_interpretation),
             grepl("descriptive", disposition$silhouette_interpretation),
             disposition$representation_interpretation == "sensitivity_no_ranking",
             disposition$method_interpretation == "sensitivity_no_ranking",
             disposition$internal_only, !disposition$labels_used,
             !disposition$outcomes_used, !disposition$biological_interpretation,
             !disposition$manuscript_claims),
  stringsAsFactors = FALSE
)
if (!all(validation$passed)) stop("MV10-O validation failed")
decision <- data.frame(
  contract_id = "mv10o_decision_v1",
  decision = disposition$decision,
  selected_H0_k = disposition$selected_H0_k,
  selected_H1_k = disposition$selected_H1_k,
  structurally_nondegenerate =
    disposition$all_selected_partitions_structurally_nondegenerate,
  next_stage = if (disposition$all_selected_partitions_structurally_nondegenerate) {
    "separate_label_evaluation_prefreeze"
  } else "stop_before_labels",
  labels_outcomes_state = "closed", biological_interpretation_state = "closed",
  manuscript_claims_state = "closed", stringsAsFactors = FALSE
)
artifacts <- list(
  "mv10o-rehash.csv" = rehash, "mv10o-validation.csv" = validation,
  "mv10o-decision.csv" = decision
)
for (name in names(artifacts)) atomic(artifacts[[name]], file.path(output, name))
writeLines(c(
  "# MV10-O clustering methodological-disposition closure", "",
  "The five MV10-N outputs independently reproduce and all 25 gates pass.",
  "The disposition remains internal, descriptive, and label closed."
), file.path(output, "MV10O_CLUSTERING_DISPOSITION_CLOSURE_2026-08-25.md"))
files <- sort(setdiff(list.files(output), "mv10o-artifact-manifest.csv"))
manifest <- data.frame(
  contract_id = "mv10o_artifact_manifest_v1", artifact = files,
  bytes = as.numeric(file.info(file.path(output, files))$size),
  sha256 = vapply(file.path(output, files), sha, character(1L)),
  stringsAsFactors = FALSE
)
atomic(manifest, file.path(output, "mv10o-artifact-manifest.csv"))
message("Closed MV10-O clustering disposition; checks=25")
