#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop(paste(
  "usage: build_mv09c_robustness_synthesis_closure.R <prefreeze>",
  "<mv08zy-public-summary> <mv09b-output> <closure-output>"
), call. = FALSE)
prefreeze <- normalizePath(args[[1L]], mustWork = TRUE)
source_path <- normalizePath(args[[2L]], mustWork = TRUE)
production <- normalizePath(args[[3L]], mustWork = TRUE)
output <- args[[4L]]
if (dir.exists(output) && length(list.files(output, all.files = TRUE,
                                            no.. = TRUE))) {
  stop("refusing to overwrite MV9-C closure")
}
if (!dir.exists(output)) dir.create(output, recursive = TRUE)
source("R/mv08z_landscape_production.R")
source("R/mv09_robustness_synthesis.R")
.mv08z_verify_manifest(prefreeze, "mv09a-artifact-manifest.csv")
.mv08z_verify_manifest(production, "mv09b-artifact-manifest.csv")
readc <- .mv08z_read_csv; sha <- .mv08z_sha256_file; atomic <- .mv08z_atomic_csv
contract <- readc(file.path(prefreeze, "mv09a-contract.csv"))
receipt <- readc(file.path(production, "mv09b-terminal-receipt.csv"))
if (sha(source_path) != contract$source_summary_sha256 ||
    receipt$completion_state != "complete") stop("MV9-C source drift")
fresh <- mv09_build_robustness_synthesis_v1(readc(source_path))
mapping <- c(
  aggregate = "mv09b-aggregate-comparisons.csv",
  plot_data = "mv09b-plot-data.csv",
  internal_seed_summary = "mv09b-internal-seed-summary.csv",
  external_singleton = "mv09b-external-singleton.csv",
  dimension_delta = "mv09b-dimension-delta.csv",
  dimension_delta_summary = "mv09b-dimension-delta-summary.csv"
)
rehash <- lapply(seq_along(mapping), function(i) {
  name <- names(mapping)[[i]]
  path <- file.path(production, mapping[[i]])
  saved <- readc(path)
  expected <- fresh[[name]]
  if (!isTRUE(all.equal(saved, expected, tolerance = 1e-14,
                        check.attributes = FALSE))) {
    stop("MV9-C reconstruction drift: ", name)
  }
  data.frame(
    contract_id = "mv09c_rehash_v1", artifact_order = i,
    artifact = mapping[[i]], rows = nrow(saved),
    bytes = as.numeric(file.info(path)$size), sha256 = sha(path),
    independently_recomputed = TRUE, stringsAsFactors = FALSE
  )
})
rehash <- do.call(rbind, rehash)
validation <- data.frame(
  contract_id = "mv09c_validation_v1",
  check_id = c("prefreeze_manifest", "production_manifest", "source_hash",
               "terminal_complete", "aggregate_40", "plot_440",
               "internal_summary_66", "external_singleton_110",
               "dimension_delta_220", "dimension_summary_88",
               "independent_recomputation", "finite_metrics",
               "no_inference", "label_outcome_firewall",
               "clustering_fusion_claim_firewall"),
  passed = c(TRUE, TRUE, sha(source_path) == contract$source_summary_sha256,
             receipt$completion_state == "complete",
             nrow(fresh$aggregate) == 40L, nrow(fresh$plot_data) == 440L,
             nrow(fresh$internal_seed_summary) == 66L,
             nrow(fresh$external_singleton) == 110L,
             nrow(fresh$dimension_delta) == 220L,
             nrow(fresh$dimension_delta_summary) == 88L,
             all(rehash$independently_recomputed),
             all(is.finite(fresh$plot_data$value)),
             !receipt$inference_performed,
             !receipt$labels_used && !receipt$outcomes_used,
             receipt$clustering_jobs == 0L && receipt$fusion_jobs == 0L &&
               !receipt$manuscript_claims),
  stringsAsFactors = FALSE
)
if (!all(validation$passed)) stop("MV9-C validation failed")
decision <- data.frame(
  contract_id = "mv09c_decision_v1",
  decision = "close_descriptive_label_closed_robustness_synthesis",
  interpretation = "sensitivity_structure_not_equivalence_or_biology",
  next_stage = "separate_scientific_review_and_figure_prefreeze",
  labels_outcomes_state = "closed", clustering_state = "closed",
  fusion_state = "closed", manuscript_claims_state = "closed",
  stringsAsFactors = FALSE
)
atomic(rehash, file.path(output, "mv09c-production-rehash.csv"))
atomic(validation, file.path(output, "mv09c-validation.csv"))
atomic(decision, file.path(output, "mv09c-decision.csv"))
writeLines(c(
  "# MV9-C aggregate robustness-synthesis closure", "",
  "All six public aggregate tables were independently reconstructed from the",
  "immutable 40-row comparison summary. Results remain descriptive and",
  "label-closed; external rows are explicitly singleton-cohort evidence."
), file.path(output, "MV09C_ROBUSTNESS_SYNTHESIS_CLOSURE_2026-08-25.md"))
files <- sort(setdiff(list.files(output), "mv09c-artifact-manifest.csv"))
manifest <- data.frame(
  contract_id = "mv09c_artifact_manifest_v1", artifact = files,
  bytes = as.numeric(file.info(file.path(output, files))$size),
  sha256 = vapply(file.path(output, files), sha, character(1L)),
  stringsAsFactors = FALSE
)
atomic(manifest, file.path(output, "mv09c-artifact-manifest.csv"))
message("Closed MV9-C robustness synthesis; checks=15")
