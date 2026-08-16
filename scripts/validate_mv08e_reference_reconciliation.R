#!/usr/bin/env Rscript

options(warn = 2)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required.")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("usage: validate_mv08e_reference_reconciliation.R EVIDENCE_DIR REPEAT_DIR OUTPUT_CSV")
}
evidence_dir <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
repeat_dir <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
output <- args[[3L]]
source("R/provenance_utils.R")
readc <- function(name) utils::read.csv(
  file.path(evidence_dir, name), stringsAsFactors = FALSE, check.names = FALSE)
sha <- function(path) digest::digest(
  file = path, algo = "sha256", serialize = FALSE)
truth <- function(value) tolower(as.character(value)) %in% c("true", "t", "1")

fingerprint <- readc("mv08e-reference-fingerprint.csv")
panel <- readc("mv08e-common475-panel.csv")
missing <- readc("mv08e-missing-biotypes.csv")
biotypes <- readc("mv08e-missing-biotype-summary.csv")
fastq <- readc("mv08e-hca-fastq-resource.csv")
fastq_summary <- readc("mv08e-hca-fastq-resource-summary.csv")
sources <- readc("mv08e-authoritative-sources.csv")
decision <- readc("mv08e-decision.csv")
manifest <- readc("mv08e-artifact-manifest.csv")
spec_path <- "docs/specifications/MV08E_REFERENCE_RECONCILIATION_AND_PANEL_SENSITIVITY_PREFREEZE_V1.md"
spec <- paste(readLines(spec_path, warn = FALSE), collapse = "\n")

expected_files <- sort(c(
  "mv08e-reference-fingerprint.csv", "mv08e-common475-panel.csv",
  "mv08e-missing-biotypes.csv", "mv08e-missing-biotype-summary.csv",
  "mv08e-hca-fastq-resource.csv", "mv08e-hca-fastq-resource-summary.csv",
  "mv08e-authoritative-sources.csv", "mv08e-decision.csv"))
manifest_paths <- file.path(evidence_dir, manifest$file)
repeat_files <- c(expected_files, "mv08e-artifact-manifest.csv")
repeat_paths <- file.path(repeat_dir, repeat_files)
public_paths <- file.path(evidence_dir, repeat_files)
repeat_validation <- data.frame(
  contract_id = "mv08e_reference_reconciliation_repeat_v1",
  file = repeat_files,
  production_bytes = as.numeric(file.info(public_paths)$size),
  repeat_bytes = as.numeric(file.info(repeat_paths)$size),
  production_sha256 = vapply(public_paths, sha, character(1L)),
  repeat_sha256 = vapply(repeat_paths, sha, character(1L)),
  byte_identical = vapply(seq_along(public_paths), function(index) {
    identical(readBin(public_paths[[index]], "raw", file.info(public_paths[[index]])$size),
              readBin(repeat_paths[[index]], "raw", file.info(repeat_paths[[index]])$size))
  }, logical(1L)), stringsAsFactors = FALSE)
write_provenance_csv(repeat_validation,
  file.path(dirname(output), "mv08e-repeat-validation.csv"))

add <- function(category, passed, detail) data.frame(
  contract_id = "mv08e_reference_reconciliation_independent_validation_v1",
  category = category, passed = isTRUE(passed), detail = detail,
  stringsAsFactors = FALSE)
checks <- list(
  add("exact_reference_fingerprint",
      nrow(fingerprint) == 1L && fingerprint$unfiltered_ensembl93_genes == 58395L &&
        fingerprint$filtered_reference_genes == 33538L &&
        fingerprint$hca_h5_feature_genes == 33538L &&
        truth(fingerprint$id_axis_byte_exact) && truth(fingerprint$name_axis_byte_exact) &&
        fingerprint$filtered_gtf_id_axis_sha256 == fingerprint$hca_h5_id_axis_sha256,
      "Ensembl93 -> Cell Ranger 3.0.0 -> HCA ID/name axes exact"),
  add("common475_axis",
      nrow(panel) == 475L && identical(as.integer(panel$panel_order_475), 1:475) &&
        all(diff(panel$panel_order_500) > 0L) &&
        length(unique(panel$common475_axis_sha256)) == 1L &&
        all(!truth(panel$hca_expression_used_for_selection)) &&
        all(!truth(panel$biological_outcomes_used_for_selection)),
      "ordered 475 subset fixed without HCA expression or outcomes"),
  add("missing25_annotation",
      nrow(missing) == 25L && length(unique(missing$ensembl_stable_id)) == 25L &&
        all(missing$cellranger_3_0_0_filter_status ==
          "excluded_by_documented_biotype_filter") &&
        all(truth(missing$current_ensembl_lookup_present)) &&
        all(!truth(missing$crosswalk_can_restore_matrix_count)),
      "25 current stable IDs excluded by historical biotype filter"),
  add("biotype_accountability",
      sum(biotypes$genes) == 25L && nrow(biotypes) == 5L &&
        sum(missing$ensembl93_gene_biotype == "processed_transcript") == 4L &&
        sum(grepl("pseudogene", missing$ensembl93_gene_biotype, fixed = TRUE)) == 21L,
      "4 processed transcripts plus 21 pseudogene-related annotations"),
  add("panel_accountability",
      fingerprint$panel_genes == 500L &&
        fingerprint$panel_in_unfiltered_ensembl93 == 500L &&
        fingerprint$panel_in_filtered_cellranger_reference == 475L &&
        fingerprint$panel_excluded_by_reference_filter == 25L,
      "500 unfiltered; 475 quantified; 25 filtered"),
  add("raw_read_inventory",
      nrow(fastq) == 8L && identical(as.integer(fastq$unit_order), 1:8) &&
        sum(fastq$fastq_files) == 48L && sum(fastq$fastq_bytes) == 85034239918 &&
        all(!truth(fastq$raw_reads_downloaded)) &&
        fastq_summary$reprocessing_status ==
          "deferred_until_500_vs_475_sensitivity_decision",
      "8 units; 48 FASTQs; 85,034,239,918 bytes; not downloaded"),
  add("authoritative_sources",
      nrow(sources) == 4L && identical(as.integer(sources$source_order), 1:4) &&
        all(grepl("^https://", sources$url)) &&
        setequal(sources$source_id, c("ensembl_release_93_gtf",
          "cellranger_3_0_0_reference_build_steps",
          "hca_bone_marrow_publication", "ensembl_rest_lookup")),
      "Ensembl, 10x, HCA publication, and current Ensembl lookup recorded"),
  add("decision_boundary",
      nrow(decision) == 1L && decision$decision ==
        "exact_reference_identified_crosswalk_cannot_restore_counts" &&
        decision$next_authorized_stage ==
          "reference_only_500_vs_475_sensitivity_prefreeze" &&
        !truth(decision$hca_fastq_download_authorized_now) &&
        truth(decision$raw_read_reprocessing_preference_recorded),
      "sensitivity next; raw-read preference retained; no FASTQ authorization"),
  add("landscape_contract",
      all(vapply(c("every consecutive active landscape level",
        "exact or error-controlled squared-L2", "no primary fixed grid",
        "H0 and H1 calculated and reported separately"),
        function(value) grepl(value, spec, fixed = TRUE), logical(1L))),
      "all levels, exact/error-controlled, grid-free, H0/H1 separate"),
  add("artifact_manifest",
      nrow(manifest) == 8L && !anyDuplicated(manifest$file) &&
        setequal(manifest$file, expected_files) &&
        all(file.exists(manifest_paths)) &&
        identical(unname(vapply(manifest_paths, sha, character(1L))),
                  unname(manifest$sha256)) &&
        all(as.numeric(file.info(manifest_paths)$size) ==
              as.numeric(manifest$bytes)),
      "8 production artifacts have exact public hashes and sizes"),
  add("repeat_identity", all(file.exists(repeat_paths)) &&
        all(repeat_validation$byte_identical),
      "9 of 9 production artifacts repeat byte-for-byte"),
  add("privacy_and_stop_boundary",
      all(!truth(manifest$contains_expression)) &&
        all(!truth(manifest$contains_cell_barcode)) &&
        all(!truth(manifest$contains_local_absolute_path)) &&
        !any(grepl("[A-Za-z]:\\\\|/mnt/", vapply(public_paths, function(path) {
          paste(readLines(path, warn = FALSE), collapse = "\n")
        }, character(1L)))) &&
        !truth(decision$biological_outcomes_computed),
      "no expression/barcode/private path; outcomes remain closed")
)
validation <- do.call(rbind, checks)
if (any(!validation$passed)) {
  stop("MV8-E validation failed: ",
       paste(validation$category[!validation$passed], collapse = ", "))
}
write_provenance_csv(validation, output)
message("MV8-E independent validation passed ", nrow(validation), "/",
        nrow(validation), ".")
