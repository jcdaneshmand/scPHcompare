#!/usr/bin/env Rscript

options(warn = 2, digits = 17, scipen = 999)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required.")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(paste("usage: validate_mv08b_read_only_dataset_admission_audit.R",
    "PRODUCTION_DIR REPEAT_DIR VALIDATION_OUTPUT EXPECTED_HEAD"), call. = FALSE)
}
production_dir <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
repeat_dir <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
validation_output <- args[[3L]]
expected_head <- tolower(trimws(args[[4L]]))
readc <- function(path) utils::read.csv(path, stringsAsFactors = FALSE,
                                         check.names = FALSE)
truth <- function(x) if (is.logical(x)) !is.na(x) & x else
  tolower(trimws(as.character(x))) == "true"
sha <- function(path) digest::digest(file = path, algo = "sha256",
                                     serialize = FALSE)
required <- c("mv08b-source-manifest.csv", "mv08b-gse120221-overlap-map.csv",
  "mv08b-gse120221-source-compatibility.csv",
  "mv08b-gse120221-overlap-summary.csv",
  "mv08b-candidate-admission-audit.csv",
  "mv08b-landscape-continuation-contract.csv", "mv08b-decision.csv",
  "mv08b-artifact-manifest.csv")
if (any(!file.exists(file.path(production_dir, required))) ||
    any(!file.exists(file.path(repeat_dir, required)))) {
  stop("MV8-B production or repeat output is incomplete.", call. = FALSE)
}
sources <- readc(file.path(production_dir, required[[1L]]))
mapping <- readc(file.path(production_dir, required[[2L]]))
compatibility <- readc(file.path(production_dir, required[[3L]]))
summary <- readc(file.path(production_dir, required[[4L]]))
candidates <- readc(file.path(production_dir, required[[5L]]))
landscape <- readc(file.path(production_dir, required[[6L]]))
decision <- readc(file.path(production_dir, required[[7L]]))
manifest <- readc(file.path(production_dir, required[[8L]]))

source_manifest_ok <- nrow(sources) == 8L && !anyDuplicated(sources$source_id) &&
  all(sources$accepted_head == expected_head) &&
  !any(truth(sources$contains_expression_values)) &&
  all(nchar(sources$sha256) == 64L) && all(sources$bytes > 0) &&
  !any(grepl("/mnt/|E:\\\\Repositories", sources$locator, ignore.case = TRUE))
mapping_ok <- nrow(mapping) == 25L && !anyDuplicated(mapping$srs) &&
  !anyDuplicated(mapping$gsm) && all(truth(mapping$exact_accession_overlap)) &&
  !any(truth(mapping$independent_external_library)) &&
  all(mapping$submission == "SRA779509") &&
  all(mapping$sra_study == "SRP162214") &&
  all(mapping$bioproject == "PRJNA492227")
source_ok <- nrow(compatibility) == 25L &&
  all(truth(compatibility$fixed_384_cell_depth_eligible)) &&
  all(compatibility$shared_cell_barcodes > 0L) &&
  all(compatibility$validation_panel_genes == 500L) &&
  all(compatibility$main_panel_genes == 500L) &&
  !any(truth(compatibility$expression_values_exported))
summary_ok <- nrow(summary) == 1L &&
  summary$exact_accession_overlaps == 25L &&
  summary$fixed_depth_eligible_libraries == 25L &&
  summary$libraries_with_shared_barcodes == 25L &&
  !truth(summary$independent_external_validation) &&
  truth(summary$technical_reprocessing_validation)
candidate_ok <- nrow(candidates) == 6L &&
  sum(candidates$disposition == "preferred_pending_download_authorization") == 1L &&
  candidates$disposition[candidates$candidate_id == "existing_gse120221_25"] ==
    "technical_reprocessing_only" &&
  all(candidates$known_accession_overlap[
    candidates$candidate_id == "existing_gse120221_25"] == 25L) &&
  !any(truth(candidates$selection_outcomes_opened)) &&
  !any(truth(candidates$download_authorized)) &&
  !any(truth(candidates$new_ph_authorized)) &&
  all(candidates$observed_ph_peak_rss_bytes > 0) &&
  all(candidates$projected_ph_jobs >= 0)
landscape_expected <- c("comparison_unit_sample", "cell_view_observations_cells",
  "gene_view_observations_fixed_500_genes", "finite_positive_intervals",
  "essential_H0_excluded", "all_consecutive_active_levels", "H0_H1_separate",
  "squared_L2_exact_or_error_controlled", "no_universal_grid",
  "no_universal_level_cap", "streamed_or_chunked_execution")
landscape_ok <- nrow(landscape) == 11L &&
  identical(landscape$item, landscape_expected) &&
  all(truth(landscape$required_state)) && all(truth(landscape$current_state)) &&
  !any(truth(landscape$external_data_may_change))
decision_ok <- nrow(decision) == 1L &&
  decision$gse120221_disposition == "technical_reprocessing_only" &&
  decision$next_authorized_action ==
    "owner_review_of_dataset_and_compute_dossier" &&
  !truth(decision$external_expression_download_authorized) &&
  !truth(decision$new_ph_authorized) &&
  !truth(decision$manuscript_claim_promoted) &&
  !truth(decision$confidential_material_published) &&
  !truth(decision$final_author_roles_decided)
manifest_ok <- nrow(manifest) == 7L && !anyDuplicated(manifest$filename) &&
  all(manifest$accepted_head == expected_head) &&
  all(vapply(seq_len(nrow(manifest)), function(i) {
    path <- file.path(production_dir, manifest$filename[[i]])
    file.exists(path) && identical(tolower(sha(path)),
      tolower(manifest$sha256[[i]])) && file.info(path)$size == manifest$bytes[[i]]
  }, logical(1L))) && !any(truth(manifest$contains_expression_values)) &&
  !any(truth(manifest$contains_private_paths)) &&
  !any(truth(manifest$contains_confidential_review_text))
repeat_ok <- all(vapply(required, function(name)
  identical(sha(file.path(production_dir, name)),
            sha(file.path(repeat_dir, name))), logical(1L)))
public_text <- paste(vapply(file.path(production_dir, required), function(path)
  paste(readLines(path, warn = FALSE), collapse = "\n"), character(1L)),
  collapse = "\n")
privacy_ok <- !grepl("/mnt/|E:\\\\Repositories|pasted-text|Dear Dr Rouchka|reviewer",
  public_text, ignore.case = TRUE) &&
  !grepl("\\.pdf", public_text, ignore.case = TRUE)

checks <- data.frame(
  contract_id = "mv08b_independent_validation_v1",
  check = c("source_manifest", "exact_accession_overlap", "local_source_compatibility",
    "gse120221_summary", "candidate_admission", "landscape_contract",
    "decision_boundary", "artifact_manifest", "byte_identical_repeat",
    "privacy_and_confidentiality"),
  passed = c(source_manifest_ok, mapping_ok, source_ok, summary_ok, candidate_ok, landscape_ok,
    decision_ok, manifest_ok, repeat_ok, privacy_ok),
  detail = c("eight frozen input sources hashed without publishing private paths",
    "25 official GSM/BioSample/SRS/SRR mappings close to SRA779509",
    "25 local source pairs share barcodes and retain the fixed panel/depth",
    "separate processing retained but external independence rejected",
    "six frozen candidates have explicit view and estimand dispositions",
    "all active levels and separate H0/H1 remain immutable",
    "stops before expression download, PH, claims, or author decisions",
    "seven production artifacts independently rehashed",
    "production and repeat bundles are byte-identical",
    "no private path, PDF, expression value, or confidential review text"),
  stringsAsFactors = FALSE)
if (any(!checks$passed)) {
  stop("MV8-B independent validation failed: ",
       paste(checks$check[!checks$passed], collapse = ", "), call. = FALSE)
}
utils::write.table(checks, validation_output, sep = ",", row.names = FALSE,
                   col.names = TRUE, quote = TRUE, na = "", qmethod = "double")
cat("MV8-B independent validation passed:", nrow(checks), "/", nrow(checks),
    "\n")
