#!/usr/bin/env Rscript

options(warn = 2)
for (package in c("digest", "pkgload", "Matrix")) {
  if (!requireNamespace(package, quietly = TRUE)) stop(package, " required.")
}
pkgload::load_all(".", quiet = TRUE, export_all = TRUE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8L) {
  stop("usage: build_mv07e_full_corpus_prefreeze.R CANDIDATE_CSV RETAINED_CSV ",
    "RECONCILIATION_CSV PANEL_CSV ACCESSION_EVIDENCE_CSV HISTORICAL_HEURISTIC_R ",
    "AUDIT_DIR EXPECTED_HEAD", call. = FALSE)
}
candidate_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
retained_path <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
reconciliation_path <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
panel_path <- normalizePath(args[[4L]], winslash = "/", mustWork = TRUE)
accession_path <- normalizePath(args[[5L]], winslash = "/", mustWork = TRUE)
heuristic_path <- normalizePath(args[[6L]], winslash = "/", mustWork = TRUE)
audit_dir <- args[[7L]]; expected_head <- tolower(trimws(args[[8L]]))
if (!grepl("^[0-9a-f]{40}$", expected_head)) stop("Full EXPECTED_HEAD required.")
head <- tolower(trimws(system2("git", c("rev-parse", "HEAD"), stdout = TRUE)))
if (!identical(head, expected_head)) stop("MV7-E prospective HEAD mismatch.")
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
if (length(list.files(audit_dir, all.files = TRUE, no.. = TRUE))) {
  stop("MV7-E audit directory must be empty.", call. = FALSE)
}
sha <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
write_once <- function(value, name) {
  path <- file.path(audit_dir, name)
  if (file.exists(path)) stop("Refusing to overwrite: ", path, call. = FALSE)
  write_provenance_csv(value, path)
}

candidate <- read.csv(candidate_path, stringsAsFactors = FALSE, check.names = FALSE)
retained <- read.csv(retained_path, stringsAsFactors = FALSE, check.names = FALSE)
reconciliation <- read.csv(reconciliation_path, stringsAsFactors = FALSE,
                           check.names = FALSE)
panel <- read.csv(panel_path, stringsAsFactors = FALSE, check.names = FALSE)
accession <- read.csv(accession_path, stringsAsFactors = FALSE, check.names = FALSE)
canonical <- mv07e_resolve_approaches_v1(reconciliation, retained, accession)
correction <- mv07e_approach_correction_v1(canonical)

candidate$sample_id <- paste(candidate$SRA, candidate[["Sample Name"]], sep = "_")
omitted <- reconciliation$sample_id[
  reconciliation$corpus_class == "single_study_tissue_descriptive_only"]
if (length(omitted) != 34L || anyDuplicated(candidate$sample_id)) {
  stop("MV7-E omitted sample axis is invalid.", call. = FALSE)
}
availability_rows <- lapply(omitted, function(sample_id) {
  source_path <- candidate[["File Path"]][candidate$sample_id == sample_id]
  if (length(source_path) != 1L || !file.exists(source_path) ||
      !grepl("/PH_ClusteringApp/data/VastlyDifferentTissues/", source_path,
             fixed = TRUE)) {
    stop("Candidate sparse source is missing or outside the frozen root: ",
         sample_id, call. = FALSE)
  }
  environment <- new.env(parent = emptyenv())
  loaded <- load(source_path, envir = environment)
  objects <- mget(loaded, envir = environment, inherits = FALSE)
  matrix_index <- which(vapply(objects, function(value) {
    inherits(value, "Matrix") || is.matrix(value)
  }, logical(1L)))
  if (length(matrix_index) != 1L) {
    stop("Expected exactly one matrix in sparse source: ", sample_id,
         call. = FALSE)
  }
  value <- objects[[matrix_index]]
  feature_ids <- rownames(value)
  if (is.null(feature_ids) || anyDuplicated(feature_ids)) {
    stop("Sparse source feature axis is invalid: ", sample_id, call. = FALSE)
  }
  row <- mv07e_panel_availability_row_v1(sample_id, feature_ids, panel)
  row$source_basename <- basename(source_path)
  row$source_bytes <- as.numeric(file.info(source_path)$size)
  row$source_sha256 <- sha(source_path)
  rm(environment, loaded, objects, value); invisible(gc())
  row
})
availability <- do.call(rbind, availability_rows)
availability <- availability[order(availability$sample_id, method = "radix"),
                             , drop = FALSE]
panel_sha <- unique(tolower(panel$panel_sha256))
if (length(panel_sha) != 1L || sha(panel_path) == panel_sha) {
  # The semantic panel hash and file hash are intentionally distinct.
  stop("Panel semantic identity is missing or unexpectedly equals file hash.",
       call. = FALSE)
}
panel_decision <- mv07e_panel_decision_v1(availability, panel_sha)
sample_ids <- canonical$sample_id
axis <- mv07e_sample_seed_axis_v1(sample_ids)
fit_scopes <- mv07e_fit_scopes_v1()
topology <- mv07e_topology_contract_v1()
pairs <- mv07e_pair_scope_v1()
landscape <- mv07e_landscape_contract_v1()
resources <- mv07e_resource_contract_v1()
decision <- mv07e_decision_v1(panel_decision, correction)

lineage <- data.frame(
  contract_id = "mv07e_approach_field_lineage_v1",
  field = c("Approach.x", "Approach.y", "canonical_approach"),
  origin = c("historical filtered Seurat metadata",
    "public candidate metadata joined by sample ID",
    "public candidate metadata with accession confirmation for conflicts"),
  generation = c("expression heuristic using feature count mitochondrial fraction and UMI count",
    "declared source metadata", "deterministic provenance policy"),
  matches_public_124 = c(108L, 124L, 124L), disagreements_public_124 = c(16L, 0L, 0L),
  scientific_status = c("prohibited", "authoritative", "authoritative"),
  stringsAsFactors = FALSE)
firewall <- data.frame(
  contract_id = "mv07e_label_firewall_v1",
  stage = c("mv07f_upstream", "panel_fit", "standardization_pca",
    "typed_view_ph", "landscape_distance"),
  allowed_identity = "sample_id;seed;cache_identity;feature_identity",
  prohibited_inputs = "tissue;study;approach;biological_outcomes;prior_topology_results",
  labels_may_select = FALSE, labels_may_stop = FALSE,
  label_access_state = c("closed", "closed", "closed", "closed", "closed"),
  stringsAsFactors = FALSE)
stage_contract <- data.frame(
  contract_id = "mv07e_stage_authorization_v1",
  stage = c("MV7-F", "post-MV7-F-panel", "MV7-G", "MV7-H", "MV7-I"),
  authorized_now = c(TRUE, TRUE, FALSE, FALSE, FALSE),
  permitted_work = c("34 raw shards and 170 SCT caches only",
    "availability-only 620-cache global-core panel fit and exact panel lock",
    "six-sample five-seed typed-view and PH sentinel after exact panel lock",
    "full topology and landscapes after MV7-G resource gate",
    "descriptive evaluation after immutable distance completion"),
  prerequisite = c("MV7-E independent validation", "MV7-F complete validation",
    "exact 124-panel independent validation", "MV7-G measured resource authorization",
    "MV7-H independent validation and label-access prefreeze"),
  stringsAsFactors = FALSE)

source_paths <- c(candidate_metadata = candidate_path,
  retained_metadata = retained_path, reconciliation = reconciliation_path,
  accepted_panel = panel_path, accession_evidence = accession_path,
  historical_heuristic = heuristic_path,
  helper = "R/mv07e_full_corpus_prefreeze.R",
  builder = "scripts/build_mv07e_full_corpus_prefreeze.R",
  validator = "scripts/validate_mv07e_full_corpus_prefreeze.R",
  repeat_validator = "scripts/validate_mv07e_prefreeze_repeat.R",
  tests = "tests/testthat/test-mv07e-full-corpus-prefreeze.R",
  specification =
    "docs/specifications/MV07E_METADATA_PROVENANCE_AND_FULL124_ESTIMAND_PREFREEZE_V1.md")
if (any(!file.exists(source_paths))) stop("MV7-E source freeze is incomplete.")
private <- names(source_paths) %in% c("retained_metadata", "historical_heuristic")
locators <- unname(source_paths)
locators[names(source_paths) == "retained_metadata"] <-
  "private:historical_joined_metadata_cellcounts.csv"
locators[names(source_paths) == "historical_heuristic"] <-
  "private:historical_Ph_Pipeline_Part_1_V5.R"
source_freeze <- data.frame(
  contract_id = "mv07e_source_freeze_v1", source_id = names(source_paths),
  artifact_locator = locators, sha256 = vapply(source_paths, sha, character(1L)),
  bytes = as.numeric(file.info(source_paths)$size), accepted_head = expected_head,
  private_source = private, stringsAsFactors = FALSE)

criteria <- data.frame(
  contract_id = "mv07e_acceptance_criteria_v1", criterion_order = 1:16,
  criterion_id = c("source_freeze", "approach_conflicts_resolved",
    "public_field_matches_124", "primary_approach_estimability_corrected",
    "panel_availability_complete", "fallback_selected_by_availability",
    "sample_seed_axis", "fit_scopes", "typed_views", "pair_scope",
    "landscape_contract", "resource_contract", "label_firewall",
    "stage_order", "no_ph_or_outcomes", "primary_90_immutable"),
  passed = c(nrow(source_freeze) == 12L,
    nrow(accession) == 16L && all(accession$canonical_approach == accession$public_approach),
    all(canonical$canonical_approach == canonical$public_metadata_approach),
    correction$primary_samples == 90L && correction$primary_mixed_approach_studies == 0L,
    nrow(availability) == 34L && sum(availability$missing_features > 0L) == 1L,
    panel_decision$branch == "fit_deterministic_global_core_over_124",
    nrow(axis) == 620L, nrow(fit_scopes) == 5L && all(fit_scopes$transductive),
    nrow(topology) == 2L, pairs$component_distance_rows == 152520,
    nrow(landscape) == 8L && !any(landscape$changed_by_mv07e),
    nrow(resources) == 5L && all(resources$atomic_write),
    nrow(firewall) == 5L && !any(firewall$labels_may_select),
    identical(stage_contract$authorized_now, c(TRUE, TRUE, FALSE, FALSE, FALSE)),
    !decision$ph_authorized && !decision$outcomes_authorized,
    !decision$primary_90_recalculation_authorized),
  evidence = c("12 exact source identities", "16 official GEO resolutions",
    "Approach.y/public agreement 124 of 124", "90 scRNA and zero mixed studies",
    "33 complete and one one-feature absence", "fixed MV6-C algorithm over 620 caches",
    "124 samples times five seeds", "one global descriptive fit per seed",
    "cell and gene topology frozen", "152520 H0/H1 component rows",
    "all-active-level H0/H1 definition unchanged", "atomic caps and staged projections",
    "labels closed through landscape production", "upstream then panel then sentinel",
    "zero PCA PH landscape clustering outcome jobs", "accepted primary result untouched"),
  stringsAsFactors = FALSE)
if (!all(criteria$passed)) stop("MV7-E acceptance criteria failed: ",
  paste(criteria$criterion_id[!criteria$passed], collapse = ", "), call. = FALSE)

outputs <- list(
  "mv07e-source-freeze.csv" = source_freeze,
  "mv07e-accession-resolution.csv" = accession,
  "mv07e-approach-field-lineage.csv" = lineage,
  "mv07e-canonical-approach.csv" = canonical,
  "mv07e-approach-correction.csv" = correction,
  "mv07e-panel-availability.csv" = availability,
  "mv07e-panel-decision.csv" = panel_decision,
  "mv07e-sample-seed-axis.csv" = axis,
  "mv07e-fit-scopes.csv" = fit_scopes,
  "mv07e-typed-topology.csv" = topology,
  "mv07e-pair-scope.csv" = pairs,
  "mv07e-landscape-contract.csv" = landscape,
  "mv07e-resource-resume-contract.csv" = resources,
  "mv07e-label-firewall.csv" = firewall,
  "mv07e-stage-authorization.csv" = stage_contract,
  "mv07e-acceptance-criteria.csv" = criteria,
  "mv07e-decision.csv" = decision)
for (name in names(outputs)) write_once(outputs[[name]], name)
message("MV7-E prefreeze complete: provenance resolved; 124-panel fallback activated")
