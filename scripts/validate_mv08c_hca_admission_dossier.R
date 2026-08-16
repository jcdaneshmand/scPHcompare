#!/usr/bin/env Rscript

options(warn = 2, digits = 17, scipen = 999)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required.")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(paste("usage: validate_mv08c_hca_admission_dossier.R",
    "PRODUCTION_DIR REPEAT_DIR VALIDATION_OUTPUT EXPECTED_HEAD"),
    call. = FALSE)
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

required <- c("mv08c-source-manifest.csv", "mv08c-project-summary.csv",
  "mv08c-component-accountability.csv", "mv08c-primary-file-manifest.csv",
  "mv08c-primary-unit-ledger.csv", "mv08c-independence-audit.csv",
  "mv08c-access-terms.csv", "mv08c-view-gates.csv",
  "mv08c-resource-plan.csv", "mv08c-resource-summary.csv",
  "mv08c-execution-plan.csv", "mv08c-stop-resume.csv",
  "mv08c-landscape-contract.csv", "mv08c-decision.csv",
  "mv08c-artifact-manifest.csv")
if (any(!file.exists(file.path(production_dir, required))) ||
    any(!file.exists(file.path(repeat_dir, required)))) {
  stop("MV8-C production or repeat dossier is incomplete.", call. = FALSE)
}

sources <- readc(file.path(production_dir, required[[1L]]))
project <- readc(file.path(production_dir, required[[2L]]))
components <- readc(file.path(production_dir, required[[3L]]))
files <- readc(file.path(production_dir, required[[4L]]))
units <- readc(file.path(production_dir, required[[5L]]))
independence <- readc(file.path(production_dir, required[[6L]]))
terms <- readc(file.path(production_dir, required[[7L]]))
views <- readc(file.path(production_dir, required[[8L]]))
resource_plan <- readc(file.path(production_dir, required[[9L]]))
resources <- readc(file.path(production_dir, required[[10L]]))
execution <- readc(file.path(production_dir, required[[11L]]))
stop_resume <- readc(file.path(production_dir, required[[12L]]))
landscape <- readc(file.path(production_dir, required[[13L]]))
decision <- readc(file.path(production_dir, required[[14L]]))
artifact_manifest <- readc(file.path(production_dir, required[[15L]]))

source_ok <- nrow(sources) == 10L && !anyDuplicated(sources$source_id) &&
  all(sources$accepted_head == expected_head) &&
  all(grepl("^[0-9a-f]{64}$", sources$sha256)) && all(sources$bytes > 0) &&
  !any(truth(sources$contains_expression_values)) &&
  sum(sources$access_class == "private_hash_only") == 1L &&
  !any(grepl("/mnt/|[A-Z]:\\\\", sources$locator))
project_ok <- nrow(project) == 1L &&
  project$project_id == "cc95ff89-2e68-4a08-a234-480eca21ce79" &&
  project$catalog == "dcp60" && project$api_version == "19.0" &&
  project$project_label == "HematopoieticImmuneCellAtlas" &&
  project$project_cell_estimate == 1453784L &&
  project$project_donor_count == 28L && project$data_use_restriction == "NRES" &&
  truth(project$accessible) && !truth(project$expression_values_opened)
component_ok <- nrow(components) == 4L &&
  identical(as.integer(components$entity_count), c(8L, 16L, 63L, 1L)) &&
  identical(truth(components$primary_eligible), c(TRUE, FALSE, FALSE, FALSE)) &&
  !any(truth(components$outcome_or_topology_used))
file_ok <- nrow(files) == 8L && !anyDuplicated(files$file_uuid) &&
  !anyDuplicated(files$file_name) && all(files$file_format == "h5") &&
  all(grepl("^[0-9a-f]{64}$", files$file_sha256)) &&
  sum(as.numeric(files$file_size_bytes)) == 202770089 &&
  all(grepl("^drs://", files$drs_uri)) &&
  all(grepl("^https://service[.]azul[.]data[.]humancellatlas[.]org/repository/files/",
            files$azul_download_url)) &&
  !any(grepl("X-Amz-|Signature=|Expires=", files$azul_download_url,
             ignore.case = TRUE)) &&
  !any(truth(files$expression_download_authorized))
unit_ok <- nrow(units) == 8L && !anyDuplicated(units$unit_id) &&
  !anyDuplicated(units$donor_document_sha256) &&
  !anyDuplicated(units$specimen_document_sha256) &&
  !anyDuplicated(units$suspension_document_sha256) &&
  !anyDuplicated(units$sequencing_process_sha256) &&
  all(grepl("^[0-9a-f]{64}$", units$donor_document_sha256)) &&
  all(truth(units$sequencing_input_matches_suspension)) &&
  all(truth(units$one_donor)) && all(truth(units$one_specimen)) &&
  all(truth(units$one_suspension)) &&
  all(truth(units$one_sequencing_process)) &&
  all(units$library_construction_approach == "10x 3' v2") &&
  all(units$metadata_cell_estimate == 7000L) &&
  !any(truth(units$estimate_is_post_qc_cell_count)) &&
  !any(truth(units$protected_donor_attributes_published))
independence_ok <- nrow(independence) == 9L &&
  !any(truth(independence$overlap_found)) &&
  all(independence$exact_overlap_count == 0L) &&
  all(independence$canonical_overlap_count == 0L) &&
  independence$evidence_state[independence$identifier_class ==
    "biosample_accessions"] == "not_represented_in_hca_compact_manifest" &&
  independence$evidence_state[independence$identifier_class ==
    "real_world_person_identity"] ==
      "not_provable_from_pseudonymized_metadata"
terms_ok <- nrow(terms) == 5L && all(grepl("^https://", terms$locator)) &&
  any(terms$state == "CC BY 4.0") && all(truth(terms$recorded_before_download))
view_ok <- nrow(views) == 10L &&
  all(c("all_500_ordered_genes", "minimum_384_usable_cells",
        "reference_transform_reuse", "landscape_definition") %in% views$gate) &&
  !truth(views$current_sprint_passed[views$gate == "all_500_ordered_genes"]) &&
  !truth(views$current_sprint_passed[views$gate == "minimum_384_usable_cells"]) &&
  truth(views$current_sprint_passed[views$gate == "landscape_definition"]) &&
  all(views$failure_action[views$gate %in% c("all_500_ordered_genes",
    "minimum_384_usable_cells")] == "block both views") &&
  !any(truth(views$expression_values_opened))
resource_ok <- nrow(resource_plan) == 6L &&
  identical(as.integer(resource_plan$unit_count), c(5L, 8L, 8L, 40L, 80L,
                                                     19840L)) &&
  sum(truth(resource_plan$current_authorized)) == 1L &&
  nrow(resources) == 1L && resources$primary_files == 8L &&
  resources$total_payload_bytes == 202770089 && truth(resources$payload_within_cap) &&
  resources$ph_jobs == 80L &&
  resources$query_reference_component_distances == 19840L &&
  resources$fallback_cap_gib == 12 && resources$query_query_distances == 0L &&
  resources$reference_recomputations == 0L &&
  resources$observed_ph_peak_rss_bytes > 0 &&
  truth(resources$projection_excludes_download_conversion_and_outcome_time)
execution_ok <- nrow(execution) == 10L && execution$current_state[[1L]] == "next" &&
  all(execution$current_state[-1L] == "not_authorized") &&
  sum(truth(execution$expression_download_required)) == 1L &&
  sum(truth(execution$outcome_access_required)) == 1L &&
  nrow(stop_resume) == 9L && !any(truth(stop_resume$scientific_exclusion_allowed))
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
  decision$disposition == "metadata_closed_download_candidate" &&
  decision$primary_biological_units == 8L &&
  !truth(decision$known_identifier_overlap) && truth(decision$metadata_gates_passed) &&
  !truth(decision$structural_gates_passed) &&
  decision$cell_view_admission == "metadata_admitted_structural_pending" &&
  decision$gene_view_admission == "metadata_admitted_structural_pending" &&
  decision$next_owner_decision == "authorize_exact_eight_file_structural_download" &&
  !truth(decision$external_expression_download_authorized) &&
  !truth(decision$new_pca_ph_landscape_distance_clustering_authorized) &&
  !truth(decision$manuscript_claim_promoted) &&
  !truth(decision$final_author_roles_decided)
artifact_ok <- nrow(artifact_manifest) == 14L &&
  !anyDuplicated(artifact_manifest$filename) &&
  all(artifact_manifest$accepted_head == expected_head) &&
  all(vapply(seq_len(nrow(artifact_manifest)), function(i) {
    path <- file.path(production_dir, artifact_manifest$filename[[i]])
    file.exists(path) && identical(tolower(sha(path)),
      tolower(artifact_manifest$sha256[[i]])) &&
      file.info(path)$size == artifact_manifest$bytes[[i]]
  }, logical(1L))) &&
  !any(truth(artifact_manifest$contains_expression_values)) &&
  !any(truth(artifact_manifest$contains_protected_donor_attributes)) &&
  !any(truth(artifact_manifest$contains_expiring_signed_urls)) &&
  !any(truth(artifact_manifest$contains_private_paths)) &&
  !any(truth(artifact_manifest$contains_confidential_review_text))
repeat_ok <- all(vapply(required, function(name)
  identical(sha(file.path(production_dir, name)),
            sha(file.path(repeat_dir, name))), logical(1L)))
public_text <- paste(vapply(file.path(production_dir, required), function(path)
  paste(readLines(path, warn = FALSE), collapse = "\n"), character(1L)),
  collapse = "\n")
privacy_ok <- !grepl("/mnt/|[A-Z]:\\\\|X-Amz-|organism_age|biological_sex|contactName|email@|pasted-text|Dear Dr Rouchka|reviewer",
  public_text, ignore.case = TRUE) &&
  !grepl("donor_organism[.]sex|donor_organism[.]organism_age",
         public_text, ignore.case = TRUE)

checks <- data.frame(
  contract_id = "mv08c_hca_independent_validation_v1",
  check = c("source_manifest", "project_identity", "component_accountability",
    "primary_file_manifest", "primary_unit_linkage", "independence_limits",
    "access_terms", "view_specific_gates", "resource_queue",
    "execution_and_stop_resume", "landscape_contract", "decision_boundary",
    "artifact_manifest", "byte_identical_repeat", "privacy_and_confidentiality"),
  passed = c(source_ok, project_ok, component_ok, file_ok, unit_ok,
    independence_ok, terms_ok, view_ok, resource_ok, execution_ok,
    landscape_ok, decision_ok, artifact_ok, repeat_ok, privacy_ok),
  detail = c(
    "ten sources hashed; private main metadata remains hash-only",
    "official dcp60 project identity, access, restriction, and update metadata close",
    "24 whole-marrow samples partition 8 individual + 16 pooled; 63 selected samples excluded",
    "eight stable H5 files total 202770089 bytes with versions and SHA-256",
    "eight unique donor/specimen/suspension/process chains; donor attributes suppressed",
    "no known exact/canonical overlap; BioSample and real-person limitations explicit",
    "CC BY 4.0, DUA, release policy, and manifest meaning recorded",
    "500-gene, 384-cell, transform, and separate-view gates remain prospective",
    "80 PH jobs and 19840 component distances reproduce within resource caps",
    "only owner review is next; exact stop and resume receipts are frozen",
    "all active levels, separate H0/H1, and grid-free exact/error-controlled L2 retained",
    "stops at metadata-closed candidate before H5 download or scientific computation",
    "fourteen public artifacts independently rehashed",
    "production and repeat fifteen-file dossiers are byte-identical",
    "no protected donor attributes, signed URLs, private paths, expression, or review text"),
  stringsAsFactors = FALSE)
if (any(!checks$passed)) {
  stop("MV8-C independent validation failed: ",
       paste(checks$check[!checks$passed], collapse = ", "), call. = FALSE)
}
utils::write.table(checks, validation_output, sep = ",", row.names = FALSE,
                   col.names = TRUE, quote = TRUE, na = "", qmethod = "double")
cat("MV8-C independent validation passed:", nrow(checks), "/", nrow(checks),
    "\n")
