#!/usr/bin/env Rscript

options(warn = 2, digits = 17, scipen = 999)

mv08c_truth <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  tolower(trimws(as.character(x))) == "true"
}

mv08c_sha <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

mv08c_hash_text <- function(x) {
  vapply(as.character(x), digest::digest, character(1L), algo = "sha256",
         serialize = FALSE)
}

mv08c_scalar <- function(x, default = NA_character_) {
  if (is.null(x) || !length(x)) return(default)
  value <- unlist(x, recursive = TRUE, use.names = FALSE)
  value <- value[!is.na(value)]
  if (!length(value)) default else as.character(value[[1L]])
}

mv08c_canonical <- function(x) {
  gsub("[^a-z0-9]", "", tolower(trimws(as.character(x))))
}

mv08c_write_csv <- function(x, output_dir, name) {
  path <- file.path(output_dir, name)
  if (file.exists(path)) stop("Refusing to overwrite: ", path, call. = FALSE)
  utils::write.table(x, path, sep = ",", row.names = FALSE, col.names = TRUE,
                     quote = TRUE, na = "", qmethod = "double")
}

mv08c_read_json <- function(path) {
  jsonlite::fromJSON(path, simplifyVector = FALSE)
}

mv08c_file_summary <- function(hit, format) {
  found <- Filter(function(x) identical(mv08c_scalar(x$format), format),
                  hit$fileTypeSummaries)
  if (length(found) != 1L) {
    stop("Expected one ", format, " file summary for sample ", hit$entryId,
         ".", call. = FALSE)
  }
  found[[1L]]
}

mv08c_protocol_values <- function(hit, field) {
  unique(unlist(lapply(hit$protocols, function(x) x[[field]]),
                recursive = TRUE, use.names = FALSE))
}

mv08c_main <- function() {
  if (!requireNamespace("digest", quietly = TRUE)) stop("digest required.")
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("jsonlite required.")
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 12L) {
    stop(paste("usage: build_mv08c_hca_admission_dossier.R",
      "PROJECT_JSON WHOLE_SAMPLES_JSON SELECTED_SAMPLES_JSON",
      "PRIMARY_MANIFEST OPENAPI_JSON QUERY_CONTRACT MAIN_METADATA",
      "PANEL_CSV PH_METRICS SOURCE_METRICS OUTPUT_DIR EXPECTED_HEAD"),
      call. = FALSE)
  }
  input_paths <- args[1:10]
  if (any(!file.exists(input_paths))) {
    stop("MV8-C source inputs are incomplete: ",
         paste(input_paths[!file.exists(input_paths)], collapse = ", "),
         call. = FALSE)
  }
  output_dir <- args[[11L]]
  expected_head <- tolower(trimws(args[[12L]]))
  observed_head <- tolower(trimws(system2("git", c("rev-parse", "HEAD"),
                                          stdout = TRUE)))
  if (!grepl("^[0-9a-f]{40}$", expected_head) ||
      !identical(observed_head, expected_head)) {
    stop("MV8-C prospective HEAD mismatch: expected ", expected_head,
         "; observed ", observed_head, ".", call. = FALSE)
  }
  if (dir.exists(output_dir) &&
      length(list.files(output_dir, all.files = TRUE, no.. = TRUE))) {
    stop("MV8-C output directory must be empty.", call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  project_json <- mv08c_read_json(input_paths[[1L]])
  whole_json <- mv08c_read_json(input_paths[[2L]])
  selected_json <- mv08c_read_json(input_paths[[3L]])
  manifest <- utils::read.delim(input_paths[[4L]], stringsAsFactors = FALSE,
                                check.names = FALSE, quote = "", sep = "\t")
  openapi <- mv08c_read_json(input_paths[[5L]])
  query_contract <- utils::read.csv(input_paths[[6L]], stringsAsFactors = FALSE,
                                     check.names = FALSE)
  main_metadata <- utils::read.csv(input_paths[[7L]], stringsAsFactors = FALSE,
                                    check.names = FALSE)
  panel <- utils::read.csv(input_paths[[8L]], stringsAsFactors = FALSE,
                           check.names = FALSE)
  ph_metrics <- utils::read.csv(input_paths[[9L]], stringsAsFactors = FALSE,
                                check.names = FALSE)
  source_metrics <- utils::read.csv(input_paths[[10L]], stringsAsFactors = FALSE,
                                    check.names = FALSE)

  expected_counts <- c(1L, 24L, 63L, 8L)
  if (nrow(query_contract) != 4L ||
      !identical(as.integer(query_contract$query_order), 1:4) ||
      !identical(as.integer(query_contract$expected_entities), expected_counts) ||
      any(mv08c_truth(query_contract$expression_download_authorized)) ||
      any(mv08c_truth(query_contract$outcome_access_authorized))) {
    stop("MV8-C frozen query contract changed or opened a forbidden payload.",
         call. = FALSE)
  }
  if (nrow(panel) != 500L || anyDuplicated(panel$panel_order) ||
      anyDuplicated(panel$gene) ||
      !identical(as.integer(panel$panel_order), 1:500)) {
    stop("MV8-C requires the exact ordered 500-gene MV7-FP panel.",
         call. = FALSE)
  }
  observed_counts <- c(project_json$pagination$count,
    whole_json$pagination$count, selected_json$pagination$count, nrow(manifest))
  if (!identical(as.integer(observed_counts), expected_counts)) {
    stop("Official HCA entity counts do not match the prefreeze.", call. = FALSE)
  }
  required_manifest_columns <- c("file_name", "file_format", "file_size",
    "file_uuid", "file_version", "file_sha256", "file_content_type",
    "file_drs_uri", "file_azul_url", "file_mirror_uri",
    "cell_suspension.provenance.document_id",
    "cell_suspension.selected_cell_type",
    "sequencing_process.provenance.document_id",
    "library_preparation_protocol.library_construction_approach",
    "project.provenance.document_id",
    "specimen_from_organism.provenance.document_id",
    "donor_organism.provenance.document_id",
    "donor_organism.genus_species", "sample.provenance.document_id",
    "sequencing_input.provenance.document_id")
  missing_manifest <- setdiff(required_manifest_columns, names(manifest))
  if (length(missing_manifest)) {
    stop("Compact manifest missing: ", paste(missing_manifest, collapse = ", "),
         call. = FALSE)
  }
  unique_columns <- c("file_uuid", "file_name",
    "cell_suspension.provenance.document_id",
    "sequencing_process.provenance.document_id",
    "specimen_from_organism.provenance.document_id",
    "donor_organism.provenance.document_id")
  if (nrow(manifest) != 8L ||
      any(vapply(manifest[unique_columns], anyDuplicated, integer(1L)) > 0L) ||
      any(manifest$file_format != "h5") ||
      any(!grepl("^[0-9a-f]{64}$", manifest$file_sha256)) ||
      any(as.numeric(manifest$file_size) <= 0) ||
      any(manifest$`library_preparation_protocol.library_construction_approach` !=
          "10x 3' v2") ||
      any(manifest$`donor_organism.genus_species` != "Homo sapiens")) {
    stop("The HCA primary manifest fails identity, checksum, or protocol gates.",
         call. = FALSE)
  }

  whole_hits <- whole_json$hits
  selected_hits <- selected_json$hits
  donor_counts <- vapply(whole_hits, function(hit)
    as.integer(hit$donorOrganisms[[1L]]$donorCount), integer(1L))
  if (sum(donor_counts == 1L) != 8L || sum(donor_counts == 8L) != 16L ||
      any(!donor_counts %in% c(1L, 8L))) {
    stop("Whole-marrow donor-count partition is not 8 + 16.", call. = FALSE)
  }
  selected_donor_counts <- vapply(selected_hits, function(hit)
    as.integer(hit$donorOrganisms[[1L]]$donorCount), integer(1L))
  if (length(selected_hits) != 63L || any(selected_donor_counts != 1L)) {
    stop("Selected marrow exclusion component changed.", call. = FALSE)
  }
  hit_ids <- vapply(whole_hits, function(hit) mv08c_scalar(hit$entryId),
                    character(1L))
  if (anyDuplicated(hit_ids)) stop("Duplicate whole-marrow sample entry ID.")
  hit_by_id <- setNames(whole_hits, hit_ids)
  sample_ids <- manifest$`sample.provenance.document_id`
  if (anyNA(match(sample_ids, names(hit_by_id)))) {
    stop("Every manifest file must map to a whole-marrow sample hit.", call. = FALSE)
  }
  primary_hits <- hit_by_id[sample_ids]
  primary_checks <- vapply(primary_hits, function(hit) {
    h5 <- mv08c_file_summary(hit, "h5")
    selected_type <- unlist(hit$cellSuspensions[[1L]]$selectedCellType,
                            use.names = FALSE)
    species <- unlist(hit$donorOrganisms[[1L]]$genusSpecies, use.names = FALSE)
    library <- mv08c_protocol_values(hit, "libraryConstructionApproach")
    isTRUE(hit$samples[[1L]]$accessible) &&
      identical(as.integer(hit$donorOrganisms[[1L]]$donorCount), 1L) &&
      "Homo sapiens" %in% species && "10x 3' v2" %in% library &&
      "mononuclear cell of bone marrow" %in% selected_type &&
      identical(as.integer(h5$count), 1L) &&
      identical(h5$isIntermediate, FALSE) &&
      "Count matrix" %in% unlist(h5$contentDescription, use.names = FALSE)
  }, logical(1L))
  if (!all(primary_checks)) stop("A primary sample fails exact metadata gates.")

  file_order <- order(as.integer(sub(".*BM([0-9]+).*", "\\1",
                                      manifest$file_name)))
  manifest <- manifest[file_order, , drop = FALSE]
  primary_hits <- hit_by_id[manifest$`sample.provenance.document_id`]
  cell_estimates <- vapply(primary_hits, function(hit)
    as.integer(hit$cellSuspensions[[1L]]$totalCells), integer(1L))
  if (any(cell_estimates != 7000L)) {
    stop("Primary HCA metadata cell estimates changed.", call. = FALSE)
  }

  project <- project_json$hits[[1L]]$projects[[1L]]
  dates <- project_json$hits[[1L]]$dates[[1L]]
  accessions <- vapply(project$accessions, function(x)
    paste0(mv08c_scalar(x$namespace), ":", mv08c_scalar(x$accession)),
    character(1L))
  if (!identical(mv08c_scalar(project$projectId),
                 "cc95ff89-2e68-4a08-a234-480eca21ce79") ||
      any(manifest$`project.provenance.document_id` !=
          "cc95ff89-2e68-4a08-a234-480eca21ce79")) {
    stop("HCA project identity differs from the prefreeze.", call. = FALSE)
  }
  project_summary <- data.frame(
    contract_id = "mv08c_hca_project_summary_v1",
    accessed_date = "2026-08-16",
    catalog = "dcp60",
    api_version = mv08c_scalar(openapi$info$version),
    project_id = mv08c_scalar(project$projectId),
    project_label = mv08c_scalar(project$projectShortname),
    project_title = mv08c_scalar(project$projectTitle),
    project_accessions = paste(sort(accessions), collapse = ";"),
    aggregate_update_date = mv08c_scalar(dates$aggregateUpdateDate),
    project_cell_estimate = as.integer(project$estimatedCellCount),
    project_donor_count = as.integer(
      project_json$hits[[1L]]$donorOrganisms[[1L]]$donorCount),
    genus_species = "Homo sapiens",
    data_use_restriction = mv08c_scalar(project$dataUseRestriction),
    accessible = isTRUE(project$accessible),
    project_url = paste0(
      "https://explore.data.humancellatlas.org/projects/",
      mv08c_scalar(project$projectId)),
    expression_values_opened = FALSE,
    stringsAsFactors = FALSE)

  component_accountability <- data.frame(
    contract_id = "mv08c_hca_component_accountability_v1",
    component_order = 1:4,
    component_id = c("primary_one_donor_whole_marrow",
      "excluded_eight_donor_whole_marrow_pools",
      "excluded_selected_marrow_hematopoietic_cells",
      "excluded_aggregate_marrow_h5ad"),
    entity_type = c("sample", "sample", "sample", "file"),
    entity_count = c(sum(donor_counts == 1L), sum(donor_counts == 8L),
      length(selected_hits), 1L),
    donor_count_per_entity = c("1", "8", "1", "multiple_or_aggregated"),
    selected_cell_scope = c("whole marrow mononuclear cells",
      "whole marrow mononuclear cells", "selected marrow hematopoietic cells",
      "mixed genetic pools and controls"),
    primary_eligible = c(TRUE, FALSE, FALSE, FALSE),
    disposition = c("metadata_admitted_structural_pending",
      "excluded_pooled_donors", "excluded_composition_restricted",
      "excluded_aggregated_matrix"),
    outcome_or_topology_used = FALSE,
    stringsAsFactors = FALSE)

  primary_file_manifest <- data.frame(
    contract_id = "mv08c_hca_primary_file_manifest_v1",
    file_order = seq_len(nrow(manifest)),
    unit_id = sprintf("HCA_BM_%03d", seq_len(nrow(manifest))),
    file_name = manifest$file_name,
    file_uuid = manifest$file_uuid,
    file_version = manifest$file_version,
    file_size_bytes = as.numeric(manifest$file_size),
    file_sha256 = manifest$file_sha256,
    file_format = manifest$file_format,
    file_content_type = manifest$file_content_type,
    drs_uri = manifest$file_drs_uri,
    azul_download_url = manifest$file_azul_url,
    mirror_uri = manifest$file_mirror_uri,
    expression_download_authorized = FALSE,
    stringsAsFactors = FALSE)

  primary_unit_ledger <- data.frame(
    contract_id = "mv08c_hca_primary_unit_ledger_v1",
    unit_order = seq_len(nrow(manifest)),
    unit_id = primary_file_manifest$unit_id,
    file_uuid = manifest$file_uuid,
    donor_document_sha256 = mv08c_hash_text(
      manifest$`donor_organism.provenance.document_id`),
    specimen_document_sha256 = mv08c_hash_text(
      manifest$`specimen_from_organism.provenance.document_id`),
    suspension_document_sha256 = mv08c_hash_text(
      manifest$`cell_suspension.provenance.document_id`),
    sequencing_process_sha256 = mv08c_hash_text(
      manifest$`sequencing_process.provenance.document_id`),
    sequencing_input_matches_suspension =
      manifest$`sequencing_input.provenance.document_id` ==
      manifest$`cell_suspension.provenance.document_id`,
    one_donor = TRUE,
    one_specimen = TRUE,
    one_suspension = TRUE,
    one_sequencing_process = TRUE,
    library_construction_approach =
      manifest$`library_preparation_protocol.library_construction_approach`,
    metadata_cell_estimate = cell_estimates,
    estimate_is_post_qc_cell_count = FALSE,
    protected_donor_attributes_published = FALSE,
    stringsAsFactors = FALSE)

  main_strings <- unique(trimws(as.character(unlist(main_metadata,
    recursive = TRUE, use.names = FALSE))))
  main_strings <- main_strings[nzchar(main_strings) & !is.na(main_strings)]
  main_canonical <- unique(mv08c_canonical(main_strings))
  hca_accession_values <- vapply(project$accessions, function(x)
    mv08c_scalar(x$accession), character(1L))
  fields <- list(
    project_uuid = mv08c_scalar(project$projectId),
    project_accessions = hca_accession_values,
    donor_codes = manifest$`donor_organism.biomaterial_core.biomaterial_id`,
    donor_documents = manifest$`donor_organism.provenance.document_id`,
    specimen_documents = manifest$`specimen_from_organism.provenance.document_id`,
    sequencing_processes = manifest$`sequencing_process.provenance.document_id`,
    file_identifiers = c(manifest$file_uuid, manifest$file_name))
  overlap_rows <- lapply(seq_along(fields), function(i) {
    values <- unique(fields[[i]])
    exact <- sum(values %in% main_strings)
    canonical <- sum(mv08c_canonical(values) %in% main_canonical)
    data.frame(check_order = i, identifier_class = names(fields)[[i]],
      hca_identifier_count = length(values), exact_overlap_count = exact,
      canonical_overlap_count = canonical,
      overlap_found = exact > 0L || canonical > 0L,
      evidence_state = "tested_against_all_frozen_main_metadata_fields",
      limitation = "exact and punctuation-insensitive pseudonym comparison only",
      stringsAsFactors = FALSE)
  })
  independence_audit <- do.call(rbind, overlap_rows)
  independence_audit <- rbind(independence_audit,
    data.frame(check_order = 8L, identifier_class = "biosample_accessions",
      hca_identifier_count = 0L, exact_overlap_count = 0L,
      canonical_overlap_count = 0L, overlap_found = FALSE,
      evidence_state = "not_represented_in_hca_compact_manifest",
      limitation = "no direct BioSample comparison is possible from the frozen compact format",
      stringsAsFactors = FALSE),
    data.frame(check_order = 9L, identifier_class = "real_world_person_identity",
      hca_identifier_count = 8L, exact_overlap_count = 0L,
      canonical_overlap_count = 0L, overlap_found = FALSE,
      evidence_state = "not_provable_from_pseudonymized_metadata",
      limitation = "absence of known metadata overlap is not proof of identity separation",
      stringsAsFactors = FALSE))
  independence_audit$contract_id <-
    "mv08c_hca_independence_audit_v1"
  independence_audit <- independence_audit[, c("contract_id", setdiff(
    names(independence_audit), "contract_id"))]
  if (any(independence_audit$overlap_found)) {
    stop("A known HCA/main-corpus identifier overlap was detected.", call. = FALSE)
  }

  access_terms <- data.frame(
    contract_id = "mv08c_hca_access_terms_v1",
    term_order = 1:5,
    item = c("project_access_state", "dataset_license", "data_use_agreement",
      "data_release_policy", "manifest_payload_definition"),
    state = c("access_granted; project metadata NRES", "CC BY 4.0",
      "recorded", "recorded", "compact TSV is metadata, not data payload"),
    locator = c(project_summary$project_url,
      "https://creativecommons.org/licenses/by/4.0/",
      "https://data.humancellatlas.org/about/data-use-agreement",
      "https://www.humancellatlas.org/data-release-policy/",
      "https://data.humancellatlas.org/guides/quick-start-guide"),
    obligation = c("honor current project access status and portal terms",
      "attribute source, link license, and indicate modifications",
      "apply the dataset license and release policy",
      "respect data-generator publication interests and local requirements",
      "do not treat or publish the compact manifest as expression data"),
    recorded_before_download = TRUE,
    stringsAsFactors = FALSE)

  view_gates <- data.frame(
    contract_id = "mv08c_hca_view_gates_v1",
    gate_order = 1:10,
    gate = c("eight_independent_units", "count_file_identity",
      "shared_frozen_500_gene_source", "stable_feature_mapping",
      "all_500_ordered_genes", "minimum_384_usable_cells",
      "reference_transform_reuse", "cell_30pc_coordinates",
      "gene_pearson_chord_coordinates", "landscape_definition"),
    cell_view_state = c("pass_metadata", "pass_metadata", "required",
      "pending_structural_download", "pending_structural_download",
      "pending_structural_download", "pending_transform_validation",
      "pending_calculation", "not_applicable", "pass_contract"),
    gene_view_state = c("pass_metadata", "pass_metadata", "required",
      "pending_structural_download", "pending_structural_download",
      "pending_structural_download", "pending_transform_validation",
      "not_applicable", "pending_calculation", "pass_contract"),
    hard_gate = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
    current_sprint_passed = c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, TRUE),
    failure_action = c("block cohort", "block file", "block both views",
      "block both views", "block both views", "block both views",
      "block external projection", "block cell view", "block gene view",
      "block both views"),
    expression_values_opened = FALSE,
    stringsAsFactors = FALSE)

  ph_elapsed <- as.numeric(ph_metrics$elapsed_seconds)
  ph_peak <- as.numeric(ph_metrics$peak_process_tree_rss_bytes)
  source_peak <- as.numeric(source_metrics$peak_process_tree_rss_bytes)
  payload_bytes <- sum(as.numeric(manifest$file_size))
  ph_jobs <- 8L * 5L * 2L
  distance_rows <- 8L * 124L * 5L * 2L * 2L
  resource_plan <- data.frame(
    contract_id = "mv08c_hca_resource_plan_v1",
    plan_order = 1:6,
    stage = c("metadata_complete", "future_exact_file_download",
      "future_structural_gate", "future_query_projection",
      "future_ph", "future_landscape_distance"),
    unit_count = c(5L, 8L, 8L, 40L, ph_jobs, distance_rows),
    unit = c("metadata inputs", "H5 files", "H5 files", "donor-seeds",
      "view-level PH jobs", "query-reference component distances"),
    planned_bytes = c(NA, payload_bytes, payload_bytes, NA, NA, NA),
    current_authorized = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE),
    resource_gate = c("complete", "2 GiB exact-payload cap",
      "sequential bounded inspection", "immutable 500-gene and 30-PC transform",
      "Ripserr primary; exact GUDHI fallback under 12 GiB",
      "streamed/chunked exact or error-controlled squared-L2"),
    stringsAsFactors = FALSE)
  resource_summary <- data.frame(
    contract_id = "mv08c_hca_resource_summary_v1",
    primary_files = 8L,
    total_payload_bytes = payload_bytes,
    total_payload_gib = payload_bytes / 1024^3,
    metadata_stage_download_cap_gib = 2,
    payload_within_cap = payload_bytes <= 2 * 1024^3,
    deterministic_seeds_per_donor = 5L,
    ph_jobs = ph_jobs,
    query_reference_component_distances = distance_rows,
    observed_mv07h_ph_median_seconds = stats::median(ph_elapsed),
    observed_mv07h_ph_p95_seconds = unname(stats::quantile(ph_elapsed, .95,
                                                           type = 7L)),
    projected_ph_worker_hours_median = ph_jobs * stats::median(ph_elapsed) / 3600,
    projected_ph_worker_hours_p95 = ph_jobs * unname(
      stats::quantile(ph_elapsed, .95, type = 7L)) / 3600,
    observed_ph_peak_rss_bytes = max(ph_peak),
    observed_source_fit_peak_rss_bytes = max(source_peak),
    fallback_cap_gib = 12,
    query_query_distances = 0L,
    reference_recomputations = 0L,
    projection_excludes_download_conversion_and_outcome_time = TRUE,
    stringsAsFactors = FALSE)

  execution_plan <- data.frame(
    contract_id = "mv08c_hca_execution_plan_v1",
    step_order = 1:10,
    step = c("owner_review", "exact_eight_file_download",
      "byte_and_sha256_verification", "structural_h5_inspection",
      "identifier_and_500_gene_gate", "384_cell_gate",
      "immutable_reference_projection", "separate_cell_gene_ph",
      "all_level_h0_h1_landscape_distances", "label_opened_retrieval_evaluation"),
    current_state = c("next", rep("not_authorized", 9L)),
    resume_receipt = c("mv08c decision row", "one receipt per exact file",
      "verified byte/checksum ledger", "shape/type/axis ledger",
      "ordered feature mapping ledger", "per-donor QC depth ledger",
      "transform and coordinate hashes", "per-view PH records",
      "chunk manifests and component distances", "prediction lock and endpoint ledger"),
    stop_condition = c("owner declines or changes scope", "file identity changed",
      "byte count or checksum mismatch", "ambiguous or invalid matrix structure",
      "unstable IDs or any frozen gene absent", "any donor below 384 usable cells",
      "reference transform cannot be reused without refit",
      "timeout, memory, engine, or repeat failure", "distance or repeat mismatch",
      "label firewall or prediction lock failure"),
    expression_download_required = c(FALSE, TRUE, rep(FALSE, 8L)),
    outcome_access_required = c(rep(FALSE, 9L), TRUE),
    stringsAsFactors = FALSE)

  stop_resume <- data.frame(
    contract_id = "mv08c_hca_stop_resume_v1",
    rule_order = 1:9,
    gate = c("manifest_identity", "download_checksum", "matrix_structure",
      "feature_mapping", "frozen_panel", "usable_cell_depth",
      "reference_transform", "ph_resource", "reproducibility"),
    stop_state = c("file set/version changes", "bytes or SHA-256 mismatch",
      "orientation/type/axes ambiguous", "IDs do not map uniquely",
      "one or more ordered genes absent", "one or more donors below 384",
      "center/scale/rotation requires refit", "12 GiB or timeout guard reached",
      "production/repeat scientific hashes differ"),
    scientific_exclusion_allowed = FALSE,
    resume_requirement = c("new prospective manifest decision",
      "fresh exact download with matching receipt", "prefrozen parser remediation",
      "prefrozen mapping remediation", "separate pre-outcome sensitivity contract",
      "separate pre-outcome depth contract", "validated immutable projection interface",
      "approved exact fallback or revised resource prefreeze",
      "root-cause remediation and clean repeat"),
    stringsAsFactors = FALSE)

  landscape_contract <- data.frame(
    contract_id = "mv08c_hca_landscape_contract_v1",
    item_order = 1:11,
    item = c("comparison_unit_sample", "cell_view_observations_cells",
      "gene_view_observations_fixed_500_genes", "finite_positive_intervals",
      "essential_H0_excluded", "all_consecutive_active_levels",
      "H0_H1_separate", "squared_L2_exact_or_error_controlled",
      "no_universal_grid", "no_universal_level_cap",
      "streamed_or_chunked_execution"),
    required_state = TRUE,
    current_state = TRUE,
    external_data_may_change = FALSE,
    stringsAsFactors = FALSE)

  decision <- data.frame(
    contract_id = "mv08c_hca_admission_decision_v1",
    project_id = mv08c_scalar(project$projectId),
    disposition = "metadata_closed_download_candidate",
    primary_biological_units = 8L,
    known_identifier_overlap = FALSE,
    independence_claim = "no_known_metadata_overlap_with_named_pseudonym_and_biosample_limits",
    metadata_gates_passed = TRUE,
    structural_gates_passed = FALSE,
    cell_view_admission = "metadata_admitted_structural_pending",
    gene_view_admission = "metadata_admitted_structural_pending",
    next_owner_decision = "authorize_exact_eight_file_structural_download",
    external_expression_download_authorized = FALSE,
    new_pca_ph_landscape_distance_clustering_authorized = FALSE,
    manuscript_claim_promoted = FALSE,
    final_author_roles_decided = FALSE,
    stringsAsFactors = FALSE)

  source_locators <- c(
    "https://service.azul.data.humancellatlas.org/index/projects",
    "https://service.azul.data.humancellatlas.org/index/samples",
    "https://service.azul.data.humancellatlas.org/index/samples",
    "official compact manifest generated from frozen MV8-C query",
    "https://service.azul.data.humancellatlas.org/openapi.json",
    "docs/specifications/mv08c-hca-query-contract-v1.csv",
    "private_local_metadata_not_published",
    "docs/audits/mv07fp-panel-evidence/mv07fp-panel.csv",
    "docs/audits/mv07h-full-ph-evidence/mv07h-ph-metrics.csv",
    "docs/audits/mv07h-full-ph-evidence/mv07h-source-metrics.csv")
  source_manifest <- data.frame(
    contract_id = "mv08c_hca_source_manifest_v1",
    source_order = seq_along(input_paths),
    source_id = c("project_json", "whole_marrow_sample_json",
      "selected_marrow_sample_json", "primary_h5_compact_manifest",
      "openapi_schema", "query_contract", "frozen_main_metadata",
      "frozen_500_gene_panel", "mv07h_ph_metrics", "mv07h_source_metrics"),
    locator = source_locators,
    sha256 = vapply(input_paths, mv08c_sha, character(1L)),
    bytes = as.numeric(file.info(input_paths)$size),
    access_class = c(rep("official_public_metadata", 5L),
      "public_repository", "private_hash_only", rep("public_repository", 3L)),
    raw_input_published = c(rep(FALSE, 5L), TRUE, FALSE, TRUE, TRUE, TRUE),
    contains_expression_values = FALSE,
    accepted_head = expected_head,
    stringsAsFactors = FALSE)

  outputs <- list(
    "mv08c-source-manifest.csv" = source_manifest,
    "mv08c-project-summary.csv" = project_summary,
    "mv08c-component-accountability.csv" = component_accountability,
    "mv08c-primary-file-manifest.csv" = primary_file_manifest,
    "mv08c-primary-unit-ledger.csv" = primary_unit_ledger,
    "mv08c-independence-audit.csv" = independence_audit,
    "mv08c-access-terms.csv" = access_terms,
    "mv08c-view-gates.csv" = view_gates,
    "mv08c-resource-plan.csv" = resource_plan,
    "mv08c-resource-summary.csv" = resource_summary,
    "mv08c-execution-plan.csv" = execution_plan,
    "mv08c-stop-resume.csv" = stop_resume,
    "mv08c-landscape-contract.csv" = landscape_contract,
    "mv08c-decision.csv" = decision)
  for (name in names(outputs)) mv08c_write_csv(outputs[[name]], output_dir, name)
  artifact_manifest <- data.frame(
    contract_id = "mv08c_hca_artifact_manifest_v1",
    artifact_order = seq_along(outputs),
    filename = names(outputs),
    sha256 = vapply(file.path(output_dir, names(outputs)), mv08c_sha,
                    character(1L)),
    bytes = as.numeric(file.info(file.path(output_dir, names(outputs)))$size),
    contains_expression_values = FALSE,
    contains_protected_donor_attributes = FALSE,
    contains_expiring_signed_urls = FALSE,
    contains_private_paths = FALSE,
    contains_confidential_review_text = FALSE,
    accepted_head = expected_head,
    stringsAsFactors = FALSE)
  mv08c_write_csv(artifact_manifest, output_dir,
                  "mv08c-artifact-manifest.csv")
  cat("MV8-C HCA metadata dossier built: 8 files,", payload_bytes,
      "planned bytes,", ph_jobs, "PH jobs,", distance_rows,
      "component distances; expression payloads remain unopened.\n")
}

if (sys.nframe() == 0L) mv08c_main()
