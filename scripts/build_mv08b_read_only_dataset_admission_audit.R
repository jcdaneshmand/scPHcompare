#!/usr/bin/env Rscript

options(warn = 2, digits = 17, scipen = 999)

mv08b_truth <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  tolower(trimws(as.character(x))) == "true"
}

mv08b_required_candidate_columns <- function() {
  c("contract_id", "candidate_order", "candidate_id", "repository_accession",
    "cohort_class", "planned_role", "required_comparison", "official_record",
    "minimum_independent_donors", "whole_bone_marrow_required",
    "cell_view_requested", "gene_view_requested", "selection_outcomes_opened",
    "external_expression_download_authorized", "new_ph_authorized")
}

mv08b_validate_registry <- function(registry) {
  missing <- setdiff(mv08b_required_candidate_columns(), names(registry))
  if (length(missing)) stop("MV8-B candidate registry missing: ",
                            paste(missing, collapse = ", "), call. = FALSE)
  if (nrow(registry) != 6L || anyDuplicated(registry$candidate_id) ||
      !identical(as.integer(registry$candidate_order), 1:6) ||
      any(mv08b_truth(registry$selection_outcomes_opened)) ||
      any(mv08b_truth(registry$external_expression_download_authorized)) ||
      any(mv08b_truth(registry$new_ph_authorized))) {
    stop("MV8-B frozen candidate registry is stale or outcome-opened.", call. = FALSE)
  }
  required_ids <- c("existing_gse120221_25",
    "hca_hematopoietic_immune_cell_atlas", "gse185381_healthy_controls",
    "gse180298_healthy_hspc", "gse212092_day0_bone_marrow",
    "tabula_sapiens_gse201333")
  if (!identical(registry$candidate_id, required_ids)) {
    stop("MV8-B candidate order or identity changed after prefreeze.", call. = FALSE)
  }
  invisible(TRUE)
}

mv08b_required_fact_columns <- function() {
  c("candidate_id", "accessed_date", "official_url", "organism",
    "specimen_scope", "independent_donors", "biological_sample_units",
    "library_records", "count_matrix_state", "donor_mapping_state",
    "feature_identifier_state", "access_terms", "known_accession_overlap",
    "donor_overlap_state", "estimated_download_gib", "cell_depth_state",
    "cell_view_admission", "gene_view_admission", "disposition",
    "supported_estimand", "unsupported_estimand", "unresolved_fields",
    "evidence_note")
}

mv08b_validate_facts <- function(facts, registry) {
  missing <- setdiff(mv08b_required_fact_columns(), names(facts))
  if (length(missing)) stop("MV8-B source facts missing: ",
                            paste(missing, collapse = ", "), call. = FALSE)
  if (nrow(facts) != nrow(registry) || anyDuplicated(facts$candidate_id) ||
      !setequal(facts$candidate_id, registry$candidate_id)) {
    stop("MV8-B facts must contain exactly one row per frozen candidate.",
         call. = FALSE)
  }
  allowed <- c("preferred_pending_download_authorization",
    "reserve_pending_metadata_resolution", "technical_reprocessing_only",
    "composition_shift_sensitivity_only", "small_feasibility_only",
    "multi_tissue_later_stage", "reject_overlap",
    "reject_unit_not_reconstructable", "reject_incompatible_input")
  if (any(!facts$disposition %in% allowed) ||
      any(!grepl("^https://", facts$official_url)) ||
      any(!nzchar(facts$supported_estimand)) ||
      any(!nzchar(facts$unsupported_estimand))) {
    stop("MV8-B facts contain an invalid disposition, URL, or estimand.",
         call. = FALSE)
  }
  invisible(TRUE)
}

mv08b_matrix_from_rdata <- function(path) {
  env <- new.env(parent = emptyenv())
  loaded <- load(path, envir = env)
  candidates <- loaded[vapply(loaded, function(name) {
    object <- env[[name]]
    length(dim(object)) == 2L && !is.null(rownames(object)) &&
      !is.null(colnames(object))
  }, logical(1L))]
  if (length(candidates) != 1L) {
    stop("Expected exactly one named matrix in source: ", basename(path),
         call. = FALSE)
  }
  env[[candidates[[1L]]]]
}

mv08b_normalize_barcodes <- function(x) sub("-1$", "", as.character(x))

mv08b_source_axis_summary <- function(validation_path, main_path, panel_genes) {
  validation <- mv08b_matrix_from_rdata(validation_path)
  main <- mv08b_matrix_from_rdata(main_path)
  validation_barcodes <- unique(mv08b_normalize_barcodes(colnames(validation)))
  main_barcodes <- unique(mv08b_normalize_barcodes(colnames(main)))
  shared_barcodes <- intersect(validation_barcodes, main_barcodes)
  validation_genes <- unique(as.character(rownames(validation)))
  main_genes <- unique(sub("[-_]ENSG[0-9]+\\..*$", "",
                           as.character(rownames(main))))
  data.frame(
    validation_source_cells = length(validation_barcodes),
    main_source_cells = length(main_barcodes),
    shared_cell_barcodes = length(shared_barcodes),
    barcode_jaccard = length(shared_barcodes) /
      length(union(validation_barcodes, main_barcodes)),
    validation_source_features = length(validation_genes),
    main_source_features = length(main_genes),
    shared_gene_symbols = length(intersect(validation_genes, main_genes)),
    validation_panel_genes = length(intersect(panel_genes, validation_genes)),
    main_panel_genes = length(intersect(panel_genes, main_genes)),
    expression_values_exported = FALSE,
    stringsAsFactors = FALSE)
}

mv08b_write_csv <- function(x, output_dir, name) {
  path <- file.path(output_dir, name)
  if (file.exists(path)) stop("Refusing to overwrite: ", path, call. = FALSE)
  utils::write.table(x, path, sep = ",", row.names = FALSE, col.names = TRUE,
                     quote = TRUE, na = "", qmethod = "double")
}

mv08b_main <- function() {
  if (!requireNamespace("digest", quietly = TRUE)) stop("digest required.")
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Matrix required.")
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 8L) {
    stop(paste("usage: build_mv08b_read_only_dataset_admission_audit.R",
      "MAIN_METADATA VALIDATION_METADATA NCBI_RUNINFO PANEL_CSV",
      "CANDIDATE_REGISTRY CANDIDATE_FACTS OUTPUT_DIR EXPECTED_HEAD"),
      call. = FALSE)
  }
  input_paths <- args[1:6]
  if (any(!file.exists(input_paths))) {
    stop("MV8-B source inputs are incomplete: ",
         paste(input_paths[!file.exists(input_paths)], collapse = ", "))
  }
  output_dir <- args[[7L]]
  expected_head <- tolower(trimws(args[[8L]]))
  observed_head <- tolower(trimws(system2("git", c("rev-parse", "HEAD"),
                                          stdout = TRUE)))
  if (!grepl("^[0-9a-f]{40}$", expected_head) ||
      !identical(observed_head, expected_head)) {
    stop("MV8-B prospective HEAD mismatch: expected ", expected_head,
         "; observed ", observed_head, ".", call. = FALSE)
  }
  if (dir.exists(output_dir) &&
      length(list.files(output_dir, all.files = TRUE, no.. = TRUE))) {
    stop("MV8-B output directory must be empty.", call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  readc <- function(path) utils::read.csv(path, stringsAsFactors = FALSE,
                                           check.names = FALSE)
  main <- readc(input_paths[[1L]])
  validation <- readc(input_paths[[2L]])
  runinfo <- readc(input_paths[[3L]])
  panel <- readc(input_paths[[4L]])
  registry <- readc(input_paths[[5L]])
  facts <- readc(input_paths[[6L]])
  mv08b_validate_registry(registry)
  mv08b_validate_facts(facts, registry)

  main_bm <- main[main$SRA == "SRA779509", , drop = FALSE]
  validation_gsm <- validation[["GEO_Accession (exp)"]]
  if (nrow(main_bm) != 25L || nrow(validation) != 25L ||
      nrow(runinfo) != 25L || anyDuplicated(runinfo$Sample) ||
      anyDuplicated(runinfo$SampleName)) {
    stop("MV8-B GSE120221 expected 25 unique main, validation, and NCBI rows.",
         call. = FALSE)
  }
  main_srs <- main_bm[["SRS Number"]]
  main_idx <- match(runinfo$Sample, main_srs)
  validation_idx <- match(runinfo$SampleName, validation_gsm)
  if (anyNA(main_idx) || anyNA(validation_idx) ||
      !all(runinfo$Submission == "SRA779509") ||
      !all(runinfo$SRAStudy == "SRP162214") ||
      !all(runinfo$BioProject == "PRJNA492227")) {
    stop("MV8-B official accession mapping does not close exactly.", call. = FALSE)
  }
  panel_genes <- unique(as.character(panel$gene))
  if (length(panel_genes) != 500L) stop("MV8-B requires the frozen 500-gene panel.")

  mapping <- data.frame(
    contract_id = "mv08b_gse120221_overlap_map_v1",
    mapping_order = seq_len(nrow(runinfo)),
    main_sample_id = main_bm$orig.ident[main_idx],
    submission = runinfo$Submission,
    sra_study = runinfo$SRAStudy,
    srs = runinfo$Sample,
    srr = runinfo$Run,
    biosample = runinfo$BioSample,
    gsm = runinfo$SampleName,
    bioproject = runinfo$BioProject,
    validation_sample = validation$orig.ident[validation_idx],
    exact_accession_overlap = TRUE,
    independent_external_library = FALSE,
    stringsAsFactors = FALSE)

  source_rows <- vector("list", nrow(mapping))
  for (i in seq_len(nrow(mapping))) {
    validation_path <- validation[["File Path"]][validation_idx[[i]]]
    main_path <- main_bm[["File Path"]][main_idx[[i]]]
    if (!file.exists(validation_path) || !file.exists(main_path)) {
      stop("MV8-B local source missing for accession ", mapping$srs[[i]])
    }
    axis <- mv08b_source_axis_summary(validation_path, main_path, panel_genes)
    source_rows[[i]] <- cbind(
      data.frame(contract_id = "mv08b_gse120221_source_compatibility_v1",
        mapping_order = i, srs = mapping$srs[[i]], gsm = mapping$gsm[[i]],
        validation_post_qc_cells = as.integer(
          validation$Number_of_Cells_After_Filtering[validation_idx[[i]]]),
        main_post_qc_cells = as.integer(
          main_bm$Number_of_Cells_After_Filtering[main_idx[[i]]]),
        fixed_384_cell_depth_eligible =
          as.integer(validation$Number_of_Cells_After_Filtering[validation_idx[[i]]]) >= 384L &&
          as.integer(main_bm$Number_of_Cells_After_Filtering[main_idx[[i]]]) >= 384L,
        stringsAsFactors = FALSE), axis)
  }
  compatibility <- do.call(rbind, source_rows)

  facts <- facts[match(registry$candidate_id, facts$candidate_id), , drop = FALSE]
  candidate_audit <- cbind(registry, facts[setdiff(names(facts), "candidate_id")])
  candidate_audit$preferred_by_read_only_audit <-
    candidate_audit$disposition == "preferred_pending_download_authorization"
  candidate_audit$all_hard_gates_resolved <-
    candidate_audit$preferred_by_read_only_audit &
    tolower(trimws(candidate_audit$unresolved_fields)) == "none"
  candidate_audit$download_authorized <- FALSE
  candidate_audit$new_ph_authorized <- FALSE
  candidate_audit$landscape_contract <-
    "all finite positive intervals; essential H0 excluded; all active levels; H0/H1 separate; exact or error-controlled squared-L2"

  resource <- readc("docs/audits/mv07h-full-ph-evidence/mv07h-ph-metrics.csv")
  source_resource <- readc(
    "docs/audits/mv07h-full-ph-evidence/mv07h-source-metrics.csv")
  elapsed <- as.numeric(resource$elapsed_seconds)
  peak <- as.numeric(resource$peak_process_tree_rss_bytes)
  candidate_audit$projected_ph_jobs <-
    2 * 5 * as.numeric(candidate_audit$biological_sample_units)
  candidate_audit$projected_ph_worker_hours_median <-
    candidate_audit$projected_ph_jobs * stats::median(elapsed) / 3600
  candidate_audit$projected_ph_worker_hours_conservative_p95 <-
    candidate_audit$projected_ph_jobs *
      unname(stats::quantile(elapsed, .95, type = 7L)) / 3600
  candidate_audit$observed_ph_peak_rss_bytes <- max(peak)
  candidate_audit$observed_124_source_fit_peak_rss_bytes <-
    max(as.numeric(source_resource$peak_process_tree_rss_bytes))
  candidate_audit$resource_projection_excludes_download_and_conversion <- TRUE

  summary <- data.frame(
    contract_id = "mv08b_gse120221_overlap_summary_v1",
    validation_libraries = nrow(validation),
    main_sra779509_libraries = nrow(main_bm),
    exact_accession_overlaps = sum(mapping$exact_accession_overlap),
    official_submission = unique(runinfo$Submission),
    official_study = unique(runinfo$SRAStudy),
    official_bioproject = unique(runinfo$BioProject),
    validation_min_post_qc_cells = min(compatibility$validation_post_qc_cells),
    main_min_post_qc_cells = min(compatibility$main_post_qc_cells),
    fixed_depth_eligible_libraries = sum(compatibility$fixed_384_cell_depth_eligible),
    libraries_with_shared_barcodes = sum(compatibility$shared_cell_barcodes > 0L),
    minimum_validation_panel_genes = min(compatibility$validation_panel_genes),
    minimum_main_panel_genes = min(compatibility$main_panel_genes),
    independent_external_validation = FALSE,
    technical_reprocessing_validation = TRUE,
    stringsAsFactors = FALSE)

  landscape <- data.frame(
    contract_id = "mv08b_landscape_continuation_contract_v1",
    item_order = 1:11,
    item = c("comparison_unit_sample", "cell_view_observations_cells",
      "gene_view_observations_fixed_500_genes", "finite_positive_intervals",
      "essential_H0_excluded", "all_consecutive_active_levels",
      "H0_H1_separate", "squared_L2_exact_or_error_controlled",
      "no_universal_grid", "no_universal_level_cap",
      "streamed_or_chunked_execution"),
    required_state = TRUE, current_state = TRUE,
    external_data_may_change = FALSE, stringsAsFactors = FALSE)

  decision <- data.frame(
    contract_id = "mv08b_dataset_admission_decision_v1",
    gse120221_disposition = "technical_reprocessing_only",
    preferred_candidate = facts$candidate_id[
      facts$disposition == "preferred_pending_download_authorization"][[1L]],
    next_authorized_action = "owner_review_of_dataset_and_compute_dossier",
    external_expression_download_authorized = FALSE,
    new_ph_authorized = FALSE,
    manuscript_claim_promoted = FALSE,
    confidential_material_published = FALSE,
    final_author_roles_decided = FALSE,
    stringsAsFactors = FALSE)

  source_paths <- c(input_paths,
    "docs/audits/mv07h-full-ph-evidence/mv07h-ph-metrics.csv",
    "docs/audits/mv07h-full-ph-evidence/mv07h-source-metrics.csv")
  source_manifest <- data.frame(
    contract_id = "mv08b_source_manifest_v1",
    source_order = seq_along(source_paths),
    source_id = c("main_metadata", "validation_metadata", "ncbi_runinfo",
      "frozen_panel", "candidate_registry", "candidate_facts",
      "mv07h_ph_metrics", "mv07h_source_metrics"),
    locator = c("private_local_metadata_not_published",
      "private_local_metadata_not_published",
      "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc=SRP162214",
      "docs/audits/mv07fp-panel-evidence/mv07fp-panel.csv",
      "docs/specifications/mv08b-dataset-candidate-registry-v1.csv",
      "docs/audits/mv08b-dataset-source-facts-v1.csv",
      "docs/audits/mv07h-full-ph-evidence/mv07h-ph-metrics.csv",
      "docs/audits/mv07h-full-ph-evidence/mv07h-source-metrics.csv"),
    sha256 = vapply(source_paths, function(path)
      digest::digest(file = path, algo = "sha256", serialize = FALSE),
      character(1L)),
    bytes = as.numeric(file.info(source_paths)$size),
    access_class = c("private_hash_only", "private_hash_only",
      "official_public_metadata", rep("public_repository", 5L)),
    contains_expression_values = FALSE,
    accepted_head = expected_head,
    stringsAsFactors = FALSE)

  outputs <- list(
    "mv08b-source-manifest.csv" = source_manifest,
    "mv08b-gse120221-overlap-map.csv" = mapping,
    "mv08b-gse120221-source-compatibility.csv" = compatibility,
    "mv08b-gse120221-overlap-summary.csv" = summary,
    "mv08b-candidate-admission-audit.csv" = candidate_audit,
    "mv08b-landscape-continuation-contract.csv" = landscape,
    "mv08b-decision.csv" = decision)
  for (name in names(outputs)) mv08b_write_csv(outputs[[name]], output_dir, name)
  manifest <- data.frame(
    contract_id = "mv08b_artifact_manifest_v1",
    filename = names(outputs),
    sha256 = vapply(file.path(output_dir, names(outputs)), function(path)
      digest::digest(file = path, algo = "sha256", serialize = FALSE),
      character(1L)),
    bytes = as.numeric(file.info(file.path(output_dir, names(outputs)))$size),
    contains_expression_values = FALSE,
    contains_private_paths = FALSE,
    contains_confidential_review_text = FALSE,
    accepted_head = expected_head,
    stringsAsFactors = FALSE)
  mv08b_write_csv(manifest, output_dir, "mv08b-artifact-manifest.csv")
}

if (sys.nframe() == 0L) mv08b_main()
