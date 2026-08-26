# Prospective contracts and source-inventory helpers for MV17.

.mv17a_sha256_file <- function(path) {
  if (!file.exists(path) || dir.exists(path)) {
    stop("MV17-A source artifact is absent: ", path, call. = FALSE)
  }
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

.mv17a_set_sha256 <- function(values) {
  values <- sort(tolower(as.character(values)), method = "radix")
  if (!length(values) || anyNA(values) ||
      any(!grepl("^[0-9a-f]{64}$", values))) {
    stop("MV17-A artifact-set hashes are invalid", call. = FALSE)
  }
  digest::digest(paste(values, collapse = "\n"), algo = "sha256",
                 serialize = FALSE)
}

mv17a_estimand_registry_v1 <- function() {
  grid <- expand.grid(
    view = c("cell", "gene"),
    homology_dimension = c("H0", "H1"),
    stringsAsFactors = FALSE
  )
  grid <- grid[order(grid$view, grid$homology_dimension,
                     method = "radix"), , drop = FALSE]
  rownames(grid) <- NULL
  grid$contract_id <- "mv17a_calibration_estimand_v1"
  grid$point_axis <- ifelse(grid$view == "cell", "selected_cells", "genes")
  grid$coordinate_axis <- ifelse(
    grid$view == "cell", "shared_30_pc_coordinates",
    "selected_384_cell_coordinates"
  )
  grid$geometry <- ifelse(
    grid$view == "cell", "euclidean_shared_pca_v1",
    "euclidean_correlation_chord_v1"
  )
  grid$filtration <- "complete_vietoris_rips_all_finite_positive_intervals"
  grid$essential_h0 <- ifelse(grid$homology_dimension == "H0",
                              "excluded", "not_applicable")
  grid$calibration_question <- ifelse(
    grid$homology_dimension == "H0",
    "does_finite_component_merger_structure_exceed_axis_matched_null_structure",
    "does_finite_loop_structure_exceed_axis_matched_null_structure"
  )
  grid$primary_object <- "persistence_diagram_dimension_separate"
  grid$null_family_selected <- FALSE
  grid$summary_selected <- FALSE
  grid$real_calibration_authorized <- FALSE
  grid[c(
    "contract_id", "view", "homology_dimension", "point_axis",
    "coordinate_axis", "geometry", "filtration", "essential_h0",
    "calibration_question", "primary_object", "null_family_selected",
    "summary_selected", "real_calibration_authorized"
  )]
}

mv17a_localization_registry_v1 <- function() {
  data.frame(
    contract_id = "mv17a_localization_output_v1",
    view = rep(c("cell", "gene"), each = 2L),
    homology_dimension = rep(c("H0", "H1"), 2L),
    target = c(
      "mst_component_merger_edges_between_cells",
      "representative_cycle_support_over_cells",
      "mst_component_merger_edges_between_genes",
      "representative_cycle_support_over_genes"
    ),
    private_atomic_output = c(
      "canonical_endpoint_tokens_and_death_scale",
      "canonical_support_tokens_scale_and_nonuniqueness",
      "canonical_endpoint_tokens_and_death_scale",
      "canonical_support_tokens_scale_and_nonuniqueness"
    ),
    public_output = c(
      "aggregate_merger_influence_only", "aggregate_cycle_stability_only",
      "aggregate_merger_influence_only", "aggregate_cycle_stability_only"
    ),
    coefficient_field = c("not_applicable", "freeze_in_MV17_D",
                          "not_applicable", "freeze_in_MV17_D"),
    tie_policy = "freeze_in_MV17_D_before_fixture_execution",
    method_selected = FALSE,
    real_localization_authorized = FALSE,
    stringsAsFactors = FALSE
  )
}

mv17a_h2_fixture_registry_v1 <- function() {
  data.frame(
    contract_id = "mv17a_h2_fixture_class_v1",
    fixture_class = c(
      "sphere", "torus", "circle", "gaussian_cloud",
      "shuffled_sphere", "shuffled_torus", "shuffled_circle"
    ),
    control_role = c(
      "H2_positive", "H1_H2_positive", "H2_negative_H1_positive",
      "H2_negative_or_null_calibrated", "shuffled_null",
      "shuffled_null", "shuffled_null"
    ),
    expected_H1 = c("not_primary", "positive", "positive", "null",
                    "destroyed_or_attenuated", "destroyed_or_attenuated",
                    "destroyed_or_attenuated"),
    expected_H2 = c("positive", "positive", "negative", "null",
                    "destroyed_or_attenuated", "destroyed_or_attenuated",
                    "negative"),
    generator_parameters = "freeze_in_MV17_E_before_fixture_execution",
    seed = "freeze_in_MV17_E_before_fixture_execution",
    tolerance = "freeze_in_MV17_E_before_fixture_execution",
    execution_authorized = FALSE,
    real_data_H2_authorized = FALSE,
    stringsAsFactors = FALSE
  )
}

mv17a_stage_gate_registry_v1 <- function() {
  data.frame(
    contract_id = "mv17a_stage_gate_v1",
    stage = c("MV17-B", "MV17-C", "MV17-D", "MV17-E", "MV17-F",
              "MV17-G", "MV17-H", "MV17-I"),
    lane = c("null", "null", "localization", "H2_fixture", "H2_resource",
             "calibration_localization", "H2_real", "incremental_evaluation"),
    prerequisite = c(
      "MV17_A_closed", "MV17_B_closed", "MV17_A_closed", "MV17_A_closed",
      "MV17_E_closed", "MV17_C_and_MV17_D_closed", "MV17_F_closed",
      "MV17_G_and_MV17_H_closed_plus_explicit_owner_authorization"
    ),
    implementation_eligible_after_MV17A = c(TRUE, FALSE, TRUE, TRUE, FALSE,
                                            FALSE, FALSE, FALSE),
    scientific_execution_authorized = FALSE,
    exact_head_prefreeze_required = TRUE,
    stringsAsFactors = FALSE
  )
}

mv17a_firewall_v1 <- function() {
  data.frame(
    contract_id = "mv17a_firewall_v1",
    surface = c(
      "labels", "outcomes", "tissues", "biology", "clustering",
      "view_ranking", "fusion", "fusion_weight_tuning", "manuscript_claims",
      "existing_PH_mutation", "existing_landscape_mutation",
      "existing_comparison_mutation"
    ),
    state = "closed",
    authorized = FALSE,
    stringsAsFactors = FALSE
  )
}

mv17a_validate_contract_v1 <- function(estimands, localization, fixtures,
                                       gates, firewall) {
  if (nrow(estimands) != 4L ||
      !setequal(estimands$view, c("cell", "gene")) ||
      !setequal(estimands$homology_dimension, c("H0", "H1")) ||
      any(estimands$null_family_selected) || any(estimands$summary_selected) ||
      any(estimands$real_calibration_authorized)) {
    stop("MV17-A calibration estimand drift", call. = FALSE)
  }
  if (nrow(localization) != 4L || any(localization$method_selected) ||
      any(localization$real_localization_authorized) ||
      !all(grepl("aggregate_", localization$public_output))) {
    stop("MV17-A localization contract drift", call. = FALSE)
  }
  if (!setequal(fixtures$fixture_class, c(
      "sphere", "torus", "circle", "gaussian_cloud", "shuffled_sphere",
      "shuffled_torus", "shuffled_circle")) ||
      any(fixtures$execution_authorized) || any(fixtures$real_data_H2_authorized)) {
    stop("MV17-A H2 fixture-class contract drift", call. = FALSE)
  }
  if (!identical(gates$stage, c("MV17-B", "MV17-C", "MV17-D", "MV17-E",
                                "MV17-F", "MV17-G", "MV17-H", "MV17-I")) ||
      any(gates$scientific_execution_authorized) ||
      !all(gates$exact_head_prefreeze_required)) {
    stop("MV17-A stage-gate drift", call. = FALSE)
  }
  if (nrow(firewall) != 12L || any(firewall$authorized) ||
      !all(firewall$state == "closed")) {
    stop("MV17-A firewall drift", call. = FALSE)
  }
  TRUE
}

mv17a_inventory_rds_groups_v1 <- function(group_paths) {
  if (!length(group_paths) || any(!file.exists(group_paths))) {
    stop("MV17-A cell group set is incomplete", call. = FALSE)
  }
  rows <- lapply(group_paths, function(path) {
    group <- readRDS(path)
    required <- c("dataset_scope", "panel_id", "seed", "model", "records")
    if (!all(required %in% names(group)) || !length(group$records) ||
        !all(c("gene_ids", "rotation", "n_components", "cache_key") %in%
             names(group$model))) {
      stop("MV17-A cell group schema drift", call. = FALSE)
    }
    record_ok <- vapply(group$records, function(record) {
      all(c("view", "result", "oracle") %in% names(record)) &&
        identical(dim(record$view$payload), c(384L, 30L)) &&
        length(record$view$point_ids) == 384L &&
        identical(colnames(record$result$diagram),
                  c("dimension", "birth", "death")) &&
        identical(as.integer(record$oracle$finite_h0_intervals), 383L) &&
        isTRUE(record$oracle$passed)
    }, logical(1L))
    if (!all(record_ok) || ncol(group$model$rotation) != 30L ||
        nrow(group$model$rotation) != length(group$model$gene_ids) ||
        group$model$n_components != 30L) {
      stop("MV17-A cell axes or diagrams drift", call. = FALSE)
    }
    data.frame(
      source_family = "cell", dataset_scope = group$dataset_scope,
      panel_id = group$panel_id, seed = as.integer(group$seed),
      units = length(group$records), feature_axis = length(group$model$gene_ids),
      pca_components = as.integer(group$model$n_components),
      selected_observations = 384L,
      dimension_records = 2L * length(group$records),
      artifact_bytes = as.numeric(file.info(path)$size),
      artifact_sha256 = .mv17a_sha256_file(path),
      schema_sha256 = digest::digest(
        paste(sort(c(names(group), names(group$model),
                     names(group$records[[1L]]$view),
                     names(group$records[[1L]]$result))), collapse = "|"),
        algo = "sha256", serialize = FALSE
      ), stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

mv17a_inventory_gene_ph_v1 <- function(bindings, source_roots) {
  required <- c(
    "group_order", "job_id", "source_role", "output_file", "output_bytes",
    "output_sha256", "diagram_sha256", "finite_h0_intervals",
    "finite_h1_intervals", "outcome_label_state", "biological_outcomes_computed"
  )
  if (!all(required %in% names(bindings)) || nrow(bindings) != 2544L ||
      !setequal(bindings$source_role, names(source_roots))) {
    stop("MV17-A gene PH bindings drift", call. = FALSE)
  }
  keys <- paste(bindings$source_role, bindings$output_file, sep = "\r")
  first <- !duplicated(keys)
  artifacts <- bindings[first, , drop = FALSE]
  paths <- mapply(file.path, source_roots[artifacts$source_role],
                  artifacts$output_file, USE.NAMES = FALSE)
  if (length(paths) != 1272L || any(!file.exists(paths))) {
    stop("MV17-A gene PH corpus is incomplete", call. = FALSE)
  }
  hashes <- vapply(paths, .mv17a_sha256_file, character(1L))
  bytes <- as.numeric(file.info(paths)$size)
  selected_axis <- character(length(paths))
  panel_axis <- character(length(paths))
  diagram_hash <- character(length(paths))
  point_count <- integer(length(paths))
  schema_hash <- character(length(paths))
  for (index in seq_along(paths)) {
    record <- readRDS(paths[[index]])
    if (!all(c("identity", "topology_result") %in% names(record)) ||
        !identical(colnames(record$topology_result$diagram),
                   c("dimension", "birth", "death"))) {
      stop("MV17-A gene PH record schema drift", call. = FALSE)
    }
    selected_axis[[index]] <- record$identity$selected_cell_sha256
    panel_axis[[index]] <- record$identity$panel_sha256
    diagram_hash[[index]] <- record$topology_result$provenance$diagram_sha256
    point_count[[index]] <- as.integer(record$identity$point_count)
    schema_hash[[index]] <- digest::digest(paste(sort(c(
      names(record), names(record$identity), names(record$topology_result),
      names(record$topology_result$provenance)
    )), collapse = "|"), algo = "sha256", serialize = FALSE)
  }
  if (!identical(unname(bytes), unname(as.numeric(artifacts$output_bytes))) ||
      !identical(unname(tolower(hashes)),
                 unname(tolower(artifacts$output_sha256))) ||
      !identical(unname(tolower(diagram_hash)),
                 unname(tolower(artifacts$diagram_sha256))) ||
      any(!grepl("^[0-9a-f]{64}$", selected_axis)) ||
      any(!grepl("^[0-9a-f]{64}$", panel_axis)) ||
      !setequal(point_count, c(475L, 500L)) ||
      length(unique(schema_hash)) != 1L ||
      any(bindings$outcome_label_state != "closed") ||
      any(as.logical(bindings$biological_outcomes_computed))) {
    stop("MV17-A gene PH identities or hashes drift", call. = FALSE)
  }
  data.frame(
    source_family = "gene_PH",
    artifacts = length(paths), dimension_records = nrow(bindings),
    selected_axis_identities = length(selected_axis),
    panel_axis_identities = length(panel_axis),
    point_count_min = min(point_count), point_count_max = max(point_count),
    bytes = sum(bytes), artifact_set_sha256 = .mv17a_set_sha256(hashes),
    diagram_set_sha256 = .mv17a_set_sha256(diagram_hash),
    selected_axis_set_sha256 = .mv17a_set_sha256(selected_axis),
    panel_axis_set_sha256 = .mv17a_set_sha256(panel_axis),
    schema_sha256 = unique(schema_hash), stringsAsFactors = FALSE
  )
}

mv17a_inventory_delimited_chunks_v1 <- function(root, completion,
                                                prefix) {
  paths <- sort(list.files(root, pattern = "^distances[.]csv$",
                           recursive = TRUE, full.names = TRUE),
                method = "radix")
  if (!length(paths) || nrow(completion) != length(paths) ||
      !all(c("pair_count", "distances_bytes", "distances_sha256") %in%
           names(completion))) {
    stop("MV17-A distance chunk inventory drift", call. = FALSE)
  }
  hashes <- vapply(paths, .mv17a_sha256_file, character(1L))
  bytes <- as.numeric(file.info(paths)$size)
  rows <- vapply(paths, function(path) length(readLines(path, warn = FALSE)) - 1L,
                 integer(1L))
  headers <- vapply(paths, function(path) readLines(path, n = 1L, warn = FALSE),
                    character(1L))
  completion <- completion[order(as.integer(completion$production_order)),
                           , drop = FALSE]
  if (!identical(unname(as.integer(rows)),
                 unname(as.integer(completion$pair_count))) ||
      !identical(unname(bytes), unname(as.numeric(completion$distances_bytes))) ||
      !identical(unname(tolower(hashes)),
                 unname(tolower(completion$distances_sha256))) ||
      length(unique(headers)) != 1L) {
    stop("MV17-A distance artifact rehash drift", call. = FALSE)
  }
  data.frame(
    source_family = prefix,
    artifacts = length(paths),
    records = sum(rows),
    bytes = sum(bytes),
    artifact_set_sha256 = .mv17a_set_sha256(hashes),
    schema_sha256 = digest::digest(headers[[1L]], algo = "sha256",
                                   serialize = FALSE),
    stringsAsFactors = FALSE
  )
}
