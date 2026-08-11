# MV5-O prospective label-closed complete-matrix production prefreeze.

.mv05o_digest <- function(value) {
  digest::digest(value, algo = "sha256", serialize = TRUE)
}

.mv05o_hash_ok <- function(value) {
  is.character(value) && all(grepl("^[0-9a-f]{64}$", value))
}

.mv05o_assert_label_closed <- function(value, label = "input") {
  prohibited <- c("tissue", "approach", "class", "label", "outcome",
                  "biological_label", "technical_label")
  found <- intersect(tolower(names(value)), prohibited)
  if (length(found)) {
    stop(label, " contains prohibited outcome columns: ",
         paste(found, collapse = ", "), ".", call. = FALSE)
  }
  if ("outcome_label_state" %in% names(value) &&
      any(value$outcome_label_state != "closed")) {
    stop(label, " is not label-closed.", call. = FALSE)
  }
  if ("biological_outcomes_computed" %in% names(value) &&
      any(as.logical(value$biological_outcomes_computed))) {
    stop(label, " reports biological outcomes.", call. = FALSE)
  }
  invisible(TRUE)
}

mv05o_validate_mv05n_inventories_v1 <- function(groups, chunks) {
  group_required <- c(
    "source_group_id", "group_order", "fold_id", "seed", "representation",
    "training_samples", "unordered_training_pairs", "h0_h1_request_rows",
    "chunk_count", "request_identity_set_sha256", "pair_scope",
    "landscape_definition_id", "outcome_label_state",
    "biological_outcomes_computed"
  )
  chunk_required <- c(
    "chunk_id", "source_group_id", "group_order", "fold_id", "seed",
    "representation", "homology_dimension", "request_rows",
    "request_identity_set_sha256", "first_pair_request_id",
    "last_pair_request_id", "outcome_label_state",
    "biological_outcomes_computed"
  )
  if (!is.data.frame(groups) || !is.data.frame(chunks) ||
      !all(group_required %in% names(groups)) ||
      !all(chunk_required %in% names(chunks))) {
    stop("MV5-O requires the complete MV5-N group/chunk inventories.",
         call. = FALSE)
  }
  .mv05o_assert_label_closed(groups, "groups")
  .mv05o_assert_label_closed(chunks, "chunks")
  representations <- c("sct_whole", "inductive_integrated")
  if (nrow(groups) != 150L || nrow(chunks) != 4340L ||
      anyDuplicated(groups[c("source_group_id", "representation")]) ||
      anyDuplicated(chunks$chunk_id) ||
      !setequal(groups$representation, representations) ||
      !setequal(chunks$representation, representations) ||
      any(table(groups$representation) != 75L) ||
      any(table(chunks$representation) != 2170L) ||
      any(tapply(groups$unordered_training_pairs, groups$representation, sum) !=
          262675L) ||
      any(tapply(groups$h0_h1_request_rows, groups$representation, sum) !=
          525350L) ||
      any(tapply(chunks$request_rows, chunks$representation, sum) != 525350L) ||
      any(groups$pair_scope != "training_training_unordered") ||
      any(groups$landscape_definition_id !=
          "all_active_exact_critical_pairs_v1") ||
      !setequal(chunks$homology_dimension, c("H0", "H1")) ||
      any(chunks$request_rows < 1L | chunks$request_rows > 250L) ||
      !.mv05o_hash_ok(groups$request_identity_set_sha256) ||
      !.mv05o_hash_ok(chunks$request_identity_set_sha256)) {
    stop("MV5-N inventories violate the exact MV5-O production scope.",
         call. = FALSE)
  }
  key <- paste(groups$source_group_id, groups$representation, sep = "\r")
  chunk_key <- paste(chunks$source_group_id, chunks$representation, sep = "\r")
  if (!all(chunk_key %in% key)) {
    stop("MV5-O chunks do not map to frozen production groups.", call. = FALSE)
  }
  chunk_counts <- table(chunk_key)
  row_counts <- tapply(chunks$request_rows, chunk_key, sum)
  if (any(unname(chunk_counts[key]) != groups$chunk_count) ||
      any(unname(row_counts[key]) != groups$h0_h1_request_rows)) {
    stop("MV5-O chunk aggregation does not reproduce group inventories.",
         call. = FALSE)
  }
  invisible(TRUE)
}

mv05o_source_freeze_v1 <- function(paths, roles, base_revision) {
  if (!length(paths) || length(paths) != length(roles) ||
      any(!file.exists(paths)) ||
      !grepl("^[0-9a-f]{40}$", base_revision)) {
    stop("MV5-O source freeze inputs are incomplete.", call. = FALSE)
  }
  normalized <- gsub("\\\\", "/", paths)
  if (any(grepl("(?i)(docs/private|reviewer|[.]pdf$|example_run[.]r$)",
                normalized, perl = TRUE))) {
    stop("MV5-O source freeze contains prohibited public paths.", call. = FALSE)
  }
  result <- data.frame(
    contract_id = "mv05o_source_freeze_v1",
    base_revision = base_revision,
    artifact_path = normalized,
    artifact_role = roles,
    size_bytes = unname(file.info(paths)$size),
    sha256 = vapply(paths, function(path) {
      digest::digest(file = path, algo = "sha256", serialize = FALSE)
    }, character(1L)),
    outcome_label_state = "closed",
    biological_outcomes_computed = FALSE,
    public_safe = TRUE,
    stringsAsFactors = FALSE
  )
  result <- result[order(result$artifact_path, method = "radix"), ]
  result$source_freeze_sha256 <- .mv05o_digest(
    paste(result$artifact_path, result$sha256, sep = "\r")
  )
  rownames(result) <- NULL
  result
}

mv05o_build_queues_v1 <- function(groups, chunks, source_freeze_sha256,
                                   implementation_hashes) {
  mv05o_validate_mv05n_inventories_v1(groups, chunks)
  if (!.mv05o_hash_ok(source_freeze_sha256) ||
      !is.character(implementation_hashes) ||
      !identical(sort(names(implementation_hashes)),
                 sort(c("prefreeze", "stager", "landscape", "baseline"))) ||
      !.mv05o_hash_ok(unname(implementation_hashes)) ||
      length(source_freeze_sha256) != 1L ||
      length(implementation_hashes) != 4L) {
    stop("MV5-O production hashes are invalid.", call. = FALSE)
  }
  implementation_set_sha256 <- .mv05o_digest(
    paste(names(implementation_hashes), implementation_hashes, sep = "\r")
  )
  groups <- groups[order(groups$representation, groups$fold_id, groups$seed,
                         method = "radix"), ]
  groups$production_group_id <- vapply(seq_len(nrow(groups)), function(index) {
    paste0("mv05o_group_v1:", .mv05o_digest(list(
      source_group_id = groups$source_group_id[[index]],
      request_root = groups$request_identity_set_sha256[[index]],
      source_freeze = source_freeze_sha256,
      implementation_set = implementation_set_sha256
    )))
  }, character(1L))
  groups$source_freeze_sha256 <- source_freeze_sha256
  groups$production_implementation_set_sha256 <- implementation_set_sha256
  groups$stager_implementation_sha256 <- implementation_hashes[["stager"]]
  groups$execution_order <- seq_len(nrow(groups))
  groups$landscape_request_rows <- groups$h0_h1_request_rows
  groups$energy_pair_rows <- groups$unordered_training_pairs
  groups$shared_pseudobulk_pair_rows <- ifelse(
    groups$representation == "sct_whole", groups$unordered_training_pairs, 0L
  )
  groups$shared_pseudobulk_disposition <- ifelse(
    groups$representation == "sct_whole", "compute_once",
    "reuse_fold_seed_identity_from_sct_whole"
  )
  groups$max_parallel_workers <- 2L
  groups$per_unit_elapsed_cap_seconds <- 900L
  groups$per_process_rss_cap_bytes <- 4 * 1024^3
  groups$stage_worker_hour_cap <- 21.6
  groups$stage_private_storage_cap_bytes <- 10 * 1024^3
  groups$status_schema_id <- "mv05o_atomic_status_v1"
  groups$output_schema_id <- "mv05o_distance_output_v1"
  groups$outcome_label_state <- "closed"
  groups$biological_outcomes_computed <- FALSE
  groups$clustering_jobs_executed <- 0L
  groups$production_executed <- FALSE

  group_map <- groups[c("source_group_id", "representation",
                        "production_group_id", "execution_order")]
  landscape <- merge(chunks, group_map,
                     by = c("source_group_id", "representation"), sort = FALSE)
  landscape <- landscape[order(landscape$execution_order,
                               landscape$homology_dimension,
                               landscape$first_pair_request_id,
                               method = "radix"), ]
  landscape$production_chunk_id <- vapply(seq_len(nrow(landscape)), function(index) {
    paste0("mv05o_landscape_chunk_v1:", .mv05o_digest(list(
      chunk_id = landscape$chunk_id[[index]],
      request_root = landscape$request_identity_set_sha256[[index]],
      implementation = implementation_hashes[["landscape"]]
    )))
  }, character(1L))
  landscape$source_freeze_sha256 <- source_freeze_sha256
  landscape$implementation_sha256 <- implementation_hashes[["landscape"]]
  landscape$unit_order <- seq_len(nrow(landscape))
  landscape$engine_id <- "persim_exact_critical_pairs_all_active_v1"
  landscape$per_unit_elapsed_cap_seconds <- 900L
  landscape$per_process_rss_cap_bytes <- 4 * 1024^3
  landscape$atomic_write_required <- TRUE
  landscape$resume_requires_hash_validation <- TRUE
  landscape$execution_authorized_after_prefreeze_commit <- TRUE
  landscape$production_executed <- FALSE
  landscape$outcome_label_state <- "closed"
  landscape$biological_outcomes_computed <- FALSE

  energy <- groups
  energy$baseline_method <- "cell_distribution_energy_v1"
  energy$baseline_unit_id <- vapply(seq_len(nrow(energy)), function(index) {
    paste0("mv05o_energy_group_v1:", .mv05o_digest(list(
      production_group_id = energy$production_group_id[[index]],
      pair_count = energy$unordered_training_pairs[[index]],
      implementation = implementation_hashes[["baseline"]]
    )))
  }, character(1L))
  energy$pair_rows <- energy$unordered_training_pairs
  energy$compute_disposition <- "compute_representation_specific"
  pseudo <- groups[groups$representation == "sct_whole", ]
  pseudo$baseline_method <- "pseudobulk_training_standardized_panel_v1"
  pseudo$baseline_unit_id <- vapply(seq_len(nrow(pseudo)), function(index) {
    paste0("mv05o_pseudobulk_group_v1:", .mv05o_digest(list(
      fold_id = pseudo$fold_id[[index]], seed = pseudo$seed[[index]],
      pair_count = pseudo$unordered_training_pairs[[index]],
      implementation = implementation_hashes[["baseline"]]
    )))
  }, character(1L))
  pseudo$pair_rows <- pseudo$unordered_training_pairs
  pseudo$compute_disposition <- "compute_once_reuse_across_representations"
  baseline_names <- c(
    "baseline_method", "baseline_unit_id", "fold_id", "seed",
    "representation", "production_group_id", "pair_rows",
    "compute_disposition", "per_unit_elapsed_cap_seconds",
    "per_process_rss_cap_bytes", "status_schema_id", "output_schema_id",
    "outcome_label_state", "biological_outcomes_computed",
    "clustering_jobs_executed", "production_executed"
  )
  baseline <- rbind(energy[baseline_names], pseudo[baseline_names])
  baseline$source_freeze_sha256 <- source_freeze_sha256
  baseline$implementation_sha256 <- implementation_hashes[["baseline"]]
  baseline$unit_order <- seq_len(nrow(baseline))
  baseline$atomic_write_required <- TRUE
  baseline$resume_requires_hash_validation <- TRUE

  group_keep <- c(
    "production_group_id", "execution_order", "source_group_id", "group_order",
    "fold_id", "seed", "representation", "training_samples",
    "unordered_training_pairs", "request_identity_set_sha256", "chunk_count",
    "landscape_request_rows", "energy_pair_rows", "shared_pseudobulk_pair_rows",
    "shared_pseudobulk_disposition", "max_parallel_workers",
    "per_unit_elapsed_cap_seconds", "per_process_rss_cap_bytes",
    "stage_worker_hour_cap", "stage_private_storage_cap_bytes",
    "status_schema_id", "output_schema_id", "outcome_label_state",
    "biological_outcomes_computed", "clustering_jobs_executed",
    "production_executed", "source_freeze_sha256",
    "production_implementation_set_sha256", "stager_implementation_sha256"
  )
  list(groups = groups[group_keep], landscape = landscape, baseline = baseline)
}

mv05o_build_validation_plan_v1 <- function(group_queue, landscape_queue,
                                            baseline_queue) {
  .mv05o_assert_label_closed(group_queue, "group_queue")
  profiles <- do.call(rbind, lapply(split(group_queue, group_queue$representation),
                                    function(current) {
    current <- current[order(current$training_samples, current$fold_id,
                             current$seed, method = "radix"), ]
    sizes <- sort(unique(current$training_samples))
    targets <- c(min(sizes), stats::median(sizes), max(sizes))
    names(targets) <- c("minimum", "representative", "maximum")
    do.call(rbind, lapply(names(targets), function(profile) {
      eligible <- current[current$training_samples == targets[[profile]], ]
      eligible <- eligible[order(eligible$fold_id, eligible$seed,
                                 method = "radix"), ]
      data.frame(profile = profile,
                 representation = eligible$representation[[1L]],
                 production_group_id = eligible$production_group_id[[1L]],
                 stringsAsFactors = FALSE)
    }))
  }))
  oracle <- do.call(rbind, lapply(seq_len(nrow(profiles)), function(index) {
    current <- landscape_queue[
      landscape_queue$production_group_id == profiles$production_group_id[[index]], ]
    do.call(rbind, lapply(c("H0", "H1"), function(dimension) {
      eligible <- current[current$homology_dimension == dimension, ]
      eligible <- eligible[order(eligible$first_pair_request_id,
                                 eligible$production_chunk_id,
                                 method = "radix"), ]
      data.frame(
        contract_id = "mv05o_validation_plan_v1",
        validation_id = paste("exact_r_oracle", profiles$profile[[index]],
                              profiles$representation[[index]], dimension,
                              sep = ":"),
        validation_type = "independent_exact_r_oracle_first_request_v1",
        profile = profiles$profile[[index]],
        representation = profiles$representation[[index]],
        homology_dimension = dimension,
        production_group_id = profiles$production_group_id[[index]],
        production_unit_id = eligible$production_chunk_id[[1L]],
        pair_request_id = eligible$first_pair_request_id[[1L]],
        required_count = 1L, tolerance = 1e-10,
        outcome_label_state = "closed", biological_outcomes_computed = FALSE,
        stringsAsFactors = FALSE
      )
    }))
  }))
  maximum <- profiles[profiles$profile == "maximum", ]
  repeat_rows <- do.call(rbind, lapply(seq_len(nrow(maximum)), function(index) {
    data.frame(
      contract_id = "mv05o_validation_plan_v1",
      validation_id = paste0("maximum_group_clean_repeat:",
                             maximum$representation[[index]]),
      validation_type = "all_landscape_and_energy_outputs_byte_repeat_v1",
      profile = "maximum", representation = maximum$representation[[index]],
      homology_dimension = "H0_and_H1",
      production_group_id = maximum$production_group_id[[index]],
      production_unit_id = "all_group_units", pair_request_id = "all_group_pairs",
      required_count = sum(
        landscape_queue$production_group_id == maximum$production_group_id[[index]]
      ) + 1L,
      tolerance = 0, outcome_label_state = "closed",
      biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
    )
  }))
  resume <- data.frame(
    contract_id = "mv05o_validation_plan_v1",
    validation_id = "immutable_resume_all_units",
    validation_type = "hash_size_timestamp_unchanged_zero_rebuild_v1",
    profile = "all", representation = "all", homology_dimension = "all",
    production_group_id = "all_150_groups", production_unit_id = "all_4565_units",
    pair_request_id = "all_pair_roots", required_count =
      nrow(landscape_queue) + nrow(baseline_queue), tolerance = 0,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
  result <- rbind(oracle, repeat_rows, resume)
  rownames(result) <- NULL
  result
}

mv05o_abort_rules_v1 <- function() {
  data.frame(
    contract_id = "mv05o_abort_rules_v1",
    rule_id = sprintf("MV5O-ABORT-%02d", 1:10),
    trigger = c(
      "source_or_implementation_hash_mismatch",
      "pair_group_or_chunk_identity_root_mismatch",
      "missing_duplicate_or_query_query_request",
      "outcome_or_label_column_or_nonzero_downstream_counter",
      "unit_elapsed_seconds_greater_than_900",
      "process_tree_rss_bytes_greater_than_4294967296",
      "projected_or_observed_worker_hours_greater_than_21.6",
      "projected_or_observed_private_storage_bytes_greater_than_10737418240",
      "partial_or_stale_output_status_pair",
      "oracle_repeat_resume_or_baseline_formula_failure"
    ),
    disposition = c(
      rep("abort_before_new_unit_and_preserve_completed_immutable_artifacts", 8L),
      "quarantine_partial_artifact_abort_and_require_root_cause_review",
      "abort_stage_revoke_production_authorization_and_open_corrective_sprint"
    ),
    automatic_retry = FALSE,
    outcome_label_state = "closed",
    biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
}

mv05o_validate_prefreeze_v1 <- function(source_freeze, queues,
                                        validation_plan, abort_rules) {
  reconstructed_groups <- queues$groups
  reconstructed_groups$contract_id <-
    "mv05n_training_pair_group_inventory_v1"
  reconstructed_groups$h0_h1_request_rows <-
    reconstructed_groups$landscape_request_rows
  reconstructed_groups$pair_scope <- "training_training_unordered"
  reconstructed_groups$landscape_definition_id <-
    "all_active_exact_critical_pairs_v1"
  reconstructed_chunks <- queues$landscape
  reconstructed_chunks$contract_id <-
    "mv05n_training_pair_chunk_inventory_v1"
  mv05o_validate_mv05n_inventories_v1(
    groups = reconstructed_groups,
    chunks = reconstructed_chunks
  )
  lapply(c(queues, list(source_freeze = source_freeze,
                        validation_plan = validation_plan,
                        abort_rules = abort_rules)),
         .mv05o_assert_label_closed)
  if (nrow(source_freeze) < 12L || nrow(queues$groups) != 150L ||
      nrow(queues$landscape) != 4340L || nrow(queues$baseline) != 225L ||
      nrow(validation_plan) != 15L || nrow(abort_rules) != 10L ||
      sum(queues$landscape$request_rows) != 1050700L ||
      sum(queues$baseline$pair_rows[
        queues$baseline$baseline_method == "cell_distribution_energy_v1"]) !=
          525350L ||
      sum(queues$baseline$pair_rows[
        queues$baseline$baseline_method ==
          "pseudobulk_training_standardized_panel_v1"]) != 262675L ||
      any(queues$groups$production_executed) ||
      any(queues$groups$clustering_jobs_executed != 0L) ||
      any(abort_rules$automatic_retry) ||
      sum(validation_plan$validation_type ==
            "independent_exact_r_oracle_first_request_v1") != 12L ||
      sum(validation_plan$validation_type ==
            "all_landscape_and_energy_outputs_byte_repeat_v1") != 2L) {
    stop("MV5-O prefreeze validation failed.", call. = FALSE)
  }
  invisible(TRUE)
}
