# Internal MV6-F label-closed matched dual-view production contracts.

.mv06f_sha256 <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

.mv06f_required_hash <- function(value, name) {
  value <- tolower(as.character(value))
  if (length(value) != 1L || is.na(value) ||
      !grepl("^[0-9a-f]{64}$", value)) {
    stop("MV6-F requires one SHA-256 for ", name, ".", call. = FALSE)
  }
  value
}

.mv06f_forbidden_fields <- c(
  "tissue", "approach", "endpoint", "outcome", "class", "label",
  "biological_label", "technical_label", "cluster", "ari", "nmi"
)

mv06f_validate_prefreeze_inputs_v1 <- function(candidate, folds, resources,
                                                panel) {
  required_candidate <- c(
    "sample_id", "study", "outcome_label_state",
    "biological_outcomes_computed"
  )
  required_folds <- c(
    "fold_id", "fit_scope_id", "held_out_study", "seed",
    "training_samples", "held_out_samples", "outcome_label_state",
    "biological_outcomes_computed"
  )
  required_resources <- c(
    "sample_id", "seed", "selected_cell_sha256",
    "normalization_cache_key", "private_cache_file",
    "private_cache_size_bytes", "private_cache_sha256", "disposition",
    "outcome_label_state", "biological_outcomes_computed"
  )
  required_panel <- c(
    "panel_sha256", "panel_order", "feature_id", "gene",
    "selected_all_cache_core"
  )
  frames <- list(candidate = candidate, folds = folds, resources = resources,
                 panel = panel)
  if (!all(vapply(frames, is.data.frame, logical(1L))) ||
      nrow(candidate) != 90L || nrow(folds) != 75L ||
      nrow(resources) != 450L || nrow(panel) != 500L ||
      !all(required_candidate %in% names(candidate)) ||
      !all(required_folds %in% names(folds)) ||
      !all(required_resources %in% names(resources)) ||
      !all(required_panel %in% names(panel))) {
    stop("MV6-F input schemas or cardinalities differ from the freeze.",
         call. = FALSE)
  }
  if (any(vapply(frames, function(value) {
    any(tolower(names(value)) %in% .mv06f_forbidden_fields)
  }, logical(1L)))) {
    stop("MV6-F inputs cross the label firewall.", call. = FALSE)
  }
  closed <- function(value) {
    all(value$outcome_label_state == "closed") &&
      !any(as.logical(value$biological_outcomes_computed))
  }
  if (!closed(candidate) || !closed(folds) || !closed(resources) ||
      anyDuplicated(candidate$sample_id) ||
      length(unique(candidate$study)) != 15L ||
      anyDuplicated(paste(folds$fold_id, folds$seed, sep = "\r")) ||
      anyDuplicated(paste(resources$sample_id, resources$seed, sep = "\r")) ||
      !identical(sort(unique(as.integer(folds$seed))), 20260805:20260809) ||
      any(folds$training_samples + folds$held_out_samples != 90L) ||
      any(resources$disposition != "built_atomic") ||
      any(as.numeric(resources$private_cache_size_bytes) <= 0)) {
    stop("MV6-F input identities or closed-label state are invalid.",
         call. = FALSE)
  }
  candidate_ids <- sort(candidate$sample_id, method = "radix")
  for (seed in 20260805:20260809) {
    observed <- sort(resources$sample_id[resources$seed == seed],
                     method = "radix")
    if (!identical(observed, candidate_ids)) {
      stop("MV6-F cache axis is incomplete for seed ", seed, ".",
           call. = FALSE)
    }
  }
  fold_studies <- sort(unique(folds$held_out_study), method = "radix")
  if (!identical(fold_studies,
                 sort(unique(candidate$study), method = "radix"))) {
    stop("MV6-F fold studies differ from the candidate axis.", call. = FALSE)
  }
  for (study in fold_studies) {
    observed <- sum(candidate$study == study)
    rows <- folds[folds$held_out_study == study, , drop = FALSE]
    if (nrow(rows) != 5L || any(rows$held_out_samples != observed) ||
        any(rows$training_samples != 90L - observed)) {
      stop("MV6-F fold counts do not reconstruct for ", study, ".",
           call. = FALSE)
    }
  }
  panel <- panel[order(panel$panel_order), , drop = FALSE]
  panel_hash <- unique(tolower(panel$panel_sha256))
  if (!identical(as.integer(panel$panel_order), seq_len(500L)) ||
      length(panel_hash) != 1L ||
      panel_hash !=
        "7be22cdb9056fed427c78d58be2b19e258c7c6807e6f7ac3900bd1bfa1d19eb8" ||
      any(!as.logical(panel$selected_all_cache_core)) ||
      anyDuplicated(panel$feature_id)) {
    stop("MV6-F panel differs from the accepted global core.", call. = FALSE)
  }
  invisible(TRUE)
}

mv06f_build_group_queue_v1 <- function(candidate, folds, resources, panel,
                                       input_hashes) {
  mv06f_validate_prefreeze_inputs_v1(candidate, folds, resources, panel)
  required_hashes <- c("candidate", "folds", "resources", "panel")
  if (!all(required_hashes %in% names(input_hashes))) {
    stop("MV6-F queue requires all four input hashes.", call. = FALSE)
  }
  input_hashes <- vapply(input_hashes[required_hashes], .mv06f_required_hash,
                         character(1L), name = "input")
  folds <- folds[order(folds$fold_id, as.integer(folds$seed),
                       method = "radix"), , drop = FALSE]
  rows <- lapply(seq_len(nrow(folds)), function(index) {
    fold <- folds[index, , drop = FALSE]
    identity <- list(
      contract_id = "mv06f_group_identity_v1", fold_id = fold$fold_id,
      fit_scope_id = fold$fit_scope_id,
      held_out_study = fold$held_out_study,
      seed = as.integer(fold$seed),
      training_samples = as.integer(fold$training_samples),
      held_out_samples = as.integer(fold$held_out_samples),
      candidate_sha256 = input_hashes[["candidate"]],
      fold_plan_sha256 = input_hashes[["folds"]],
      resource_sha256 = input_hashes[["resources"]],
      panel_file_sha256 = input_hashes[["panel"]],
      panel_scientific_sha256 =
        "7be22cdb9056fed427c78d58be2b19e258c7c6807e6f7ac3900bd1bfa1d19eb8"
    )
    biological_pairs <- identity$training_samples * identity$held_out_samples
    data.frame(
      contract_id = "mv06f_group_queue_v1",
      group_id = paste0("mv06f_group_v1:", digest::digest(
        identity, algo = "sha256", serialize = TRUE
      )), fold_id = identity$fold_id, fit_scope_id = identity$fit_scope_id,
      held_out_study = identity$held_out_study, seed = identity$seed,
      training_samples = identity$training_samples,
      held_out_samples = identity$held_out_samples,
      cell_ph_jobs = 90L, gene_ph_jobs = 90L,
      diagram_dimension_records = 360L,
      biological_pairs = biological_pairs,
      landscape_component_rows = 4L * biological_pairs,
      outcome_label_state = "closed", biological_outcomes_computed = FALSE,
      fusion_jobs = 0L, clustering_jobs = 0L, outcome_jobs = 0L,
      stringsAsFactors = FALSE
    )
  })
  queue <- do.call(rbind, rows)
  maximum <- max(queue$biological_pairs)
  maximum_folds <- sort(unique(queue$fold_id[
    queue$biological_pairs == maximum
  ]), method = "radix")
  stage_fold <- maximum_folds[[1L]]
  stage_rows <- which(queue$fold_id == stage_fold &
                        queue$biological_pairs == maximum)
  preferred <- stage_rows[queue$seed[stage_rows] == 20260807L]
  stage_index <- if (length(preferred) == 1L) preferred else stage_rows[[1L]]
  queue$stage <- ifelse(seq_len(nrow(queue)) == stage_index,
                        "stage_1_maximum", "stage_2")
  remaining <- setdiff(seq_len(nrow(queue)), stage_index)
  remaining <- remaining[order(-queue$biological_pairs[remaining],
                               queue$fold_id[remaining], queue$seed[remaining],
                               method = "radix")]
  queue$execution_order <- NA_integer_
  queue$execution_order[stage_index] <- 1L
  queue$execution_order[remaining] <- seq_along(remaining) + 1L
  queue <- queue[order(queue$execution_order), , drop = FALSE]
  rownames(queue) <- NULL
  if (nrow(queue) != 75L || sum(queue$cell_ph_jobs) != 6750L ||
      sum(queue$gene_ph_jobs) != 6750L ||
      sum(queue$diagram_dimension_records) != 27000L ||
      sum(queue$biological_pairs) != 35350L ||
      sum(queue$landscape_component_rows) != 141400L ||
      sum(queue$stage == "stage_1_maximum") != 1L ||
      anyDuplicated(queue$group_id)) {
    stop("MV6-F queue totals differ from the prospective contract.",
         call. = FALSE)
  }
  root_fields <- queue[, c(
    "group_id", "fold_id", "fit_scope_id", "held_out_study", "seed",
    "training_samples", "held_out_samples", "biological_pairs",
    "landscape_component_rows", "stage", "execution_order"
  )]
  attr(queue, "queue_root_sha256") <- digest::digest(
    root_fields, algo = "sha256", serialize = TRUE
  )
  queue
}

mv06f_validate_queue_v1 <- function(queue) {
  required <- c(
    "contract_id", "group_id", "fold_id", "fit_scope_id",
    "held_out_study", "seed", "training_samples", "held_out_samples",
    "cell_ph_jobs", "gene_ph_jobs", "diagram_dimension_records",
    "biological_pairs", "landscape_component_rows", "stage",
    "execution_order", "outcome_label_state",
    "biological_outcomes_computed", "fusion_jobs", "clustering_jobs",
    "outcome_jobs"
  )
  if (!is.data.frame(queue) || nrow(queue) != 75L ||
      !all(required %in% names(queue)) ||
      any(queue$contract_id != "mv06f_group_queue_v1") ||
      anyDuplicated(queue$group_id) ||
      !identical(as.integer(queue$execution_order), seq_len(75L)) ||
      sum(queue$stage == "stage_1_maximum") != 1L ||
      sum(queue$cell_ph_jobs) != 6750L || sum(queue$gene_ph_jobs) != 6750L ||
      sum(queue$diagram_dimension_records) != 27000L ||
      sum(queue$biological_pairs) != 35350L ||
      sum(queue$landscape_component_rows) != 141400L ||
      any(queue$outcome_label_state != "closed") ||
      any(as.logical(queue$biological_outcomes_computed)) ||
      any(queue$fusion_jobs != 0L) || any(queue$clustering_jobs != 0L) ||
      any(queue$outcome_jobs != 0L)) {
    stop("MV6-F public queue is stale or crosses its stop boundary.",
         call. = FALSE)
  }
  invisible(queue)
}

mv06f_queue_root_v1 <- function(queue) {
  mv06f_validate_queue_v1(queue)
  root_fields <- queue[, c(
    "group_id", "fold_id", "fit_scope_id", "held_out_study", "seed",
    "training_samples", "held_out_samples", "biological_pairs",
    "landscape_component_rows", "stage", "execution_order"
  )]
  digest::digest(root_fields, algo = "sha256", serialize = TRUE)
}

mv06f_new_ph_record_v1 <- function(group_key, sample_id, role, view, result) {
  oracle <- mv06d_validate_ph_result_v1(result, view)
  identity <- list(
    contract_id = "mv06f_matched_ph_identity_v1",
    group_key = as.character(group_key), sample_id = as.character(sample_id),
    role = match.arg(role, c("held_out", "training")),
    seed = as.integer(view$subsample_seed), view_id = view$view_id,
    view_cache_key = view$cache_key, view_payload_sha256 = view$payload_sha256,
    max_dim = 1L, threshold = -1, field = 2L,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE
  )
  payload <- list(identity = identity, topology_result = result,
                  h0_mst_oracle = oracle)
  payload_sha <- digest::digest(payload, algo = "sha256", serialize = TRUE)
  structure(list(
    contract_id = "mv06f_matched_ph_record_v1", identity = identity,
    topology_result = result, h0_mst_oracle = oracle,
    payload_sha256 = payload_sha,
    cache_key = paste0("mv06f_matched_ph_record_v1:", payload_sha),
    downstream_execution = list(
      landscape_pairs = 0L, fusion_jobs = 0L, clustering_jobs = 0L,
      outcome_jobs = 0L
    )
  ), class = c("scph_mv06f_ph_record_v1", "list"))
}

mv06f_validate_ph_record_v1 <- function(record) {
  payload <- list(identity = record$identity,
                  topology_result = record$topology_result,
                  h0_mst_oracle = record$h0_mst_oracle)
  if (!inherits(record, "scph_mv06f_ph_record_v1") ||
      !identical(record$contract_id, "mv06f_matched_ph_record_v1") ||
      !identical(record$payload_sha256, digest::digest(
        payload, algo = "sha256", serialize = TRUE
      )) || !identical(record$cache_key,
                       paste0("mv06f_matched_ph_record_v1:",
                              record$payload_sha256)) ||
      record$identity$outcome_label_state != "closed" ||
      isTRUE(record$identity$biological_outcomes_computed) ||
      any(unlist(record$downstream_execution, use.names = FALSE) != 0)) {
    stop("MV6-F PH record is stale or crossed its stop boundary.",
         call. = FALSE)
  }
  invisible(record)
}

mv06f_finite_intervals_v1 <- function(record, dimension) {
  mv06f_validate_ph_record_v1(record)
  dimension <- match.arg(dimension, c("H0", "H1"))
  value <- as.integer(sub("H", "", dimension, fixed = TRUE))
  diagram <- record$topology_result$diagram
  selected <- diagram[
    diagram[, "dimension"] == value & is.finite(diagram[, "death"]) &
      diagram[, "death"] > diagram[, "birth"], c("birth", "death"),
    drop = FALSE
  ]
  storage.mode(selected) <- "double"
  selected[order(selected[, "birth"], -selected[, "death"]), , drop = FALSE]
}

mv06f_pair_id_v1 <- function(group_id, query_sample_id, training_sample_id,
                             view_id, homology_dimension) {
  identity <- list(
    contract_id = "mv06f_landscape_pair_identity_v1",
    group_id = as.character(group_id),
    query_sample_id = as.character(query_sample_id),
    training_sample_id = as.character(training_sample_id),
    view_id = match.arg(view_id, c("cell_topology_v1", "gene_topology_v1")),
    homology_dimension = match.arg(homology_dimension, c("H0", "H1"))
  )
  paste0("mv06f_pair_v1:", digest::digest(
    identity, algo = "sha256", serialize = TRUE
  ))
}

mv06f_validate_group_directory_v1 <- function(path, queue_row,
                                               queue_root_sha256,
                                               implementation_root_sha256,
                                               rust_library_sha256) {
  if (!dir.exists(path) || !is.data.frame(queue_row) || nrow(queue_row) != 1L) {
    stop("MV6-F completed group directory or queue row is absent.",
         call. = FALSE)
  }
  required_files <- c("diagrams.rds", "diagram-manifest.csv",
                      "distances.csv", "metrics.csv", "status.csv")
  paths <- file.path(path, required_files)
  if (!all(file.exists(paths))) {
    stop("MV6-F completed group has a one-sided artifact set.",
         call. = FALSE)
  }
  status <- utils::read.csv(paths[[5L]], stringsAsFactors = FALSE,
                            check.names = FALSE)
  expected_status <- c(
    "contract_id", "group_id", "queue_root_sha256",
    "implementation_root_sha256", "rust_library_sha256", "completion_state",
    "diagrams_sha256", "diagrams_bytes", "diagram_manifest_sha256",
    "diagram_manifest_bytes", "distances_sha256", "distances_bytes",
    "metrics_sha256", "metrics_bytes", "ph_jobs",
    "diagram_dimension_records", "biological_pairs",
    "landscape_component_rows", "outcome_label_state",
    "biological_outcomes_computed", "fusion_jobs", "clustering_jobs",
    "outcome_jobs"
  )
  if (nrow(status) != 1L || !all(expected_status %in% names(status)) ||
      status$contract_id != "mv06f_group_status_v1" ||
      status$group_id != queue_row$group_id ||
      status$queue_root_sha256 != queue_root_sha256 ||
      status$implementation_root_sha256 != implementation_root_sha256 ||
      status$rust_library_sha256 != rust_library_sha256 ||
      status$completion_state != "complete" ||
      status$outcome_label_state != "closed" ||
      as.logical(status$biological_outcomes_computed) ||
      any(unlist(status[c("fusion_jobs", "clustering_jobs", "outcome_jobs")]) !=
            0)) {
    stop("MV6-F completed group status is stale.", call. = FALSE)
  }
  hashes <- vapply(paths[1:4], .mv06f_sha256, character(1L))
  bytes <- unname(file.info(paths[1:4])$size)
  expected_hashes <- as.character(unlist(status[c(
    "diagrams_sha256", "diagram_manifest_sha256", "distances_sha256",
    "metrics_sha256"
  )], use.names = FALSE))
  if (!identical(unname(hashes), expected_hashes) ||
      !identical(as.numeric(bytes), as.numeric(unlist(status[c(
    "diagrams_bytes", "diagram_manifest_bytes", "distances_bytes",
    "metrics_bytes"
  )], use.names = FALSE)))) {
    stop("MV6-F completed group artifact hashes or sizes are stale.",
         call. = FALSE)
  }
  records <- readRDS(paths[[1L]])
  if (!is.list(records) || length(records) != 180L) {
    stop("MV6-F completed group does not contain 180 PH records.",
         call. = FALSE)
  }
  invisible(lapply(records, mv06f_validate_ph_record_v1))
  manifest <- utils::read.csv(paths[[2L]], stringsAsFactors = FALSE,
                              check.names = FALSE)
  distances <- utils::read.csv(paths[[3L]], stringsAsFactors = FALSE,
                               check.names = FALSE)
  if (nrow(manifest) != 180L || nrow(distances) !=
      queue_row$landscape_component_rows || anyDuplicated(manifest$ph_cache_key) ||
      anyDuplicated(distances$pair_id) ||
      any(!is.finite(distances$squared_distance)) ||
      any(distances$squared_distance < 0) ||
      any(!as.logical(distances$exact)) ||
      any(!as.logical(distances$all_active_levels)) ||
      any(as.logical(distances$level_cap_applied)) ||
      any(distances$rust_status != 0L) ||
      any(manifest$outcome_label_state != "closed") ||
      any(as.logical(manifest$biological_outcomes_computed)) ||
      any(distances$outcome_label_state != "closed") ||
      any(as.logical(distances$biological_outcomes_computed)) ||
      status$ph_jobs != 180L || status$diagram_dimension_records != 360L ||
      status$biological_pairs != queue_row$biological_pairs ||
      status$landscape_component_rows != queue_row$landscape_component_rows) {
    stop("MV6-F completed group scientific artifacts are invalid.",
         call. = FALSE)
  }
  invisible(status)
}
