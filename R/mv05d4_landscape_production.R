# MV5-D4 label-closed exact all-level cell-landscape contracts.

.mv05d4_hash <- function(value, label) {
  value <- as.character(value)
  if (length(value) != 1L || is.na(value) ||
      !grepl("^[0-9a-f]{64}$", value)) {
    stop(label, " must be one lowercase SHA-256 string.", call. = FALSE)
  }
  value
}

.mv05d4_pair_identity <- function(row) {
  list(
    contract_id = "mv05d4_cell_landscape_pair_v1",
    group_id = row$group_id,
    fold_id = row$fold_id,
    seed = as.integer(row$seed),
    homology_dimension = row$homology_dimension,
    query_job_id = row$query_job_id,
    training_job_id = row$training_job_id,
    query_record_cache_key = row$query_record_cache_key,
    training_record_cache_key = row$training_record_cache_key,
    query_diagram_sha256 = row$query_diagram_sha256,
    training_diagram_sha256 = row$training_diagram_sha256,
    query_result_file_sha256 = row$query_result_file_sha256,
    training_result_file_sha256 = row$training_result_file_sha256,
    landscape_definition_id = "all_active_exact_critical_pairs_v1"
  )
}

mv05d4_build_pair_manifest_v1 <- function(ph_manifest, view_metrics,
                                           max_pairs = 250L) {
  mv05d3_validate_full_manifest_v1(ph_manifest)
  mv05d3_validate_view_metrics_v1(view_metrics)
  if (!setequal(ph_manifest$job_id, view_metrics$job_id)) {
    stop("MV5-D4 inputs do not identify the same 6,750 views.", call. = FALSE)
  }
  metrics <- view_metrics[match(ph_manifest$job_id, view_metrics$job_id), ]
  rows <- list()
  groups <- split(seq_len(nrow(ph_manifest)), ph_manifest$group_id)
  for (group_id in names(groups)) {
    indices <- groups[[group_id]]
    group <- ph_manifest[indices, , drop = FALSE]
    group_metrics <- metrics[indices, , drop = FALSE]
    query <- which(group$execution_role == "held_out")
    training <- which(group$execution_role == "training")
    if (!length(query) || !length(training)) {
      stop("MV5-D4 group lacks held-out or training views.", call. = FALSE)
    }
    for (query_index in query) for (training_index in training) {
      for (dimension in c("H0", "H1")) {
        row <- data.frame(
          contract_id = "mv05d4_cell_landscape_pair_v1",
          group_id = group_id, group_order = group$group_order[[1L]],
          fold_id = group$fold_id[[1L]], seed = group$seed[[1L]],
          homology_dimension = dimension,
          query_sample_id = group$sample_id[[query_index]],
          training_sample_id = group$sample_id[[training_index]],
          query_job_id = group$job_id[[query_index]],
          training_job_id = group$job_id[[training_index]],
          query_record_cache_key =
            group_metrics$record_cache_key[[query_index]],
          training_record_cache_key =
            group_metrics$record_cache_key[[training_index]],
          query_diagram_sha256 = group_metrics$diagram_sha256[[query_index]],
          training_diagram_sha256 =
            group_metrics$diagram_sha256[[training_index]],
          query_result_file = group_metrics$result_file[[query_index]],
          training_result_file = group_metrics$result_file[[training_index]],
          query_result_file_sha256 =
            group_metrics$result_file_sha256[[query_index]],
          training_result_file_sha256 =
            group_metrics$result_file_sha256[[training_index]],
          representation = "sct_whole", view_id = "cell_topology_v1",
          pair_scope = "held_out_query_to_training_reference",
          landscape_definition_id = "all_active_exact_critical_pairs_v1",
          essential_h0_policy = "exclude",
          level_policy = "all_consecutive_active_levels",
          integration_policy = "exact_linear_critical_pair_segments",
          supports_primary_retrieval = TRUE,
          supports_full_matrix_clustering = FALSE,
          supports_within_study_pair_contrasts = FALSE,
          outcome_label_state = "closed",
          biological_outcomes_computed = FALSE,
          stringsAsFactors = FALSE
        )
        row$pair_request_id <- paste0(
          "mv05d4_pair_v1:", digest::digest(
            .mv05d4_pair_identity(row), algo = "sha256", serialize = TRUE
          )
        )
        rows[[length(rows) + 1L]] <- row
      }
    }
  }
  result <- do.call(rbind, rows)
  result <- result[order(
    result$group_order, result$homology_dimension, result$pair_request_id,
    method = "radix"
  ), , drop = FALSE]
  locality <- paste(result$group_id, result$homology_dimension, sep = "\r")
  within <- stats::ave(seq_len(nrow(result)), locality, FUN = seq_along)
  result$local_chunk_index <- (within - 1L) %/% as.integer(max_pairs) + 1L
  chunk_group <- paste(locality, result$local_chunk_index, sep = "\r")
  result$chunk_index <- match(chunk_group, unique(chunk_group))
  ids <- split(result$pair_request_id, result$chunk_index)
  chunk_ids <- vapply(ids, function(value) paste0(
    "mv05d4_chunk_v1:",
    digest::digest(value, algo = "sha256", serialize = TRUE)
  ), character(1L))
  result$chunk_id <- unname(chunk_ids[as.character(result$chunk_index)])
  result$chunk_offset <- stats::ave(
    seq_len(nrow(result)), result$chunk_index, FUN = seq_along
  )
  rownames(result) <- NULL
  mv05d4_validate_pair_manifest_v1(result, max_pairs = max_pairs)
  result
}

mv05d4_validate_pair_manifest_v1 <- function(manifest, max_pairs = 250L) {
  required <- c(
    "contract_id", "pair_request_id", "group_id", "group_order", "fold_id",
    "seed", "homology_dimension", "query_sample_id", "training_sample_id",
    "query_job_id", "training_job_id", "query_record_cache_key",
    "training_record_cache_key", "query_diagram_sha256",
    "training_diagram_sha256", "query_result_file_sha256",
    "training_result_file_sha256", "representation", "view_id", "pair_scope",
    "landscape_definition_id", "essential_h0_policy", "level_policy",
    "integration_policy", "chunk_index", "chunk_id", "chunk_offset",
    "outcome_label_state", "biological_outcomes_computed"
  )
  hash_fields <- c(
    "query_diagram_sha256", "training_diagram_sha256",
    "query_result_file_sha256", "training_result_file_sha256"
  )
  if (!is.data.frame(manifest) || !all(required %in% names(manifest)) ||
      nrow(manifest) != 70700L ||
      any(manifest$contract_id != "mv05d4_cell_landscape_pair_v1") ||
      anyDuplicated(manifest$pair_request_id) ||
      length(unique(manifest$group_id)) != 75L ||
      !setequal(manifest$homology_dimension, c("H0", "H1")) ||
      any(manifest$query_sample_id == manifest$training_sample_id) ||
      any(manifest$representation != "sct_whole") ||
      any(manifest$view_id != "cell_topology_v1") ||
      any(manifest$pair_scope !=
            "held_out_query_to_training_reference") ||
      any(manifest$landscape_definition_id !=
            "all_active_exact_critical_pairs_v1") ||
      any(manifest$essential_h0_policy != "exclude") ||
      any(manifest$level_policy != "all_consecutive_active_levels") ||
      any(manifest$integration_policy !=
            "exact_linear_critical_pair_segments") ||
      any(manifest$outcome_label_state != "closed") ||
      any(as.logical(manifest$biological_outcomes_computed)) ||
      any(c("tissue", "approach") %in% names(manifest)) ||
      any(!vapply(unlist(manifest[hash_fields], use.names = FALSE),
                  function(x) grepl("^[0-9a-f]{64}$", x), logical(1L))) ||
      any(table(manifest$chunk_id) > as.integer(max_pairs))) {
    stop("MV5-D4 pair manifest violates its frozen contract.", call. = FALSE)
  }
  pair_key <- paste(
    manifest$group_id, manifest$query_job_id, manifest$training_job_id,
    sep = "\r"
  )
  dimension_counts <- table(pair_key)
  if (length(dimension_counts) != 35350L || any(dimension_counts != 2L) ||
      any(table(manifest$group_id, manifest$homology_dimension)[, "H0"] !=
          table(manifest$group_id, manifest$homology_dimension)[, "H1"])) {
    stop("MV5-D4 does not contain exactly two dimensions per eligible pair.",
         call. = FALSE)
  }
  invisible(manifest)
}

mv05d4_chunk_summary_v1 <- function(manifest) {
  mv05d4_validate_pair_manifest_v1(manifest)
  groups <- split(manifest, manifest$chunk_index)
  result <- do.call(rbind, lapply(groups, function(chunk) data.frame(
    contract_id = "mv05d4_chunk_manifest_v1",
    chunk_index = chunk$chunk_index[[1L]],
    chunk_id = chunk$chunk_id[[1L]],
    group_id = chunk$group_id[[1L]], group_order = chunk$group_order[[1L]],
    fold_id = chunk$fold_id[[1L]], seed = chunk$seed[[1L]],
    homology_dimension = chunk$homology_dimension[[1L]],
    pair_count = nrow(chunk),
    first_pair_request_id = chunk$pair_request_id[[1L]],
    last_pair_request_id = chunk$pair_request_id[[nrow(chunk)]],
    execution_disposition = "pending_or_resume_by_chunk_identity",
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )))
  rownames(result) <- NULL
  result
}
