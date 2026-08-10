# Locked paired SCT-versus-integrated retrieval comparison for MV5-L.

.mv05l_contract_id <- "mv05l_locked_representation_comparison_v1"
.mv05l_freeze_commit <- "b3f7e280181f4aba4de26cfdfad5e637605e31e1"
.mv05l_seeds <- 20260805:20260809
.mv05l_endpoints <- c(
  "cross_study_tissue_mrr_v1",
  "cross_study_tissue_1nn_balanced_accuracy_v1"
)
.mv05l_estimands <- c(
  "did_h0_topology_minus_energy",
  "did_h1_topology_minus_energy",
  "direct_h0_integrated_minus_sct",
  "direct_h1_integrated_minus_sct",
  "direct_raw_composite_integrated_minus_sct",
  "direct_energy_integrated_minus_sct",
  "direct_pseudobulk_integrated_minus_sct"
)
.mv05l_estimand_roles <- c(
  did_h0_topology_minus_energy = "primary_difference_in_differences",
  did_h1_topology_minus_energy = "primary_difference_in_differences",
  direct_h0_integrated_minus_sct = "secondary_direct_representation",
  direct_h1_integrated_minus_sct = "secondary_direct_representation",
  direct_raw_composite_integrated_minus_sct = "descriptive_direct_representation",
  direct_energy_integrated_minus_sct = "secondary_direct_representation",
  direct_pseudobulk_integrated_minus_sct = "pipeline_identity_control"
)
.mv05l_method_map <- data.frame(
  family_id = c("h0", "h1", "raw_composite", "energy", "pseudobulk"),
  sct_method_id = c(
    "cell_landscape_h0_v1",
    "cell_landscape_h1_v1",
    "cell_landscape_h0_h1_raw_euclidean_v1",
    "cell_distribution_energy_shared_pca_v1",
    "pseudobulk_shared_panel_euclidean_v1"
  ),
  integrated_method_id = c(
    "integrated_cell_landscape_h0_v1",
    "integrated_cell_landscape_h1_v1",
    "integrated_cell_landscape_h0_h1_raw_euclidean_v1",
    "integrated_cell_distribution_energy_v1",
    "pseudobulk_training_standardized_panel_v1"
  ),
  analysis_role = c(
    "primary_topology_component", "primary_topology_component",
    "descriptive_only", "matched_within_representation_baseline",
    "shared_context_and_identity_control"
  ),
  stringsAsFactors = FALSE
)

.mv05l_expected_hashes <- c(
  mv05e_query_endpoints = "b65dfd0bae4889fd6297372e7ea180ec8f8465487f68e3c6a5bafd4b9920565f",
  mv05e_prediction_lock = "ed2ff4e8d5c52ddabb527853ba5ff60f81dc6017f79bbea068bdcd07939cd053",
  mv05e_independent_validation = "0f8173d57e1bef27311d8e655de867e2b694c6e39c05e7253a52361c3eb2561b",
  mv05e_deterministic_repeat = "8974196f2e28f712b3b726181b7ebdcc056b7ef8cf779bc4c53ac362e9e7bd07",
  mv05e_artifact_manifest = "2a2401ac30a52067530139fcfb4c6d642d681da3e4ada5a8270db13edf97329a",
  mv05k_query_endpoints = "8fa9f70a64f35b87c5e408d32bd7fc4440d5df611c87cd201ff68968b959df2b",
  mv05k_prediction_lock = "e19334a5adfc2fed09dfd0297958f8a837ddbf5ff6753d97ba8e4ea60696ae56",
  mv05k_independent_validation = "bd80254875753b4488b7488545f348ac7ccf478a1b9f2b89a4cc87aefeffbfd2",
  mv05k_deterministic_repeat = "296ed828e30436782b6053c8ca5a35b71e844d1bd6387ac9ef02c0b10d205f49",
  mv05k_artifact_manifest = "983b6be20b1a5858d14b4d08959f382fb09453b64f08c70ed3931e8d10c51a06"
)

.mv05l_file_sha256 <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

.mv05l_with_seed <- function(seed, expression) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  force(expression)
}

.mv05l_boundary_exceedances <- function(null, observed) {
  tolerance <- 64 * .Machine$double.eps *
    pmax(1, abs(null), abs(observed))
  sum(abs(null) + tolerance >= abs(observed))
}

mv05l_verify_input_lock_v1 <- function(paths, source_commit,
                                       expected_hashes = .mv05l_expected_hashes,
                                       expected_commit = .mv05l_freeze_commit) {
  if (!identical(as.character(source_commit), expected_commit)) {
    stop("MV5-L must execute from the frozen pre-join commit.", call. = FALSE)
  }
  if (!is.character(paths) || is.null(names(paths)) ||
      !identical(sort(names(paths)), sort(names(expected_hashes)))) {
    stop("MV5-L input paths do not match the frozen artifact registry.",
         call. = FALSE)
  }
  paths <- paths[names(expected_hashes)]
  if (any(!file.exists(paths))) {
    stop("One or more locked MV5-L inputs are absent.", call. = FALSE)
  }
  observed <- vapply(paths, .mv05l_file_sha256, character(1L))
  if (!identical(unname(observed), unname(expected_hashes))) {
    bad <- names(expected_hashes)[observed != expected_hashes]
    stop("MV5-L immutable input hash mismatch: ", paste(bad, collapse = ", "),
         call. = FALSE)
  }
  data.frame(
    contract_id = "mv05l_input_lock_v1",
    comparison_freeze_commit = expected_commit,
    input_id = names(expected_hashes),
    input_file = basename(unname(paths)),
    expected_sha256 = unname(expected_hashes),
    observed_sha256 = unname(observed),
    lock_passed_before_endpoint_read = TRUE,
    marginal_aggregate_outcomes_known_at_specification = TRUE,
    joint_sample_contrasts_known_at_specification = FALSE,
    stringsAsFactors = FALSE
  )
}

.mv05l_require_endpoint_schema <- function(x, representation) {
  required <- c(
    "contract_id", "fold_id", "held_out_study", "seed", "method_id",
    "query_sample_id", "query_tissue", "training_samples",
    "training_studies", "first_same_tissue_sample_id",
    "first_same_tissue_rank", "reciprocal_rank", "nearest_sample_id",
    "nearest_tissue", "one_nn_correct", "nearest_distance_tied",
    "endpoint_status", "labels_opened_after_prediction_lock",
    "upstream_refit", "reranked_after_label_open"
  )
  if (!is.data.frame(x) || !all(required %in% names(x))) {
    stop("Incomplete ", representation, " endpoint schema.", call. = FALSE)
  }
  if (any(x$endpoint_status != "estimable") || any(x$upstream_refit) ||
      any(x$reranked_after_label_open) ||
      any(!x$labels_opened_after_prediction_lock)) {
    stop("The ", representation, " endpoint artifact is not accepted.",
         call. = FALSE)
  }
  invisible(TRUE)
}

mv05l_pair_locked_endpoints_v1 <- function(
    sct, integrated, expected_rows_per_family = 450L,
    require_pseudobulk_identity = TRUE) {
  .mv05l_require_endpoint_schema(sct, "SCT")
  .mv05l_require_endpoint_schema(integrated, "integrated")
  keys <- c("fold_id", "held_out_study", "seed", "query_sample_id",
            "query_tissue")
  common <- c(
    "training_samples", "training_studies", "endpoint_status",
    "labels_opened_after_prediction_lock", "upstream_refit",
    "reranked_after_label_open"
  )
  pieces <- lapply(seq_len(nrow(.mv05l_method_map)), function(index) {
    map <- .mv05l_method_map[index, , drop = FALSE]
    left <- sct[sct$method_id == map$sct_method_id, , drop = FALSE]
    right <- integrated[
      integrated$method_id == map$integrated_method_id, , drop = FALSE
    ]
    if (nrow(left) != expected_rows_per_family ||
        nrow(right) != expected_rows_per_family) {
      stop("Incomplete endpoint family: ", map$family_id, call. = FALSE)
    }
    left <- left[order(left$fold_id, left$seed, left$query_sample_id,
                       method = "radix"), , drop = FALSE]
    right <- right[order(right$fold_id, right$seed, right$query_sample_id,
                         method = "radix"), , drop = FALSE]
    if (!all(vapply(keys, function(key) {
      identical(left[[key]], right[[key]])
    }, logical(1L)))) {
      stop("Cross-representation endpoint identities do not align for ",
           map$family_id, ".", call. = FALSE)
    }
    for (field in common) {
      if (!identical(left[[field]], right[[field]])) {
        stop("Cross-representation field mismatch for ", map$family_id,
             ": ", field, call. = FALSE)
      }
    }
    data.frame(
      contract_id = "mv05l_paired_query_endpoint_v1",
      family_id = map$family_id,
      analysis_role = map$analysis_role,
      fold_id = left$fold_id,
      held_out_study = left$held_out_study,
      seed = as.integer(left$seed),
      query_sample_id = left$query_sample_id,
      query_tissue = left$query_tissue,
      training_samples = as.integer(left$training_samples),
      training_studies = as.integer(left$training_studies),
      sct_method_id = left$method_id,
      integrated_method_id = right$method_id,
      sct_first_same_tissue_sample_id = left$first_same_tissue_sample_id,
      integrated_first_same_tissue_sample_id =
        right$first_same_tissue_sample_id,
      sct_first_same_tissue_rank = as.integer(left$first_same_tissue_rank),
      integrated_first_same_tissue_rank =
        as.integer(right$first_same_tissue_rank),
      sct_reciprocal_rank = as.numeric(left$reciprocal_rank),
      integrated_reciprocal_rank = as.numeric(right$reciprocal_rank),
      direct_reciprocal_rank_difference =
        as.numeric(right$reciprocal_rank - left$reciprocal_rank),
      sct_nearest_sample_id = left$nearest_sample_id,
      integrated_nearest_sample_id = right$nearest_sample_id,
      sct_nearest_tissue = left$nearest_tissue,
      integrated_nearest_tissue = right$nearest_tissue,
      sct_one_nn_correct = as.logical(left$one_nn_correct),
      integrated_one_nn_correct = as.logical(right$one_nn_correct),
      direct_one_nn_difference = as.numeric(right$one_nn_correct) -
        as.numeric(left$one_nn_correct),
      sct_nearest_distance_tied = as.logical(left$nearest_distance_tied),
      integrated_nearest_distance_tied =
        as.logical(right$nearest_distance_tied),
      endpoint_status = left$endpoint_status,
      stringsAsFactors = FALSE
    )
  })
  paired <- do.call(rbind, pieces)
  paired$family_id <- factor(paired$family_id,
                             levels = .mv05l_method_map$family_id)
  paired <- paired[order(
    paired$family_id, paired$fold_id, paired$seed, paired$query_sample_id,
    method = "radix"
  ), , drop = FALSE]
  paired$family_id <- as.character(paired$family_id)
  rownames(paired) <- NULL
  pseudo <- paired[paired$family_id == "pseudobulk", , drop = FALSE]
  pseudo_identical <- all(
    pseudo$sct_first_same_tissue_sample_id ==
      pseudo$integrated_first_same_tissue_sample_id,
    pseudo$sct_first_same_tissue_rank ==
      pseudo$integrated_first_same_tissue_rank,
    pseudo$sct_reciprocal_rank == pseudo$integrated_reciprocal_rank,
    pseudo$sct_nearest_sample_id == pseudo$integrated_nearest_sample_id,
    pseudo$sct_nearest_tissue == pseudo$integrated_nearest_tissue,
    pseudo$sct_one_nn_correct == pseudo$integrated_one_nn_correct,
    pseudo$sct_nearest_distance_tied ==
      pseudo$integrated_nearest_distance_tied
  )
  if (require_pseudobulk_identity && !pseudo_identical) {
    stop("The required shared-pseudobulk identity control failed.",
         call. = FALSE)
  }
  attr(paired, "pseudobulk_identical") <- pseudo_identical
  paired
}

mv05l_build_compatibility_v1 <- function(paired) {
  do.call(rbind, lapply(.mv05l_method_map$family_id, function(family) {
    part <- paired[paired$family_id == family, , drop = FALSE]
    data.frame(
      contract_id = "mv05l_endpoint_compatibility_v1",
      family_id = family,
      paired_query_seed_rows = nrow(part),
      samples = length(unique(part$query_sample_id)),
      held_out_studies = length(unique(part$held_out_study)),
      tissues = length(unique(part$query_tissue)),
      seeds = length(unique(part$seed)),
      identity_mismatches = 0L,
      denominator_mismatches = 0L,
      nonestimable_pairs = sum(part$endpoint_status != "estimable"),
      pairing_status = "complete",
      stringsAsFactors = FALSE
    )
  }))
}

mv05l_build_identity_control_v1 <- function(paired) {
  part <- paired[paired$family_id == "pseudobulk", , drop = FALSE]
  checks <- c(
    first_same_tissue_sample_id = all(
      part$sct_first_same_tissue_sample_id ==
        part$integrated_first_same_tissue_sample_id),
    first_same_tissue_rank = all(
      part$sct_first_same_tissue_rank ==
        part$integrated_first_same_tissue_rank),
    reciprocal_rank = all(
      part$sct_reciprocal_rank == part$integrated_reciprocal_rank),
    nearest_sample_id = all(
      part$sct_nearest_sample_id == part$integrated_nearest_sample_id),
    nearest_tissue = all(
      part$sct_nearest_tissue == part$integrated_nearest_tissue),
    one_nn_correct = all(
      part$sct_one_nn_correct == part$integrated_one_nn_correct),
    nearest_distance_tied = all(
      part$sct_nearest_distance_tied ==
        part$integrated_nearest_distance_tied)
  )
  data.frame(
    contract_id = "mv05l_pseudobulk_identity_control_v1",
    check_id = names(checks),
    rows_checked = nrow(part),
    mismatches = vapply(names(checks), function(name) {
      switch(name,
        first_same_tissue_sample_id = sum(
          part$sct_first_same_tissue_sample_id !=
            part$integrated_first_same_tissue_sample_id),
        first_same_tissue_rank = sum(
          part$sct_first_same_tissue_rank !=
            part$integrated_first_same_tissue_rank),
        reciprocal_rank = sum(
          part$sct_reciprocal_rank != part$integrated_reciprocal_rank),
        nearest_sample_id = sum(
          part$sct_nearest_sample_id != part$integrated_nearest_sample_id),
        nearest_tissue = sum(
          part$sct_nearest_tissue != part$integrated_nearest_tissue),
        one_nn_correct = sum(
          part$sct_one_nn_correct != part$integrated_one_nn_correct),
        nearest_distance_tied = sum(
          part$sct_nearest_distance_tied !=
            part$integrated_nearest_distance_tied)
      )
    }, integer(1L)),
    exact_identity_passed = unname(checks),
    expected_relationship = "shared_pseudobulk_endpoint_identity",
    stringsAsFactors = FALSE
  )
}

.mv05l_query_estimands <- function(paired) {
  keys <- c("fold_id", "held_out_study", "seed", "query_sample_id",
            "query_tissue")
  direct_map <- c(
    h0 = "direct_h0_integrated_minus_sct",
    h1 = "direct_h1_integrated_minus_sct",
    raw_composite = "direct_raw_composite_integrated_minus_sct",
    energy = "direct_energy_integrated_minus_sct",
    pseudobulk = "direct_pseudobulk_integrated_minus_sct"
  )
  direct <- do.call(rbind, lapply(names(direct_map), function(family) {
    part <- paired[paired$family_id == family, , drop = FALSE]
    rbind(
      data.frame(part[keys], endpoint_id = .mv05l_endpoints[[1L]],
                 estimand_id = direct_map[[family]],
                 value = part$direct_reciprocal_rank_difference,
                 stringsAsFactors = FALSE),
      data.frame(part[keys], endpoint_id = .mv05l_endpoints[[2L]],
                 estimand_id = direct_map[[family]],
                 value = part$direct_one_nn_difference,
                 stringsAsFactors = FALSE)
    )
  }))
  energy <- paired[paired$family_id == "energy", , drop = FALSE]
  did <- do.call(rbind, lapply(c(h0 = "did_h0_topology_minus_energy",
                                 h1 = "did_h1_topology_minus_energy"),
    function(estimand) {
      family <- if (grepl("h0", estimand, fixed = TRUE)) "h0" else "h1"
      topology <- paired[paired$family_id == family, , drop = FALSE]
      topology <- topology[order(topology$fold_id, topology$seed,
                                 topology$query_sample_id, method = "radix"), ]
      energy_ordered <- energy[order(energy$fold_id, energy$seed,
                                     energy$query_sample_id, method = "radix"), ]
      if (!all(vapply(keys, function(key) {
        identical(topology[[key]], energy_ordered[[key]])
      }, logical(1L)))) {
        stop("Topology and energy query identities do not align.", call. = FALSE)
      }
      rbind(
        data.frame(topology[keys], endpoint_id = .mv05l_endpoints[[1L]],
          estimand_id = estimand,
          value = topology$direct_reciprocal_rank_difference -
            energy_ordered$direct_reciprocal_rank_difference,
          stringsAsFactors = FALSE),
        data.frame(topology[keys], endpoint_id = .mv05l_endpoints[[2L]],
          estimand_id = estimand,
          value = topology$direct_one_nn_difference -
            energy_ordered$direct_one_nn_difference,
          stringsAsFactors = FALSE)
      )
    }))
  result <- rbind(did, direct)
  result$endpoint_id <- factor(result$endpoint_id, levels = .mv05l_endpoints)
  result$estimand_id <- factor(result$estimand_id, levels = .mv05l_estimands)
  result <- result[order(result$endpoint_id, result$estimand_id,
                         result$query_tissue, result$held_out_study,
                         result$query_sample_id, result$seed,
                         method = "radix"), , drop = FALSE]
  result$endpoint_id <- as.character(result$endpoint_id)
  result$estimand_id <- as.character(result$estimand_id)
  rownames(result) <- NULL
  result
}

mv05l_summarize_estimands_v1 <- function(paired) {
  query <- .mv05l_query_estimands(paired)
  sample <- stats::aggregate(
    value ~ endpoint_id + estimand_id + query_sample_id + query_tissue +
      held_out_study,
    data = query, FUN = mean
  )
  counts <- stats::aggregate(
    seed ~ endpoint_id + estimand_id + query_sample_id + query_tissue +
      held_out_study,
    data = query, FUN = length
  )
  names(counts)[names(counts) == "seed"] <- "completed_seeds"
  sample <- merge(sample, counts, sort = FALSE)
  sample$contract_id <- "mv05l_sample_estimand_v1"
  sample$estimand_role <- unname(.mv05l_estimand_roles[sample$estimand_id])
  names(sample)[names(sample) == "value"] <- "estimate"
  sample <- sample[, c(
    "contract_id", "endpoint_id", "estimand_id", "estimand_role",
    "query_sample_id", "query_tissue", "held_out_study", "estimate",
    "completed_seeds"
  )]
  tissue <- stats::aggregate(
    estimate ~ endpoint_id + estimand_id + estimand_role + query_tissue,
    data = sample, FUN = mean
  )
  tissue_n <- stats::aggregate(
    query_sample_id ~ endpoint_id + estimand_id + estimand_role + query_tissue,
    data = sample, FUN = length
  )
  names(tissue_n)[names(tissue_n) == "query_sample_id"] <- "samples"
  tissue <- merge(tissue, tissue_n, sort = FALSE)
  tissue$contract_id <- "mv05l_tissue_estimand_v1"
  tissue$study_blocks <- vapply(tissue$query_tissue, function(tissue_id) {
    length(unique(sample$held_out_study[sample$query_tissue == tissue_id]))
  }, integer(1L))
  tissue <- tissue[, c(
    "contract_id", "endpoint_id", "estimand_id", "estimand_role",
    "query_tissue", "estimate", "samples", "study_blocks"
  )]
  macro <- stats::aggregate(
    estimate ~ endpoint_id + estimand_id + estimand_role,
    data = tissue, FUN = mean
  )
  macro$contract_id <- "mv05l_macro_estimand_v1"
  macro$tissues <- 5L
  macro$samples <- 90L
  macro$study_blocks <- 15L
  macro$completed_seeds <- 5L
  macro <- macro[, c(
    "contract_id", "endpoint_id", "estimand_id", "estimand_role",
    "estimate", "tissues", "samples", "study_blocks", "completed_seeds"
  )]
  order_rows <- function(x) {
    x$endpoint_id <- factor(x$endpoint_id, levels = .mv05l_endpoints)
    x$estimand_id <- factor(x$estimand_id, levels = .mv05l_estimands)
    extra <- intersect(c("query_tissue", "held_out_study", "query_sample_id"),
                       names(x))
    ordering <- c(list(x$endpoint_id, x$estimand_id), x[extra],
                  list(method = "radix"))
    x <- x[do.call(order, ordering), , drop = FALSE]
    x$endpoint_id <- as.character(x$endpoint_id)
    x$estimand_id <- as.character(x$estimand_id)
    rownames(x) <- NULL
    x
  }
  list(sample = order_rows(sample), tissue = order_rows(tissue),
       macro = order_rows(macro), query = query)
}

.mv05l_macro_from_rows <- function(base, values) {
  tissue <- tapply(values, base$query_tissue, mean)
  mean(tissue[sort(names(tissue), method = "radix")])
}

mv05l_block_inference_v1 <- function(
    sample_summary, bootstrap_replicates = 2000L,
    bootstrap_seed = 20260810L, randomization_replicates = 9999L,
    randomization_seed = 20260811L) {
  required <- c("endpoint_id", "estimand_id", "query_sample_id",
                "query_tissue", "held_out_study", "estimate",
                "completed_seeds")
  if (!is.data.frame(sample_summary) || !all(required %in% names(sample_summary)) ||
      any(sample_summary$completed_seeds != 5L)) {
    stop("Incomplete MV5-L sample estimands.", call. = FALSE)
  }
  sample_ids <- sort(unique(sample_summary$query_sample_id), method = "radix")
  base_part <- sample_summary[
    sample_summary$endpoint_id == .mv05l_endpoints[[1L]] &
      sample_summary$estimand_id == .mv05l_estimands[[1L]], , drop = FALSE
  ]
  base <- base_part[match(sample_ids, base_part$query_sample_id),
                    c("query_sample_id", "query_tissue", "held_out_study")]
  if (anyNA(base$query_sample_id)) {
    stop("MV5-L sample identity base is incomplete.", call. = FALSE)
  }
  arrays <- lapply(.mv05l_endpoints, function(endpoint) {
    matrix(vapply(.mv05l_estimands, function(estimand) {
      part <- sample_summary[
        sample_summary$endpoint_id == endpoint &
          sample_summary$estimand_id == estimand, , drop = FALSE
      ]
      part$estimate[match(sample_ids, part$query_sample_id)]
    }, numeric(length(sample_ids))), nrow = length(sample_ids),
    dimnames = list(sample_ids, .mv05l_estimands))
  })
  names(arrays) <- .mv05l_endpoints
  if (any(!is.finite(unlist(arrays)))) {
    stop("MV5-L estimand arrays contain non-finite values.", call. = FALSE)
  }
  strata <- split(unique(base[, c("query_tissue", "held_out_study")]),
                  unique(base[, c("query_tissue", "held_out_study")])$query_tissue)
  strata <- strata[sort(names(strata), method = "radix")]
  bootstrap <- .mv05l_with_seed(bootstrap_seed, {
    result <- vector("list", length(arrays))
    for (endpoint_index in seq_along(arrays)) {
      values <- arrays[[endpoint_index]]
      estimates <- matrix(NA_real_, bootstrap_replicates,
                          length(.mv05l_estimands),
                          dimnames = list(NULL, .mv05l_estimands))
      for (replicate_id in seq_len(bootstrap_replicates)) {
        tissue_values <- matrix(NA_real_, length(strata),
                                length(.mv05l_estimands))
        for (tissue_index in seq_along(strata)) {
          studies <- sort(strata[[tissue_index]]$held_out_study,
                          method = "radix")
          drawn <- sample(studies, length(studies), replace = TRUE)
          rows <- unlist(lapply(drawn, function(study) {
            which(base$query_tissue == names(strata)[[tissue_index]] &
                    base$held_out_study == study)
          }), use.names = FALSE)
          tissue_values[tissue_index, ] <- colMeans(values[rows, , drop = FALSE])
        }
        estimates[replicate_id, ] <- colMeans(tissue_values)
      }
      result[[endpoint_index]] <- estimates
    }
    result
  })
  names(bootstrap) <- .mv05l_endpoints
  intervals <- do.call(rbind, lapply(.mv05l_endpoints, function(endpoint) {
    values <- arrays[[endpoint]]
    do.call(rbind, lapply(seq_along(.mv05l_estimands), function(index) {
      interval <- stats::quantile(bootstrap[[endpoint]][, index],
                                  c(0.025, 0.975), names = FALSE, type = 7)
      data.frame(
        contract_id = "mv05l_estimand_interval_v1",
        endpoint_id = endpoint,
        estimand_id = .mv05l_estimands[[index]],
        estimand_role = unname(.mv05l_estimand_roles[index]),
        estimate = .mv05l_macro_from_rows(base, values[, index]),
        ci_lower = interval[[1L]], ci_upper = interval[[2L]],
        bootstrap_replicates = bootstrap_replicates,
        bootstrap_seed = bootstrap_seed, tissues = length(strata),
        study_blocks = length(unique(base$held_out_study)),
        samples = nrow(base), completed_seeds = 5L,
        inference_status = "estimable", stringsAsFactors = FALSE
      )
    }))
  }))
  primary_ids <- .mv05l_estimands[1:2]
  primary <- do.call(rbind, lapply(.mv05l_endpoints, function(endpoint) {
    do.call(rbind, lapply(primary_ids, function(estimand) {
      row <- intervals[intervals$endpoint_id == endpoint &
                         intervals$estimand_id == estimand, , drop = FALSE]
      data.frame(
        contract_id = "mv05l_primary_contrast_v1",
        family_id = if (endpoint == .mv05l_endpoints[[1L]])
          "F1_representation_did_mrr" else "none_supportive",
        endpoint_id = endpoint, estimand_id = estimand,
        contrast = "integrated_topology_minus_energy_minus_sct_topology_minus_energy",
        estimate = row$estimate, ci_lower = row$ci_lower,
        ci_upper = row$ci_upper, favorable_direction = "positive",
        bootstrap_replicates = bootstrap_replicates,
        bootstrap_seed = bootstrap_seed,
        randomization_replicates = if (endpoint == .mv05l_endpoints[[1L]])
          randomization_replicates else NA_integer_,
        randomization_seed = if (endpoint == .mv05l_endpoints[[1L]])
          randomization_seed else NA_integer_,
        tissues = length(strata), study_blocks = length(unique(base$held_out_study)),
        samples = nrow(base), completed_seeds = 5L,
        raw_p_value = NA_real_, holm_p_value = NA_real_,
        inference_status = "estimable", stringsAsFactors = FALSE
      )
    }))
  }))
  primary_rows <- which(primary$family_id == "F1_representation_did_mrr")
  study_ids <- sort(unique(base$held_out_study), method = "radix")
  signs <- .mv05l_with_seed(randomization_seed, matrix(
    sample(c(-1, 1), randomization_replicates * length(study_ids),
           replace = TRUE),
    nrow = randomization_replicates, ncol = length(study_ids),
    dimnames = list(NULL, study_ids)
  ))
  null_hash <- character(length(primary_rows))
  exceedances <- integer(length(primary_rows))
  for (position in seq_along(primary_rows)) {
    row_index <- primary_rows[[position]]
    estimand <- primary$estimand_id[[row_index]]
    differences <- arrays[[.mv05l_endpoints[[1L]]]][, estimand]
    null <- vapply(seq_len(randomization_replicates), function(replicate_id) {
      signed <- differences * signs[
        replicate_id, match(base$held_out_study, study_ids)
      ]
      .mv05l_macro_from_rows(base, signed)
    }, numeric(1L))
    exceedances[[position]] <- .mv05l_boundary_exceedances(
      null, primary$estimate[[row_index]]
    )
    primary$raw_p_value[[row_index]] <-
      (exceedances[[position]] + 1) / (randomization_replicates + 1)
    null_hash[[position]] <- digest::digest(
      null, algo = "sha256", serialize = TRUE
    )
  }
  primary$holm_p_value[primary_rows] <- stats::p.adjust(
    primary$raw_p_value[primary_rows], method = "holm"
  )
  bootstrap_audit <- data.frame(
    contract_id = "mv05l_bootstrap_audit_v1",
    endpoint_id = .mv05l_endpoints,
    bootstrap_replicates = bootstrap_replicates,
    bootstrap_seed = bootstrap_seed,
    tissues = length(strata), study_blocks = length(study_ids),
    estimands = length(.mv05l_estimands),
    replicate_matrix_sha256 = vapply(
      bootstrap, digest::digest, character(1L), algo = "sha256",
      serialize = TRUE
    ),
    resampling_unit = "tissue_stratified_held_out_study_block",
    paired_across_representations_methods_and_estimands = TRUE,
    seeds_treated_as_independent = FALSE,
    stringsAsFactors = FALSE
  )
  randomization_audit <- data.frame(
    contract_id = "mv05l_randomization_audit_v1",
    family_id = "F1_representation_did_mrr",
    estimand_id = primary$estimand_id[primary_rows],
    randomization_replicates = randomization_replicates,
    randomization_seed = randomization_seed,
    study_blocks = length(study_ids),
    sign_matrix_sha256 = digest::digest(
      signs, algo = "sha256", serialize = TRUE
    ),
    null_distribution_sha256 = null_hash,
    exceedance_count = exceedances,
    boundary_policy = "absolute_two_sided_64eps_ties_count_as_exceedance_v1",
    epsilon_multiplier = 64L, two_sided = TRUE,
    exact_zero_p_value_possible = FALSE,
    stringsAsFactors = FALSE
  )
  list(intervals = intervals, primary_contrasts = primary,
       bootstrap_audit = bootstrap_audit,
       randomization_audit = randomization_audit)
}
