# Prediction-locked MV6-H blocked fusion-outcome helpers.

.mv06h_methods <- c(
  "cell_H0", "cell_H1", "gene_H0", "gene_H1", "cell_composite",
  "fusion_gene_weight_025", "fusion_gene_weight_050",
  "fusion_gene_weight_075", "gene_composite"
)
.mv06h_seeds <- 20260805:20260809
.mv06h_eligible_tissues <- c(
  "bone marrow", "colon", "liver", "pbmc", "testis"
)
.mv06h_primary <- "fusion_gene_weight_050"
.mv06h_comparators <- c("cell_composite", "gene_composite")

.mv06h_sha256 <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

.mv06h_digest <- function(value) {
  digest::digest(value, algo = "sha256", serialize = TRUE)
}

.mv06h_is_true <- function(value) {
  tolower(as.character(value)) == "true"
}

.mv06h_safe_group <- function(value) {
  gsub("[^A-Za-z0-9._-]", "_", value)
}

.mv06h_with_seed <- function(seed, expression) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  old_kind <- RNGkind()
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(as.integer(seed))
  force(expression)
}

mv06h_validate_group_manifest_v1 <- function(manifest) {
  required <- c(
    "group_id", "fold_id", "held_out_study", "seed", "execution_order",
    "group_root_kind", "group_locator", "rankings_sha256",
    "rankings_bytes", "ranking_rows", "scales_sha256", "status_sha256",
    "production_implementation_root_sha256", "outcome_label_state",
    "biological_outcomes_computed", "fusion_evaluations", "outcome_jobs"
  )
  if (!is.data.frame(manifest) || nrow(manifest) != 75L ||
      !all(required %in% names(manifest)) || anyDuplicated(manifest$group_id) ||
      !identical(sort(as.integer(manifest$execution_order)), 1:75) ||
      !identical(sort(unique(as.integer(manifest$seed))), .mv06h_seeds) ||
      sum(manifest$group_root_kind == "accepted_stage1_sentinel") != 1L ||
      sum(manifest$group_root_kind == "corrected_serial_completion") != 74L ||
      sum(as.integer(manifest$ranking_rows)) != 318150L ||
      any(manifest$outcome_label_state != "closed") ||
      any(.mv06h_is_true(manifest$biological_outcomes_computed)) ||
      any(as.integer(manifest$fusion_evaluations) != 0L) ||
      any(as.integer(manifest$outcome_jobs) != 0L) ||
      any(!grepl("^[0-9a-f]{64}$", manifest$rankings_sha256)) ||
      any(!grepl("^[0-9a-f]{64}$", manifest$scales_sha256)) ||
      any(!grepl("^[0-9a-f]{64}$", manifest$status_sha256))) {
    stop("MV6-H prediction group manifest is incomplete or outcome-open.",
         call. = FALSE)
  }
  invisible(manifest)
}

mv06h_verify_prediction_lock_v1 <- function(lock_dir) {
  files <- c(
    groups = "mv06h-prediction-group-manifest.csv",
    sources = "mv06h-source-manifest.csv",
    implementation = "mv06h-implementation-manifest.csv",
    methods = "mv06h-method-registry.csv",
    endpoints = "mv06h-endpoint-registry.csv",
    contrasts = "mv06h-contrast-registry.csv",
    inference = "mv06h-inference-registry.csv",
    firewall = "mv06h-label-firewall.csv"
  )
  paths <- stats::setNames(file.path(lock_dir, unname(files)), names(files))
  lock_path <- file.path(lock_dir, "mv06h-prediction-lock.csv")
  if (any(!file.exists(paths)) || !file.exists(lock_path)) {
    stop("MV6-H prediction-lock evidence is incomplete.", call. = FALSE)
  }
  lock <- utils::read.csv(lock_path, stringsAsFactors = FALSE,
                          check.names = FALSE)
  if (nrow(lock) != 1L || !.mv06h_is_true(lock$label_open_authorized) ||
      lock$prediction_lock_status != "passed_before_label_open" ||
      .mv06h_is_true(lock$labels_opened) ||
      .mv06h_is_true(lock$outcomes_computed)) {
    stop("MV6-H prediction lock does not authorize label opening.", call. = FALSE)
  }
  actual <- vapply(paths, .mv06h_sha256, character(1L))
  expected <- unlist(lock[paste0(names(files), "_sha256")], use.names = TRUE)
  names(expected) <- sub("_sha256$", "", names(expected))
  if (!identical(unname(actual), unname(as.character(expected)))) {
    stop("MV6-H prediction-lock ledger hash drift.", call. = FALSE)
  }
  groups <- utils::read.csv(paths[["groups"]], stringsAsFactors = FALSE,
                            check.names = FALSE)
  mv06h_validate_group_manifest_v1(groups)
  list(lock = lock, groups = groups)
}

mv06h_open_frozen_labels_v1 <- function(metadata_path, prediction_lock,
    expected_metadata_sha256 =
      "e2b6f1e7d7dc05ba8e80627aad4ea1677c3b718a844397acd06596d7e22a18d0") {
  if (!is.data.frame(prediction_lock) || nrow(prediction_lock) != 1L ||
      !.mv06h_is_true(prediction_lock$label_open_authorized) ||
      prediction_lock$prediction_lock_status != "passed_before_label_open") {
    stop("MV6-H labels may open only after a passing prediction lock.",
         call. = FALSE)
  }
  metadata_sha <- .mv06h_sha256(metadata_path)
  if (!identical(metadata_sha, expected_metadata_sha256)) {
    stop("MV6-H metadata SHA-256 drifted.", call. = FALSE)
  }
  raw <- utils::read.csv(metadata_path, stringsAsFactors = FALSE,
                         check.names = TRUE)
  required <- c("orig.ident", "SRA", "Tissue.x", "Approach.x",
                "Number_of_Cells_After_Filtering")
  if (!all(required %in% names(raw))) {
    stop("MV6-H frozen metadata schema drifted.", call. = FALSE)
  }
  labels <- data.frame(
    sample_id = trimws(as.character(raw$orig.ident)),
    study = trimws(as.character(raw$SRA)),
    tissue = tolower(trimws(as.character(raw$Tissue.x))),
    approach = trimws(as.character(raw$Approach.x)),
    filtered_cells = as.integer(raw$Number_of_Cells_After_Filtering),
    stringsAsFactors = FALSE
  )
  if (nrow(labels) != 124L || anyNA(labels) || anyDuplicated(labels$sample_id) ||
      length(unique(labels$study)) != 18L) {
    stop("MV6-H frozen metadata identity drifted.", call. = FALSE)
  }
  candidates <- labels[labels$tissue %in% .mv06h_eligible_tissues, , drop = FALSE]
  candidates <- candidates[order(candidates$sample_id, method = "radix"), ]
  counts <- c("bone marrow" = 31L, colon = 13L, liver = 6L,
              pbmc = 12L, testis = 28L)
  observed <- table(factor(candidates$tissue, levels = names(counts)))
  if (nrow(candidates) != 90L || length(unique(candidates$study)) != 15L ||
      !identical(as.integer(observed), as.integer(counts)) ||
      any(tapply(candidates$tissue, candidates$study,
                 function(x) length(unique(x))) != 1L)) {
    stop("MV6-H candidate label axis drifted.", call. = FALSE)
  }
  list(labels = candidates, provenance = data.frame(
    contract_id = "mv06h_frozen_label_source_v1",
    source_file = basename(metadata_path), source_sha256 = metadata_sha,
    source_samples = nrow(labels), source_studies = length(unique(labels$study)),
    candidate_samples = nrow(candidates),
    candidate_studies = length(unique(candidates$study)),
    eligible_tissues = length(unique(candidates$tissue)),
    label_open_boundary = "after_committed_mv06h_prediction_lock",
    labels_used_for_upstream_fit = FALSE, stringsAsFactors = FALSE))
}

mv06h_read_locked_rankings_v1 <- function(stage1_root, completion_root,
                                           group_manifest) {
  mv06h_validate_group_manifest_v1(group_manifest)
  pieces <- vector("list", nrow(group_manifest))
  for (index in seq_len(nrow(group_manifest))) {
    row <- group_manifest[index, , drop = FALSE]
    root <- if (row$group_root_kind == "accepted_stage1_sentinel")
      stage1_root else completion_root
    path <- file.path(root, .mv06h_safe_group(row$group_id), "rankings.csv")
    if (!file.exists(path) || .mv06h_sha256(path) != row$rankings_sha256) {
      stop("MV6-H locked ranking artifact hash drift: ", row$group_id,
           call. = FALSE)
    }
    part <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    if (nrow(part) != as.integer(row$ranking_rows)) {
      stop("MV6-H locked ranking row count drift.", call. = FALSE)
    }
    pieces[[index]] <- part
  }
  result <- do.call(rbind, pieces)
  rownames(result) <- NULL
  result
}

mv06h_validate_rankings_v1 <- function(rankings, labels,
    expected_rows = 318150L, expected_samples = 90L,
    expected_methods = .mv06h_methods, expected_seeds = .mv06h_seeds) {
  required <- c(
    "group_id", "fold_id", "seed", "query_sample_id", "training_sample_id",
    "method_id", "normalized_distance", "tie_break", "rank",
    "outcome_label_state", "biological_outcomes_computed",
    "fusion_evaluations", "outcome_jobs"
  )
  if (!is.data.frame(rankings) || !all(required %in% names(rankings)) ||
      nrow(rankings) != expected_rows ||
      anyNA(rankings[c("group_id", "fold_id", "seed", "query_sample_id",
                       "training_sample_id", "method_id",
                       "normalized_distance", "tie_break", "rank")]) ||
      any(!is.finite(rankings$normalized_distance)) ||
      any(rankings$normalized_distance < 0) ||
      any(rankings$outcome_label_state != "closed") ||
      any(.mv06h_is_true(rankings$biological_outcomes_computed)) ||
      any(as.integer(rankings$fusion_evaluations) != 0L) ||
      any(as.integer(rankings$outcome_jobs) != 0L) ||
      any(rankings$tie_break != "canonical_training_sample_id") ||
      !setequal(rankings$method_id, expected_methods) ||
      !identical(sort(unique(as.integer(rankings$seed))),
                 sort(as.integer(expected_seeds)))) {
    stop("MV6-H ranking corpus is malformed or outcome-open.", call. = FALSE)
  }
  ids <- unique(c(rankings$query_sample_id, rankings$training_sample_id))
  if (length(ids) != expected_samples || !setequal(ids, labels$sample_id)) {
    stop("MV6-H ranking and label sample axes differ.", call. = FALSE)
  }
  groups <- split(seq_len(nrow(rankings)), interaction(
    rankings$group_id, rankings$method_id, rankings$query_sample_id,
    drop = TRUE, lex.order = TRUE))
  valid <- vapply(groups, function(ix) {
    part <- rankings[ix, , drop = FALSE]
    part <- part[order(as.integer(part$rank)), , drop = FALSE]
    identical(as.integer(part$rank), seq_len(nrow(part))) &&
      !anyDuplicated(part$training_sample_id) &&
      identical(part$training_sample_id,
                part$training_sample_id[order(
                  part$normalized_distance, part$training_sample_id,
                  method = "radix")])
  }, logical(1L))
  if (!all(valid)) stop("MV6-H locked rankings are not canonical.", call. = FALSE)
  invisible(TRUE)
}

mv06h_evaluate_rankings_v1 <- function(rankings, labels,
    expected_rows = 318150L, expected_observations = 4050L,
    expected_samples = 90L, expected_methods = .mv06h_methods,
    expected_seeds = .mv06h_seeds) {
  mv06h_validate_rankings_v1(rankings, labels, expected_rows, expected_samples,
                              expected_methods, expected_seeds)
  train_index <- match(rankings$training_sample_id, labels$sample_id)
  query_index <- match(rankings$query_sample_id, labels$sample_id)
  heldout <- sub("^large_loso_v1:", "", rankings$fold_id)
  if (any(labels$study[query_index] != heldout) ||
      any(labels$study[train_index] == heldout)) {
    stop("MV6-H held-out-study isolation failed.", call. = FALSE)
  }
  rankings$.training_tissue <- labels$tissue[train_index]
  groups <- split(seq_len(nrow(rankings)), interaction(
    rankings$group_id, rankings$method_id, rankings$query_sample_id,
    drop = TRUE, lex.order = TRUE))
  result <- do.call(rbind, lapply(groups, function(ix) {
    part <- rankings[ix, , drop = FALSE]
    part <- part[order(as.integer(part$rank)), , drop = FALSE]
    qi <- match(part$query_sample_id[[1L]], labels$sample_id)
    q_tissue <- labels$tissue[[qi]]; q_study <- labels$study[[qi]]
    same <- which(part$.training_tissue == q_tissue)
    status <- if (!q_tissue %in% .mv06h_eligible_tissues)
      "single_study_tissue_not_estimable" else if (!length(same))
      "training_tissue_absent_not_estimable" else "estimable"
    first <- if (length(same)) same[[1L]] else NA_integer_
    data.frame(
      contract_id = "mv06h_query_method_outcome_v1",
      group_id = part$group_id[[1L]], fold_id = part$fold_id[[1L]],
      held_out_study = q_study, seed = as.integer(part$seed[[1L]]),
      method_id = part$method_id[[1L]],
      query_sample_id = part$query_sample_id[[1L]], query_tissue = q_tissue,
      training_samples = nrow(part),
      first_same_tissue_sample_id = if (is.na(first)) NA_character_ else
        part$training_sample_id[[first]],
      first_same_tissue_rank = if (is.na(first)) NA_integer_ else
        as.integer(part$rank[[first]]),
      reciprocal_rank = if (is.na(first)) NA_real_ else 1 / as.numeric(part$rank[[first]]),
      nearest_sample_id = part$training_sample_id[[1L]],
      nearest_tissue = part$.training_tissue[[1L]],
      one_nn_correct = if (status == "estimable")
        part$.training_tissue[[1L]] == q_tissue else NA,
      endpoint_status = status, labels_opened_after_prediction_lock = TRUE,
      upstream_refit = FALSE, reranked_after_label_open = FALSE,
      method_selection_executed = FALSE, outcomes_computed = TRUE,
      stringsAsFactors = FALSE)
  }))
  result <- result[order(match(result$method_id, expected_methods),
                          result$seed, result$query_tissue,
                          result$held_out_study, result$query_sample_id,
                          method = "radix"), ]
  rownames(result) <- NULL
  if (nrow(result) != expected_observations ||
      any(result$endpoint_status != "estimable")) {
    stop("MV6-H query-outcome completion/estimability failed.", call. = FALSE)
  }
  result
}

mv06h_summarize_outcomes_v1 <- function(observations,
                                         expected_seeds = .mv06h_seeds) {
  groups <- split(seq_len(nrow(observations)), interaction(
    observations$method_id, observations$query_sample_id,
    drop = TRUE, lex.order = TRUE))
  sample <- do.call(rbind, lapply(groups, function(ix) {
    part <- observations[ix, , drop = FALSE]
    data.frame(
      contract_id = "mv06h_sample_method_summary_v1",
      method_id = part$method_id[[1L]],
      query_sample_id = part$query_sample_id[[1L]],
      query_tissue = part$query_tissue[[1L]],
      held_out_study = part$held_out_study[[1L]],
      mean_reciprocal_rank = mean(part$reciprocal_rank),
      one_nn_balanced_accuracy = mean(as.numeric(part$one_nn_correct)),
      seeds = length(unique(part$seed)), stringsAsFactors = FALSE)
  }))
  if (any(sample$seeds != length(expected_seeds))) {
    stop("MV6-H technical seeds were not complete within sample.", call. = FALSE)
  }
  tissue_groups <- split(seq_len(nrow(sample)), interaction(
    sample$method_id, sample$query_tissue, drop = TRUE, lex.order = TRUE))
  tissue <- do.call(rbind, lapply(tissue_groups, function(ix) {
    part <- sample[ix, , drop = FALSE]
    data.frame(contract_id = "mv06h_tissue_method_summary_v1",
      method_id = part$method_id[[1L]], query_tissue = part$query_tissue[[1L]],
      mean_reciprocal_rank = mean(part$mean_reciprocal_rank),
      one_nn_balanced_accuracy = mean(part$one_nn_balanced_accuracy),
      samples = nrow(part), studies = length(unique(part$held_out_study)),
      stringsAsFactors = FALSE)
  }))
  method_groups <- split(seq_len(nrow(tissue)), tissue$method_id)
  method <- do.call(rbind, lapply(method_groups, function(ix) {
    part <- tissue[ix, , drop = FALSE]
    data.frame(contract_id = "mv06h_method_macro_summary_v1",
      method_id = part$method_id[[1L]],
      macro_mean_reciprocal_rank = mean(part$mean_reciprocal_rank),
      macro_one_nn_balanced_accuracy = mean(part$one_nn_balanced_accuracy),
      tissues = nrow(part), samples = sum(part$samples),
      stringsAsFactors = FALSE)
  }))
  order_method <- function(x) x[order(match(x$method_id, .mv06h_methods),
                                       method = "radix"), , drop = FALSE]
  list(sample = order_method(sample), tissue = order_method(tissue),
       method = order_method(method))
}

.mv06h_macro <- function(values, tissue) {
  means <- tapply(values, tissue, mean)
  mean(means[sort(names(means), method = "radix")])
}

mv06h_block_inference_v1 <- function(sample_summary,
    bootstrap_replicates = 2000L, bootstrap_seed = 20260815L,
    randomization_replicates = 9999L, randomization_seed = 20260816L,
    expected_samples = 90L, methods = .mv06h_methods,
    expected_seeds = .mv06h_seeds) {
  required <- c("method_id", "query_sample_id", "query_tissue",
                "held_out_study", "mean_reciprocal_rank",
                "one_nn_balanced_accuracy", "seeds")
  if (!is.data.frame(sample_summary) || !all(required %in% names(sample_summary)) ||
      nrow(sample_summary) != expected_samples * length(methods) ||
      any(sample_summary$seeds != length(expected_seeds))) {
    stop("MV6-H sample summaries are incomplete for blocked inference.",
         call. = FALSE)
  }
  sample_ids <- sort(unique(sample_summary$query_sample_id), method = "radix")
  base <- sample_summary[sample_summary$method_id == methods[[1L]],
    c("query_sample_id", "query_tissue", "held_out_study")]
  base <- base[match(sample_ids, base$query_sample_id), , drop = FALSE]
  endpoints <- c(
    cross_study_tissue_mrr_v1 = "mean_reciprocal_rank",
    cross_study_tissue_1nn_balanced_accuracy_v1 = "one_nn_balanced_accuracy")
  arrays <- lapply(endpoints, function(metric) matrix(vapply(methods,
    function(method) {
      part <- sample_summary[sample_summary$method_id == method, ]
      part[[metric]][match(sample_ids, part$query_sample_id)]
    }, numeric(length(sample_ids))), nrow = length(sample_ids),
    dimnames = list(sample_ids, methods)))
  strata_table <- unique(base[c("query_tissue", "held_out_study")])
  strata <- split(strata_table, strata_table$query_tissue)
  strata <- strata[sort(names(strata), method = "radix")]
  bootstrap <- .mv06h_with_seed(bootstrap_seed, lapply(arrays, function(values) {
    estimates <- matrix(NA_real_, bootstrap_replicates, length(methods),
                        dimnames = list(NULL, methods))
    for (b in seq_len(bootstrap_replicates)) {
      tissue_values <- matrix(NA_real_, length(strata), length(methods))
      for (t in seq_along(strata)) {
        studies <- sort(strata[[t]]$held_out_study, method = "radix")
        drawn <- sample(studies, length(studies), replace = TRUE)
        ix <- unlist(lapply(drawn, function(study) which(
          base$query_tissue == names(strata)[[t]] &
            base$held_out_study == study)), use.names = FALSE)
        tissue_values[t, ] <- colMeans(values[ix, , drop = FALSE])
      }
      estimates[b, ] <- colMeans(tissue_values)
    }
    estimates
  }))
  intervals <- do.call(rbind, lapply(names(arrays), function(endpoint) {
    values <- arrays[[endpoint]]
    do.call(rbind, lapply(seq_along(methods), function(j) {
      ci <- stats::quantile(bootstrap[[endpoint]][, j], c(.025, .975),
                            names = FALSE, type = 7)
      data.frame(contract_id = "mv06h_method_interval_v1",
        endpoint_id = endpoint, method_id = methods[[j]],
        estimate = .mv06h_macro(values[, j], base$query_tissue),
        ci_lower = ci[[1L]], ci_upper = ci[[2L]],
        bootstrap_replicates = bootstrap_replicates,
        bootstrap_seed = bootstrap_seed, tissues = length(strata),
        study_blocks = length(unique(base$held_out_study)),
        samples = nrow(base), completed_seeds = length(expected_seeds),
        stringsAsFactors = FALSE)
    }))
  }))
  primary_values <- arrays[["cross_study_tissue_mrr_v1"]]
  fusion_index <- match(.mv06h_primary, methods)
  contrasts <- do.call(rbind, lapply(.mv06h_comparators, function(comparator) {
    comparator_index <- match(comparator, methods)
    differences <- primary_values[, fusion_index] -
      primary_values[, comparator_index]
    boot <- bootstrap[["cross_study_tissue_mrr_v1"]][, fusion_index] -
      bootstrap[["cross_study_tissue_mrr_v1"]][, comparator_index]
    ci <- stats::quantile(boot, c(.025, .975), names = FALSE, type = 7)
    data.frame(contract_id = "mv06h_primary_paired_contrast_v1",
      family_id = "mv06g_F1_primary_mrr",
      endpoint_id = "cross_study_tissue_mrr_v1",
      method_id = .mv06h_primary, comparator_id = comparator,
      estimate = .mv06h_macro(differences, base$query_tissue),
      ci_lower = ci[[1L]], ci_upper = ci[[2L]],
      favorable_direction = "positive", raw_p_value = NA_real_,
      holm_p_value = NA_real_, bootstrap_replicates = bootstrap_replicates,
      bootstrap_seed = bootstrap_seed,
      randomization_replicates = randomization_replicates,
      randomization_seed = randomization_seed,
      tissues = length(strata), study_blocks = length(unique(base$held_out_study)),
      samples = nrow(base), completed_seeds = length(expected_seeds),
      stringsAsFactors = FALSE)
  }))
  study_ids <- sort(unique(base$held_out_study), method = "radix")
  signs <- .mv06h_with_seed(randomization_seed, matrix(sample(c(-1, 1),
    randomization_replicates * length(study_ids), replace = TRUE),
    nrow = randomization_replicates, ncol = length(study_ids),
    dimnames = list(NULL, study_ids)))
  null_hash <- character(nrow(contrasts))
  for (i in seq_len(nrow(contrasts))) {
    comparator <- contrasts$comparator_id[[i]]
    differences <- primary_values[, fusion_index] -
      primary_values[, match(comparator, methods)]
    null <- vapply(seq_len(randomization_replicates), function(b) {
      .mv06h_macro(differences * signs[b, match(base$held_out_study, study_ids)],
                   base$query_tissue)
    }, numeric(1L))
    contrasts$raw_p_value[[i]] <-
      (sum(abs(null) >= abs(contrasts$estimate[[i]])) + 1) /
      (randomization_replicates + 1)
    null_hash[[i]] <- .mv06h_digest(null)
  }
  contrasts$holm_p_value <- stats::p.adjust(contrasts$raw_p_value, "holm")
  contrasts$fusion_benefit_both_positive <- all(contrasts$estimate > 0)
  bootstrap_audit <- data.frame(
    contract_id = "mv06h_bootstrap_audit_v1", endpoint_id = names(bootstrap),
    bootstrap_replicates = bootstrap_replicates,
    bootstrap_seed = bootstrap_seed, tissues = length(strata),
    study_blocks = length(study_ids),
    replicate_matrix_sha256 = vapply(bootstrap, .mv06h_digest, character(1L)),
    resampling_unit = "tissue_stratified_held_out_study_block",
    seeds_treated_as_independent = FALSE, stringsAsFactors = FALSE)
  randomization_audit <- data.frame(
    contract_id = "mv06h_randomization_audit_v1",
    family_id = "mv06g_F1_primary_mrr", method_id = .mv06h_primary,
    comparator_id = contrasts$comparator_id,
    randomization_replicates = randomization_replicates,
    randomization_seed = randomization_seed, study_blocks = length(study_ids),
    sign_matrix_sha256 = .mv06h_digest(signs),
    null_distribution_sha256 = null_hash, two_sided = TRUE,
    exact_zero_p_value_possible = FALSE, stringsAsFactors = FALSE)
  list(method_intervals = intervals, contrasts = contrasts,
       bootstrap_audit = bootstrap_audit,
       randomization_audit = randomization_audit)
}
