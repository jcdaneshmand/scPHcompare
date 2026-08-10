# Prediction-locked MV5-K integrated cell retrieval evaluation helpers.

.mv05k_methods <- c(
  "integrated_cell_landscape_h0_v1",
  "integrated_cell_landscape_h1_v1",
  "integrated_cell_landscape_h0_h1_raw_euclidean_v1",
  "integrated_cell_distribution_energy_v1",
  "pseudobulk_training_standardized_panel_v1"
)

.mv05k_method_roles <- c(
  integrated_cell_landscape_h0_v1 = "confirmatory",
  integrated_cell_landscape_h1_v1 = "confirmatory",
  integrated_cell_landscape_h0_h1_raw_euclidean_v1 = "descriptive_secondary",
  integrated_cell_distribution_energy_v1 = "matched_integrated_cell_baseline",
  pseudobulk_training_standardized_panel_v1 = "context_baseline"
)

.mv05k_eligible_tissues <- c(
  "bone marrow", "colon", "liver", "pbmc", "testis"
)

.mv05k_seeds <- 20260805:20260809

.mv05k_file_sha256 <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

.mv05k_restore_rng <- function(had_seed, old_seed) {
  if (had_seed) {
    assign(".Random.seed", old_seed, envir = .GlobalEnv)
  } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    rm(".Random.seed", envir = .GlobalEnv)
  }
}

.mv05k_with_seed <- function(seed, expression) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit(.mv05k_restore_rng(had_seed, if (had_seed) old_seed else NULL),
          add = TRUE)
  set.seed(as.integer(seed))
  force(expression)
}

.mv05k_randomization_exceedance_count_v1 <- function(null, observed,
                                                       epsilon_multiplier = 64) {
  null <- as.numeric(null)
  observed <- as.numeric(observed)
  if (!length(null) || length(observed) != 1L || !is.finite(observed) ||
      any(!is.finite(null)) || !is.finite(epsilon_multiplier) ||
      epsilon_multiplier < 1) {
    stop("Randomization exceedance inputs are invalid.", call. = FALSE)
  }
  boundary <- abs(observed)
  tolerance <- epsilon_multiplier * .Machine$double.eps *
    pmax(1, abs(null), boundary)
  sum(abs(null) + tolerance >= boundary)
}

mv05k_verify_prediction_lock_v1 <- function(
    ranking_path, completion_path, group_index_path, method_registry_path,
    scale_disposition_path, assembly_summary_path,
    source_commit,
    expected_ranking_sha256 =
      "4588902bce89a04cae0c7676b4f21f81e83013a29120ca2a4b39f3ffacfb677e",
    expected_completion_sha256 =
      "971e6902c1b6c2b8d39c74a9d07abc9c6d05c99c7b47bb8bd10e4a2209a08ed9",
    expected_group_index_sha256 =
      "0fd3e05a312079cebd4e49f4d0f7326cda037fd460772176385d6ea9e0e023ab",
    expected_method_registry_sha256 =
      "becb8451948e421916915fa231a50a83adb22fbd62cc10cadb3a6d934fd7de5d",
    expected_scale_disposition_sha256 =
      "9f6f93e73c072088e9917d395dd2f98005bfeb9d4ec6ec6dec33c0d50333e8bf",
    expected_assembly_sha256 =
      "aabfdb0338826451d08efb27ebacea4eb6b56cad5fc30bd6f9159f60baeff68c") {
  paths <- c(ranking_path, completion_path, group_index_path,
             method_registry_path, scale_disposition_path,
             assembly_summary_path)
  if (any(!file.exists(paths))) {
    stop("MV5-J prediction-lock input is missing.", call. = FALSE)
  }
  source_commit <- as.character(source_commit)
  expected_source_commit <-
    "6836d1cec28bf23f2b920c3e88bffb2e12385e24"
  if (length(source_commit) != 1L ||
      !identical(source_commit, expected_source_commit)) {
    stop("Git source commit does not match the frozen MV5-K base.",
         call. = FALSE)
  }
  ranking_sha <- .mv05k_file_sha256(ranking_path)
  completion_sha <- .mv05k_file_sha256(completion_path)
  group_sha <- .mv05k_file_sha256(group_index_path)
  registry_sha <- .mv05k_file_sha256(method_registry_path)
  scale_sha <- .mv05k_file_sha256(scale_disposition_path)
  assembly_sha <- .mv05k_file_sha256(assembly_summary_path)
  if (!identical(ranking_sha, expected_ranking_sha256) ||
      !identical(completion_sha, expected_completion_sha256) ||
      !identical(group_sha, expected_group_index_sha256) ||
      !identical(registry_sha, expected_method_registry_sha256) ||
      !identical(scale_sha, expected_scale_disposition_sha256) ||
      !identical(assembly_sha, expected_assembly_sha256)) {
    stop("An MV5-J SHA-256 does not match the frozen prediction lock.",
         call. = FALSE)
  }
  completion <- utils::read.csv(
    completion_path, stringsAsFactors = FALSE, check.names = FALSE
  )
  groups <- utils::read.csv(
    group_index_path, stringsAsFactors = FALSE, check.names = FALSE
  )
  registry <- utils::read.csv(
    method_registry_path, stringsAsFactors = FALSE, check.names = FALSE
  )
  scale <- utils::read.csv(
    scale_disposition_path, stringsAsFactors = FALSE, check.names = FALSE
  )
  assembly <- utils::read.csv(
    assembly_summary_path, stringsAsFactors = FALSE, check.names = FALSE
  )
  required_completion <- c(
    "group_id", "seed", "method_id", "status", "outcome_label_state",
    "biological_outcomes_computed"
  )
  required_groups <- c(
    "group_id", "seed", "completed_methods", "failed_methods",
    "outcome_label_state", "biological_outcomes_computed"
  )
  if (!all(required_completion %in% names(completion)) ||
      !all(required_groups %in% names(groups)) || nrow(assembly) != 1L ||
      !all(c("method_id", "role", "primary_retrieval") %in% names(registry)) ||
      !all(c("method_id", "held_out_query_pairs_used_for_scale",
             "additional_topology_jobs_executed") %in% names(scale))) {
    stop("MV5-J prediction-lock evidence has an unexpected schema.",
         call. = FALSE)
  }
  completion_keys <- paste(completion$group_id, completion$method_id, sep = "|")
  valid <- nrow(completion) == 375L && nrow(groups) == 75L &&
    !anyDuplicated(completion_keys) && !anyDuplicated(groups$group_id) &&
    nrow(registry) == 5L && !anyDuplicated(registry$method_id) &&
    identical(as.character(registry$method_id), .mv05k_methods) &&
    nrow(scale) == 225L &&
    all(scale$held_out_query_pairs_used_for_scale == 0L) &&
    all(scale$additional_topology_jobs_executed == 0L) &&
    identical(sort(unique(as.character(completion$method_id)), method = "radix"),
              sort(.mv05k_methods, method = "radix")) &&
    identical(sort(unique(as.integer(completion$seed))), .mv05k_seeds) &&
    all(completion$status == "completed") &&
    all(completion$outcome_label_state == "closed") &&
    !any(as.logical(completion$biological_outcomes_computed)) &&
    all(groups$completed_methods == 5L) && all(groups$failed_methods == 0L) &&
    all(groups$outcome_label_state == "closed") &&
    !any(as.logical(groups$biological_outcomes_computed)) &&
    identical(as.character(assembly$ranking_file_sha256), ranking_sha) &&
    identical(as.character(assembly$completion_file_sha256), completion_sha) &&
    identical(as.character(assembly$group_index_file_sha256), group_sha) &&
    identical(as.character(assembly$method_registry_file_sha256), registry_sha) &&
    identical(as.character(assembly$scale_disposition_file_sha256), scale_sha) &&
    assembly$groups == 75L && assembly$methods == 5L &&
    assembly$retrieval_rows == 176750L &&
    assembly$method_completion_rows == 375L &&
    assembly$failed_method_groups == 0L &&
    identical(as.character(assembly$outcome_label_state), "closed") &&
    !as.logical(assembly$biological_outcomes_computed)
  if (!isTRUE(valid)) {
    stop("MV5-J prediction-lock evidence failed identity/completion checks.",
         call. = FALSE)
  }
  data.frame(
    contract_id = "mv05k_prediction_lock_v1",
    source_commit = source_commit,
    ranking_sha256 = ranking_sha,
    completion_sha256 = completion_sha,
    group_index_sha256 = group_sha,
    method_registry_sha256 = registry_sha,
    scale_disposition_sha256 = scale_sha,
    assembly_summary_sha256 = assembly_sha,
    ranking_rows = 176750L,
    bundle_groups = nrow(groups),
    completion_groups = nrow(completion),
    completed_method_groups = sum(completion$status == "completed"),
    failed_method_groups = sum(completion$status != "completed"),
    methods = length(unique(completion$method_id)),
    seeds = length(unique(completion$seed)),
    label_open_authorized = TRUE,
    prediction_lock_status = "passed_before_label_open",
    stringsAsFactors = FALSE
  )
}

mv05k_open_frozen_labels_v1 <- function(metadata_path, prediction_lock,
                                         expected_metadata_sha256 =
    "e2b6f1e7d7dc05ba8e80627aad4ea1677c3b718a844397acd06596d7e22a18d0") {
  if (!is.data.frame(prediction_lock) || nrow(prediction_lock) != 1L ||
      !isTRUE(prediction_lock$label_open_authorized[[1L]])) {
    stop("Frozen labels may open only after a passing prediction lock.",
         call. = FALSE)
  }
  metadata_sha <- .mv05k_file_sha256(metadata_path)
  if (!identical(metadata_sha, expected_metadata_sha256)) {
    stop("Metadata SHA-256 does not match the frozen MV-05 label source.",
         call. = FALSE)
  }
  raw <- utils::read.csv(
    metadata_path, stringsAsFactors = FALSE, check.names = TRUE
  )
  required <- c(
    "orig.ident", "SRA", "Tissue.x", "Approach.x",
    "Number_of_Cells_After_Filtering"
  )
  if (!all(required %in% names(raw))) {
    stop("Frozen MV-05 metadata lacks required fields.", call. = FALSE)
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
    stop("Frozen MV-05 metadata identity/count checks failed.", call. = FALSE)
  }
  studies_by_tissue <- tapply(labels$study, labels$tissue,
                              function(x) length(unique(x)))
  eligible <- sort(names(studies_by_tissue)[studies_by_tissue >= 2L],
                   method = "radix")
  if (!identical(eligible, sort(.mv05k_eligible_tissues, method = "radix"))) {
    stop("Frozen cross-study tissue eligibility has changed.", call. = FALSE)
  }
  candidates <- labels[labels$tissue %in% .mv05k_eligible_tissues, , drop = FALSE]
  candidates <- candidates[order(candidates$sample_id, method = "radix"), , drop = FALSE]
  if (nrow(candidates) != 90L || length(unique(candidates$study)) != 15L ||
      any(tapply(candidates$tissue, candidates$study,
                 function(x) length(unique(x))) != 1L)) {
    stop("Frozen candidate sample/study structure has changed.", call. = FALSE)
  }
  expected_counts <- c(
    "bone marrow" = 31L, colon = 13L, liver = 6L, pbmc = 12L, testis = 28L
  )
  if (!identical(as.integer(table(factor(candidates$tissue,
                                         levels = names(expected_counts)))),
                 as.integer(expected_counts))) {
    stop("Frozen candidate tissue counts have changed.", call. = FALSE)
  }
  list(
    labels = candidates,
    provenance = data.frame(
      contract_id = "mv05k_frozen_label_source_v1",
      source_file = basename(metadata_path),
      source_sha256 = metadata_sha,
      source_samples = nrow(labels), source_studies = length(unique(labels$study)),
      candidate_samples = nrow(candidates),
      candidate_studies = length(unique(candidates$study)),
      eligible_tissues = length(unique(candidates$tissue)),
      label_open_boundary = "after_mv05j_prediction_lock_passed",
      labels_used_for_upstream_fit = FALSE,
      stringsAsFactors = FALSE
    )
  )
}

.mv05k_validate_rankings <- function(rankings, labels,
                                      expected_rows = 176750L,
                                      expected_methods = .mv05k_methods,
                                      expected_seeds = .mv05k_seeds) {
  required <- c(
    "fold_id", "seed", "method_id", "method_role", "query_sample_id",
    "training_sample_id", "distance", "pair_id", "neighbor_rank",
    "distance_tie_size", "distance_tied", "tie_break_policy",
    "outcome_label_state", "biological_outcomes_computed"
  )
  if (!is.data.frame(rankings) || !all(required %in% names(rankings)) ||
      nrow(rankings) != expected_rows) {
    stop("MV5-J rankings have an invalid schema or row count.", call. = FALSE)
  }
  if (anyNA(rankings[, required]) || anyDuplicated(rankings$pair_id) ||
      any(!is.finite(rankings$distance)) || any(rankings$distance < 0) ||
      any(rankings$outcome_label_state != "closed") ||
      any(as.logical(rankings$biological_outcomes_computed)) ||
      any(rankings$tie_break_policy !=
            "exact_distance_then_canonical_sample_id_v1") ||
      !identical(sort(unique(rankings$method_id), method = "radix"),
                 sort(expected_methods, method = "radix")) ||
      !identical(sort(unique(as.integer(rankings$seed))),
                 sort(as.integer(expected_seeds)))) {
    stop("MV5-J rankings failed immutable-content validation.", call. = FALSE)
  }
  ids <- unique(c(rankings$query_sample_id, rankings$training_sample_id))
  if (!setequal(ids, labels$sample_id)) {
    stop("Ranking sample IDs do not equal the frozen 90-sample label set.",
         call. = FALSE)
  }
  keys <- interaction(
    rankings$fold_id, rankings$seed, rankings$method_id,
    rankings$query_sample_id, drop = TRUE, lex.order = TRUE
  )
  groups <- split(seq_len(nrow(rankings)), keys)
  valid_group <- vapply(groups, function(index) {
    part <- rankings[index, , drop = FALSE]
    ord <- order(part$neighbor_rank)
    part <- part[ord, , drop = FALSE]
    identical(as.integer(part$neighbor_rank), seq_len(nrow(part))) &&
      !anyDuplicated(part$training_sample_id) &&
      all(diff(part$distance) >= 0) &&
      identical(
        part$training_sample_id,
        part$training_sample_id[order(part$distance, part$training_sample_id,
                                      method = "radix")]
      )
  }, logical(1L))
  if (!all(valid_group)) {
    stop("MV5-J rankings are not canonical immutable rank sequences.",
         call. = FALSE)
  }
  invisible(TRUE)
}

mv05k_evaluate_retrieval_v1 <- function(
    rankings, labels, expected_rows = 176750L,
    expected_methods = .mv05k_methods, expected_seeds = .mv05k_seeds,
    expected_observations = 2250L, require_all_estimable = TRUE) {
  .mv05k_validate_rankings(
    rankings, labels, expected_rows, expected_methods, expected_seeds
  )
  label_index <- match(rankings$training_sample_id, labels$sample_id)
  query_index <- match(rankings$query_sample_id, labels$sample_id)
  training_tissue <- labels$tissue[label_index]
  training_study <- labels$study[label_index]
  query_tissue <- labels$tissue[query_index]
  query_study <- labels$study[query_index]
  held_out <- sub("^large_loso_v1:", "", rankings$fold_id)
  if (any(query_study != held_out) || any(training_study == held_out)) {
    stop("Ranking fold membership violates held-out study isolation.",
         call. = FALSE)
  }
  rankings$.training_tissue <- training_tissue
  keys <- interaction(
    rankings$fold_id, rankings$seed, rankings$method_id,
    rankings$query_sample_id, drop = TRUE, lex.order = TRUE
  )
  groups <- split(seq_len(nrow(rankings)), keys)
  observations <- do.call(rbind, lapply(groups, function(index) {
    part <- rankings[index, , drop = FALSE]
    part <- part[order(part$neighbor_rank), , drop = FALSE]
    q_index <- match(part$query_sample_id[[1L]], labels$sample_id)
    q_tissue <- labels$tissue[[q_index]]
    q_study <- labels$study[[q_index]]
    same <- which(part$.training_tissue == q_tissue)
    eligible <- q_tissue %in% .mv05k_eligible_tissues
    status <- if (!eligible) {
      "single_study_tissue_not_estimable"
    } else if (!length(same)) {
      "training_tissue_absent_not_estimable"
    } else {
      "estimable"
    }
    first <- if (length(same)) same[[1L]] else NA_integer_
    data.frame(
      contract_id = "mv05k_query_endpoint_v1",
      fold_id = part$fold_id[[1L]],
      held_out_study = q_study,
      seed = as.integer(part$seed[[1L]]),
      method_id = part$method_id[[1L]],
      method_role = unname(.mv05k_method_roles[part$method_id[[1L]]]),
      query_sample_id = part$query_sample_id[[1L]],
      query_tissue = q_tissue,
      training_samples = nrow(part),
      training_studies = length(unique(training_study[index])),
      first_same_tissue_sample_id = if (is.na(first)) NA_character_ else
        part$training_sample_id[[first]],
      first_same_tissue_rank = if (is.na(first)) NA_integer_ else
        as.integer(part$neighbor_rank[[first]]),
      reciprocal_rank = if (is.na(first)) NA_real_ else
        1 / as.numeric(part$neighbor_rank[[first]]),
      nearest_sample_id = part$training_sample_id[[1L]],
      nearest_tissue = part$.training_tissue[[1L]],
      one_nn_correct = if (identical(status, "estimable"))
        part$.training_tissue[[1L]] == q_tissue else NA,
      nearest_distance_tied = as.logical(part$distance_tied[[1L]]),
      endpoint_status = status,
      labels_opened_after_prediction_lock = TRUE,
      upstream_refit = FALSE,
      reranked_after_label_open = FALSE,
      stringsAsFactors = FALSE
    )
  }))
  rownames(observations) <- NULL
  observations <- observations[order(
    match(observations$method_id, .mv05k_methods), observations$seed,
    observations$query_tissue, observations$held_out_study,
    observations$query_sample_id, method = "radix"
  ), , drop = FALSE]
  if (nrow(observations) != expected_observations ||
      (require_all_estimable &&
       sum(observations$endpoint_status == "estimable") !=
         expected_observations)) {
    stop("Unexpected MV5-K endpoint estimability/completion.", call. = FALSE)
  }
  observations
}

.mv05k_group_summary <- function(observations, group_fields, contract_id) {
  parts <- split(
    observations,
    interaction(observations[, group_fields, drop = FALSE], drop = TRUE,
                lex.order = TRUE)
  )
  result <- do.call(rbind, lapply(parts, function(part) {
    values <- part[1L, group_fields, drop = FALSE]
    cbind(
      data.frame(contract_id = contract_id, stringsAsFactors = FALSE), values,
      data.frame(
        mean_reciprocal_rank = mean(part$reciprocal_rank),
        one_nn_balanced_accuracy = mean(as.numeric(part$one_nn_correct)),
        observations = nrow(part),
        samples = length(unique(part$query_sample_id)),
        studies = length(unique(part$held_out_study)),
        seeds = length(unique(part$seed)),
        failures = sum(part$endpoint_status != "estimable"),
        stringsAsFactors = FALSE
      )
    )
  }))
  rownames(result) <- NULL
  result
}

mv05k_summarize_retrieval_v1 <- function(observations) {
  estimable <- observations[observations$endpoint_status == "estimable", , drop = FALSE]
  tissue_seed <- .mv05k_group_summary(
    estimable, c("method_id", "method_role", "seed", "query_tissue"),
    "mv05k_tissue_seed_summary_v1"
  )
  seed_macro <- do.call(rbind, lapply(split(
    tissue_seed,
    interaction(tissue_seed$method_id, tissue_seed$seed,
                drop = TRUE, lex.order = TRUE)
  ), function(part) {
    data.frame(
      contract_id = "mv05k_seed_macro_summary_v1",
      method_id = part$method_id[[1L]], method_role = part$method_role[[1L]],
      seed = part$seed[[1L]],
      macro_mean_reciprocal_rank = mean(part$mean_reciprocal_rank),
      macro_one_nn_balanced_accuracy = mean(part$one_nn_balanced_accuracy),
      tissues = nrow(part), samples = sum(part$samples),
      studies = length(unique(estimable$held_out_study[
        estimable$method_id == part$method_id[[1L]] &
          estimable$seed == part$seed[[1L]]
      ])),
      failures = sum(part$failures), stringsAsFactors = FALSE
    )
  }))
  sample_summary <- .mv05k_group_summary(
    estimable,
    c("method_id", "method_role", "query_sample_id", "query_tissue",
      "held_out_study"),
    "mv05k_sample_summary_v1"
  )
  tissue_summary <- .mv05k_group_summary(
    estimable, c("method_id", "method_role", "query_tissue"),
    "mv05k_tissue_summary_v1"
  )
  method_summary <- do.call(rbind, lapply(split(
    tissue_summary, tissue_summary$method_id
  ), function(part) {
    data.frame(
      contract_id = "mv05k_method_macro_summary_v1",
      method_id = part$method_id[[1L]], method_role = part$method_role[[1L]],
      macro_mean_reciprocal_rank = mean(part$mean_reciprocal_rank),
      macro_one_nn_balanced_accuracy = mean(part$one_nn_balanced_accuracy),
      tissues = nrow(part), samples = sum(part$samples),
      studies = length(unique(estimable$held_out_study[
        estimable$method_id == part$method_id[[1L]]
      ])),
      completed_seeds = length(unique(estimable$seed[
        estimable$method_id == part$method_id[[1L]]
      ])),
      failures = sum(part$failures), stringsAsFactors = FALSE
    )
  }))
  order_method <- function(x) {
    x[order(match(x$method_id, .mv05k_methods), method = "radix"), , drop = FALSE]
  }
  list(
    tissue_seed = tissue_seed[order(
      match(tissue_seed$method_id, .mv05k_methods), tissue_seed$seed,
      tissue_seed$query_tissue, method = "radix"), , drop = FALSE],
    seed_macro = seed_macro[order(
      match(seed_macro$method_id, .mv05k_methods), seed_macro$seed,
      method = "radix"), , drop = FALSE],
    sample = sample_summary[order(
      match(sample_summary$method_id, .mv05k_methods),
      sample_summary$query_tissue, sample_summary$held_out_study,
      sample_summary$query_sample_id, method = "radix"), , drop = FALSE],
    tissue = tissue_summary[order(
      match(tissue_summary$method_id, .mv05k_methods),
      tissue_summary$query_tissue, method = "radix"), , drop = FALSE],
    method = order_method(method_summary)
  )
}

.mv05k_metric_column <- function(endpoint_id) {
  switch(endpoint_id,
         cross_study_tissue_mrr_v1 = "mean_reciprocal_rank",
         cross_study_tissue_1nn_balanced_accuracy_v1 =
           "one_nn_balanced_accuracy",
         stop("Unknown MV5-K endpoint ID.", call. = FALSE))
}

.mv05k_macro_from_sample_rows <- function(rows, metric) {
  tissue <- tapply(rows[[metric]], rows$query_tissue, mean)
  mean(tissue[sort(names(tissue), method = "radix")])
}

mv05k_block_inference_v1 <- function(
    sample_summary, bootstrap_replicates = 2000L,
    bootstrap_seed = 20260808L, randomization_replicates = 9999L,
    randomization_seed = 20260809L) {
  required <- c(
    "method_id", "query_sample_id", "query_tissue", "held_out_study",
    "mean_reciprocal_rank", "one_nn_balanced_accuracy", "seeds"
  )
  if (!is.data.frame(sample_summary) ||
      !all(required %in% names(sample_summary)) ||
      any(sample_summary$seeds != length(.mv05k_seeds))) {
    stop("Sample summaries are incomplete for MV5-K inference.", call. = FALSE)
  }
  methods <- .mv05k_methods
  sample_ids <- sort(unique(sample_summary$query_sample_id), method = "radix")
  if (length(sample_ids) != 90L ||
      nrow(sample_summary) != length(sample_ids) * length(methods)) {
    stop("MV5-K inference requires a complete 90-sample by five-method table.",
         call. = FALSE)
  }
  base <- sample_summary[sample_summary$method_id == methods[[1L]],
                         c("query_sample_id", "query_tissue", "held_out_study")]
  base <- base[match(sample_ids, base$query_sample_id), , drop = FALSE]
  metric_arrays <- lapply(c(
    "cross_study_tissue_mrr_v1",
    "cross_study_tissue_1nn_balanced_accuracy_v1"
  ), function(endpoint) {
    metric <- .mv05k_metric_column(endpoint)
    matrix(vapply(methods, function(method) {
      part <- sample_summary[sample_summary$method_id == method, , drop = FALSE]
      part[[metric]][match(sample_ids, part$query_sample_id)]
    }, numeric(length(sample_ids))), nrow = length(sample_ids),
    dimnames = list(sample_ids, methods))
  })
  names(metric_arrays) <- c(
    "cross_study_tissue_mrr_v1",
    "cross_study_tissue_1nn_balanced_accuracy_v1"
  )
  strata <- split(unique(base[, c("query_tissue", "held_out_study")]),
                  unique(base[, c("query_tissue", "held_out_study")])$query_tissue)
  strata <- strata[sort(names(strata), method = "radix")]
  bootstrap <- .mv05k_with_seed(bootstrap_seed, {
    result <- vector("list", length(metric_arrays))
    for (endpoint_index in seq_along(metric_arrays)) {
      values <- metric_arrays[[endpoint_index]]
      estimates <- matrix(NA_real_, nrow = bootstrap_replicates,
                          ncol = length(methods), dimnames = list(NULL, methods))
      for (replicate_id in seq_len(bootstrap_replicates)) {
        tissue_values <- matrix(NA_real_, nrow = length(strata),
                                ncol = length(methods))
        for (tissue_index in seq_along(strata)) {
          study_ids <- sort(strata[[tissue_index]]$held_out_study,
                            method = "radix")
          drawn <- sample(study_ids, length(study_ids), replace = TRUE)
          row_index <- unlist(lapply(drawn, function(study) {
            which(base$query_tissue == names(strata)[[tissue_index]] &
                    base$held_out_study == study)
          }), use.names = FALSE)
          tissue_values[tissue_index, ] <- colMeans(values[row_index, , drop = FALSE])
        }
        estimates[replicate_id, ] <- colMeans(tissue_values)
      }
      result[[endpoint_index]] <- estimates
    }
    result
  })
  names(bootstrap) <- names(metric_arrays)
  method_intervals <- do.call(rbind, lapply(names(metric_arrays), function(endpoint) {
    values <- metric_arrays[[endpoint]]
    point <- vapply(seq_along(methods), function(index) {
      rows <- base
      rows$value <- values[, index]
      .mv05k_macro_from_sample_rows(rows, "value")
    }, numeric(1L))
    do.call(rbind, lapply(seq_along(methods), function(index) {
      interval <- stats::quantile(bootstrap[[endpoint]][, index],
                                  c(0.025, 0.975), names = FALSE, type = 7)
      data.frame(
        contract_id = "mv05k_method_interval_v1", endpoint_id = endpoint,
        method_id = methods[[index]], method_role = unname(.mv05k_method_roles[index]),
        estimate = point[[index]], ci_lower = interval[[1L]],
        ci_upper = interval[[2L]], bootstrap_replicates = bootstrap_replicates,
        bootstrap_seed = bootstrap_seed, tissues = length(strata),
        study_blocks = length(unique(base$held_out_study)), samples = nrow(base),
        completed_seeds = length(.mv05k_seeds),
        inference_status = "estimable", stringsAsFactors = FALSE
      )
    }))
  }))
  contrast_methods <- c("integrated_cell_landscape_h0_v1", "integrated_cell_landscape_h1_v1")
  energy <- "integrated_cell_distribution_energy_v1"
  contrasts <- do.call(rbind, lapply(names(metric_arrays), function(endpoint) {
    values <- metric_arrays[[endpoint]]
    do.call(rbind, lapply(contrast_methods, function(method) {
      method_index <- match(method, methods)
      energy_index <- match(energy, methods)
      sample_difference <- values[, method_index] - values[, energy_index]
      rows <- base
      rows$value <- sample_difference
      point <- .mv05k_macro_from_sample_rows(rows, "value")
      boot_difference <- bootstrap[[endpoint]][, method_index] -
        bootstrap[[endpoint]][, energy_index]
      interval <- stats::quantile(boot_difference, c(0.025, 0.975),
                                  names = FALSE, type = 7)
      data.frame(
        contract_id = "mv05k_paired_contrast_v1",
        family_id = if (endpoint == "cross_study_tissue_mrr_v1")
          "F1_primary_retrieval" else "none_supportive",
        endpoint_id = endpoint, method_id = method, baseline_method_id = energy,
        contrast = "topology_minus_matched_energy",
        estimate = point, ci_lower = interval[[1L]], ci_upper = interval[[2L]],
        favorable_direction = "positive",
        bootstrap_replicates = bootstrap_replicates,
        bootstrap_seed = bootstrap_seed,
        randomization_replicates = if (endpoint ==
          "cross_study_tissue_mrr_v1") randomization_replicates else NA_integer_,
        randomization_seed = if (endpoint ==
          "cross_study_tissue_mrr_v1") randomization_seed else NA_integer_,
        tissues = length(strata),
        study_blocks = length(unique(base$held_out_study)), samples = nrow(base),
        completed_seeds = length(.mv05k_seeds),
        raw_p_value = NA_real_, holm_p_value = NA_real_,
        inference_status = if (length(unique(base$held_out_study)) >= 4L)
          "estimable" else "insufficient_independent_study_blocks",
        stringsAsFactors = FALSE
      )
    }))
  }))
  primary_rows <- which(contrasts$family_id == "F1_primary_retrieval" &
                          contrasts$inference_status == "estimable")
  sign_audit <- data.frame()
  if (length(primary_rows)) {
    study_ids <- sort(unique(base$held_out_study), method = "radix")
    signs <- .mv05k_with_seed(randomization_seed, matrix(
      sample(c(-1, 1), randomization_replicates * length(study_ids),
             replace = TRUE),
      nrow = randomization_replicates, ncol = length(study_ids),
      dimnames = list(NULL, study_ids)
    ))
    null_digest <- character(length(primary_rows))
    exceedance_count <- integer(length(primary_rows))
    for (position in seq_along(primary_rows)) {
      row_index <- primary_rows[[position]]
      method <- contrasts$method_id[[row_index]]
      values <- metric_arrays[["cross_study_tissue_mrr_v1"]]
      differences <- values[, match(method, methods)] -
        values[, match(energy, methods)]
      null <- vapply(seq_len(randomization_replicates), function(replicate_id) {
        rows <- base
        rows$value <- differences * signs[
          replicate_id, match(rows$held_out_study, study_ids)
        ]
        .mv05k_macro_from_sample_rows(rows, "value")
      }, numeric(1L))
      exceedances <- .mv05k_randomization_exceedance_count_v1(
        null, contrasts$estimate[[row_index]]
      )
      exceedance_count[[position]] <- exceedances
      contrasts$raw_p_value[[row_index]] <-
        mv05_monte_carlo_p_v1(exceedances, randomization_replicates)
      null_digest[[position]] <- digest::digest(
        null, algo = "sha256", serialize = TRUE
      )
    }
    contrasts$holm_p_value[primary_rows] <- stats::p.adjust(
      contrasts$raw_p_value[primary_rows], method = "holm"
    )
    sign_audit <- data.frame(
      contract_id = "mv05k_randomization_audit_v1",
      family_id = "F1_primary_retrieval",
      method_id = contrasts$method_id[primary_rows],
      baseline_method_id = energy,
      randomization_replicates = randomization_replicates,
      randomization_seed = randomization_seed,
      study_blocks = length(study_ids),
      sign_matrix_sha256 = digest::digest(signs, algo = "sha256", serialize = TRUE),
      null_distribution_sha256 = null_digest,
      exceedance_count = exceedance_count,
      boundary_policy = "absolute_two_sided_64eps_ties_count_as_exceedance_v1",
      epsilon_multiplier = 64,
      two_sided = TRUE,
      exact_zero_p_value_possible = FALSE,
      stringsAsFactors = FALSE
    )
  }
  bootstrap_audit <- data.frame(
    contract_id = "mv05k_bootstrap_audit_v1",
    endpoint_id = names(bootstrap),
    bootstrap_replicates = bootstrap_replicates,
    bootstrap_seed = bootstrap_seed,
    tissues = length(strata),
    study_blocks = length(unique(base$held_out_study)),
    replicate_matrix_sha256 = vapply(bootstrap, digest::digest, character(1L),
                                      algo = "sha256", serialize = TRUE),
    resampling_unit = "tissue_stratified_held_out_study_block",
    seeds_treated_as_independent = FALSE,
    stringsAsFactors = FALSE
  )
  list(method_intervals = method_intervals, contrasts = contrasts,
       bootstrap_audit = bootstrap_audit,
       randomization_audit = sign_audit)
}
