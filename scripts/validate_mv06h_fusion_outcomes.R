#!/usr/bin/env Rscript

options(warn = 2)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest is required.")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) stop(
  paste("usage: validate_mv06h_fusion_outcomes.R METADATA STAGE1_ROOT",
        "COMPLETION_ROOT LOCK_DIR OUTCOME_DIR LOCK_COMMIT OUTPUT"),
  call. = FALSE)
metadata <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
stage1 <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
completion <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
lock_dir <- normalizePath(args[[4L]], winslash = "/", mustWork = TRUE)
outcome_dir <- normalizePath(args[[5L]], winslash = "/", mustWork = TRUE)
lock_commit <- args[[6L]]; output <- args[[7L]]
if (file.exists(output)) stop("Refusing to overwrite MV6-H validation.",
                              call. = FALSE)
readc <- function(path) read.csv(path, stringsAsFactors = FALSE,
                                 check.names = FALSE)
sha <- function(path) digest::digest(file = path, algo = "sha256",
                                     serialize = FALSE)
truth <- function(x) tolower(as.character(x)) == "true"
safe <- function(x) gsub("[^A-Za-z0-9._-]", "_", x)
close <- function(a, b) {
  if (length(a) != length(b)) return(Inf)
  if (!length(a)) return(0)
  max(abs(as.numeric(a) - as.numeric(b)))
}
align <- function(expected, observed, fields) {
  key <- function(x) do.call(paste, c(x[fields], sep = "\r"))
  observed[match(key(expected), key(observed)), , drop = FALSE]
}
methods <- c("cell_H0", "cell_H1", "gene_H0", "gene_H1",
  "cell_composite", "fusion_gene_weight_025", "fusion_gene_weight_050",
  "fusion_gene_weight_075", "gene_composite")
primary <- "fusion_gene_weight_050"
comparators <- c("cell_composite", "gene_composite")

lock <- readc(file.path(lock_dir, "mv06h-prediction-lock.csv"))
groups <- readc(file.path(lock_dir, "mv06h-prediction-group-manifest.csv"))
receipt <- readc(file.path(outcome_dir, "mv06h-label-open-receipt.csv"))
if (nrow(lock) != 1L || nrow(groups) != 75L || nrow(receipt) != 1L ||
    receipt$lock_commit != lock_commit ||
    receipt$prediction_lock_sha256 != sha(file.path(
      lock_dir, "mv06h-prediction-lock.csv")) ||
    receipt$prediction_manifest_root_sha256 !=
      lock$prediction_manifest_root_sha256 ||
    !truth(receipt$predictions_loaded_and_hash_verified) ||
    !truth(receipt$receipt_written_before_metadata_hash_or_read) ||
    truth(receipt$labels_opened) || truth(receipt$outcomes_computed)) {
  stop("MV6-H durable label-open receipt failed.", call. = FALSE)
}
if (sha(metadata) != lock$metadata_expected_sha256) {
  stop("MV6-H authoritative metadata hash drift.", call. = FALSE)
}
raw <- read.csv(metadata, stringsAsFactors = FALSE, check.names = TRUE)
required <- c("orig.ident", "SRA", "Tissue.x", "Approach.x",
              "Number_of_Cells_After_Filtering")
if (!all(required %in% names(raw))) stop("Metadata schema drift.", call. = FALSE)
labels <- data.frame(
  sample_id = trimws(as.character(raw$orig.ident)),
  study = trimws(as.character(raw$SRA)),
  tissue = tolower(trimws(as.character(raw$Tissue.x))),
  approach = trimws(as.character(raw$Approach.x)), stringsAsFactors = FALSE)
eligible <- c("bone marrow", "colon", "liver", "pbmc", "testis")
labels <- labels[labels$tissue %in% eligible, , drop = FALSE]
labels <- labels[order(labels$sample_id, method = "radix"), ]
if (nrow(labels) != 90L || length(unique(labels$study)) != 15L ||
    anyDuplicated(labels$sample_id)) stop("Candidate label axis drift.", call. = FALSE)

pieces <- vector("list", nrow(groups)); post_rows <- vector("list", nrow(groups))
for (i in seq_len(nrow(groups))) {
  row <- groups[i, , drop = FALSE]
  root <- if (row$group_root_kind == "accepted_stage1_sentinel") stage1 else completion
  path <- file.path(root, safe(row$group_id), "rankings.csv")
  actual <- sha(path)
  if (actual != row$rankings_sha256) stop("Prediction hash drift.", call. = FALSE)
  rankings <- readc(path)
  if (nrow(rankings) != as.integer(row$ranking_rows) ||
      !setequal(rankings$method_id, methods) ||
      any(rankings$outcome_label_state != "closed") ||
      any(truth(rankings$biological_outcomes_computed))) {
    stop("Prediction group schema/firewall drift.", call. = FALSE)
  }
  train_index <- match(rankings$training_sample_id, labels$sample_id)
  query_index <- match(rankings$query_sample_id, labels$sample_id)
  heldout <- sub("^large_loso_v1:", "", rankings$fold_id)
  if (anyNA(train_index) || anyNA(query_index) ||
      any(labels$study[query_index] != heldout) ||
      any(labels$study[train_index] == heldout)) {
    stop("Independent LOSO label join failed.", call. = FALSE)
  }
  rankings$.training_tissue <- labels$tissue[train_index]
  parts <- split(seq_len(nrow(rankings)), interaction(
    rankings$method_id, rankings$query_sample_id,
    drop = TRUE, lex.order = TRUE))
  pieces[[i]] <- do.call(rbind, lapply(parts, function(ix) {
    part <- rankings[ix, , drop = FALSE]
    part <- part[order(as.integer(part$rank)), , drop = FALSE]
    q <- match(part$query_sample_id[[1L]], labels$sample_id)
    same <- which(part$.training_tissue == labels$tissue[[q]])
    if (!length(same)) stop("Unexpected nonestimable query.", call. = FALSE)
    first <- same[[1L]]
    data.frame(
      group_id = part$group_id[[1L]], method_id = part$method_id[[1L]],
      query_sample_id = part$query_sample_id[[1L]],
      held_out_study = labels$study[[q]], query_tissue = labels$tissue[[q]],
      first_same_tissue_sample_id = part$training_sample_id[[first]],
      first_same_tissue_rank = as.integer(part$rank[[first]]),
      reciprocal_rank = 1 / as.numeric(part$rank[[first]]),
      nearest_sample_id = part$training_sample_id[[1L]],
      nearest_tissue = part$.training_tissue[[1L]],
      one_nn_correct = part$.training_tissue[[1L]] == labels$tissue[[q]],
      stringsAsFactors = FALSE)
  }))
  post_rows[[i]] <- data.frame(group_id = row$group_id,
    expected_sha256 = row$rankings_sha256, actual_sha256 = actual,
    matched = TRUE, stringsAsFactors = FALSE)
}
independent_query <- do.call(rbind, pieces); rownames(independent_query) <- NULL
query <- readc(file.path(outcome_dir, "mv06h-query-method-outcomes.csv"))
keys <- c("group_id", "method_id", "query_sample_id")
query <- align(independent_query, query, keys)
if (nrow(independent_query) != 4050L || anyNA(query$query_sample_id) ||
    !identical(query$first_same_tissue_sample_id,
               independent_query$first_same_tissue_sample_id) ||
    !identical(as.integer(query$first_same_tissue_rank),
               independent_query$first_same_tissue_rank) ||
    !identical(query$nearest_sample_id, independent_query$nearest_sample_id) ||
    !identical(query$nearest_tissue, independent_query$nearest_tissue) ||
    !identical(as.logical(query$one_nn_correct),
               independent_query$one_nn_correct) ||
    close(query$reciprocal_rank, independent_query$reciprocal_rank) > 1e-14 ||
    any(query$endpoint_status != "estimable") || any(query$upstream_refit) ||
    any(query$reranked_after_label_open) || any(query$method_selection_executed)) {
  stop("Independent MV6-H query endpoint reconstruction failed.", call. = FALSE)
}

sample_parts <- split(seq_len(nrow(independent_query)), interaction(
  independent_query$method_id, independent_query$query_sample_id,
  drop = TRUE, lex.order = TRUE))
independent_sample <- do.call(rbind, lapply(sample_parts, function(ix) {
  part <- independent_query[ix, , drop = FALSE]
  data.frame(method_id = part$method_id[[1L]],
    query_sample_id = part$query_sample_id[[1L]],
    query_tissue = part$query_tissue[[1L]],
    held_out_study = part$held_out_study[[1L]],
    mean_reciprocal_rank = mean(part$reciprocal_rank),
    one_nn_balanced_accuracy = mean(as.numeric(part$one_nn_correct)),
    seeds = nrow(part), stringsAsFactors = FALSE)
}))
sample_observed <- readc(file.path(outcome_dir, "mv06h-sample-method-summaries.csv"))
sample_observed <- align(independent_sample, sample_observed,
                         c("method_id", "query_sample_id"))
sample_diff <- max(close(sample_observed$mean_reciprocal_rank,
                         independent_sample$mean_reciprocal_rank),
                   close(sample_observed$one_nn_balanced_accuracy,
                         independent_sample$one_nn_balanced_accuracy))
if (nrow(independent_sample) != 810L || any(independent_sample$seeds != 5L) ||
    sample_diff > 1e-14) stop("Sample aggregation reconstruction failed.",
                              call. = FALSE)

tissue_parts <- split(seq_len(nrow(independent_sample)), interaction(
  independent_sample$method_id, independent_sample$query_tissue,
  drop = TRUE, lex.order = TRUE))
independent_tissue <- do.call(rbind, lapply(tissue_parts, function(ix) {
  part <- independent_sample[ix, , drop = FALSE]
  data.frame(method_id = part$method_id[[1L]],
    query_tissue = part$query_tissue[[1L]],
    mean_reciprocal_rank = mean(part$mean_reciprocal_rank),
    one_nn_balanced_accuracy = mean(part$one_nn_balanced_accuracy),
    samples = nrow(part), studies = length(unique(part$held_out_study)),
    stringsAsFactors = FALSE)
}))
tissue_observed <- readc(file.path(outcome_dir, "mv06h-tissue-method-summaries.csv"))
tissue_observed <- align(independent_tissue, tissue_observed,
                         c("method_id", "query_tissue"))
tissue_diff <- max(close(tissue_observed$mean_reciprocal_rank,
                         independent_tissue$mean_reciprocal_rank),
                   close(tissue_observed$one_nn_balanced_accuracy,
                         independent_tissue$one_nn_balanced_accuracy))
method_independent <- do.call(rbind, lapply(methods, function(method) {
  part <- independent_tissue[independent_tissue$method_id == method, ]
  data.frame(method_id = method,
    macro_mean_reciprocal_rank = mean(part$mean_reciprocal_rank),
    macro_one_nn_balanced_accuracy = mean(part$one_nn_balanced_accuracy),
    stringsAsFactors = FALSE)
}))
method_observed <- readc(file.path(outcome_dir, "mv06h-method-summaries.csv"))
method_observed <- align(method_independent, method_observed, "method_id")
method_diff <- max(close(method_observed$macro_mean_reciprocal_rank,
                         method_independent$macro_mean_reciprocal_rank),
                   close(method_observed$macro_one_nn_balanced_accuracy,
                         method_independent$macro_one_nn_balanced_accuracy))
if (tissue_diff > 1e-14 || method_diff > 1e-14) {
  stop("Tissue/method macro aggregation reconstruction failed.", call. = FALSE)
}

sample_ids <- sort(unique(independent_sample$query_sample_id), method = "radix")
base <- independent_sample[independent_sample$method_id == methods[[1L]],
  c("query_sample_id", "query_tissue", "held_out_study")]
base <- base[match(sample_ids, base$query_sample_id), , drop = FALSE]
arrays <- list(
  cross_study_tissue_mrr_v1 = matrix(vapply(methods, function(method) {
    part <- independent_sample[independent_sample$method_id == method, ]
    part$mean_reciprocal_rank[match(sample_ids, part$query_sample_id)]
  }, numeric(90L)), nrow = 90L, dimnames = list(sample_ids, methods)),
  cross_study_tissue_1nn_balanced_accuracy_v1 = matrix(vapply(methods,
    function(method) {
      part <- independent_sample[independent_sample$method_id == method, ]
      part$one_nn_balanced_accuracy[match(sample_ids, part$query_sample_id)]
    }, numeric(90L)), nrow = 90L, dimnames = list(sample_ids, methods)))
strata_table <- unique(base[c("query_tissue", "held_out_study")])
strata <- split(strata_table, strata_table$query_tissue)
strata <- strata[sort(names(strata), method = "radix")]
RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(20260815L)
bootstrap <- lapply(arrays, function(values) {
  estimates <- matrix(NA_real_, 2000L, length(methods),
                      dimnames = list(NULL, methods))
  for (b in 1:2000) {
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
})
bootstrap_audit <- readc(file.path(outcome_dir, "mv06h-bootstrap-audit.csv"))
expected_boot_hash <- vapply(bootstrap, digest::digest, character(1L),
                             algo = "sha256", serialize = TRUE)
bootstrap_audit <- bootstrap_audit[match(names(bootstrap),
  bootstrap_audit$endpoint_id), , drop = FALSE]
if (!identical(expected_boot_hash,
    stats::setNames(bootstrap_audit$replicate_matrix_sha256,
                    bootstrap_audit$endpoint_id)) ||
    any(bootstrap_audit$bootstrap_replicates != 2000L) ||
    any(bootstrap_audit$bootstrap_seed != 20260815L) ||
    any(truth(bootstrap_audit$seeds_treated_as_independent))) {
  stop("Blocked bootstrap reconstruction failed.", call. = FALSE)
}

intervals <- readc(file.path(outcome_dir, "mv06h-method-intervals.csv"))
interval_diff <- 0
for (endpoint in names(arrays)) for (method in methods) {
  values <- arrays[[endpoint]][, method]
  point <- mean(tapply(values, base$query_tissue, mean))
  ci <- quantile(bootstrap[[endpoint]][, method], c(.025, .975),
                 names = FALSE, type = 7)
  observed <- intervals[intervals$endpoint_id == endpoint &
                          intervals$method_id == method, ]
  interval_diff <- max(interval_diff, close(observed$estimate, point),
                       close(c(observed$ci_lower, observed$ci_upper), ci))
}
if (nrow(intervals) != 18L || interval_diff > 1e-14) {
  stop("Method interval reconstruction failed.", call. = FALSE)
}

contrasts <- readc(file.path(outcome_dir, "mv06h-primary-contrasts.csv"))
fusion_index <- match(primary, methods); contrast_diff <- 0
study_ids <- sort(unique(base$held_out_study), method = "radix")
RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(20260816L)
signs <- matrix(sample(c(-1, 1), 9999L * length(study_ids), replace = TRUE),
  nrow = 9999L, ncol = length(study_ids), dimnames = list(NULL, study_ids))
raw_p <- numeric(2L); null_hash <- character(2L)
for (i in seq_along(comparators)) {
  comparator <- comparators[[i]]
  differences <- arrays$cross_study_tissue_mrr_v1[, fusion_index] -
    arrays$cross_study_tissue_mrr_v1[, match(comparator, methods)]
  point <- mean(tapply(differences, base$query_tissue, mean))
  boot <- bootstrap$cross_study_tissue_mrr_v1[, fusion_index] -
    bootstrap$cross_study_tissue_mrr_v1[, match(comparator, methods)]
  ci <- quantile(boot, c(.025, .975), names = FALSE, type = 7)
  observed <- contrasts[contrasts$comparator_id == comparator, ]
  contrast_diff <- max(contrast_diff, close(observed$estimate, point),
                       close(c(observed$ci_lower, observed$ci_upper), ci))
  null <- vapply(1:9999, function(b) mean(tapply(
    differences * signs[b, match(base$held_out_study, study_ids)],
    base$query_tissue, mean)), numeric(1L))
  raw_p[[i]] <- (sum(abs(null) >= abs(point)) + 1) / 10000
  null_hash[[i]] <- digest::digest(null, algo = "sha256", serialize = TRUE)
}
contrasts <- contrasts[match(comparators, contrasts$comparator_id), ]
p_diff <- max(close(contrasts$raw_p_value, raw_p),
              close(contrasts$holm_p_value, p.adjust(raw_p, "holm")))
randomization <- readc(file.path(outcome_dir, "mv06h-randomization-audit.csv"))
randomization <- randomization[match(comparators,
  randomization$comparator_id), , drop = FALSE]
if (nrow(contrasts) != 2L || contrast_diff > 1e-14 || p_diff > 1e-14 ||
    any(randomization$sign_matrix_sha256 != digest::digest(
      signs, algo = "sha256", serialize = TRUE)) ||
    !identical(randomization$null_distribution_sha256, null_hash)) {
  stop("Contrast/randomization/Holm reconstruction failed.", call. = FALSE)
}

post <- readc(file.path(outcome_dir, "mv06h-prediction-postcheck.csv"))
expected_post <- do.call(rbind, post_rows)
post <- align(expected_post, post, "group_id")
production <- readc(file.path(outcome_dir, "mv06h-production-summary.csv"))
decision <- readc(file.path(outcome_dir, "mv06h-decision.csv"))
expected_disposition <- if (!all(contrasts$estimate > 0))
  "no_reliable_outperformance_report_views_separately" else if
    (all(contrasts$holm_p_value <= .05))
  "blocked_primary_supports_fusion_pending_robustness_synthesis" else
  "directionally_positive_but_not_confirmatory_report_views_separately"
if (!all(truth(post$matched)) ||
    any(post$actual_sha256 != expected_post$actual_sha256) ||
    production$ranking_rows != 318150L || production$query_method_outcomes != 4050L ||
    production$upstream_refits != 0L || production$reranking_operations != 0L ||
    production$method_selection_operations != 0L ||
    production$advanced_fusion_jobs != 0L || production$clustering_jobs != 0L ||
    decision$disposition != expected_disposition ||
    truth(decision$advanced_fusion_authorized) || truth(decision$clustering_authorized) ||
    truth(decision$outcome_driven_selection_executed)) {
  stop("Post-label firewall or disposition reconstruction failed.", call. = FALSE)
}

validation <- data.frame(
  contract_id = "mv06h_outcome_independent_validation_v1",
  validation_id = c(
    "durable_prediction_lock_receipt", "authoritative_metadata_hash_schema",
    "complete_prediction_hash_postcheck", "independent_LOSO_label_join",
    "query_MRR_and_1NN_reconstruction", "five_seed_sample_aggregation",
    "equal_tissue_macro_aggregation", "blocked_bootstrap_matrix",
    "all_method_intervals", "two_primary_MRR_contrasts",
    "paired_sign_matrix_and_nulls", "Holm_two_test_family",
    "no_refit_rerank_or_selection", "advanced_fusion_and_clustering_closed",
    "prespecified_negative_disposition"),
  passed = TRUE,
  rows_checked = c(1L, nrow(raw), 75L, 318150L, 4050L, 810L, 45L,
                   4000L, 18L, 2L, 19998L, 2L, 1L, 1L, 1L),
  maximum_absolute_difference = c(0, 0, 0, 0, 0, sample_diff,
    max(tissue_diff, method_diff), 0, interval_diff, contrast_diff, 0,
    p_diff, 0, 0, 0),
  evidence = c(
    paste0("lock_commit=", lock_commit), "124_rows_90_candidates_15_studies",
    "75_of_75_ranking_hashes", "query_heldout_training_exclusion_exact",
    "4050_of_4050", "5_seeds_within_each_of_810_sample_methods",
    "45_tissue_methods_9_equal_five_tissue_macros",
    "two_2000_replicate_matrices_reproduced_by_hash",
    "18_of_18", "F050_minus_cell_and_F050_minus_gene",
    "9999_by_15_sign_matrix_and_two_null_hashes",
    "two_raw_and_Holm_adjusted_p_values", "all_zero", "both_false",
    expected_disposition),
  independent_production_helper_called = FALSE,
  stringsAsFactors = FALSE)
write.csv(validation, output, row.names = FALSE, na = "")
message("MV6-H independent outcome validation passed: 15/15 disposition=",
        expected_disposition)
