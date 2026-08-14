#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(
    "usage: validate_mv05k_retrieval_evaluation.R FROZEN_METADATA_CSV ",
    "RESULT_DIR OUTPUT_CSV", call. = FALSE
  )
}
metadata_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
result_dir <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
output_path <- args[[3L]]
source("R/provenance_utils.R")

file_sha <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}
read_result <- function(stem) {
  utils::read.csv(
    file.path(result_dir, paste0(stem, "-2026-08-10.csv")),
    stringsAsFactors = FALSE, check.names = FALSE
  )
}
max_difference <- function(x, y) {
  if (!length(x) && !length(y)) return(0)
  max(abs(as.numeric(x) - as.numeric(y)))
}
assert_close <- function(x, y, label, tolerance = 1e-12) {
  difference <- max_difference(x, y)
  if (!is.finite(difference) || difference > tolerance) {
    stop(label, " differs by ", format(difference, digits = 17), ".",
         call. = FALSE)
  }
  difference
}
macro <- function(values, tissues) {
  tissue_values <- tapply(values, tissues, mean)
  mean(tissue_values[sort(names(tissue_values), method = "radix")])
}

methods <- c(
  "integrated_cell_landscape_h0_v1", "integrated_cell_landscape_h1_v1",
  "integrated_cell_landscape_h0_h1_raw_euclidean_v1",
  "integrated_cell_distribution_energy_v1",
  "pseudobulk_training_standardized_panel_v1"
)
seeds <- 20260805:20260809
expected_ranking_sha <-
  "4588902bce89a04cae0c7676b4f21f81e83013a29120ca2a4b39f3ffacfb677e"
expected_metadata_sha <-
  "e2b6f1e7d7dc05ba8e80627aad4ea1677c3b718a844397acd06596d7e22a18d0"
expected_completion_sha <-
  "971e6902c1b6c2b8d39c74a9d07abc9c6d05c99c7b47bb8bd10e4a2209a08ed9"
expected_group_sha <-
  "0fd3e05a312079cebd4e49f4d0f7326cda037fd460772176385d6ea9e0e023ab"
expected_registry_sha <-
  "becb8451948e421916915fa231a50a83adb22fbd62cc10cadb3a6d934fd7de5d"
expected_scale_sha <-
  "9f6f93e73c072088e9917d395dd2f98005bfeb9d4ec6ec6dec33c0d50333e8bf"
expected_assembly_sha <-
  "aabfdb0338826451d08efb27ebacea4eb6b56cad5fc30bd6f9159f60baeff68c"
ranking_path <-
  "docs/audits/mv05j-integrated-cell-retrieval-rankings-2026-08-09.csv.gz"
completion_path <- "docs/audits/mv05j-method-completion-2026-08-09.csv"
group_path <- "docs/audits/mv05j-group-bundle-index-2026-08-09.csv"
registry_path <- "docs/audits/mv05j-method-registry-2026-08-09.csv"
scale_path <- "docs/audits/mv05j-component-scale-disposition-2026-08-09.csv"
assembly_path <- "docs/audits/mv05j-public-assembly-summary-2026-08-09.csv"

if (!identical(file_sha(ranking_path), expected_ranking_sha) ||
    !identical(file_sha(metadata_path), expected_metadata_sha) ||
    !identical(file_sha(completion_path), expected_completion_sha) ||
    !identical(file_sha(group_path), expected_group_sha) ||
    !identical(file_sha(registry_path), expected_registry_sha) ||
    !identical(file_sha(scale_path), expected_scale_sha) ||
    !identical(file_sha(assembly_path), expected_assembly_sha)) {
  stop("Independent source-hash validation failed.", call. = FALSE)
}
completion <- utils::read.csv(
  completion_path, stringsAsFactors = FALSE, check.names = FALSE
)
groups <- utils::read.csv(group_path, stringsAsFactors = FALSE,
                          check.names = FALSE)
registry <- utils::read.csv(registry_path, stringsAsFactors = FALSE,
                            check.names = FALSE)
scale <- utils::read.csv(scale_path, stringsAsFactors = FALSE,
                         check.names = FALSE)
if (nrow(completion) != 375L || nrow(groups) != 75L ||
    nrow(registry) != 5L || !identical(registry$method_id, methods) ||
    nrow(scale) != 225L ||
    any(scale$held_out_query_pairs_used_for_scale != 0L) ||
    any(scale$additional_topology_jobs_executed != 0L) ||
    any(completion$status != "completed") ||
    any(groups$completed_methods != 5L) || any(groups$failed_methods != 0L)) {
  stop("Independent MV5-J completion validation failed.", call. = FALSE)
}

raw_labels <- utils::read.csv(
  metadata_path, stringsAsFactors = FALSE, check.names = TRUE
)
labels <- data.frame(
  sample_id = as.character(raw_labels$orig.ident),
  study = as.character(raw_labels$SRA),
  tissue = tolower(as.character(raw_labels$Tissue.x)),
  stringsAsFactors = FALSE
)
study_counts <- tapply(labels$study, labels$tissue,
                       function(x) length(unique(x)))
eligible_tissues <- sort(names(study_counts)[study_counts >= 2L],
                         method = "radix")
labels <- labels[labels$tissue %in% eligible_tissues, , drop = FALSE]
labels <- labels[order(labels$sample_id, method = "radix"), , drop = FALSE]
if (nrow(labels) != 90L || length(unique(labels$study)) != 15L ||
    !identical(eligible_tissues,
               c("bone marrow", "colon", "liver", "pbmc", "testis")) ||
    any(tapply(labels$tissue, labels$study,
               function(x) length(unique(x))) != 1L)) {
  stop("Independent label eligibility reconstruction failed.", call. = FALSE)
}

rankings <- utils::read.csv(
  ranking_path, stringsAsFactors = FALSE, check.names = FALSE
)
if (nrow(rankings) != 176750L || anyDuplicated(rankings$pair_id) ||
    any(c("tissue", "approach") %in% names(rankings)) ||
    any(rankings$outcome_label_state != "closed") ||
    any(as.logical(rankings$biological_outcomes_computed))) {
  stop("Independent immutable-ranking validation failed.", call. = FALSE)
}
rankings$training_tissue <- labels$tissue[
  match(rankings$training_sample_id, labels$sample_id)
]
rankings$query_tissue <- labels$tissue[
  match(rankings$query_sample_id, labels$sample_id)
]
rankings$query_study <- labels$study[
  match(rankings$query_sample_id, labels$sample_id)
]
rankings$training_study <- labels$study[
  match(rankings$training_sample_id, labels$sample_id)
]
if (anyNA(rankings[, c("training_tissue", "query_tissue", "query_study",
                       "training_study")]) ||
    any(rankings$query_study != sub("^large_loso_v1:", "", rankings$fold_id)) ||
    any(rankings$query_study == rankings$training_study)) {
  stop("Independent label join or held-out-study isolation failed.",
       call. = FALSE)
}

key_fields <- c("fold_id", "seed", "method_id", "query_sample_id")
rank_groups <- split(
  rankings, interaction(rankings[, key_fields], drop = TRUE, lex.order = TRUE)
)
independent_query <- do.call(rbind, lapply(rank_groups, function(part) {
  part <- part[order(part$neighbor_rank), , drop = FALSE]
  if (!identical(as.integer(part$neighbor_rank), seq_len(nrow(part))) ||
      any(diff(part$distance) < 0)) {
    stop("An immutable rank sequence is invalid.", call. = FALSE)
  }
  same <- which(part$training_tissue == part$query_tissue[[1L]])
  first <- same[[1L]]
  data.frame(
    fold_id = part$fold_id[[1L]], seed = as.integer(part$seed[[1L]]),
    method_id = part$method_id[[1L]],
    query_sample_id = part$query_sample_id[[1L]],
    query_tissue = part$query_tissue[[1L]],
    held_out_study = part$query_study[[1L]],
    first_same_tissue_sample_id = part$training_sample_id[[first]],
    first_same_tissue_rank = as.integer(part$neighbor_rank[[first]]),
    reciprocal_rank = 1 / as.numeric(part$neighbor_rank[[first]]),
    nearest_sample_id = part$training_sample_id[[1L]],
    nearest_tissue = part$training_tissue[[1L]],
    one_nn_correct = part$training_tissue[[1L]] == part$query_tissue[[1L]],
    stringsAsFactors = FALSE
  )
}))
rownames(independent_query) <- NULL
query <- read_result("mv05k-query-endpoints")
endpoint_key <- function(x) paste(
  x$fold_id, x$seed, x$method_id, x$query_sample_id, sep = "\r"
)
query <- query[match(endpoint_key(independent_query), endpoint_key(query)),
               , drop = FALSE]
if (anyNA(query$query_sample_id) ||
    !identical(query$first_same_tissue_sample_id,
               independent_query$first_same_tissue_sample_id) ||
    !identical(as.integer(query$first_same_tissue_rank),
               independent_query$first_same_tissue_rank) ||
    !identical(query$nearest_sample_id, independent_query$nearest_sample_id) ||
    !identical(query$nearest_tissue, independent_query$nearest_tissue) ||
    !identical(as.logical(query$one_nn_correct),
               independent_query$one_nn_correct) ||
    any(query$endpoint_status != "estimable")) {
  stop("Independent query-endpoint reconstruction failed.", call. = FALSE)
}
query_difference <- assert_close(
  query$reciprocal_rank, independent_query$reciprocal_rank,
  "query reciprocal ranks"
)

sample_summary <- read_result("mv05k-sample-endpoint-summaries")
sample_groups <- split(
  independent_query,
  interaction(independent_query$method_id, independent_query$query_sample_id,
              drop = TRUE, lex.order = TRUE)
)
independent_sample <- do.call(rbind, lapply(sample_groups, function(part) {
  data.frame(
    method_id = part$method_id[[1L]],
    query_sample_id = part$query_sample_id[[1L]],
    query_tissue = part$query_tissue[[1L]],
    held_out_study = part$held_out_study[[1L]],
    mean_reciprocal_rank = mean(part$reciprocal_rank),
    one_nn_balanced_accuracy = mean(as.numeric(part$one_nn_correct)),
    seeds = length(unique(part$seed)), stringsAsFactors = FALSE
  )
}))
sample_key <- function(x) paste(x$method_id, x$query_sample_id, sep = "\r")
sample_summary <- sample_summary[
  match(sample_key(independent_sample), sample_key(sample_summary)),,
  drop = FALSE
]
sample_mrr_difference <- assert_close(
  sample_summary$mean_reciprocal_rank,
  independent_sample$mean_reciprocal_rank, "sample MRR summaries"
)
sample_1nn_difference <- assert_close(
  sample_summary$one_nn_balanced_accuracy,
  independent_sample$one_nn_balanced_accuracy, "sample 1-NN summaries"
)
if (any(sample_summary$seeds != 5L) ||
    any(independent_sample$seeds != 5L)) {
  stop("Seeds were not treated as five repeated technical realizations.",
       call. = FALSE)
}

production <- read_result("mv05k-production-summary")
independent_macro <- do.call(rbind, lapply(methods, function(method) {
  part <- independent_sample[independent_sample$method_id == method, ]
  data.frame(
    method_id = method,
    mrr = macro(part$mean_reciprocal_rank, part$query_tissue),
    one_nn = macro(part$one_nn_balanced_accuracy, part$query_tissue),
    stringsAsFactors = FALSE
  )
}))
production <- production[match(methods, production$method_id), , drop = FALSE]
method_mrr_difference <- assert_close(
  production$macro_mean_reciprocal_rank, independent_macro$mrr,
  "method macro MRR"
)
method_1nn_difference <- assert_close(
  production$macro_one_nn_balanced_accuracy, independent_macro$one_nn,
  "method macro 1-NN"
)
if (any(production$upstream_refits != 0L) ||
    any(production$reranking_operations != 0L) ||
    any(production$clustering_jobs_executed != 0L) ||
    any(production$integration_jobs_executed != 0L) ||
    any(production$gene_view_jobs_executed != 0L) ||
    any(production$fusion_jobs_executed != 0L) ||
    any(production$new_data_jobs_executed != 0L) ||
    any(production$sct_outcome_files_read != 0L) ||
    any(production$nonestimable_observations != 0L)) {
  stop("A forbidden or undeclared MV5-K operation is recorded.", call. = FALSE)
}

# Independently reconstruct the exact frozen block-resampling sequence.
sample_ids <- sort(unique(independent_sample$query_sample_id), method = "radix")
base <- independent_sample[independent_sample$method_id == methods[[1L]],
                           c("query_sample_id", "query_tissue",
                             "held_out_study")]
base <- base[match(sample_ids, base$query_sample_id), , drop = FALSE]
arrays <- list(
  cross_study_tissue_mrr_v1 = matrix(vapply(methods, function(method) {
    part <- independent_sample[independent_sample$method_id == method, ]
    part$mean_reciprocal_rank[match(sample_ids, part$query_sample_id)]
  }, numeric(length(sample_ids))), nrow = length(sample_ids),
  dimnames = list(sample_ids, methods)),
  cross_study_tissue_1nn_balanced_accuracy_v1 = matrix(
    vapply(methods, function(method) {
      part <- independent_sample[independent_sample$method_id == method, ]
      part$one_nn_balanced_accuracy[match(sample_ids, part$query_sample_id)]
    }, numeric(length(sample_ids))), nrow = length(sample_ids),
    dimnames = list(sample_ids, methods))
)
study_table <- unique(base[, c("query_tissue", "held_out_study")])
strata <- split(study_table, study_table$query_tissue)
strata <- strata[sort(names(strata), method = "radix")]
set.seed(20260808L)
bootstrap <- vector("list", length(arrays))
names(bootstrap) <- names(arrays)
for (endpoint_index in seq_along(arrays)) {
  values <- arrays[[endpoint_index]]
  estimates <- matrix(NA_real_, nrow = 2000L, ncol = length(methods),
                      dimnames = list(NULL, methods))
  for (replicate_id in seq_len(2000L)) {
    tissue_values <- matrix(NA_real_, nrow = length(strata),
                            ncol = length(methods))
    for (tissue_index in seq_along(strata)) {
      study_ids <- sort(strata[[tissue_index]]$held_out_study,
                        method = "radix")
      drawn <- sample(study_ids, length(study_ids), replace = TRUE)
      rows <- unlist(lapply(drawn, function(study) {
        which(base$query_tissue == names(strata)[[tissue_index]] &
                base$held_out_study == study)
      }), use.names = FALSE)
      tissue_values[tissue_index, ] <- colMeans(values[rows, , drop = FALSE])
    }
    estimates[replicate_id, ] <- colMeans(tissue_values)
  }
  bootstrap[[endpoint_index]] <- estimates
}
bootstrap_audit <- read_result("mv05k-bootstrap-audit")
observed_bootstrap_hash <- vapply(
  bootstrap, digest::digest, character(1L), algo = "sha256", serialize = TRUE
)
bootstrap_audit <- bootstrap_audit[
  match(names(bootstrap), bootstrap_audit$endpoint_id), , drop = FALSE
]
if (!identical(observed_bootstrap_hash,
               stats::setNames(bootstrap_audit$replicate_matrix_sha256,
                               bootstrap_audit$endpoint_id)) ||
    any(bootstrap_audit$study_blocks != 15L) ||
    any(bootstrap_audit$seeds_treated_as_independent)) {
  stop("Independent tissue-stratified study-block bootstrap failed.",
       call. = FALSE)
}

intervals <- read_result("mv05k-method-intervals")
interval_difference <- 0
for (endpoint in names(arrays)) {
  for (method in methods) {
    expected <- stats::quantile(
      bootstrap[[endpoint]][, method], c(0.025, 0.975),
      names = FALSE, type = 7
    )
    observed <- intervals[
      intervals$endpoint_id == endpoint & intervals$method_id == method,
      c("ci_lower", "ci_upper")
    ]
    interval_difference <- max(
      interval_difference,
      assert_close(unlist(observed), expected, "method interval")
    )
  }
}

contrasts <- read_result("mv05k-paired-contrasts")
energy <- "integrated_cell_distribution_energy_v1"
topology <- c("integrated_cell_landscape_h0_v1", "integrated_cell_landscape_h1_v1")
contrast_difference <- 0
for (endpoint in names(arrays)) {
  for (method in topology) {
    values <- arrays[[endpoint]][, method] - arrays[[endpoint]][, energy]
    expected_point <- macro(values, base$query_tissue)
    expected_interval <- stats::quantile(
      bootstrap[[endpoint]][, method] - bootstrap[[endpoint]][, energy],
      c(0.025, 0.975), names = FALSE, type = 7
    )
    observed <- contrasts[
      contrasts$endpoint_id == endpoint & contrasts$method_id == method,
      , drop = FALSE
    ]
    contrast_difference <- max(
      contrast_difference,
      assert_close(observed$estimate, expected_point, "contrast point"),
      assert_close(c(observed$ci_lower, observed$ci_upper), expected_interval,
                   "contrast interval")
    )
  }
}

study_ids <- sort(unique(base$held_out_study), method = "radix")
set.seed(20260809L)
signs <- matrix(
  sample(c(-1, 1), 9999L * length(study_ids), replace = TRUE),
  nrow = 9999L, ncol = length(study_ids),
  dimnames = list(NULL, study_ids)
)
randomization_audit <- read_result("mv05k-randomization-audit")
if (any(randomization_audit$sign_matrix_sha256 != digest::digest(
  signs, algo = "sha256", serialize = TRUE
))) {
  stop("Independent randomization sign-matrix validation failed.",
       call. = FALSE)
}
raw_p <- numeric(length(topology))
null_hash <- character(length(topology))
exceedance_count <- integer(length(topology))
for (index in seq_along(topology)) {
  method <- topology[[index]]
  difference <- arrays$cross_study_tissue_mrr_v1[, method] -
    arrays$cross_study_tissue_mrr_v1[, energy]
  null <- vapply(seq_len(9999L), function(replicate_id) {
    macro(
      difference * signs[replicate_id, match(base$held_out_study, study_ids)],
      base$query_tissue
    )
  }, numeric(1L))
  observed <- contrasts[
    contrasts$endpoint_id == "cross_study_tissue_mrr_v1" &
      contrasts$method_id == method, , drop = FALSE
  ]
  boundary <- abs(observed$estimate[[1L]])
  tolerance <- 64 * .Machine$double.eps * pmax(1, abs(null), boundary)
  exceedance_count[[index]] <- sum(abs(null) + tolerance >= boundary)
  raw_p[[index]] <- (exceedance_count[[index]] + 1) / 10000
  null_hash[[index]] <- digest::digest(
    null, algo = "sha256", serialize = TRUE
  )
}
randomization_audit <- randomization_audit[
  match(topology, randomization_audit$method_id), , drop = FALSE
]
if (!identical(null_hash, randomization_audit$null_distribution_sha256) ||
    !identical(exceedance_count,
               as.integer(randomization_audit$exceedance_count)) ||
    any(randomization_audit$boundary_policy !=
          "absolute_two_sided_64eps_ties_count_as_exceedance_v1") ||
    any(randomization_audit$epsilon_multiplier != 64)) {
  stop("Independent randomization null validation failed.", call. = FALSE)
}
primary_contrasts <- contrasts[
  contrasts$endpoint_id == "cross_study_tissue_mrr_v1",
]
primary_contrasts <- primary_contrasts[
  match(topology, primary_contrasts$method_id), , drop = FALSE
]
randomization_difference <- max(
  assert_close(primary_contrasts$raw_p_value, raw_p,
               "randomization p-values"),
  assert_close(primary_contrasts$holm_p_value,
               stats::p.adjust(raw_p, method = "holm"), "Holm p-values")
)

validation <- data.frame(
  contract_id = "mv05k_independent_validation_v1",
  validation = c(
    "prediction_lock_and_completion", "frozen_label_join_and_loso",
    "query_endpoint_formula", "sample_seed_aggregation",
    "tissue_macro_denominator", "blocked_bootstrap_reconstruction",
    "method_interval_reconstruction", "paired_contrast_reconstruction",
    "randomization_and_holm", "forbidden_operation_counters",
    "source_hashes_unchanged_after_evaluation"
  ),
  status = "passed",
  rows_checked = c(
    375L, nrow(rankings), nrow(query), nrow(sample_summary), nrow(production),
    2000L * length(arrays), nrow(intervals), nrow(contrasts), 9999L * 2L,
    nrow(production), 2L
  ),
  max_absolute_difference = c(
    0, 0, query_difference,
    max(sample_mrr_difference, sample_1nn_difference),
    max(method_mrr_difference, method_1nn_difference), 0,
    interval_difference, contrast_difference, randomization_difference, 0, 0
  ),
  evidence = c(
    "75 bundles; 375/375 method groups; frozen ranking SHA-256",
    "90 samples; 15 study blocks; five eligible tissues; no study overlap",
    "first same-tissue immutable rank and fixed rank-1 prediction",
    "five seeds averaged within sample; never independent n",
    "sample mean within tissue then equal five-tissue macro mean",
    "2000 paired tissue-stratified held-out-study resamples reproduced by hash",
    "all method MRR and 1-NN percentile intervals reproduced",
    "H0/H1 minus energy point estimates and intervals reproduced",
    "9999 paired study sign flips and two-test F1 Holm adjustment reproduced",
    "zero refit/rerank/clustering/integration/gene/fusion counters",
    "ranking and frozen metadata hashes unchanged"
  ),
  stringsAsFactors = FALSE
)

if (!identical(file_sha(ranking_path), expected_ranking_sha) ||
    !identical(file_sha(metadata_path), expected_metadata_sha) ||
    !identical(file_sha(completion_path), expected_completion_sha) ||
    !identical(file_sha(group_path), expected_group_sha) ||
    !identical(file_sha(registry_path), expected_registry_sha) ||
    !identical(file_sha(scale_path), expected_scale_sha) ||
    !identical(file_sha(assembly_path), expected_assembly_sha)) {
  stop("An immutable source changed during validation.", call. = FALSE)
}
write_provenance_csv(validation, output_path)
message("MV5-K independent retrieval evaluation validation passed.")
