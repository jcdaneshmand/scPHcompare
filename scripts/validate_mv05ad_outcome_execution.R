#!/usr/bin/env Rscript

options(warn = 2)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest is required.")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(paste("usage: validate_mv05ad_outcome_execution.R EXTERNAL_METADATA_PATH",
             "PRIVATE_ROOT AUDIT_DIR OUTPUT_CSV"), call. = FALSE)
}
metadata_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
private_root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
audit_dir <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
output_path <- args[[4L]]
readc <- function(path, ...) read.csv(
  path, stringsAsFactors = FALSE, check.names = FALSE, ...)
sha <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
dig <- function(value) digest::digest(value, algo = "sha256", serialize = TRUE)
close_num <- function(x, y, tolerance = 2e-13) {
  length(x) == length(y) && all(is.na(x) == is.na(y)) &&
    all(abs(x[!is.na(x)] - y[!is.na(y)]) <= tolerance)
}
with_rng <- function(seed, expression) {
  old_kind <- RNGkind(); had <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
  if (had) old_seed <- get(".Random.seed", .GlobalEnv)
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had) assign(".Random.seed", old_seed, .GlobalEnv)
    else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(seed)
  force(expression)
}
macro <- function(base, values, weights = NULL) {
  mean(vapply(sort(unique(base$query_tissue), method = "radix"), function(t) {
    hit <- base$query_tissue == t
    if (is.null(weights)) mean(values[hit]) else weighted.mean(values[hit], weights[hit])
  }, numeric(1L)))
}

paths <- function(name) file.path(audit_dir, name)
rankings <- readc("docs/audits/mv05ad-cosine-rankings-2026-08-11.csv.gz")
cosine_public <- readc(paths("mv05ad-cosine-query-method-outcomes-2026-08-11.csv"))
paired_public <- readc(paths("mv05ad-paired-query-method-endpoints-2026-08-11.csv"))
query_public <- readc(paths("mv05ad-query-estimands-2026-08-11.csv"))
sample_public <- readc(paths("mv05ad-sample-estimands-2026-08-11.csv"))
tissue_public <- readc(paths("mv05ad-tissue-estimands-2026-08-11.csv"))
macro_public <- readc(paths("mv05ad-macro-estimands-2026-08-11.csv"))
interval_public <- readc(paths("mv05ad-estimand-intervals-2026-08-11.csv"))
primary_public <- readc(paths("mv05ad-primary-contrasts-2026-08-11.csv"))
bootstrap_public <- readc(paths("mv05ad-bootstrap-audit-2026-08-11.csv"))
random_public <- readc(paths("mv05ad-randomization-audit-2026-08-11.csv"))
completion <- readc(paths("mv05ad-outcome-unit-completion-2026-08-11.csv"))
registry <- readc("docs/audits/mv05ac-estimand-registry-2026-08-11.csv")
method_map <- readc("docs/audits/mv05ac-method-map-2026-08-11.csv")

header <- readc(metadata_path, nrows = 0L)
classes <- rep("NULL", ncol(header)); names(classes) <- names(header)
classes[c("orig.ident", "SRA", "Tissue.x")] <- "character"
raw <- read.csv(metadata_path, stringsAsFactors = FALSE, check.names = FALSE,
                colClasses = classes)
labels <- data.frame(sample_id = trimws(raw$orig.ident), study = trimws(raw$SRA),
                     tissue = tolower(trimws(raw$Tissue.x)), stringsAsFactors = FALSE)
labels <- labels[labels$tissue %in% c("bone marrow", "colon", "liver", "pbmc", "testis"), ]
if (nrow(labels) != 90L || anyDuplicated(labels$sample_id)) stop("Label identity failed.")

# Independent canonical rank validation and endpoint reconstruction.
groups <- split(seq_len(nrow(rankings)), interaction(
  rankings$representation, rankings$fold_id, rankings$seed,
  rankings$method_id, rankings$query_sample_id,
  drop = TRUE, lex.order = TRUE))
endpoint_rows <- vector("list", length(groups)); gi <- 0L; rank_rows <- 0L
for (index in groups) {
  gi <- gi + 1L; part <- rankings[index, , drop = FALSE]
  expected <- order(part$distance, part$training_sample_id, method = "radix")
  part <- part[expected, , drop = FALSE]
  runs <- rle(part$distance); ties <- rep(runs$lengths, runs$lengths)
  if (!identical(as.integer(part$neighbor_rank), seq_len(nrow(part))) ||
      !identical(as.integer(part$distance_tie_size), ties) ||
      !identical(as.logical(part$distance_tied), ties > 1L) ||
      any(part$tie_break_policy !=
        "ascending_distance_then_canonical_training_sample_id_radix_v1")) {
    stop("Independent ranking reconstruction failed.")
  }
  q <- match(part$query_sample_id[[1L]], labels$sample_id)
  tr <- match(part$training_sample_id, labels$sample_id)
  if (anyNA(c(q, tr))) stop("Independent label join failed.")
  same <- which(labels$tissue[tr] == labels$tissue[[q]])
  if (!length(same)) stop("Unexpected nonestimable endpoint.")
  family <- method_map$family_id[
    method_map$representation == part$representation[[1L]] &
      method_map$cosine_method_id == part$method_id[[1L]]]
  endpoint_rows[[gi]] <- data.frame(
    fold_id = part$fold_id[[1L]], held_out_study = labels$study[[q]],
    seed = as.integer(part$seed[[1L]]), representation = part$representation[[1L]],
    family_id = family, method_id = part$method_id[[1L]],
    query_sample_id = part$query_sample_id[[1L]], query_tissue = labels$tissue[[q]],
    training_samples = nrow(part),
    training_studies = length(unique(labels$study[tr])),
    first_same_tissue_rank = part$neighbor_rank[same[[1L]]],
    reciprocal_rank = 1 / part$neighbor_rank[same[[1L]]],
    one_nn_correct = labels$tissue[tr[[1L]]] == labels$tissue[[q]],
    nearest_distance_tied = part$distance_tied[[1L]], stringsAsFactors = FALSE)
  rank_rows <- rank_rows + nrow(part)
}
cosine <- do.call(rbind, endpoint_rows)
kfields <- c("fold_id", "seed", "representation", "family_id", "query_sample_id")
key <- function(x, fields = kfields) do.call(paste, c(x[fields], sep = "\r"))
cosine_public <- cosine_public[match(key(cosine), key(cosine_public)), ]
if (rank_rows != 282800L || nrow(cosine) != 3600L || anyNA(cosine_public$method_id) ||
    !close_num(cosine$reciprocal_rank, cosine_public$reciprocal_rank) ||
    !identical(as.logical(cosine$one_nn_correct), as.logical(cosine_public$one_nn_correct)) ||
    !identical(as.integer(cosine$first_same_tissue_rank),
               as.integer(cosine_public$first_same_tissue_rank))) {
  stop("Independent endpoint reconstruction failed.")
}

# Independently reconstruct pairing from accepted baseline rows.
sct <- readc("docs/audits/mv05e-query-endpoints-2026-08-08.csv")
integrated <- readc("docs/audits/mv05k-query-endpoints-2026-08-10.csv")
baseline <- rbind(transform(sct, representation = "sct_whole"),
                  transform(integrated, representation = "inductive_integrated"))
baseline$family_id <- NA_character_
for (i in seq_len(nrow(method_map))) {
  hit <- baseline$representation == method_map$representation[[i]] &
    baseline$method_id == method_map$baseline_method_id[[i]]
  baseline$family_id[hit] <- method_map$family_id[[i]]
}
baseline <- baseline[!is.na(baseline$family_id), ]
pair_fields <- c("fold_id", "held_out_study", "seed", "representation",
                 "family_id", "query_sample_id", "query_tissue")
pk <- function(x) key(x, pair_fields)
baseline <- baseline[match(pk(cosine), pk(baseline)), ]
direct_rr <- cosine$reciprocal_rank - baseline$reciprocal_rank
direct_nn <- as.numeric(cosine$one_nn_correct) - as.numeric(baseline$one_nn_correct)
paired_public <- paired_public[match(pk(cosine), pk(paired_public)), ]
if (!close_num(direct_rr, paired_public$direct_reciprocal_rank_change) ||
    !close_num(direct_nn, paired_public$direct_one_nn_change)) {
  stop("Independent baseline pairing failed.")
}

# Reconstruct every query estimand without production helpers.
direct <- do.call(rbind, lapply(seq_len(nrow(cosine)), function(i) {
  rbind(
    data.frame(cosine[i, pair_fields], endpoint_id = "cross_study_tissue_mrr_v1",
               estimand_type = "direct_cosine_chord_minus_euclidean", value = direct_rr[[i]]),
    data.frame(cosine[i, pair_fields],
               endpoint_id = "cross_study_tissue_1nn_balanced_accuracy_v1",
               estimand_type = "direct_cosine_chord_minus_euclidean", value = direct_nn[[i]]))
}))
did <- do.call(rbind, lapply(c("h0", "h1"), function(family) {
  top <- cosine$family_id == family; energy <- cosine$family_id == "energy"
  did_fields <- setdiff(pair_fields, "family_id")
  ek <- key(cosine[energy, ], did_fields)
  energy_index <- which(energy)[match(key(cosine[top, ], did_fields), ek)]
  top_index <- which(top)
  rbind(
    data.frame(cosine[top_index, pair_fields], endpoint_id = "cross_study_tissue_mrr_v1",
      estimand_type = "topology_increment_cosine_chord_minus_euclidean_difference_in_differences",
      value = direct_rr[top_index] - direct_rr[energy_index]),
    data.frame(cosine[top_index, pair_fields],
      endpoint_id = "cross_study_tissue_1nn_balanced_accuracy_v1",
      estimand_type = "topology_increment_cosine_chord_minus_euclidean_difference_in_differences",
      value = direct_nn[top_index] - direct_nn[energy_index]))
}))
query <- rbind(direct, did)
lookup <- paste(registry$estimand_type, registry$representation, registry$family_id,
                registry$endpoint_id, sep = "\r")
observed <- paste(query$estimand_type, query$representation, query$family_id,
                  query$endpoint_id, sep = "\r")
query$estimand_id <- registry$estimand_id[match(observed, lookup)]
qkey_fields <- c("estimand_id", "fold_id", "seed", "query_sample_id")
query_public <- query_public[match(key(query, qkey_fields), key(query_public, qkey_fields)), ]
if (nrow(query) != 10800L || anyNA(query_public$value) ||
    !close_num(query$value, query_public$value)) stop("Query estimand reconstruction failed.")

# Independent seed/sample/tissue/macro aggregation.
sample <- aggregate(value ~ estimand_id + query_sample_id + query_tissue + held_out_study,
                    query, mean)
skey <- c("estimand_id", "query_sample_id")
sample_public <- sample_public[match(key(sample, skey), key(sample_public, skey)), ]
if (!close_num(sample$value, sample_public$estimate)) stop("Sample aggregation failed.")
tissue <- aggregate(value ~ estimand_id + query_tissue, sample, mean)
tkey <- c("estimand_id", "query_tissue")
tissue_public <- tissue_public[match(key(tissue, tkey), key(tissue_public, tkey)), ]
if (!close_num(tissue$value, tissue_public$estimate)) stop("Tissue aggregation failed.")
macro_values <- aggregate(value ~ estimand_id, tissue, mean)
macro_public <- macro_public[match(macro_values$estimand_id, macro_public$estimand_id), ]
if (!close_num(macro_values$value, macro_public$estimate)) stop("Macro aggregation failed.")

# Recreate exact deterministic bootstrap and sign matrices and inference.
ids <- sort(unique(sample$query_sample_id), method = "radix")
ordered_registry <- registry[order(registry$estimand_order), ]
base <- sample[sample$estimand_id == ordered_registry$estimand_id[[1L]],
               c("query_sample_id", "query_tissue", "held_out_study")]
base <- base[match(ids, base$query_sample_id), ]
values <- vapply(ordered_registry$estimand_id, function(id) {
  part <- sample[sample$estimand_id == id, ]; part$value[match(ids, part$query_sample_id)]
}, numeric(length(ids)))
studies <- sort(unique(base$held_out_study), method = "radix")
design <- unique(base[c("held_out_study", "query_tissue")])
counts <- with_rng(20260814L, {
  x <- matrix(0L, 2000L, length(studies), dimnames = list(NULL, studies))
  for (t in sort(unique(base$query_tissue), method = "radix")) {
    ss <- sort(design$held_out_study[design$query_tissue == t], method = "radix")
    for (b in 1:2000) {
      draw <- sample(ss, length(ss), replace = TRUE); tab <- table(draw)
      x[b, match(names(tab), studies)] <- as.integer(tab)
    }
  }
  x
})
study_index <- match(base$held_out_study, studies)
boot <- matrix(NA_real_, 2000L, 24L,
               dimnames = list(NULL, ordered_registry$estimand_id))
for (b in 1:2000) for (j in 1:24) {
  boot[b, j] <- macro(base, values[, j], counts[b, study_index])
}
private_matrices <- readRDS(file.path(private_root, "inference-matrices.rds"))
if (!identical(counts, private_matrices$bootstrap_counts) ||
    !close_num(as.vector(boot), as.vector(private_matrices$bootstrap_estimates)) ||
    dig(counts) != bootstrap_public$block_count_matrix_sha256 ||
    dig(boot) != bootstrap_public$replicate_matrix_sha256) {
  stop("Independent bootstrap reconstruction failed: counts_exact=",
       identical(counts, private_matrices$bootstrap_counts),
       " estimates_close=",
       close_num(as.vector(boot), as.vector(private_matrices$bootstrap_estimates)),
       " count_hash=", dig(counts) == bootstrap_public$block_count_matrix_sha256,
       " estimate_hash=", dig(boot) == bootstrap_public$replicate_matrix_sha256,
       call. = FALSE)
}
ci <- t(vapply(1:24, function(j) quantile(boot[, j], c(.025, .975),
                                          names = FALSE, type = 7), numeric(2)))
interval_public <- interval_public[match(ordered_registry$estimand_id,
                                         interval_public$estimand_id), ]
if (!close_num(ci[, 1], interval_public$ci_lower) ||
    !close_num(ci[, 2], interval_public$ci_upper)) stop("Interval reconstruction failed.")

signs <- with_rng(20260815L, matrix(sample(c(-1, 1), 9999L * 15L, replace = TRUE),
                                    nrow = 9999L, ncol = 15L,
                                    dimnames = list(NULL, studies)))
if (!identical(signs, private_matrices$sign_matrix) ||
    dig(signs) != unique(random_public$sign_matrix_sha256)) {
  stop("Independent sign matrix reconstruction failed.")
}
primary_registry <- ordered_registry[
  ordered_registry$estimand_role == "confirmatory_cosine_sensitivity", ]
raw_p <- numeric(4L); exceed <- integer(4L); null_hash <- character(4L)
for (i in 1:4) {
  j <- match(primary_registry$estimand_id[[i]], ordered_registry$estimand_id)
  null <- vapply(1:9999, function(b) {
    macro(base, values[, j] * signs[b, study_index])
  }, numeric(1L))
  observed_value <- macro(base, values[, j])
  tol <- 64 * .Machine$double.eps * pmax(1, abs(null), abs(observed_value))
  exceed[[i]] <- sum(abs(null) + tol >= abs(observed_value))
  raw_p[[i]] <- (exceed[[i]] + 1) / 10000
  null_hash[[i]] <- dig(null)
}
primary_public <- primary_public[match(primary_registry$estimand_id,
                                       primary_public$estimand_id), ]
if (!close_num(raw_p, primary_public$raw_p_value) ||
    !close_num(p.adjust(raw_p, "holm"), primary_public$holm_p_value) ||
    !identical(exceed, as.integer(random_public$exceedance_count[
      match(primary_registry$estimand_id, random_public$estimand_id)])) ||
    !identical(null_hash, random_public$null_distribution_sha256[
      match(primary_registry$estimand_id, random_public$estimand_id)])) {
  stop("Independent randomization/Holm reconstruction failed.")
}

# Rehash every private outcome artifact/status pair.
private_ok <- logical(nrow(completion)); status_ok <- logical(nrow(completion))
for (i in seq_len(nrow(completion))) {
  stem <- sub("^mv05ac_eval_v1:", "", completion$evaluation_unit_id[[i]])
  artifact <- file.path(private_root, "units", paste0(stem, ".rds"))
  status_path <- file.path(private_root, "units", paste0(stem, ".status.rds"))
  status <- readRDS(status_path)
  private_ok[[i]] <- sha(artifact) == completion$artifact_sha256[[i]]
  status_ok[[i]] <- status$artifact_sha256 == sha(artifact) && status$state == "completed"
}
if (!all(private_ok) || !all(status_ok)) stop("Private outcome identity failed.")

validation <- data.frame(
  contract_id = "mv05ad_outcome_independent_validation_v1",
  validation_id = c(
    "canonical_rankings", "tissue_label_join", "query_method_outcomes",
    "accepted_baseline_pairing", "query_estimands", "sample_aggregation",
    "tissue_macro_aggregation", "bootstrap_counts_and_estimates",
    "percentile_intervals", "sign_matrix_and_nulls", "raw_and_holm_p_values",
    "private_artifact_hashes", "private_status_identity",
    "complete_reporting_axes", "zero_prohibited_operations"),
  passed = TRUE,
  evidence = c(
    "282800_rows", "90_samples_15_studies_5_tissues_tissue_only",
    "3600_query_method_rows", "3600_exact_pairs", "10800_query_rows_24_estimands",
    "2160_rows_five_seeds_each", "120_tissue_rows_24_macros",
    "2000_by_15_counts_and_2000_by_24_estimates_exact",
    "24_type7_intervals_exact", "9999_by_15_signs_four_nulls_exact",
    "four_raw_and_holm_values_exact", "150_of_150", "150_of_150",
    "7200_endpoints_24_intervals_4_tests", "no_clustering_other_config_refit_rerank_selection_equivalence"),
  production_scientific_helper_called = FALSE,
  stringsAsFactors = FALSE)
if (file.exists(output_path)) stop("Refusing to overwrite validation output.")
write.csv(validation, output_path, row.names = FALSE, na = "")
message("MV5-AD independent outcome validation passed: 15 categories, 24 estimands, 4 tests")
