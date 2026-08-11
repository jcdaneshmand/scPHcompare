#!/usr/bin/env Rscript

options(warn = 2)
for (package in c("digest", "pkgload")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(package, " is required for MV5-Q prefreeze.", call. = FALSE)
  }
}
pkgload::load_all(".", quiet = TRUE, export_all = TRUE)

args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args)) args[[1L]] else "docs/audits"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
read_public <- function(path) utils::read.csv(path, stringsAsFactors = FALSE,
                                               check.names = FALSE)
file_sha <- function(path) digest::digest(file = path, algo = "sha256",
                                          serialize = FALSE)
write_public <- function(value, path) write_provenance_csv(value, path)

expected_head <- "e6fcc7b"
head <- trimws(system2("git", c("rev-parse", "--short=7", "HEAD"), stdout = TRUE))
if (!identical(head, expected_head)) {
  stop("MV5-Q prefreeze must be created from accepted HEAD ", expected_head,
       "; observed ", head, ".", call. = FALSE)
}

paths <- c(
  mv05n_spec = "docs/specifications/MV05N_LABEL_CLOSED_CLUSTERING_RESOURCE_GATE_SPECIFICATION_V1.md",
  mv05n_code = "R/mv05n_clustering_gate.R",
  mv05q_code = "R/mv05q_clustering_artifacts.R",
  mv05q_builder = "scripts/build_mv05q_prefreeze.R",
  mv05q_runner = "scripts/run_mv05q_production.R",
  mv05p_completion = "docs/audits/mv05p-production-completion-2026-08-10.csv",
  mv05p_units = "docs/audits/mv05p-unit-manifest-2026-08-10.csv",
  mv05p_matrix_validation = "docs/audits/mv05p-matrix-validation-2026-08-10.csv",
  mv05o_group_queue = "docs/audits/mv05o-production-group-queue-2026-08-10.csv",
  sct_query = "docs/audits/mv05d5-cell-retrieval-rankings-2026-08-08.csv.gz",
  integrated_query = "docs/audits/mv05j-integrated-cell-retrieval-rankings-2026-08-09.csv.gz"
)
if (any(!file.exists(paths))) stop("MV5-Q prefreeze source is missing.", call. = FALSE)
source_freeze <- data.frame(
  contract_id = "mv05q_source_freeze_v1", source_id = names(paths),
  artifact_path = unname(paths), sha256 = vapply(paths, file_sha, character(1L)),
  bytes = as.numeric(file.info(paths)$size), accepted_head = expected_head,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
source_freeze$source_freeze_sha256 <- .mv05q_digest(
  paste(source_freeze$artifact_path, source_freeze$sha256, sep = "\r")
)

known <- c(
  mv05p_completion = "9b858e35960cb96d533b853fc60d230b2a9cee71f3b796a9209f11069043a03b",
  mv05p_units = "8e23cf63312db1dab243aaad545449a2921aafac555ba42ea1dce5da88478740",
  mv05p_matrix_validation = "979a780ed6690a57d1fabdf7fbe693edc76dd30626e0a38e519bd0c0aea45817",
  sct_query = "ee925bb073847567dffae61e13cfa688e193267e4708ba3f9edd4a52ceb0c599",
  integrated_query = "4588902bce89a04cae0c7676b4f21f81e83013a29120ca2a4b39f3ffacfb677e"
)
observed <- stats::setNames(source_freeze$sha256, source_freeze$source_id)
if (!all(observed[names(known)] == known)) {
  stop("MV5-Q accepted input hashes drifted.", call. = FALSE)
}

aliases <- mv05q_method_alias_registry_v1()
mv05q_validate_alias_registry_v1(aliases)
groups <- read_public(paths[["mv05o_group_queue"]])
matrix_validation <- read_public(paths[["mv05p_matrix_validation"]])
completion <- read_public(paths[["mv05p_completion"]])
lapply(list(groups, matrix_validation, completion), .mv05q_assert_label_closed)
if (nrow(groups) != 150L || nrow(matrix_validation) != 525L ||
    nrow(completion) != 1L || completion$production_groups != 150L ||
    completion$total_units != 4565L || completion$total_values != 1838725L ||
    !.mv05q_is_true(completion$distance_production_executed) ||
    .mv05q_is_true(completion$production_clustering_executed) ||
    any(!.mv05q_is_true(matrix_validation$complete)) ||
    any(as.integer(matrix_validation$clustering_jobs_executed) != 0L)) {
  stop("MV5-Q accepted MV5-P completion axes are invalid.", call. = FALSE)
}

group_base <- unique(groups[c("fold_id", "representation", "training_samples")])
if (nrow(group_base) != 30L ||
    any(table(groups$fold_id, groups$representation) != 5L) ||
    !identical(sort(unique(as.integer(groups$seed))), .mv05q_required_seeds)) {
  stop("MV5-Q fold/representation/seed axes drifted.", call. = FALSE)
}

queries <- list(
  mv05d5_sct_query_bundle_v1 = read_public(paths[["sct_query"]]),
  mv05j_integrated_query_bundle_v1 = read_public(paths[["integrated_query"]])
)
lapply(queries, .mv05q_assert_label_closed)
if (any(vapply(queries, nrow, integer(1L)) != 176750L) ||
    any(vapply(queries, function(x) length(unique(x$fold_id)), integer(1L)) != 15L) ||
    any(vapply(queries, function(x) length(unique(x$query_sample_id)), integer(1L)) != 90L) ||
    any(vapply(queries, function(x) !identical(sort(unique(as.integer(x$seed))),
                                              .mv05q_required_seeds), logical(1L)))) {
  stop("MV5-Q accepted query bundles have incomplete axes.", call. = FALSE)
}

alias_query_audit <- do.call(rbind, lapply(seq_len(nrow(aliases)), function(index) {
  alias <- aliases[index, , drop = FALSE]
  bundle <- queries[[alias$query_bundle_id]]
  selected <- bundle[
    bundle$representation == alias$query_representation &
      bundle$method_id == alias$query_method_id, , drop = FALSE]
  folds <- length(unique(selected$fold_id))
  seeds <- sort(unique(as.integer(selected$seed)))
  if (nrow(selected) != 35350L || folds != 15L ||
      !identical(seeds, .mv05q_required_seeds) ||
      anyDuplicated(selected[c("fold_id", "seed", "query_sample_id",
                               "training_sample_id")])) {
    stop("MV5-Q query alias lacks its complete accepted fold/seed/pair axis: ",
         alias$representation, " / ", alias$distance_id, call. = FALSE)
  }
  data.frame(
    representation = alias$representation, distance_id = alias$distance_id,
    query_bundle_id = alias$query_bundle_id,
    query_representation = alias$query_representation,
    query_method_id = alias$query_method_id, query_rows = nrow(selected),
    folds = folds, seeds = length(seeds), complete = TRUE,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
}))
.mv05q_assert_label_closed(alias_query_audit, "alias query audit")

sct_pb <- queries[[1L]][queries[[1L]]$method_id ==
  "pseudobulk_shared_panel_euclidean_v1", ]
int_pb <- queries[[2L]][queries[[2L]]$method_id ==
  "pseudobulk_training_standardized_panel_v1", ]
pb_key <- c("fold_id", "seed", "query_sample_id", "training_sample_id")
sct_pb <- sct_pb[do.call(order, c(sct_pb[pb_key], list(method = "radix"))), ]
int_pb <- int_pb[do.call(order, c(int_pb[pb_key], list(method = "radix"))), ]
if (nrow(sct_pb) != 35350L || nrow(int_pb) != 35350L ||
    !identical(sct_pb[pb_key], int_pb[pb_key]) ||
    max(abs(sct_pb$distance - int_pb$distance)) != 0) {
  stop("MV5-Q shared pseudobulk query identity is not exact.", call. = FALSE)
}

queue <- merge(group_base, aliases, by = "representation", sort = FALSE)
queue <- queue[order(queue$fold_id, queue$representation, queue$distance_id,
                     method = "radix"), ]
queue$contract_id <- "mv05q_analysis_queue_v1"
queue$execution_order <- seq_len(nrow(queue))
queue$analysis_group_id <- vapply(seq_len(nrow(queue)), function(index) {
  paste0("mv05q_analysis_v1:", .mv05q_digest(list(
    contract_id = "mv05q_analysis_queue_identity_v1",
    fold_id = queue$fold_id[[index]],
    representation = queue$representation[[index]],
    distance_id = queue$distance_id[[index]],
    seeds = .mv05q_required_seeds,
    training_samples = queue$training_samples[[index]],
    source_freeze_sha256 = source_freeze$source_freeze_sha256[[1L]]
  )))
}, character(1L))
queue$candidate_k_min <- 2L
queue$candidate_k_max <- pmin(10L, queue$training_samples - 1L)
queue$candidate_k_count <- queue$candidate_k_max - 1L
queue$seed_count <- 5L
queue$pam_candidate_fit_count <- 45L
queue$selected_partition_fit_count <- 10L
queue$heldout_assignment_algorithms <- 2L
queue$source_freeze_sha256 <- source_freeze$source_freeze_sha256[[1L]]
queue$outcome_label_state <- "closed"
queue$biological_outcomes_computed <- FALSE
queue$production_executed <- FALSE
queue$labels_opened <- FALSE
queue <- queue[c("contract_id", "analysis_group_id", "execution_order",
                 "fold_id", "representation", "distance_id", "training_samples",
                 "candidate_k_min", "candidate_k_max", "candidate_k_count",
                 "seed_count", "pam_candidate_fit_count",
                 "selected_partition_fit_count", "heldout_assignment_algorithms",
                 "training_representation", "training_method_id",
                 "training_component", "training_source_policy", "query_bundle_id",
                 "query_representation", "query_method_id",
                 "source_freeze_sha256", "outcome_label_state",
                 "biological_outcomes_computed", "production_executed",
                 "labels_opened")]
if (nrow(queue) != 150L || anyDuplicated(queue$analysis_group_id) ||
    any(queue$candidate_k_count != 9L) ||
    sum(queue$pam_candidate_fit_count) != 6750L ||
    any(queue$production_executed) || any(queue$labels_opened)) {
  stop("MV5-Q immutable analysis queue is invalid.", call. = FALSE)
}

validation <- data.frame(
  contract_id = "mv05q_validation_plan_v1",
  validation_id = c(
    "source_hashes_and_head", "complete_training_matrix_axes",
    "complete_query_training_axes", "shared_pseudobulk_exact_identity",
    "candidate_grid_and_five_seed_stability", "canonical_cluster_signatures",
    "pam_exact_tie_rule", "average_exact_tie_rule",
    "maximum_fold_clean_repeat_all_methods", "immutable_resume_all_groups",
    "resource_caps", "public_label_safety"
  ),
  validation_type = c(
    "sha256_and_starting_head_v1", "all_750_matrices_complete_symmetric_v1",
    "all_750_query_slices_complete_v1", "distance_keys_and_values_identical_v1",
    "all_150_groups_45_pam_fits_one_se_v1", "sorted_member_signature_v1",
    "lexicographic_medoid_id_v1", "canonical_cluster_signature_v1",
    "byte_identical_selected_public_outputs_v1",
    "hash_size_timestamp_unchanged_zero_rebuild_v1",
    "two_workers_900_seconds_4gib_per_group_v1",
    "closed_labels_no_outcome_columns_or_metrics_v1"
  ),
  required = TRUE, tolerance = c(rep(0, 8L), 0, 0, 0, 0),
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
abort <- data.frame(
  contract_id = "mv05q_abort_rules_v1",
  rule_id = sprintf("MV5Q-ABORT-%02d", 1:10),
  trigger = c(
    "accepted_source_or_implementation_hash_mismatch",
    "fold_seed_representation_distance_or_sample_axis_mismatch",
    "missing_duplicate_nonfinite_negative_or_asymmetric_distance",
    "outcome_label_column_label_open_or_biological_metric_detected",
    "candidate_k_or_five_seed_partition_grid_incomplete",
    "noncanonical_cluster_or_medoid_identity",
    "heldout_sample_influences_training_fit_selection_or_naming",
    "group_elapsed_seconds_above_900_or_process_tree_rss_above_4gib",
    "partial_stale_or_hash_invalid_output_status_pair",
    "repeat_resume_tie_or_public_safety_validation_failure"
  ),
  disposition = c(rep("abort_before_new_group_preserve_completed_immutable_artifacts", 8L),
                  "quarantine_partial_group_abort_and_review",
                  "abort_stage_revoke_mv05q_production_authorization"),
  automatic_retry = FALSE, outcome_label_state = "closed",
  biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
)
lapply(list(source_freeze, aliases, queue, validation, abort),
       .mv05q_assert_label_closed)

write_public(source_freeze, file.path(output_dir,
  "mv05q-source-freeze-2026-08-10.csv"))
write_public(aliases, file.path(output_dir,
  "mv05q-method-alias-registry-2026-08-10.csv"))
write_public(alias_query_audit, file.path(output_dir,
  "mv05q-query-alias-audit-2026-08-10.csv"))
write_public(queue, file.path(output_dir,
  "mv05q-analysis-queue-2026-08-10.csv"))
write_public(validation, file.path(output_dir,
  "mv05q-validation-plan-2026-08-10.csv"))
write_public(abort, file.path(output_dir,
  "mv05q-abort-rules-2026-08-10.csv"))
message("MV5-Q prefreeze passed: sources=", nrow(source_freeze),
        " query_aliases=", nrow(alias_query_audit),
        " analyses=", nrow(queue), " candidate_pam_fits=",
        sum(queue$pam_candidate_fit_count), " query_rows=", sum(vapply(queries, nrow,
        integer(1L))), " labels=closed production_executed=0")
