#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop(
    "usage: validate_mv05g_coordinate_production.R MANIFEST RESOURCE_CSV ",
    "RESULT_ROOT DETAIL_OUTPUT SUMMARY_OUTPUT", call. = FALSE
  )
}
manifest_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
resource_path <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
result_root <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
detail_path <- args[[4L]]
summary_path <- args[[5L]]
dir.create(dirname(detail_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(summary_path), recursive = TRUE, showWarnings = FALSE)

source("R/provenance_utils.R")
source("R/mv05_benchmark_execution.R")
source("R/mv05_inductive_mapping.R")
source("R/mv05f_integration_gate.R")
source("R/mv05g_coordinate_production.R")
manifest <- utils::read.csv(
  manifest_path, stringsAsFactors = FALSE, check.names = FALSE
)
metrics <- utils::read.csv(
  resource_path, stringsAsFactors = FALSE, check.names = FALSE
)
mv05g_validate_full_manifest_v1(manifest)
mv05g_validate_group_metrics_v1(metrics)

checks <- vector("list", nrow(manifest))
for (index in seq_len(nrow(manifest))) {
  job <- manifest[index, , drop = FALSE]
  metric <- metrics[metrics$group_id == job$group_id, , drop = FALSE]
  stem <- gsub("[^A-Za-z0-9_.-]", "_", job$group_id)
  path <- file.path(result_root, stem, paste0(stem, ".rds"))
  if (nrow(metric) != 1L || !file.exists(path)) {
    stop("MV5-G group resource or private result is missing: ", job$group_id,
         call. = FALSE)
  }
  record <- readRDS(path)
  identity <- record$identity
  payload <- record$payload
  identity_without_key <- identity
  identity_without_key$cache_key <- NULL
  payload_hash_pass <- identical(
    record$payload_sha256,
    digest::digest(payload, algo = "sha256", serialize = TRUE)
  )
  cache_identity_pass <- identical(
    record$cache_key,
    paste0(
      "mv05f_mapping_group_v1:",
      digest::digest(identity_without_key, algo = "sha256", serialize = TRUE)
    )
  )
  inherited_axes_pass <- identity$fold_id == job$fold_id &&
    identity$fit_scope_id == job$fit_scope_id &&
    identity$held_out_study == job$held_out_study &&
    identity$seed == as.integer(job$seed) &&
    length(identity$training_ids) == job$training_samples &&
    length(identity$query_ids) == job$held_out_samples &&
    !length(intersect(identity$training_ids, identity$query_ids)) &&
    length(c(identity$training_ids, identity$query_ids)) == 90L &&
    nrow(identity$panel) == 500L &&
    identical(identity$dimensions, 1:30) &&
    identity$cells_per_sample == 384L &&
    identity$d1_fold_cache_key == job$fold_cache_key &&
    identical(identity$outcome_label_state, "closed") &&
    identical(identity$biological_outcomes_computed, FALSE)
  ids <- sort(c(identity$training_ids, identity$query_ids), method = "radix")
  coordinates <- payload$coordinates
  coordinate_axis_pass <- is.list(coordinates) &&
    identical(sort(names(coordinates), method = "radix"), ids) &&
    all(vapply(names(coordinates), function(sample_id) {
      value <- coordinates[[sample_id]]
      is.matrix(value) && is.numeric(value) && nrow(value) == 384L &&
        ncol(value) == 30L && !anyNA(value) && all(is.finite(value)) &&
        identical(colnames(value), paste0("PC", 1:30)) &&
        length(unique(rownames(value))) == 384L &&
        all(startsWith(rownames(value), paste0(sample_id, "__")))
    }, logical(1L))) &&
    identical(
      payload$coordinate_set_sha256,
      digest::digest(coordinates, algo = "sha256", serialize = TRUE)
    )
  mappings <- payload$query_mappings
  active <- payload$active_features
  mapping_axis_pass <- is.list(mappings) && is.list(active) &&
    identical(sort(names(mappings), method = "radix"),
              sort(identity$query_ids, method = "radix")) &&
    identical(sort(names(active), method = "radix"),
              sort(identity$query_ids, method = "radix")) &&
    all(vapply(names(mappings), function(sample_id) {
      value <- mappings[[sample_id]]
      map_identity <- list(
        contract_id = value$contract_id, engine_id = value$engine_id,
        fold_id = value$fold_id, fit_scope_id = value$fit_scope_id,
        reference_sample_ids = value$reference_sample_ids,
        held_out_sample_id = value$held_out_sample_id,
        features = value$features, dimensions = value$dimensions,
        seed = value$seed, outcome_label_state = value$outcome_label_state,
        biological_outcomes_computed = value$biological_outcomes_computed,
        reference_identity_sha256 = value$reference_identity_sha256,
        anchor_count = value$anchor_count,
        query_embedding_sha256 = value$query_embedding_sha256
      )
      identical(value$cache_key, paste0(
        "mv05_inductive_mapping_v1:",
        digest::digest(map_identity, algo = "sha256", serialize = TRUE)
      )) && value$held_out_sample_id == sample_id &&
        !sample_id %in% value$reference_sample_ids &&
        identical(value$reference_sample_ids,
                  sort(identity$training_ids, method = "radix")) &&
        value$fold_id == identity$fold_id &&
        value$fit_scope_id == identity$fit_scope_id &&
        value$seed == identity$seed &&
        identical(value$dimensions, identity$dimensions) &&
        identical(value$features, active[[sample_id]]) &&
        length(active[[sample_id]]) >= 31L &&
        all(active[[sample_id]] %in% identity$panel$gene) &&
        !anyDuplicated(active[[sample_id]]) &&
        identical(value$query_embeddings, coordinates[[sample_id]]) &&
        identical(value$query_embedding_sha256, digest::digest(
          value$query_embeddings, algo = "sha256", serialize = TRUE
        )) && value$anchor_count >= 1L
    }, logical(1L)))
  reference_pass <- identical(
    payload$reference_identity_sha256_before,
    payload$reference_identity_sha256_after
  ) && all(vapply(mappings, function(value) {
    grepl("^[0-9a-f]{64}$", value$reference_identity_sha256)
  }, logical(1L)))
  forbidden <- c(
    "diagrams", "landscapes", "distances", "rankings", "retrieval",
    "clustering", "gene_views", "fusion", "outcomes", "labels", "tissue"
  )
  scope_pass <- !any(forbidden %in% names(payload)) &&
    identical(payload$downstream_execution,
      list(ph_jobs = 0L, landscape_jobs = 0L, distance_jobs = 0L,
           clustering_jobs = 0L, gene_view_jobs = 0L, fusion_jobs = 0L,
           biological_outcome_jobs = 0L)) &&
    metric$label_transfer_jobs_executed == 0L &&
    metric$ph_jobs_executed == 0L &&
    metric$landscape_jobs_executed == 0L &&
    metric$distance_jobs_executed == 0L &&
    metric$retrieval_jobs_executed == 0L &&
    metric$clustering_jobs_executed == 0L &&
    metric$gene_view_jobs_executed == 0L &&
    metric$fusion_jobs_executed == 0L &&
    metric$new_data_jobs_executed == 0L &&
    identical(record$outcome_label_state, "closed") &&
    identical(record$biological_outcomes_computed, FALSE)
  file_hash <- digest::digest(file = path, algo = "sha256", serialize = FALSE)
  file_hash_pass <- identical(file_hash, metric$private_result_sha256)
  completion_pass <- metric$disposition == "completed" &&
    metric$exit_status == 0L &&
    metric$completed_coordinate_views == 90L &&
    metric$completed_query_mappings == length(identity$query_ids) &&
    metric$elapsed_seconds <= metric$group_timeout_seconds &&
    metric$peak_process_tree_rss_bytes <= metric$rss_cap_bytes
  checks[[index]] <- data.frame(
    contract_id = "mv05g_independent_group_validation_v1",
    group_id = job$group_id, group_order = job$group_order,
    held_out_study = job$held_out_study, seed = job$seed,
    payload_hash_pass = payload_hash_pass,
    cache_identity_pass = cache_identity_pass,
    inherited_axes_pass = inherited_axes_pass,
    coordinate_axis_pass = coordinate_axis_pass,
    held_out_mapping_pass = mapping_axis_pass,
    reference_immutability_pass = reference_pass,
    scope_pass = scope_pass, file_hash_pass = file_hash_pass,
    completion_resource_pass = completion_pass,
    coordinate_views = length(coordinates), query_mappings = length(mappings),
    minimum_active_features = min(vapply(active, length, integer(1L))),
    maximum_dropped_features = max(500L - vapply(active, length, integer(1L))),
    minimum_anchor_count = min(vapply(mappings, `[[`, integer(1L),
                                       "anchor_count")),
    maximum_anchor_count = max(vapply(mappings, `[[`, integer(1L),
                                       "anchor_count")),
    result_file_sha256 = file_hash, failure_code = "",
    biological_outcomes_computed = FALSE, outcome_label_state = "closed",
    stringsAsFactors = FALSE
  )
}
detail <- do.call(rbind, checks)
detail <- detail[order(detail$group_order), , drop = FALSE]
pass_columns <- grep("_pass$", names(detail), value = TRUE)
all_pass <- all(as.matrix(detail[pass_columns]))
summary <- data.frame(
  contract_id = "mv05g_independent_validation_summary_v1",
  groups_validated = nrow(detail), studies = length(unique(detail$held_out_study)),
  seeds = length(unique(detail$seed)),
  coordinate_views = sum(detail$coordinate_views),
  query_mappings = sum(detail$query_mappings),
  validation_categories = length(pass_columns),
  validation_checks = nrow(detail) * length(pass_columns),
  all_validation_checks_pass = all_pass,
  minimum_active_features = min(detail$minimum_active_features),
  maximum_dropped_features = max(detail$maximum_dropped_features),
  minimum_anchor_count = min(detail$minimum_anchor_count),
  maximum_anchor_count = max(detail$maximum_anchor_count),
  failure_codes = sum(nzchar(detail$failure_code)),
  label_transfer_jobs_executed = 0L, ph_jobs_executed = 0L,
  landscape_jobs_executed = 0L, distance_jobs_executed = 0L,
  retrieval_jobs_executed = 0L, clustering_jobs_executed = 0L,
  gene_view_jobs_executed = 0L, fusion_jobs_executed = 0L,
  new_data_jobs_executed = 0L, biological_outcomes_computed = FALSE,
  outcome_label_state = "closed", stringsAsFactors = FALSE
)
if (!all_pass || nrow(detail) != 75L || sum(detail$coordinate_views) != 6750L ||
    sum(detail$query_mappings) != 450L) {
  stop("Independent MV5-G validation failed.", call. = FALSE)
}
write_provenance_csv(detail, detail_path)
write_provenance_csv(summary, summary_path)
message("Independently validated all 75 MV5-G coordinate groups.")
