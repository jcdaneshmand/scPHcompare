#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop(
    "usage: validate_mv05f_mapping_pilot.R MANIFEST RESOURCE_CSV RESULT_ROOT ",
    "DETAIL_OUTPUT SUMMARY_OUTPUT", call. = FALSE
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
manifest <- utils::read.csv(
  manifest_path, stringsAsFactors = FALSE, check.names = FALSE
)
metrics <- utils::read.csv(
  resource_path, stringsAsFactors = FALSE, check.names = FALSE
)
mv05f_validate_pilot_manifest_v1(manifest)
mv05f_validate_resource_metrics_v1(metrics)

checks <- vector("list", nrow(manifest))
for (index in seq_len(nrow(manifest))) {
  job <- manifest[index, , drop = FALSE]
  metric <- metrics[metrics$group_id == job$group_id, , drop = FALSE]
  stem <- gsub("[^A-Za-z0-9_.-]", "_", job$group_id)
  path <- file.path(result_root, stem, paste0(stem, ".rds"))
  record <- readRDS(path)
  payload_hash_pass <- identical(
    record$payload_sha256,
    digest::digest(record$payload, algo = "sha256", serialize = TRUE)
  )
  identity_without_key <- record$identity
  identity_without_key$cache_key <- NULL
  identity_pass <- identical(
    record$cache_key,
    paste0(
      "mv05f_mapping_group_v1:",
      digest::digest(identity_without_key, algo = "sha256", serialize = TRUE)
    )
  )
  ids <- sort(c(
    record$identity$training_ids, record$identity$query_ids
  ), method = "radix")
  coordinates <- record$payload$coordinates
  coordinate_axis_pass <- is.list(coordinates) &&
    identical(sort(names(coordinates), method = "radix"), ids) &&
    all(vapply(names(coordinates), function(sample_id) {
      value <- coordinates[[sample_id]]
      is.matrix(value) && is.numeric(value) && nrow(value) == 384L &&
        ncol(value) == 30L && !anyNA(value) && all(is.finite(value)) &&
        identical(colnames(value), paste0("PC", 1:30)) &&
        all(startsWith(rownames(value), paste0(sample_id, "__")))
    }, logical(1L)))
  mappings <- record$payload$query_mappings
  active <- record$payload$active_features
  mapping_pass <- identical(
    sort(names(mappings), method = "radix"),
    sort(record$identity$query_ids, method = "radix")
  ) && identical(sort(names(active), method = "radix"),
                 sort(record$identity$query_ids, method = "radix")) &&
    all(vapply(names(mappings), function(sample_id) {
      value <- mappings[[sample_id]]
      mapping_identity <- list(
        contract_id = value$contract_id, engine_id = value$engine_id,
        fold_id = value$fold_id, fit_scope_id = value$fit_scope_id,
        reference_sample_ids = value$reference_sample_ids,
        held_out_sample_id = value$held_out_sample_id,
        features = value$features, dimensions = value$dimensions,
        seed = value$seed,
        outcome_label_state = value$outcome_label_state,
        biological_outcomes_computed = value$biological_outcomes_computed,
        reference_identity_sha256 = value$reference_identity_sha256,
        anchor_count = value$anchor_count,
        query_embedding_sha256 = value$query_embedding_sha256
      )
      identical(value$cache_key, paste0(
        "mv05_inductive_mapping_v1:",
        digest::digest(mapping_identity, algo = "sha256", serialize = TRUE)
      )) && value$held_out_sample_id == sample_id &&
        identical(value$reference_sample_ids,
                  sort(record$identity$training_ids, method = "radix")) &&
        identical(value$features, active[[sample_id]]) &&
        all(active[[sample_id]] %in% record$identity$panel$gene) &&
        identical(value$query_embeddings, coordinates[[sample_id]]) &&
        identical(value$query_embedding_sha256, digest::digest(
          value$query_embeddings, algo = "sha256", serialize = TRUE
        ))
    }, logical(1L)))
  reference_pass <- identical(
    record$payload$reference_identity_sha256_before,
    record$payload$reference_identity_sha256_after
  )
  forbidden <- c(
    "diagrams", "landscapes", "distances", "rankings", "clustering",
    "gene_views", "fusion", "outcomes", "labels", "tissue"
  )
  scope_pass <- !any(forbidden %in% names(record$payload)) &&
    identical(record$payload$downstream_execution,
      list(ph_jobs = 0L, landscape_jobs = 0L, distance_jobs = 0L,
           clustering_jobs = 0L, gene_view_jobs = 0L, fusion_jobs = 0L,
           biological_outcome_jobs = 0L)) &&
    identical(record$outcome_label_state, "closed") &&
    identical(record$biological_outcomes_computed, FALSE)
  file_hash_pass <- identical(
    digest::digest(file = path, algo = "sha256", serialize = FALSE),
    metric$private_result_sha256
  )
  checks[[index]] <- data.frame(
    contract_id = "mv05f_independent_group_validation_v1",
    group_id = job$group_id, held_out_study = job$held_out_study,
    seed = job$seed, payload_hash_pass = payload_hash_pass,
    identity_pass = identity_pass, coordinate_axis_pass = coordinate_axis_pass,
    mapping_pass = mapping_pass, reference_immutability_pass = reference_pass,
    file_hash_pass = file_hash_pass, scope_pass = scope_pass,
    coordinate_views = length(coordinates), query_mappings = length(mappings),
    minimum_active_features = min(vapply(active, length, integer(1L))),
    maximum_dropped_features = max(500L - vapply(active, length, integer(1L))),
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
}
detail <- do.call(rbind, checks)
pass_columns <- grep("_pass$", names(detail), value = TRUE)
all_pass <- all(as.matrix(detail[pass_columns]))
summary <- data.frame(
  contract_id = "mv05f_independent_validation_summary_v1",
  groups_validated = nrow(detail), coordinate_views = sum(detail$coordinate_views),
  query_mappings = sum(detail$query_mappings),
  validation_categories = length(pass_columns),
  all_validation_categories_pass = all_pass,
  label_transfer_jobs_executed = 0L, ph_jobs_executed = 0L,
  landscape_jobs_executed = 0L, distance_jobs_executed = 0L,
  clustering_jobs_executed = 0L, gene_view_jobs_executed = 0L,
  fusion_jobs_executed = 0L, new_data_jobs_executed = 0L,
  full_integrated_execution_jobs = 0L,
  biological_outcomes_computed = FALSE, outcome_label_state = "closed",
  stringsAsFactors = FALSE
)
if (!all_pass) stop("Independent MV5-F validation failed.", call. = FALSE)
write_provenance_csv(detail, detail_path)
write_provenance_csv(summary, summary_path)
message("Independently validated all MV5-F pilot mapping bundles.")
