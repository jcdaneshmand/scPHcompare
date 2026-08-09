# Internal MV5-G full label-closed integrated-coordinate production helpers.

mv05g_build_full_manifest_v1 <- function(resources) {
  required <- c(
    "fold_id", "fit_scope_id", "held_out_study", "seed",
    "training_samples", "held_out_samples", "panel_genes",
    "cells_per_view", "pca_components", "private_cache_file",
    "private_cache_sha256", "fold_cache_key", "outcome_label_state",
    "biological_outcomes_computed"
  )
  if (!is.data.frame(resources) || !all(required %in% names(resources))) {
    stop("MV5-G requires the accepted MV5-D1 resource schema.", call. = FALSE)
  }
  .mv05f_assert_label_closed(resources, "MV5-D1 resources")
  resources <- resources[order(
    resources$held_out_study, resources$seed, method = "radix"
  ), , drop = FALSE]
  resources$contract_id <- "mv05g_integrated_coordinate_manifest_v1"
  resources$group_order <- seq_len(nrow(resources))
  resources$group_id <- paste0(
    "mv05g_group__", resources$held_out_study, "__", resources$seed
  )
  resources$production_role <- "full_label_closed_integrated_coordinate_cache"
  resources$panel_genes <- 500L
  resources$cells_per_sample <- 384L
  resources$pca_components <- 30L
  resources$group_timeout_seconds <- 1800
  resources$rss_cap_bytes <- 8 * 1024^3
  resources$stage_cap_seconds <- 43200
  resources$storage_cap_bytes <- 10 * 1000^3
  resources$maximum_heavy_workers <- 1L
  resources$label_transfer_jobs_executed <- 0L
  resources$ph_jobs_executed <- 0L
  resources$landscape_jobs_executed <- 0L
  resources$distance_jobs_executed <- 0L
  resources$retrieval_jobs_executed <- 0L
  resources$clustering_jobs_executed <- 0L
  resources$gene_view_jobs_executed <- 0L
  resources$fusion_jobs_executed <- 0L
  resources$new_data_jobs_executed <- 0L
  resources$biological_outcomes_computed <- FALSE
  resources$outcome_label_state <- "closed"
  rownames(resources) <- NULL
  mv05g_validate_full_manifest_v1(resources)
  resources
}

mv05g_validate_full_manifest_v1 <- function(manifest) {
  required <- c(
    "contract_id", "group_id", "group_order", "production_role", "fold_id",
    "fit_scope_id", "held_out_study", "seed", "training_samples",
    "held_out_samples", "panel_genes", "cells_per_sample", "pca_components",
    "private_cache_file", "private_cache_sha256", "fold_cache_key",
    "group_timeout_seconds", "rss_cap_bytes", "stage_cap_seconds",
    "storage_cap_bytes", "maximum_heavy_workers", "outcome_label_state",
    "biological_outcomes_computed", "label_transfer_jobs_executed",
    "ph_jobs_executed", "landscape_jobs_executed", "distance_jobs_executed",
    "retrieval_jobs_executed", "clustering_jobs_executed",
    "gene_view_jobs_executed", "fusion_jobs_executed",
    "new_data_jobs_executed"
  )
  zero <- c(
    "label_transfer_jobs_executed", "ph_jobs_executed",
    "landscape_jobs_executed", "distance_jobs_executed",
    "retrieval_jobs_executed", "clustering_jobs_executed",
    "gene_view_jobs_executed", "fusion_jobs_executed",
    "new_data_jobs_executed"
  )
  if (!is.data.frame(manifest) || !all(required %in% names(manifest)) ||
      nrow(manifest) != 75L ||
      !identical(as.integer(manifest$group_order), 1:75) ||
      any(manifest$contract_id != "mv05g_integrated_coordinate_manifest_v1") ||
      any(manifest$production_role !=
          "full_label_closed_integrated_coordinate_cache") ||
      anyDuplicated(manifest$group_id) ||
      anyDuplicated(paste(manifest$held_out_study, manifest$seed, sep = "\r")) ||
      length(unique(manifest$held_out_study)) != 15L ||
      length(unique(manifest$seed)) != 5L ||
      !identical(sort(unique(as.integer(manifest$seed))), 20260805:20260809) ||
      any(manifest$training_samples + manifest$held_out_samples != 90L) ||
      any(manifest$panel_genes != 500L) ||
      any(manifest$cells_per_sample != 384L) ||
      any(manifest$pca_components != 30L) ||
      any(manifest$maximum_heavy_workers != 1L) ||
      any(manifest$group_timeout_seconds != 1800) ||
      any(manifest$rss_cap_bytes != 8 * 1024^3) ||
      any(manifest$stage_cap_seconds != 43200) ||
      any(manifest$storage_cap_bytes != 10 * 1000^3) ||
      any(as.matrix(manifest[zero]) != 0) ||
      any(!grepl("^[0-9a-f]{64}$", manifest$private_cache_sha256)) ||
      any(!grepl("^mv05d1_sct_cell_fold_v1:[0-9a-f]{64}$",
                 manifest$fold_cache_key))) {
    stop("MV5-G full manifest violates its frozen 75-group contract.",
         call. = FALSE)
  }
  tab <- table(manifest$held_out_study)
  if (any(tab != 5L)) {
    stop("Every MV5-G study must appear at all five seeds.", call. = FALSE)
  }
  .mv05f_assert_label_closed(manifest, "MV5-G full manifest")
  invisible(manifest)
}

mv05g_validate_group_metrics_v1 <- function(
    metrics, expected_groups = 75L, elapsed_cap_seconds = 1800,
    rss_cap_bytes = 8 * 1024^3, storage_cap_bytes = 10 * 1000^3) {
  required <- c(
    "group_id", "held_out_study", "seed", "disposition", "exit_status",
    "completed_query_mappings", "completed_coordinate_views",
    "elapsed_seconds", "peak_process_tree_rss_bytes", "private_result_bytes",
    "reference_immutable", "label_transfer_jobs_executed",
    "ph_jobs_executed", "landscape_jobs_executed", "distance_jobs_executed",
    "retrieval_jobs_executed", "clustering_jobs_executed",
    "gene_view_jobs_executed", "fusion_jobs_executed",
    "new_data_jobs_executed", "biological_outcomes_computed",
    "outcome_label_state"
  )
  zero <- c(
    "label_transfer_jobs_executed", "ph_jobs_executed",
    "landscape_jobs_executed", "distance_jobs_executed",
    "retrieval_jobs_executed", "clustering_jobs_executed",
    "gene_view_jobs_executed", "fusion_jobs_executed",
    "new_data_jobs_executed"
  )
  if (!is.data.frame(metrics) || !all(required %in% names(metrics)) ||
      nrow(metrics) != as.integer(expected_groups) ||
      anyDuplicated(metrics$group_id) ||
      anyDuplicated(paste(metrics$held_out_study, metrics$seed, sep = "\r")) ||
      any(metrics$disposition != "completed") || any(metrics$exit_status != 0L) ||
      any(metrics$completed_coordinate_views != 90L) ||
      any(metrics$completed_query_mappings != metrics$held_out_samples) ||
      any(metrics$elapsed_seconds > elapsed_cap_seconds) ||
      any(metrics$peak_process_tree_rss_bytes > rss_cap_bytes) ||
      sum(metrics$private_result_bytes) > storage_cap_bytes ||
      any(!as.logical(metrics$reference_immutable)) ||
      any(as.matrix(metrics[zero]) != 0) ||
      any(metrics$outcome_label_state != "closed") ||
      any(as.logical(metrics$biological_outcomes_computed))) {
    stop("MV5-G metrics violate completion, scope, or resource gates.",
         call. = FALSE)
  }
  invisible(metrics)
}
