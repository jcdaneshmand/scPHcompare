mv05f_resource_fixture <- function() {
  studies <- sprintf("study_%02d", 1:15)
  query_n <- c(1L, 1L, 2L, 2L, 4L, 4L, 4L, 4L, 5L, 5L, 6L, 8L, 9L,
               10L, 25L)
  expand <- expand.grid(
    held_out_study = studies, seed = 20260805:20260809,
    stringsAsFactors = FALSE
  )
  expand$held_out_samples <- query_n[match(expand$held_out_study, studies)]
  expand$training_samples <- 90L - expand$held_out_samples
  expand$fold_id <- paste0("large_loso_v1:", expand$held_out_study)
  expand$fit_scope_id <- paste0(expand$fold_id, ":training")
  expand$held_out_missing_feature_instances <- 0L
  expand$maximum_missing_features_per_view <- 0L
  expand$held_out_missing_feature_instances[
    expand$held_out_study == "study_13"
  ] <- 3L
  expand$maximum_missing_features_per_view[
    expand$held_out_study == "study_13"
  ] <- 2L
  expand$held_out_missing_feature_instances[
    expand$held_out_study == "study_07"
  ] <- 1L
  expand$maximum_missing_features_per_view[
    expand$held_out_study == "study_07"
  ] <- 1L
  expand$private_cache_file <- paste0(
    expand$held_out_study, "__", expand$seed, ".rds"
  )
  expand$private_cache_sha256 <- paste(rep("a", 64L), collapse = "")
  expand$fold_cache_key <- paste0(
    "mv05d1_sct_cell_fold_v1:", paste(rep("b", 64L), collapse = "")
  )
  expand$outcome_label_state <- "closed"
  expand$biological_outcomes_computed <- FALSE
  expand
}

mv05f_mapping_fixture <- function() {
  sample_ids <- sprintf("sample_%03d", 1:90)
  query_id <- sample_ids[[90L]]
  training_ids <- sample_ids[-90L]
  dimensions <- 1:30
  coordinates <- lapply(sample_ids, function(sample_id) {
    matrix(
      seq_len(384L * 30L) / 1000,
      nrow = 384L, ncol = 30L,
      dimnames = list(
        paste0(sample_id, "__cell_", seq_len(384L)), paste0("PC", dimensions)
      )
    )
  })
  names(coordinates) <- sample_ids
  features <- paste0("gene_", 1:50)
  mapping <- .mv05_new_mapping_result(
    fold_id = "large_loso_v1:study", fit_scope_id = "fit",
    reference_sample_ids = training_ids, held_out_sample_id = query_id,
    features = features, dimensions = dimensions, seed = 20260805L,
    reference_identity_sha256 = paste(rep("c", 64L), collapse = ""),
    anchor_count = 5L, query_embeddings = coordinates[[query_id]]
  )
  identity <- list(
    contract_id = "mv05f_mapping_group_identity_v1",
    fold_id = "large_loso_v1:study", fit_scope_id = "fit",
    held_out_study = "study", seed = 20260805L,
    training_ids = training_ids, query_ids = query_id,
    panel = data.frame(
      feature_id = paste0("id_", 1:500), gene = paste0("gene_", 1:500),
      stringsAsFactors = FALSE
    ),
    dimensions = dimensions, cells_per_sample = 384L,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE
  )
  identity$cache_key <- paste0(
    "mv05f_mapping_group_v1:", .mv05f_sha256(identity)
  )
  mv05f_new_group_record_v1(
    identity, coordinates, stats::setNames(list(mapping), query_id),
    stats::setNames(list(features), query_id), "same", "same"
  )
}

test_that("MV5-F pilot selection is structural, fixed, and label closed", {
  manifest <- mv05f_select_pilot_groups_v1(mv05f_resource_fixture())
  expect_invisible(mv05f_validate_pilot_manifest_v1(manifest))
  expect_equal(manifest$held_out_study,
               c("study_01", "study_15", "study_13", "study_07"))
  expect_equal(manifest$held_out_samples, c(1L, 25L, 9L, 4L))
  leaking <- mv05f_resource_fixture()
  leaking$tissue <- "hidden"
  expect_error(mv05f_select_pilot_groups_v1(leaking), "label-closed")
})

test_that("MV5-F group records enforce immutable mapping-only coordinates", {
  record <- mv05f_mapping_fixture()
  expect_invisible(mv05f_validate_group_record_v1(record))
  changed <- record
  changed$payload$query_mappings[[1L]]$held_out_sample_id <- "wrong"
  changed$payload_sha256 <- .mv05f_sha256(changed$payload)
  expect_error(mv05f_validate_group_record_v1(changed), "incomplete")
  leaking <- record
  leaking$payload$outcomes <- data.frame(value = 1)
  leaking$payload_sha256 <- .mv05f_sha256(leaking$payload)
  expect_error(mv05f_validate_group_record_v1(leaking), "identity")
})

test_that("MV5-F resource and projection gates retain the stop boundary", {
  metrics <- data.frame(
    group_id = paste0("g", 1:4), held_out_study = paste0("s", 1:4),
    seed = 20260805L, disposition = "completed", exit_status = 0L,
    training_samples = c(89L, 65L, 82L, 86L),
    held_out_samples = c(1L, 25L, 8L, 4L),
    completed_query_mappings = 1:4, completed_coordinate_views = 90L,
    elapsed_seconds = c(100, 110, 120, 130),
    peak_process_tree_rss_bytes = 1000,
    private_result_bytes = 1000, input_seconds = 10,
    reference_sct_pca_seconds = c(8, 6, 7, 7.5),
    query_sct_seconds = c(1, 25, 8, 4),
    mapping_seconds = c(2, 50, 16, 8), assembly_seconds = 0.1,
    reference_immutable = TRUE,
    label_transfer_jobs_executed = 0L, ph_jobs_executed = 0L,
    landscape_jobs_executed = 0L, distance_jobs_executed = 0L,
    clustering_jobs_executed = 0L, gene_view_jobs_executed = 0L,
    fusion_jobs_executed = 0L, new_data_jobs_executed = 0L,
    biological_outcomes_computed = FALSE, outcome_label_state = "closed",
    stringsAsFactors = FALSE
  )
  expect_invisible(mv05f_validate_resource_metrics_v1(metrics))
  downstream <- data.frame(
    elapsed_seconds = c(10, 20), biological_outcomes_computed = FALSE,
    outcome_label_state = "closed"
  )
  projection <- mv05f_project_full_workload_v1(
    metrics, downstream, downstream, downstream
  )
  expect_equal(projection$projected_worker_seconds[[1L]], 115 * 75)
  expect_true(all(projection$outcome_label_state == "closed"))
  complete_folds <- mv05f_resource_fixture()
  downstream$private_result_bytes <- 100
  status <- data.frame(
    job_id = "j", representation = c("sct_fold", "inductive_integrated"),
    status = "completed", view_id = "cell_topology_v1",
    h0_finite_intervals = c(100, 90), h1_finite_intervals = c(50, 45),
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
  full <- mv05f_project_full_workload_v2(
    metrics, complete_folds, downstream, downstream, downstream, status
  )
  expect_equal(full$summary$conservative_geometry_multiplier, 1)
  expect_false(full$summary$full_execution_authorized)
  metrics$clustering_jobs_executed[[1L]] <- 1L
  expect_error(mv05f_validate_resource_metrics_v1(metrics), "scope")
})
