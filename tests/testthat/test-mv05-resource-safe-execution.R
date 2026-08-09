make_mv05c2_fold_table <- function() {
  samples <- data.frame(
    sample_id = c("a1", "a2", "b1", "c1"),
    study = c("a", "a", "b", "c"), stringsAsFactors = FALSE
  )
  do.call(rbind, lapply(sort(unique(samples$study)), function(held_out) {
    data.frame(
      fold_id = paste0("fold:", held_out),
      fit_scope_id = paste0("fold:", held_out, ":training"),
      sample_id = samples$sample_id,
      execution_role = ifelse(
        samples$study == held_out, "held_out_query", "training_reference"
      ),
      outcome_label_state = "closed", stringsAsFactors = FALSE
    )
  }))
}

test_that("sample-seed normalization identities are immutable and label-closed", {
  identity <- mv05c2_normalization_cache_identity_v1(
    sample_id = "sample_a", seed = 20260805L,
    selected_cell_sha256 = paste(rep("a", 64), collapse = ""),
    source_cache_sha256 = paste(rep("b", 64), collapse = ""),
    seurat_version = "5.3.0"
  )
  record <- mv05c2_new_normalization_cache_record_v1(
    identity, list(matrix = matrix(1:6, nrow = 2L))
  )
  expect_s3_class(record, "scph_mv05c2_normalization_cache_v1")
  expect_false(record$biological_outcomes_computed)
  expect_invisible(mv05c2_validate_normalization_cache_record_v1(record))

  tampered <- record
  tampered$payload$matrix[1, 1] <- 99
  expect_error(
    mv05c2_validate_normalization_cache_record_v1(tampered), "stale"
  )
  changed <- mv05c2_normalization_cache_identity_v1(
    sample_id = "sample_a", seed = 20260806L,
    selected_cell_sha256 = paste(rep("a", 64), collapse = ""),
    source_cache_sha256 = paste(rep("b", 64), collapse = ""),
    seurat_version = "5.3.0"
  )
  expect_false(identical(identity$cache_key, changed$cache_key))
})

test_that("normalization cache resumability refuses stale files", {
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  identity <- mv05c2_normalization_cache_identity_v1(
    "sample_a", 1L, paste(rep("a", 64), collapse = ""),
    paste(rep("b", 64), collapse = ""), "5.3.0"
  )
  expect_identical(
    mv05c2_cache_disposition_v1(path, identity$cache_key), "build_missing"
  )
  record <- mv05c2_new_normalization_cache_record_v1(identity, list(x = 1))
  saveRDS(record, path)
  expect_identical(
    mv05c2_cache_disposition_v1(path, identity$cache_key), "reuse_validated"
  )
  expect_error(
    mv05c2_cache_disposition_v1(
      path, paste0(identity$cache_key, "_stale")
    ),
    "stale identity"
  )
})

test_that("MV5-D0 SCT matrix payloads retain legacy cache semantics", {
  identity <- mv05c2_normalization_cache_identity_v1(
    "sample_a", 1L, paste(rep("a", 64), collapse = ""),
    paste(rep("b", 64), collapse = ""), "5.3.0"
  )
  value <- Matrix::Matrix(matrix(seq_len(768L), nrow = 2L), sparse = TRUE)
  rownames(value) <- c("g1", "g2")
  colnames(value) <- paste0("cell", seq_len(384L))
  record <- mv05c2_new_normalization_cache_record_v1(
    identity,
    list(payload_contract_id = "mv05d0_sct_data_matrix_v1",
         sct_data = value, selected_cell_ids = colnames(value))
  )
  expect_identical(mv05d0_sct_matrix_from_cache_v1(record), value)
  broken <- record
  broken$payload$sct_data <- broken$payload$sct_data[, -1L]
  broken$payload_sha256 <- digest::digest(
    broken$payload, algo = "sha256", serialize = TRUE
  )
  expect_error(mv05d0_sct_matrix_from_cache_v1(broken), "invalid")
})

test_that("MV5-D0 runtime-complete identities cannot collide with legacy keys", {
  runtime <- list(
    contract_id = "mv05d0_normalization_runtime_v1",
    r_version = "R test", rng_kind = c("Mersenne-Twister", "Inversion", "Rejection"),
    seurat_version = "5.3.0", seurat_object_version = "5.1.0",
    sctransform_version = "0.4.1", matrix_version = "1.7-3",
    blas = "test-blas", lapack = "test-lapack",
    future_plan = "sequential", omp_num_threads = "1",
    openblas_num_threads = "1", mkl_num_threads = "1"
  )
  sha_a <- paste(rep("a", 64), collapse = "")
  sha_b <- paste(rep("b", 64), collapse = "")
  identity <- mv05d0_normalization_cache_identity_v2(
    "sample_a", 1L, sha_a, sha_b, runtime
  )
  value <- Matrix::Matrix(matrix(seq_len(768L), nrow = 2L), sparse = TRUE)
  rownames(value) <- c("g1", "g2")
  colnames(value) <- paste0("cell", seq_len(384L))
  record <- mv05d0_new_normalization_cache_record_v2(
    identity, list(payload_contract_id = "mv05d0_sct_data_matrix_v1",
                   sct_data = value, selected_cell_ids = colnames(value))
  )
  expect_invisible(mv05d0_validate_normalization_cache_record_v2(record))
  expect_identical(mv05d0_sct_matrix_from_cache_v1(record), value)
  legacy <- mv05c2_normalization_cache_identity_v1(
    "sample_a", 1L, sha_a, sha_b, "5.3.0"
  )
  expect_false(identical(identity$cache_key, legacy$cache_key))
  changed <- runtime
  changed$sctransform_version <- "0.4.2"
  expect_false(identical(
    identity$cache_key,
    mv05d0_normalization_cache_identity_v2(
      "sample_a", 1L, sha_a, sha_b, changed
    )$cache_key
  ))
})

test_that("MV5-D0 matched-cell summaries are deterministic and label-closed", {
  first <- matrix(seq_len(2L * 400L), nrow = 2L)
  second <- first + 1
  rownames(first) <- rownames(second) <- c("g1", "g2")
  colnames(first) <- paste0("a", seq_len(400L))
  colnames(second) <- paste0("b", seq_len(400L))
  samples <- list(b = second, a = first)
  observed <- mv05d0_build_selection_summary_v1(
    samples, c("b", "a"), seeds = c(2L, 1L), n = 384L
  )
  repeated <- mv05d0_build_selection_summary_v1(
    samples, c("a", "b"), seeds = c(1L, 2L), n = 384L
  )
  expect_identical(observed, repeated)
  expect_equal(nrow(observed), 4L)
  expect_true(all(observed$outcome_label_state == "closed"))
  expect_false(any(observed$biological_outcomes_computed))
})

test_that("MV5-D0 completion rejects missing entries and cap breaches", {
  metrics <- data.frame(
    sample_id = c("a", "b"), seed = 1L,
    disposition = "built_atomic", elapsed_seconds = c(10, 11),
    peak_process_tree_rss_bytes = c(100, 200),
    private_cache_size_bytes = c(1000, 2000),
    outcome_label_state = "closed",
    biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
  expect_invisible(mv05d0_validate_resource_metrics_v1(
    metrics, expected_entries = 2L, elapsed_cap_seconds = 20,
    rss_cap_bytes = 1000, storage_cap_bytes = 10000
  ))
  expect_error(mv05d0_validate_resource_metrics_v1(
    metrics, expected_entries = 3L, elapsed_cap_seconds = 20,
    rss_cap_bytes = 1000, storage_cap_bytes = 10000
  ), "resource gates")
  metrics$peak_process_tree_rss_bytes[[2L]] <- 2000
  expect_error(mv05d0_validate_resource_metrics_v1(
    metrics, expected_entries = 2L, elapsed_cap_seconds = 20,
    rss_cap_bytes = 1000, storage_cap_bytes = 10000
  ), "resource gates")
})

test_that("MV5-D0 raw shard v2 hashes numerical count content canonically", {
  counts <- matrix(c(0, 1, 2, 0, 3, 0), nrow = 2L)
  rownames(counts) <- c("g1", "g2")
  colnames(counts) <- c("c1", "c2", "c3")
  sparse <- Matrix::Matrix(counts, sparse = TRUE)
  expect_identical(
    mv05d0_count_matrix_sha256_v1(counts),
    mv05d0_count_matrix_sha256_v1(sparse)
  )
  record <- mv05d0_new_raw_sample_cache_v2(
    "sample_a", sparse, paste(rep("a", 64), collapse = ""),
    paste(rep("b", 64), collapse = "")
  )
  expect_invisible(mv05d0_validate_raw_sample_cache_v2(record))
  changed <- record
  changed$counts[1, 1] <- 9
  expect_error(mv05d0_validate_raw_sample_cache_v2(changed), "stale")
})

test_that("MV5-D0 post-cache reprojection replaces only measured normalization", {
  previous <- data.frame(
    contract_id = "old",
    scenario = c("naive_full_mv05d", "resource_safe_sct_cell_primary"),
    normalization_worker_hours = c(NA, 10),
    cached_sct_fold_worker_hours = c(NA, 4),
    landscape_worker_hours = c(90, 3),
    integrated_reference_mapping_worker_hours = c(NA, 0),
    projected_lower_bound_worker_hours = c(240, 17),
    nominal_cap_hours = 24,
    planning_cap_with_10_percent_reserve_hours = 21.6,
    cap_passes = c(FALSE, TRUE), disposition = c("old", "old"),
    outcome_label_state = "closed",
    biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
  observed <- mv05d0_reproject_scenarios_v1(previous, 2.5)
  expect_true(is.na(observed$normalization_worker_hours[[1L]]))
  expect_equal(observed$normalization_worker_hours[[2L]], 2.5)
  expect_equal(observed$projected_lower_bound_worker_hours[[2L]], 9.5)
  expect_true(observed$cap_passes[[2L]])
  expect_match(observed$disposition[[2L]], "future_label_closed")
})

test_that("query-training scope is deterministic and excludes unsupported pairs", {
  strata <- data.frame(
    representation = c("sct_fold", "sct_fold"),
    view_id = c("cell_topology_v1", "cell_topology_v1"),
    homology_dimension = c("H0", "H1"), stringsAsFactors = FALSE
  )
  pairs <- mv05c2_build_query_training_pairs_v1(
    make_mv05c2_fold_table(), seeds = c(2L, 1L), strata = strata
  )
  # Per seed and stratum: 2*2 + 1*3 + 1*3 = 10 query-training pairs.
  expect_equal(nrow(pairs), 40L)
  expect_true(all(pairs$query_sample_id != pairs$training_sample_id))
  expect_true(all(pairs$supports_primary_retrieval))
  expect_false(any(pairs$supports_full_matrix_clustering))
  expect_false(any(pairs$supports_within_study_pair_contrasts))
  repeated <- mv05c2_build_query_training_pairs_v1(
    make_mv05c2_fold_table(), seeds = c(1L, 2L), strata = strata
  )
  expect_identical(pairs, repeated)

  leaking <- make_mv05c2_fold_table()
  leaking$tissue <- "hidden_label"
  expect_error(
    mv05c2_build_query_training_pairs_v1(leaking, 1L, strata),
    "label-closed"
  )
})

test_that("pair chunks are bounded, immutable, and deterministic", {
  strata <- data.frame(
    representation = "sct_fold", view_id = "cell_topology_v1",
    homology_dimension = "H0", stringsAsFactors = FALSE
  )
  pairs <- mv05c2_build_query_training_pairs_v1(
    make_mv05c2_fold_table(), seeds = 1:2, strata = strata
  )
  chunks <- mv05c2_assign_pair_chunks_v1(pairs, max_pairs = 6L)
  summary <- mv05c2_pair_chunk_summary_v1(chunks)
  expect_true(all(summary$pair_count <= 6L))
  expect_equal(sum(summary$pair_count), nrow(pairs))
  expect_true(all(grepl("^mv05c2_pair_chunk_v1:", summary$chunk_id)))
  expect_identical(
    chunks, mv05c2_assign_pair_chunks_v1(pairs[nrow(pairs):1, ], 6L)
  )

  invalid <- pairs
  invalid$exact[[1L]] <- FALSE
  expect_error(mv05c2_assign_pair_chunks_v1(invalid), "chunk contract")
})

test_that("MV5-D1 feature eligibility is strictly training-derived", {
  training_a <- rbind(
    GeneA = c(0, 1, 2, 3), GeneB = c(0, 2, 4, 6),
    GeneC = c(0, 4, 8, 12), GeneD = c(0, 8, 16, 24)
  )
  training_b <- sweep(training_a, 1L, c(0, 1, 0, 2), "+")
  held_out_full <- training_a + 3
  held_out_missing <- held_out_full[c("GeneA", "GeneB"), , drop = FALSE]
  full <- list(a = training_a, b = training_b, q = held_out_full)
  missing <- list(a = training_a, b = training_b, q = held_out_missing)
  first <- mv05c2_select_training_panel_v1(full, c("a", "b"), 2L)
  second <- mv05c2_select_training_panel_v1(missing, c("a", "b"), 2L)
  expect_identical(first, second)
  expect_true(any(!first$feature_id %in% rownames(held_out_missing)))
})

test_that("MV5-D1 standardization is numerically independent of held-out values", {
  genes <- paste0("Gene", sprintf("%03d", seq_len(500L)))
  cells <- paste0("Cell", sprintf("%03d", seq_len(384L)))
  base <- outer(seq_len(500L), seq_len(384L), function(gene, cell) {
    sin(gene / 17 + cell / 11) + cos(gene / 23 - cell / 7)
  })
  rownames(base) <- genes
  colnames(base) <- cells
  matrices_a <- list(train_a = base, train_b = base * 1.1 + 0.2,
                     query = base * 0.7 - 0.1)
  matrices_b <- matrices_a
  matrices_a$query <- matrices_a$query[-1L, , drop = FALSE]
  matrices_b$query <- matrices_b$query[-1L, , drop = FALSE] + 100
  panel <- data.frame(
    feature_id = genes, gene = genes, stringsAsFactors = FALSE
  )
  keys <- stats::setNames(
    paste0("mv05d0_sample_seed_sct_v2:", c(
      paste(rep("a", 64L), collapse = ""),
      paste(rep("b", 64L), collapse = ""),
      paste(rep("c", 64L), collapse = "")
    )), names(matrices_a)
  )
  first <- mv05d1_prepare_cell_sources_v1(
    matrices_a, panel, c("train_a", "train_b"), "fold:q",
    "fold:q:training", 1L, keys, cohort = "test"
  )
  second <- mv05d1_prepare_cell_sources_v1(
    matrices_b, panel, c("train_a", "train_b"), "fold:q",
    "fold:q:training", 1L, keys, cohort = "test"
  )
  expect_identical(first$center, second$center)
  expect_identical(first$scale, second$scale)
  expect_identical(first$standardization_id, second$standardization_id)
  expect_identical(first$missing_feature_counts[["query"]], 1L)
  expect_true(all(first$cell_sources$query$matrix[1L, ] == 0))
  expect_identical(
    first$training_sources$train_a$input_sha256,
    second$training_sources$train_a$input_sha256
  )
  expect_false(identical(
    first$cell_sources$query$input_sha256,
    second$cell_sources$query$input_sha256
  ))
})

test_that("MV5-D1 fold identities and resource gates are immutable and closed", {
  keys <- stats::setNames(
    paste0("mv05d0_sample_seed_sct_v2:", c(
      paste(rep("a", 64L), collapse = ""),
      paste(rep("b", 64L), collapse = ""),
      paste(rep("c", 64L), collapse = "")
    )), c("a", "b", "q")
  )
  runtime <- list(
    contract_id = "mv05d1_fold_runtime_v1", r_version = "R test",
    rng_kind = c("Mersenne-Twister", "Inversion", "Rejection"),
    matrix_version = "1", digest_version = "1", blas = "test",
    lapack = "test", omp_num_threads = "1", openblas_num_threads = "1",
    mkl_num_threads = "1"
  )
  hash <- paste(rep("d", 64L), collapse = "")
  identity <- mv05d1_cell_fold_identity_v1(
    "fold:q", "fold:q:training", "q", 1L, c("b", "a"), "q", keys,
    hash, hash, hash, runtime, panel_size = 500L, n_components = 30L
  )
  expect_identical(identity$training_ids, c("a", "b"))
  expect_identical(identity$outcome_label_state, "closed")
  changed <- identity
  changed$seed <- 2L
  expect_false(identical(
    identity$cache_key,
    mv05d1_cell_fold_identity_v1(
      "fold:q", "fold:q:training", "q", 2L, c("a", "b"), "q", keys,
      hash, hash, hash, runtime, panel_size = 500L, n_components = 30L
    )$cache_key
  ))

  metrics <- data.frame(
    fold_id = "fold:q", seed = 1L, disposition = "built_atomic",
    elapsed_seconds = 10, peak_process_tree_rss_bytes = 100,
    private_cache_size_bytes = 1000, outcome_label_state = "closed",
    biological_outcomes_computed = FALSE, ph_jobs_executed = 0L,
    landscape_jobs_executed = 0L, distance_jobs_executed = 0L,
    clustering_jobs_executed = 0L, integration_jobs_executed = 0L,
    gene_view_jobs_executed = 0L, stringsAsFactors = FALSE
  )
  expect_invisible(mv05d1_validate_resource_metrics_v1(
    metrics, 1L, 20, 1000, 10000
  ))
  metrics$ph_jobs_executed <- 1L
  expect_error(mv05d1_validate_resource_metrics_v1(
    metrics, 1L, 20, 1000, 10000
  ), "scope")
})

test_that("MV5-D1 projection keeps unmeasured PH explicit", {
  previous <- data.frame(
    contract_id = "old",
    scenario = c("naive_full_mv05d", "resource_safe_sct_cell_primary"),
    normalization_worker_hours = c(NA, 2.5),
    cached_sct_fold_worker_hours = c(NA, 4),
    landscape_worker_hours = c(90, 3.5),
    integrated_reference_mapping_worker_hours = c(NA, 0),
    projected_lower_bound_worker_hours = c(240, 10),
    nominal_cap_hours = 24,
    planning_cap_with_10_percent_reserve_hours = 21.6,
    cap_passes = c(FALSE, TRUE), disposition = c("old", "old"),
    outcome_label_state = "closed",
    biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
  )
  observed <- mv05d1_post_fold_projection_v2(previous, 1.25)
  expect_true(is.na(observed$measured_cell_coordinate_worker_hours[[1L]]))
  expect_equal(observed$measured_cell_coordinate_worker_hours[[2L]], 1.25)
  expect_equal(observed$known_components_lower_bound_worker_hours[[2L]], 7.25)
  expect_identical(observed$unmeasured_components[[2L]], "cell_ph")
  expect_false(observed$feasibility_complete[[2L]])
  expect_true(is.na(observed$cap_passes[[2L]]))
  expect_match(observed$disposition[[2L]], "pending_measured_production_cell_ph")
})
