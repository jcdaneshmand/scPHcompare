mv05d2_test_view <- function() {
  genes <- paste0("Gene", sprintf("%03d", seq_len(500L)))
  cells <- paste0("Cell", sprintf("%03d", seq_len(384L)))
  matrix <- outer(seq_len(500L), seq_len(384L), function(gene, cell) {
    sin(gene / 19 + cell / 13) + cos(gene / 31 - cell / 17)
  })
  dimnames(matrix) <- list(genes, cells)
  source <- new_cell_projection_source(
    matrix, "sample", "test", "sct_whole", "fold:test:training",
    20260805L, "training_only_standardization"
  )
  coordinates <- outer(seq_len(384L), seq_len(30L), function(cell, pc) {
    sin(cell / (pc + 3)) + cos(cell * pc / 401) + cell / 10000
  })
  dimnames(coordinates) <- list(cells, paste0("PC", seq_len(30L)))
  construct_frozen_cell_topology_view(
    source, coordinates, "mv05d1_training_fitted_pca_v1",
    "cell_pca_model_v1:test"
  )
}

mv05d2_test_identity <- function(view) {
  hash <- paste(rep("a", 64L), collapse = "")
  runtime <- list(
    contract_id = "mv05d2_ph_runtime_v1", r_version = "test",
    platform = "test", digest_version = "test", matrix_version = "test",
    ripserr_version = "test", omp_num_threads = "1",
    openblas_num_threads = "1", mkl_num_threads = "1"
  )
  mv05d2_ph_identity_v1(
    "job", "fold:test", "fold:test:training", "test", 20260805L,
    "sample", "held_out", 2L, "mv05d1:test", hash,
    view$cache_key, view$payload_sha256, hash, hash, runtime
  )
}

test_that("MV5-D2 H0 oracle equals known Euclidean MST edges", {
  points <- rbind(c(0, 0), c(1, 0), c(1, 1), c(3, 1))
  expect_equal(
    mv05d2_h0_mst_deaths_v1(points), c(1, 1, 2), tolerance = 1e-14
  )
})

test_that("MV5-D2 typed PH records are deterministic and stop at diagrams", {
  view <- mv05d2_test_view()
  identity <- mv05d2_test_identity(view)
  result <- run_topology_view_ph(view)
  record <- mv05d2_new_ph_record_v1(identity, view, result)
  repeated <- mv05d2_new_ph_record_v1(identity, view, result)

  expect_s3_class(record, "scph_mv05d2_cell_ph_record_v1")
  expect_invisible(mv05d2_validate_ph_record_v1(record, view))
  expect_invisible(mv05d3_validate_record_static_v1(record))
  expect_identical(record, repeated)
  expect_identical(record$h0_mst_oracle$finite_h0_intervals, 383L)
  expect_true(record$h0_mst_oracle$passed)
  expect_true(all(unlist(record$downstream_execution) == 0L))

  first <- tempfile(fileext = ".rds")
  second <- tempfile(fileext = ".rds")
  on.exit(unlink(c(first, second)), add = TRUE)
  saveRDS(record, first, compress = FALSE, version = 3)
  saveRDS(repeated, second, compress = FALSE, version = 3)
  expect_identical(
    digest::digest(file = first, algo = "sha256", serialize = FALSE),
    digest::digest(file = second, algo = "sha256", serialize = FALSE)
  )

  changed <- record
  changed$identity$seed <- 20260806L
  expect_error(mv05d2_validate_ph_record_v1(changed), "scientific contract")
  changed <- record
  changed$downstream_execution$landscape_jobs <- 1L
  expect_error(mv05d2_validate_ph_record_v1(changed), "stop boundary")
  changed <- record
  changed$topology_result$diagram[1L, "death"] <-
    changed$topology_result$diagram[1L, "death"] + 1
  changed$payload_sha256 <- digest::digest(
    list(identity = changed$identity,
         topology_result = changed$topology_result,
         h0_mst_oracle = changed$h0_mst_oracle),
    algo = "sha256", serialize = TRUE
  )
  changed$cache_key <- paste0(
    "mv05d2_cell_ph_record_v1:", changed$payload_sha256
  )
  expect_error(mv05d3_validate_record_static_v1(changed), "diagram")
})

test_that("MV5-D2 cache resume rejects stale identities", {
  view <- mv05d2_test_view()
  identity <- mv05d2_test_identity(view)
  record <- mv05d2_new_ph_record_v1(
    identity, view, run_topology_view_ph(view)
  )
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(record, path, compress = FALSE, version = 3)
  expect_identical(
    mv05d2_ph_cache_disposition_v1(path, identity$cache_key),
    "reuse_validated"
  )
  expect_error(
    mv05d2_ph_cache_disposition_v1(path, paste0(identity$cache_key, "x")),
    "stale identity"
  )
})

test_that("MV5-D2 resource gates and full projection remain label closed", {
  metrics <- data.frame(
    job_id = c("a", "b"), fold_id = c("f1", "f2"), seed = 1:2,
    execution_role = c("held_out", "training"),
    disposition = "completed", elapsed_seconds = c(1, 2),
    peak_process_tree_rss_bytes = c(100, 200),
    result_size_bytes = c(1000, 2000), h0_intervals = 384L,
    h1_intervals = c(5L, 7L), h0_mst_oracle_passed = TRUE,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    landscape_jobs_executed = 0L, distance_jobs_executed = 0L,
    clustering_jobs_executed = 0L, integration_jobs_executed = 0L,
    gene_view_jobs_executed = 0L, stringsAsFactors = FALSE
  )
  expect_invisible(mv05d2_validate_resource_metrics_v1(
    metrics, 2L, 10, 1000, 100, 10000
  ))
  projection <- mv05d2_project_full_ph_v1(metrics)
  expect_equal(nrow(projection), 3L)
  expect_true(all(projection$projected_jobs == 6750L))
  expect_true(all(projection$outcome_label_state == "closed"))
  expect_false(any(projection$biological_outcomes_computed))

  previous <- data.frame(
    scenario = "resource_safe_sct_cell_primary",
    normalization_worker_hours = 2.5,
    measured_cell_coordinate_worker_hours = 2.4,
    landscape_worker_hours = 3.6,
    known_components_lower_bound_worker_hours = 8.5,
    unmeasured_components = "cell_ph",
    planning_cap_with_10_percent_reserve_hours = 21.6,
    outcome_label_state = "closed",
    biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
  combined <- mv05d2_combine_primary_projection_v1(previous, projection)
  expect_equal(
    combined$projected_total_worker_hours,
    8.5 + projection$projected_worker_hours
  )
  expect_true(all(combined$cap_passes))
  expect_true(all(combined$full_production_jobs_launched == 0L))

  metrics$landscape_jobs_executed[[1L]] <- 1L
  expect_error(mv05d2_validate_resource_metrics_v1(
    metrics, 2L, 10, 1000, 100, 10000
  ), "scope")
})
