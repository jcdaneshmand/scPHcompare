test_that("MV5-I pair identity is dimension- and source-bound", {
  row <- data.frame(
    group_id = "group", fold_id = "fold", seed = 20260805L,
    homology_dimension = "H0", query_job_id = "query",
    training_job_id = "training", query_record_cache_key = "q-record",
    training_record_cache_key = "t-record",
    query_diagram_sha256 = paste(rep("a", 64L), collapse = ""),
    training_diagram_sha256 = paste(rep("b", 64L), collapse = ""),
    query_result_file_sha256 = paste(rep("c", 64L), collapse = ""),
    training_result_file_sha256 = paste(rep("d", 64L), collapse = ""),
    stringsAsFactors = FALSE
  )
  first <- .mv05i_pair_identity(row)
  changed <- row
  changed$homology_dimension <- "H1"
  expect_false(identical(first, .mv05i_pair_identity(changed)))
  changed <- row
  changed$query_result_file_sha256 <- paste(rep("e", 64L), collapse = "")
  expect_false(identical(first, .mv05i_pair_identity(changed)))
})

test_that("MV5-I rejects incomplete public pair manifests", {
  expect_error(
    mv05i_validate_pair_manifest_v1(data.frame()),
    "frozen contract"
  )
})

test_that("MV5-I measured projection replaces only the landscape term", {
  previous <- data.frame(
    measured_coordinate_worker_hours = 3,
    measured_ph_worker_hours = 1,
    measured_coordinate_ph_storage_bytes = 100,
    worker_hour_cap = 10, storage_cap_bytes = 1000,
    measured_peak_rss_bytes = 100, rss_cap_bytes = 1000,
    measured_max_group_seconds = 10, group_cap_seconds = 100,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE
  )
  completion <- data.frame(
    group_worker_seconds = 1800,
    peak_process_tree_rss_bytes = 200,
    private_total_landscape_stage_bytes = 100,
    maximum_group_seconds = 20,
    full_cell_landscape_distances_complete = TRUE,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE
  )
  retrieval <- data.frame(
    elapsed_seconds = 360, private_result_bytes = 40,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE
  )
  result <- mv05i_measured_primary_projection_v1(
    previous, completion, retrieval
  )
  expect_equal(result$measured_landscape_worker_hours, 0.5)
  expect_equal(result$projected_retrieval_worker_hours, 0.125)
  expect_equal(result$projected_total_worker_hours, 4.625)
  expect_equal(result$projected_total_storage_bytes, 250)
  expect_true(result$next_stage_authorized)
})
