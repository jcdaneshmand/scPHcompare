test_that("MV5-D4 pair identity is dimension- and source-bound", {
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
  first <- .mv05d4_pair_identity(row)
  changed <- row
  changed$homology_dimension <- "H1"
  expect_false(identical(first, .mv05d4_pair_identity(changed)))
  changed <- row
  changed$query_result_file_sha256 <- paste(rep("e", 64L), collapse = "")
  expect_false(identical(first, .mv05d4_pair_identity(changed)))
})

test_that("MV5-D4 rejects incomplete public pair manifests", {
  expect_error(
    mv05d4_validate_pair_manifest_v1(data.frame()),
    "frozen contract"
  )
})

test_that("MV5-D4 measured projection replaces only the landscape term", {
  previous <- data.frame(
    normalization_worker_hours = 2,
    measured_cell_coordinate_worker_hours = 3,
    measured_cell_ph_worker_hours = 1,
    projected_landscape_worker_hours = 4,
    planning_cap_with_10_percent_reserve_hours = 10,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE
  )
  completion <- data.frame(
    group_worker_seconds = 1800,
    full_cell_landscape_distances_complete = TRUE,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE
  )
  result <- mv05d4_measured_primary_projection_v1(previous, completion)
  expect_equal(result$measured_landscape_worker_hours, 0.5)
  expect_equal(result$measured_total_worker_hours, 6.5)
  expect_equal(result$planning_cap_margin_hours, 3.5)
  expect_true(result$cap_passes)
})
