test_that("minimum-cell boundary is recorded exactly", {
  flow <- new_sample_flow(
    sample_ids = c("below", "equal", "above"),
    cohort = "fixture",
    loaded_features = 10L,
    loaded_cells = c(249L, 250L, 251L)
  )
  flow <- record_pre_qc_sample_flow(flow, flow$sample_id, c(249L, 250L, 251L))
  flow <- record_post_qc_sample_flow(flow, flow$sample_id, c(249L, 250L, 251L), 250L)

  expect_equal(flow$ph_eligible, c(FALSE, TRUE, TRUE))
  expect_equal(
    flow$disposition,
    c("excluded_before_ph", "eligible_for_ph", "eligible_for_ph")
  )
  expect_equal(
    flow$reason_code,
    c("excluded_post_qc_min_cells", "", "")
  )
})

test_that("sample sets reconcile exactly", {
  expect_invisible(assert_sample_reconciliation(
    input_ids = c("a", "b", "c"),
    excluded_ids = "a",
    eligible_ids = c("b", "c"),
    completed_ids = "b",
    failed_ids = "c",
    strict_completion = FALSE
  ))

  expect_error(
    assert_sample_reconciliation(
      input_ids = c("a", "b"),
      excluded_ids = "a",
      eligible_ids = c("a", "b"),
      strict_completion = FALSE
    ),
    "overlap"
  )
  expect_error(
    assert_sample_reconciliation(
      input_ids = c("a", "b"),
      excluded_ids = character(),
      eligible_ids = "a",
      strict_completion = FALSE
    ),
    "do not reconcile"
  )
  expect_error(
    assert_sample_reconciliation(
      input_ids = c("a", "a"),
      excluded_ids = character(),
      eligible_ids = "a",
      strict_completion = FALSE
    ),
    "duplicated"
  )
})

test_that("retry history retains attempts and summarizes final success", {
  started <- as.POSIXct("2026-08-04 12:00:00", tz = "UTC")
  failed <- new_ph_attempt(
    "fixture", "raw", "sample-a", 1L, 1L, -1, c(20L, 10L),
    started, started + 1, exit_status = 1L, status = "ph_child_error",
    error_message = "fixture failure"
  )
  completed <- new_ph_attempt(
    "fixture", "raw", "sample-a", 1L, 2L, 10, c(20L, 10L),
    started + 2, started + 4, exit_status = 0L, pd_written = TRUE,
    status = "completed"
  )
  attempts <- rbind(failed, completed)
  final <- summarize_ph_attempts(attempts)

  expect_equal(nrow(attempts), 2L)
  expect_equal(nrow(final), 1L)
  expect_equal(final$attempt, 2L)
  expect_equal(final$status, "completed")
})

test_that("attempt numbers continue from the legacy progress history", {
  progress_path <- tempfile(fileext = ".csv")
  utils::write.csv(
    data.frame(
      job_id = c(1L, 1L, 2L),
      status = c("error", "failed", "completed")
    ),
    progress_path,
    row.names = FALSE
  )

  expect_equal(next_ph_attempt_number(progress_path, 1L), 3L)
  expect_equal(next_ph_attempt_number(progress_path, 2L), 2L)
  expect_equal(next_ph_attempt_number(tempfile(), 1L), 1L)
})

test_that("terminal failure remains distinct from completion", {
  started <- as.POSIXct("2026-08-04 12:00:00", tz = "UTC")
  attempts <- new_ph_attempt(
    "fixture", "raw", "sample-a", 1L, 1L, -1, c(20L, 10L),
    started, started + 1, exit_status = 1L, status = "ph_child_error"
  )
  final <- summarize_ph_attempts(attempts)

  expect_equal(final$status, "ph_child_error")
  expect_error(
    assert_sample_reconciliation(
      input_ids = "sample-a",
      excluded_ids = character(),
      eligible_ids = "sample-a",
      completed_ids = character(),
      failed_ids = "sample-a",
      strict_completion = TRUE
    ),
    "did not complete"
  )
})

test_that("nearest-neighbor k is bounded by observations", {
  expect_equal(bounded_knn_k(101L, 100L), 100L)
  expect_equal(bounded_knn_k(50L, 100L), 49L)
  expect_error(bounded_knn_k(1L, 100L), "At least two observations")
  expect_error(bounded_knn_k(10L, 0L), "positive integer")
})

test_that("provenance CSV round trip preserves schema and attempt order", {
  flow_path <- tempfile(fileext = ".csv")
  attempt_path <- tempfile(fileext = ".csv")
  flow <- new_sample_flow("sample-a", "fixture", 12, 20L, 10L)
  flow <- record_post_qc_sample_flow(flow, "sample-a", 10L, 5L)
  write_provenance_csv(flow, flow_path)

  started <- as.POSIXct("2026-08-04 12:00:00", tz = "UTC")
  attempt_one <- new_ph_attempt(
    "fixture", "raw", "sample-a", 1L, 1L, -1, c(20L, 10L),
    started, started + 1, status = "ph_empty_output"
  )
  attempt_two <- new_ph_attempt(
    "fixture", "raw", "sample-a", 1L, 2L, 5, c(20L, 10L),
    started + 2, started + 3, exit_status = 0L, pd_written = TRUE,
    status = "completed"
  )
  append_ph_attempts(attempt_one, attempt_path)
  append_ph_attempts(attempt_two, attempt_path)

  round_trip_flow <- utils::read.csv(flow_path, stringsAsFactors = FALSE)
  round_trip_attempts <- utils::read.csv(attempt_path, stringsAsFactors = FALSE)
  expect_equal(round_trip_flow$sample_id, "sample-a")
  expect_equal(round_trip_flow$loaded_features, 20L)
  expect_equal(round_trip_attempts$attempt, c(1L, 2L))
  expect_equal(round_trip_attempts$status, c("ph_empty_output", "completed"))
})

test_that("PH attempt logs accept additive observability columns", {
  attempt_path <- tempfile(fileext = ".csv")
  legacy <- new_ph_attempt(
    "fixture", "raw", "sample-a", 1L, 1L, -1, c(4L, 2L),
    Sys.time(), Sys.time(), status = "completed"
  )
  legacy <- legacy[, setdiff(
    names(legacy),
    c(
      "poll_interval_seconds", "memory_samples", "monitor_peak_rss_bytes",
      "child_peak_rss_bytes", "descendant_peak_rss_bytes",
      "process_tree_peak_rss_bytes", "process_tree_peak_count"
    )
  )]
  write_provenance_csv(legacy, attempt_path)

  observed <- new_ph_attempt(
    "fixture", "raw", "sample-b", 2L, 1L, -1, c(4L, 2L),
    Sys.time(), Sys.time(), status = "completed",
    poll_interval_seconds = 0.25, memory_samples = 2L,
    monitor_peak_rss_bytes = 100, child_peak_rss_bytes = 80,
    descendant_peak_rss_bytes = 20, process_tree_peak_rss_bytes = 100,
    process_tree_peak_count = 2L
  )
  expect_silent(append_ph_attempts(observed, attempt_path))
  combined <- utils::read.csv(attempt_path, stringsAsFactors = FALSE)
  expect_equal(combined$sample_id, c("sample-a", "sample-b"))
  expect_true(is.na(combined$process_tree_peak_rss_bytes[[1]]))
  expect_equal(combined$process_tree_peak_rss_bytes[[2]], 100)
})
