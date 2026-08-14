test_that("real PH subprocesses produce auditable success and failure artifacts", {
  skip_on_cran()
  skip_if_not_installed("processx")
  skip_if_not_installed("ripserr")

  fixture_dir <- tempfile("scphcompare-subprocess-")
  dir.create(fixture_dir, recursive = TRUE)
  on.exit(unlink(fixture_dir, recursive = TRUE, force = TRUE), add = TRUE)

  results_file <- file.path(fixture_dir, "intermediate_results.rds")
  progress_file <- file.path(fixture_dir, "progress.csv")
  attempt_file <- file.path(fixture_dir, "ph_attempt_log.csv")
  sample_flow_file <- file.path(fixture_dir, "sample_flow.csv")
  subprocess_dir <- file.path(fixture_dir, "subprocess")
  quiet_log <- function(...) invisible(NULL)

  success_input <- matrix(
    c(0, 0,
      1, 0,
      0, 1,
      1, 1),
    ncol = 2,
    byrow = TRUE
  )
  failure_input <- matrix(letters[seq_len(8)], ncol = 2)

  success <- process_and_monitor(
    expr_matrix = success_input,
    i = 1L,
    DIM = 1L,
    log_message = quiet_log,
    memory_threshold = 0.25,
    max_time_per_iteration = 180,
    poll_interval = 0.05,
    results_file = results_file,
    log_file = progress_file,
    temp_dataset_dir = subprocess_dir,
    sample_id = "fixture-success",
    cohort = "fixture",
    representation = "raw"
  )
  failure <- process_and_monitor(
    expr_matrix = failure_input,
    i = 2L,
    DIM = 1L,
    log_message = quiet_log,
    memory_threshold = 0.25,
    max_time_per_iteration = 180,
    poll_interval = 0.05,
    results_file = results_file,
    log_file = progress_file,
    temp_dataset_dir = subprocess_dir,
    sample_id = "fixture-failure",
    cohort = "fixture",
    representation = "raw"
  )

  attempts <- rbind(success$attempt, failure$attempt)
  append_ph_attempts(attempts, attempt_file)

  flow <- new_sample_flow(
    sample_ids = c("fixture-success", "fixture-failure"),
    cohort = "fixture",
    loaded_features = c(nrow(success_input), nrow(failure_input)),
    loaded_cells = c(ncol(success_input), ncol(failure_input))
  )
  flow <- record_pre_qc_sample_flow(flow, flow$sample_id, c(4L, 4L))
  flow <- record_post_qc_sample_flow(flow, flow$sample_id, c(4L, 4L), min_cells = 2L)
  write_provenance_csv(flow, sample_flow_file)

  expect_silent(assert_sample_reconciliation(
    input_ids = flow$sample_id,
    excluded_ids = character(),
    eligible_ids = flow$sample_id,
    completed_ids = "fixture-success",
    failed_ids = "fixture-failure",
    strict_completion = FALSE
  ))

  written_attempts <- utils::read.csv(attempt_file, stringsAsFactors = FALSE)
  written_flow <- utils::read.csv(sample_flow_file, stringsAsFactors = FALSE)

  expect_false(is.null(success$PD))
  expect_gt(NROW(success$PD), 0L)
  expect_equal(sum(success$PD[, 1] == 0), 3L)
  expect_equal(sum(success$PD[, 1] == 1), 1L)
  expect_equal(success$attempt$status, "completed")
  expect_equal(success$attempt$exit_status, 0L)
  expect_true(success$attempt$pd_written)
  expect_equal(success$attempt$poll_interval_seconds, 0.05)
  expect_gte(success$attempt$memory_samples, 1L)
  expect_true(is.finite(success$attempt$monitor_peak_rss_bytes))
  expect_true(is.finite(success$attempt$child_peak_rss_bytes))
  expect_true(is.finite(success$attempt$process_tree_peak_rss_bytes))
  expect_gte(success$attempt$process_tree_peak_rss_bytes,
             success$attempt$child_peak_rss_bytes)
  expect_gte(success$attempt$process_tree_peak_count, 1L)

  expect_null(failure$PD)
  expect_equal(failure$attempt$status, "ph_child_error")
  expect_gt(failure$attempt$exit_status, 0L)
  expect_false(failure$attempt$pd_written)
  expect_true(nzchar(failure$attempt$error_message))
  expect_gte(failure$attempt$memory_samples, 1L)

  expect_equal(written_attempts$sample_id, c("fixture-success", "fixture-failure"))
  expect_equal(written_attempts$status, c("completed", "ph_child_error"))
  expect_equal(written_attempts$input_rows, c(4L, 4L))
  expect_equal(written_attempts$input_columns, c(2L, 2L))
  expect_equal(written_flow$disposition, c("eligible_for_ph", "eligible_for_ph"))
  expect_true(all(file.exists(c(attempt_file, sample_flow_file, results_file, progress_file))))
})
