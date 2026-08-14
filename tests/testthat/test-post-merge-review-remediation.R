test_that("frozen file identity rejects a same-size substitution", {
  original <- tempfile(fileext = ".bin")
  substituted <- tempfile(fileext = ".bin")
  on.exit(unlink(c(original, substituted)), add = TRUE)
  writeBin(charToRaw("ABCD"), original)
  writeBin(charToRaw("WXYZ"), substituted)
  expected_sha256 <- digest::digest(
    file = original, algo = "sha256", serialize = FALSE
  )
  expected_size <- unname(file.info(original)$size)

  expect_invisible(verify_frozen_file_identity(
    original, expected_size, expected_sha256, "Fixture cache"
  ))
  expect_identical(unname(file.info(substituted)$size), expected_size)
  expect_error(
    verify_frozen_file_identity(
      substituted, expected_size, expected_sha256, "Fixture cache"
    ),
    "SHA-256 does not match"
  )
})

fake_live_process <- function() {
  state <- new.env(parent = emptyenv())
  state$alive <- TRUE
  state$kill_calls <- 0L
  state$is_alive <- function() state$alive
  state$kill_tree <- function() {
    state$kill_calls <- state$kill_calls + 1L
    state$alive <- FALSE
    TRUE
  }
  state
}

test_that("live SCT workers are killed at elapsed and RSS caps", {
  elapsed_process <- fake_live_process()
  elapsed <- mv05d0_enforce_live_process_caps_v1(
    elapsed_process, elapsed_seconds = 11, current_rss_bytes = 100,
    peak_rss_bytes = 90, elapsed_cap_seconds = 10, rss_cap_bytes = 1000
  )
  expect_identical(elapsed$disposition, "elapsed_cap_exceeded")
  expect_true(elapsed$kill_requested)
  expect_false(elapsed_process$alive)
  expect_identical(elapsed_process$kill_calls, 1L)

  rss_process <- fake_live_process()
  rss <- mv05d0_enforce_live_process_caps_v1(
    rss_process, elapsed_seconds = 1, current_rss_bytes = 1001,
    peak_rss_bytes = 900, elapsed_cap_seconds = 10, rss_cap_bytes = 1000
  )
  expect_identical(rss$disposition, "rss_cap_exceeded")
  expect_identical(rss$peak_rss_bytes, 1001)
  expect_true(rss$kill_requested)
  expect_false(rss_process$alive)

  healthy_process <- fake_live_process()
  healthy <- mv05d0_enforce_live_process_caps_v1(
    healthy_process, elapsed_seconds = 1, current_rss_bytes = 100,
    peak_rss_bytes = 90, elapsed_cap_seconds = 10, rss_cap_bytes = 1000
  )
  expect_true(is.na(healthy$disposition))
  expect_false(healthy$kill_requested)
  expect_true(healthy_process$alive)
})

test_that("review remediations are wired into the production scripts", {
  scripts_dir <- file.path("..", "..", "scripts")
  skip_if_not(dir.exists(scripts_dir),
              "repository-only scripts are excluded from source packages")
  cache_job <- readLines(
    file.path(scripts_dir, "run_mv05c_existing_data_job.R"), warn = FALSE
  )
  cache_verification <- grep("verify_frozen_file_identity", cache_job)
  cache_deserialization <- grep("cache <- readRDS", cache_job, fixed = TRUE)
  expect_length(cache_verification, 1L)
  expect_length(cache_deserialization, 1L)
  expect_lt(cache_verification, cache_deserialization)

  queue <- readLines(
    file.path(scripts_dir, "run_mv05d0_sct_cache_queue.R"), warn = FALSE
  )
  expect_true(any(grepl("mv05d0_enforce_live_process_caps_v1", queue,
                        fixed = TRUE)))
})
