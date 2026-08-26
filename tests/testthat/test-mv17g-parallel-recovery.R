test_that("MV17-G parallel recovery fixes eight one-thread workers and zero retries", {
  contract <- mv17g_parallel_recovery_contract_v1()
  expect_equal(contract$workers, 8L)
  expect_equal(contract$threads_per_child, 1L)
  expect_equal(contract$retries, 0L)
  expect_equal(contract$concurrent_child_RSS_cap_bytes, 8 * contract$child_RSS_cap_bytes)
  expect_equal(unname(mv17g_parallel_thread_environment_v1()), rep("1", 6L))
})

test_that("MV17-G checkpoint admission requires a complete consecutive prefix", {
  queue <- data.frame(
    job_order = 1:3,
    view = "cell",
    unit_order = 1:3,
    null_family = "observed",
    stringsAsFactors = FALSE
  )
  root <- tempfile("mv17g-checkpoint-")
  dir.create(file.path(root, "jobs"), recursive = TRUE)
  dir.create(file.path(root, "logs"), recursive = TRUE)
  for (path in mv17g_job_artifacts_v1(queue[1, ], root)) writeBin(charToRaw("x"), path)
  scan <- mv17g_checkpoint_scan_v1(queue, root)
  expect_equal(scan$state, c("complete", "absent", "absent"))
  expect_equal(mv17g_complete_prefix_v1(scan), 1L)
  expect_equal(unname(lengths(mv17g_parallel_batches_v1(1:19))), c(8L, 8L, 3L))
  writeBin(charToRaw("partial"), mv17g_job_artifacts_v1(queue[2, ], root)[["stdout"]])
  expect_error(mv17g_complete_prefix_v1(mv17g_checkpoint_scan_v1(queue, root)), "partial")
})

test_that("MV17-G parallel wave promotes distinct atomic synthetic children", {
  skip_on_os("windows")
  root <- normalizePath(testthat::test_path("..", ".."))
  private <- tempfile("mv17g-wave-")
  dir.create(file.path(private, "jobs"), recursive = TRUE)
  dir.create(file.path(private, "logs"), recursive = TRUE)
  dir.create(file.path(private, "matrices"), recursive = TRUE)
  set.seed(1717)
  saveRDS(matrix(rnorm(36), 12, 3), file.path(private, "matrices", "cell__001.rds"), version = 3)
  saveRDS(matrix(rnorm(36), 12, 3), file.path(private, "matrices", "cell__002.rds"), version = 3)
  queue <- data.frame(
    job_order = 1:2,
    view = "cell",
    unit_order = 1:2,
    null_family = "coordinate_permutation",
    seed_first = c(501L, 601L),
    replicate_count = 2L,
    stringsAsFactors = FALSE
  )
  worker <- file.path(root, "scripts", "run_mv17g_calibration_group_worker.R")
  wave <- withr::with_dir(
    root,
    mv17g_run_parallel_wave_v1(
      queue,
      private_root = private,
      matrix_root = private,
      worker_path = worker,
      contract = mv17g_parallel_recovery_contract_v1()
    )
  )
  expect_equal(wave$job_order, 1:2)
  expect_true(all(wave$exit_status == 0L))
  expect_true(all(wave$artifacts == 4L))
  scan <- mv17g_checkpoint_scan_v1(queue, private)
  expect_equal(scan$state, rep("complete", 2L))
  expect_equal(mv17g_complete_prefix_v1(scan, require_incomplete = FALSE), 2L)
  results <- lapply(seq_len(nrow(queue)), function(i) readRDS(mv17g_job_artifacts_v1(queue[i, ], private)[["result"]]))
  expect_true(all(vapply(results, function(x) x$finite && nrow(x$metrics) == 16L, logical(1L))))
  streams <- list.files(file.path(private, "logs"), pattern = "[.](stdout|stderr)[.]txt$", full.names = TRUE)
  expect_equal(length(streams), 4L)
  expect_true(all(file.info(streams)$size == 0))
})

test_that("MV17-G recovery scripts parse and bind the parallel closure path", {
  root <- testthat::test_path("..", "..")
  scripts <- file.path(root, "scripts", c(
    "build_mv17g_parallel_recovery_prefreeze.R",
    "run_mv17g_calibration_parallel_recovery.R",
    "build_mv17g_parallel_calibration_closure.R"
  ))
  expect_true(all(file.exists(scripts)))
  for (path in scripts) expect_silent(parse(path))
})

test_that("MV17-G eight-worker recovery prefreeze is exact and aggregate-only", {
  root <- testthat::test_path("..", "..")
  audit <- file.path(root, "docs", "audits", "mv17g-eight-worker-recovery-prefreeze-v1")
  skip_if_not(dir.exists(audit), "MV17-G eight-worker recovery prefreeze absent")
  read <- function(name) read.csv(file.path(audit, name), stringsAsFactors = FALSE, check.names = FALSE)
  validation <- read("mv17g-parallel-validation.csv")
  contract <- read("mv17g-parallel-contract.csv")
  summary <- read("mv17g-parallel-prefix-summary.csv")
  manifest <- read("mv17g-parallel-artifact-manifest.csv")
  expect_equal(nrow(validation), 24L)
  expect_true(all(validation$passed))
  expect_equal(contract$workers, 8L)
  expect_equal(contract$threads_per_child, 1L)
  expect_equal(contract$retries, 0L)
  expect_gt(summary$serial_prefix_children, 0L)
  expect_lt(summary$serial_prefix_children, summary$serial_prefix_children + summary$pending_children)
  expect_false(any(c("unit_id", "identity_token", "source_path", "artifact") %in% names(summary)))
  files <- file.path(audit, manifest$artifact)
  expect_equal(as.numeric(file.info(files)$size), manifest$bytes)
  expect_equal(
    unname(vapply(files, digest::digest, character(1L), algo = "sha256", file = TRUE, serialize = FALSE)),
    unname(manifest$sha256)
  )
})
