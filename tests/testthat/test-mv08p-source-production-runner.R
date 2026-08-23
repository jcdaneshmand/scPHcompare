test_that("MV8-P runner and MV8-Q closure preserve the frozen source-only contract", {
  root <- file.path("..", "..")
  worker_path <- file.path(root, "scripts", "run_mv08o_residual_source_worker.R")
  runner_path <- file.path(root, "scripts", "run_mv08p_full_source_production.R")
  closure_path <- file.path(root, "scripts", "build_mv08q_full_source_production_closure.R")
  expect_true(all(file.exists(c(worker_path, runner_path, closure_path))))
  expect_silent(parse(worker_path))
  expect_silent(parse(runner_path))
  expect_silent(parse(closure_path))

  worker <- paste(readLines(worker_path, warn = FALSE), collapse = "\n")
  runner <- paste(readLines(runner_path, warn = FALSE), collapse = "\n")
  closure <- paste(readLines(closure_path, warn = FALSE), collapse = "\n")
  expect_match(worker, "authorized_after_mv08p_commit", fixed = TRUE)
  expect_match(worker, "source hash drift", fixed = TRUE)
  expect_match(worker, "freeze_in_source_preflight", fixed = TRUE)
  expect_match(runner, "serial_one_worker_ascending_fit_cells_no_retry", fixed = TRUE)
  expect_match(runner, "cleanup_tree = TRUE", fixed = TRUE)
  expect_match(runner, "process$kill_tree()", fixed = TRUE)
  expect_match(runner, "--resume", fixed = TRUE)
  expect_match(runner, "refusing to overwrite partial MV8-P job artifacts", fixed = TRUE)
  expect_match(closure, "132/132 sources covered", fixed = TRUE)
  expect_match(closure, "topology remains closed", fixed = TRUE)

  forbidden_calls <- "ripserr::|ripsDiag\\s*\\(|landscape_distance|cluster::pam\\s*\\(|mclust::"
  expect_false(grepl(forbidden_calls, runner, perl = TRUE))
  expect_false(grepl(forbidden_calls, closure, perl = TRUE))

  queue <- read.csv(file.path(root, "docs", "audits",
    "mv08p-full-source-production-prefreeze-v2", "mv08p-remaining-source-queue.csv"),
    check.names = FALSE, stringsAsFactors = FALSE)
  expect_equal(nrow(queue), 129L)
  expect_identical(as.integer(queue$job_order), seq_len(129L))
  expect_true(all(diff(queue$fit_cells) >= 0L))
  expect_true(all(queue$authorization_state == "authorized_after_mv08p_commit"))
  expect_true(all(queue$workers == 1L & queue$retries == 0L))
  expect_false(any(queue$ph_authorized | queue$landscapes_authorized |
    queue$comparisons_authorized | queue$clustering_authorized |
    queue$fusion_authorized | queue$labels_authorized | queue$outcomes_authorized))
})
