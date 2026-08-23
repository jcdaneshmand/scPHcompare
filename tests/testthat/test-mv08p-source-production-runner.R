test_that("MV8-P runner and MV8-Q closure preserve the frozen source-only contract", {
  root <- file.path("..", "..")
  worker_path <- file.path(root, "scripts", "run_mv08o_residual_source_worker.R")
  runner_path <- file.path(root, "scripts", "run_mv08p_full_source_production.R")
  recovery_path <- file.path(root, "scripts", "run_mv08pr_source_production_recovery.R")
  second_recovery_path <- file.path(root, "scripts", "run_mv08ps_source_production_recovery.R")
  closure_path <- file.path(root, "scripts", "build_mv08q_full_source_production_closure.R")
  expect_true(all(file.exists(c(worker_path, runner_path, recovery_path,
    second_recovery_path, closure_path))))
  expect_silent(parse(worker_path))
  expect_silent(parse(runner_path))
  expect_silent(parse(recovery_path))
  expect_silent(parse(second_recovery_path))
  expect_silent(parse(closure_path))

  worker <- paste(readLines(worker_path, warn = FALSE), collapse = "\n")
  runner <- paste(readLines(runner_path, warn = FALSE), collapse = "\n")
  recovery <- paste(readLines(recovery_path, warn = FALSE), collapse = "\n")
  second_recovery <- paste(readLines(second_recovery_path, warn = FALSE), collapse = "\n")
  closure <- paste(readLines(closure_path, warn = FALSE), collapse = "\n")
  expect_match(worker, "authorized_after_mv08p_commit", fixed = TRUE)
  expect_match(worker, "source hash drift", fixed = TRUE)
  expect_match(worker, "freeze_in_source_preflight", fixed = TRUE)
  expect_match(worker, "MV08_FUTURE_GLOBALS_MAX_SIZE_BYTES", fixed = TRUE)
  expect_match(worker, "authorized bounded 2 GiB", fixed = TRUE)
  expect_match(runner, "serial_one_worker_ascending_fit_cells_no_retry", fixed = TRUE)
  expect_match(runner, "cleanup_tree = TRUE", fixed = TRUE)
  expect_match(runner, "process$kill_tree()", fixed = TRUE)
  expect_match(runner, "--resume", fixed = TRUE)
  expect_match(runner, "refusing to overwrite partial MV8-P job artifacts", fixed = TRUE)
  expect_match(recovery, "MV8-PR original stopped evidence drift", fixed = TRUE)
  expect_match(recovery, "unname(vapply(evidence_paths", fixed = TRUE)
  expect_match(recovery, "accepted-prefix private evidence drift", fixed = TRUE)
  expect_match(recovery, "future_globals_max_size_bytes = 2 * 1024^3", fixed = TRUE)
  expect_match(recovery, "attempt_number = ifelse(i == 124L, 2L, 1L)", fixed = TRUE)
  expect_match(recovery, "refusing to overwrite partial MV8-PR job artifacts", fixed = TRUE)
  expect_match(second_recovery, "verify_evidence(prior_evidence_paths", fixed = TRUE)
  expect_match(second_recovery, '"MV8-PS prior overlay"', fixed = TRUE)
  expect_match(second_recovery, "ps_contract$rss_cap_bytes != 14 * 1024^3", fixed = TRUE)
  expect_match(second_recovery, "attempt_number = ifelse(i == 125L, 2L, 1L)", fixed = TRUE)
  expect_match(second_recovery, "refusing to overwrite partial MV8-PS job artifacts", fixed = TRUE)
  expect_match(closure, "132/132 sources covered", fixed = TRUE)
  expect_match(closure, "stopped_run_evidence_preserved", fixed = TRUE)
  expect_match(closure, "unname(vapply(stopped_paths", fixed = TRUE)
  expect_match(closure, "single_explicit_retry", fixed = TRUE)
  expect_match(closure, "second_stop_evidence_preserved", fixed = TRUE)
  expect_match(closure, "explicit_retry_contract", fixed = TRUE)
  expect_match(closure, "topology remains closed", fixed = TRUE)

  forbidden_calls <- "ripserr::|ripsDiag\\s*\\(|landscape_distance|cluster::pam\\s*\\(|mclust::"
  expect_false(grepl(forbidden_calls, runner, perl = TRUE))
  expect_false(grepl(forbidden_calls, recovery, perl = TRUE))
  expect_false(grepl(forbidden_calls, second_recovery, perl = TRUE))
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

test_that("MV8-PS prefreeze binds the second stop and 14-GiB overlay", {
  root <- file.path("..", "..")
  audit <- file.path(root, "docs", "audits", "mv08ps-source-production-recovery-prefreeze-v1")
  contract <- read.csv(file.path(audit, "mv08ps-contract.csv"),
    check.names = FALSE, stringsAsFactors = FALSE)
  evidence <- read.csv(file.path(audit, "mv08ps-prior-evidence.csv"),
    check.names = FALSE, stringsAsFactors = FALSE)
  manifest <- read.csv(file.path(audit, "mv08ps-artifact-manifest.csv"),
    check.names = FALSE, stringsAsFactors = FALSE)

  expect_equal(nrow(contract), 1L)
  expect_equal(contract$prior_completed_jobs, 124L)
  expect_equal(contract$failed_job_order, 125L)
  expect_equal(contract$recovery_first_job_order, 125L)
  expect_equal(contract$recovery_last_job_order, 129L)
  expect_equal(contract$recovery_jobs, 5L)
  expect_equal(contract$explicit_retry_jobs, 1L)
  expect_equal(contract$future_globals_max_size_bytes, 2 * 1024^3)
  expect_equal(contract$workers, 1L)
  expect_equal(contract$elapsed_cap_seconds, 1800L)
  expect_equal(contract$rss_cap_bytes, 14 * 1024^3)
  expect_equal(contract$prior_evidence_state, "preserve_immutable")
  expect_equal(contract$recovery_state, "authorized_after_commit")
  expect_equal(contract$topology_execution_state, "closed")
  expect_equal(contract$outcome_label_state, "closed")
  expect_false(contract$biological_outcomes_computed)

  expect_equal(nrow(evidence), 6L)
  expect_setequal(evidence$artifact_id, c("mv08pr_stopped_resource_ledger",
    "mv08pr_stopped_progress", "failed_job125_stderr", "failed_job125_stdout",
    "accepted_job124_cache", "accepted_job124_audit"))
  expect_true(all(grepl("^[0-9a-f]{64}$", evidence$sha256)))
  expect_equal(evidence$bytes[evidence$artifact_id == "failed_job125_stdout"], 0)
  expect_equal(nrow(manifest), 3L)
  expect_setequal(manifest$artifact, c(
    "MV08PS_SOURCE_PRODUCTION_RECOVERY_PREFREEZE_2026-08-23.md",
    "mv08ps-contract.csv", "mv08ps-prior-evidence.csv"))
  for (i in seq_len(nrow(manifest))) {
    path <- file.path(audit, manifest$artifact[[i]])
    expect_equal(unname(file.info(path)$size), manifest$bytes[[i]])
    expect_equal(digest::digest(file = path, algo = "sha256", serialize = FALSE),
      manifest$sha256[[i]])
  }
})

test_that("MV8-PR prefreeze binds the preserved stop and bounded overlay", {
  root <- file.path("..", "..")
  audit <- file.path(root, "docs", "audits", "mv08pr-source-production-recovery-prefreeze-v1")
  contract <- read.csv(file.path(audit, "mv08pr-contract.csv"),
    check.names = FALSE, stringsAsFactors = FALSE)
  evidence <- read.csv(file.path(audit, "mv08pr-original-evidence.csv"),
    check.names = FALSE, stringsAsFactors = FALSE)
  manifest <- read.csv(file.path(audit, "mv08pr-artifact-manifest.csv"),
    check.names = FALSE, stringsAsFactors = FALSE)

  expect_equal(nrow(contract), 1L)
  expect_equal(contract$original_completed_jobs, 123L)
  expect_equal(contract$failed_job_order, 124L)
  expect_equal(contract$recovery_first_job_order, 124L)
  expect_equal(contract$recovery_last_job_order, 129L)
  expect_equal(contract$recovery_jobs, 6L)
  expect_equal(contract$explicit_retry_jobs, 1L)
  expect_equal(contract$future_globals_max_size_bytes, 2 * 1024^3)
  expect_equal(contract$workers, 1L)
  expect_equal(contract$elapsed_cap_seconds, 1800L)
  expect_equal(contract$rss_cap_bytes, 12 * 1024^3)
  expect_equal(contract$original_evidence_state, "preserve_immutable")
  expect_equal(contract$recovery_state, "authorized_after_commit")
  expect_equal(contract$topology_execution_state, "closed")
  expect_equal(contract$outcome_label_state, "closed")
  expect_false(contract$biological_outcomes_computed)

  expect_equal(nrow(evidence), 4L)
  expect_setequal(evidence$artifact_id, c("stopped_resource_ledger", "stopped_progress",
    "failed_job_stderr", "failed_job_stdout"))
  expect_true(all(grepl("^[0-9a-f]{64}$", evidence$sha256)))
  expect_equal(evidence$bytes[evidence$artifact_id == "failed_job_stdout"], 0)
  expect_equal(nrow(manifest), 3L)
  expect_setequal(manifest$artifact, c(
    "MV08PR_SOURCE_PRODUCTION_RECOVERY_PREFREEZE_2026-08-23.md",
    "mv08pr-contract.csv", "mv08pr-original-evidence.csv"))
  for (i in seq_len(nrow(manifest))) {
    path <- file.path(audit, manifest$artifact[[i]])
    expect_equal(unname(file.info(path)$size), manifest$bytes[[i]])
    expect_equal(digest::digest(file = path, algo = "sha256", serialize = FALSE),
      manifest$sha256[[i]])
  }
})
