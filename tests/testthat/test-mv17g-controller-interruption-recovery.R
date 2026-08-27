test_that("MV17-G controller-interruption recovery implementation parses", {
  root <- testthat::test_path("..", "..")
  scripts <- file.path(root, "scripts", c(
    "build_mv17g_controller_interruption_recovery_prefreeze.R",
    "run_mv17g_calibration_parallel_recovery.R",
    "build_mv17g_parallel_calibration_closure.R"
  ))
  expect_true(all(file.exists(scripts)))
  for (path in scripts) expect_silent(parse(path))
  recovery <- paste(readLines(scripts[[1]], warn = FALSE), collapse = "\n")
  runner <- paste(readLines(scripts[[2]], warn = FALSE), collapse = "\n")
  expect_match(recovery, "prefix != 794L", fixed = TRUE)
  expect_match(recovery, "ledger_prefix != 786L", fixed = TRUE)
  expect_match(recovery, "resume_at_job_order = prefix + 1L", fixed = TRUE)
  expect_match(runner, "pending <- scan$job_order[scan$state == \"absent\"]", fixed = TRUE)
})

test_that("MV17-G controller-interruption recovery prefreeze is exact and aggregate-only", {
  root <- testthat::test_path("..", "..")
  audit <- file.path(root, "docs", "audits", "mv17g-controller-interruption-recovery-prefreeze-v1")
  skip_if_not(dir.exists(audit), "MV17-G controller-interruption recovery prefreeze absent")
  read <- function(name) read.csv(file.path(audit, name), stringsAsFactors = FALSE, check.names = FALSE)
  validation <- read("mv17g-controller-recovery-validation.csv")
  contract <- read("mv17g-controller-recovery-contract.csv")
  state <- read("mv17g-controller-recovery-state-summary.csv")
  decision <- read("mv17g-controller-recovery-decision.csv")
  manifest <- read("mv17g-controller-recovery-artifact-manifest.csv")
  expect_equal(nrow(validation), 30L)
  expect_true(all(validation$passed))
  expect_equal(contract$completed_prefix_children, 794L)
  expect_equal(contract$durable_progress_rows, 786L)
  expect_equal(contract$adopted_unledgered_suffix_children, 8L)
  expect_equal(contract$resume_at_job_order, 795L)
  expect_equal(contract$workers, 8L)
  expect_equal(contract$threads_per_child, 1L)
  expect_equal(contract$retries, 0L)
  expect_equal(state$controller_processes, 0L)
  expect_equal(state$attempt1_GNU_time_bytes, 0)
  expect_false(decision$completed_prefix_recomputation)
  expect_false(decision$automatic_retry)
  expect_false(decision$scientific_contract_changed)
  expect_false(any(c("unit_id", "unit_order", "source_path", "artifact") %in% names(state)))
  files <- file.path(audit, manifest$artifact)
  expect_equal(as.numeric(file.info(files)$size), manifest$bytes)
  expect_equal(
    unname(vapply(files, digest::digest, character(1L), algo = "sha256", file = TRUE, serialize = FALSE)),
    unname(manifest$sha256)
  )
})
