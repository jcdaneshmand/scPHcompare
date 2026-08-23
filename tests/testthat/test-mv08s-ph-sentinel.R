test_that("MV8-S freezes only the bounded same-axis PH sentinel", {
  root <- file.path("..", "..", "docs", "audits",
                    "mv08s-ph-sentinel-prefreeze-v1")
  skip_if_not(dir.exists(root))
  read <- function(name) utils::read.csv(
    file.path(root, name), check.names = FALSE, stringsAsFactors = FALSE
  )
  contract <- read("mv08s-contract.csv")
  inputs <- read("mv08s-external-input-bindings.csv")
  sources <- read("mv08s-source-cache-bindings.csv")
  queue <- read("mv08s-ph-sentinel-queue.csv")
  cross <- read("mv08s-cross-engine-contract.csv")
  checks <- read("mv08s-validation.csv")
  manifest <- read("mv08s-artifact-manifest.csv")

  expect_equal(nrow(contract), 1L)
  expect_equal(nrow(inputs), 8L)
  expect_equal(nrow(sources), 2L)
  expect_equal(nrow(queue), 23L)
  expect_equal(sum(queue$view_kind == "cell_topology_v1"), 8L)
  expect_equal(sum(queue$view_kind == "gene_topology_v1"), 15L)
  expect_equal(sum(queue$execution_role == "source_produced_gene_ph"), 7L)
  expect_equal(nrow(cross), 4L)
  expect_setequal(unique(cross$mode), c("subset32", "full"))
  expect_true(all(queue$workers == 1L))
  expect_true(all(queue$retries == 0L))
  expect_true(all(queue$authorization_state == "authorized_after_mv08s_commit"))
  expect_equal(contract$monitored_units, 66L)
  expect_identical(contract$execution_head_state,
                   "bind_exact_committed_head_at_launch_and_publish")
  expect_equal(contract$full_ph_jobs_authorized, 0L)
  expect_equal(contract$landscape_groups_authorized, 0L)
  expect_equal(contract$label_jobs_authorized, 0L)
  expect_equal(contract$outcome_jobs_authorized, 0L)
  expect_true(all(checks$passed))
  expect_true(all(file.exists(file.path(root, manifest$artifact))))
  observed <- vapply(file.path(root, manifest$artifact), function(path) {
    digest::digest(file = path, algo = "sha256", serialize = FALSE)
  }, character(1L))
  expect_identical(unname(observed), manifest$sha256)
})

test_that("MV8-S runner retains the fallback and downstream firewalls", {
  runner_path <- file.path("..", "..", "scripts", "run_mv08s_ph_sentinel.R")
  runner <- paste(readLines(runner_path, warn = FALSE),
                  collapse = "\n")
  expect_match(runner, "primary\\$disposition != \\\"rss_cap_exceeded\\\"")
  expect_match(runner, "completed_units != 66L")
  expect_match(runner, "execution_head = expected_head", fixed = TRUE)
  expect_match(runner, "MV08S_RECOVERY_PREFREEZE", fixed = TRUE)
  baseline <- paste(readLines(file.path(
    "..", "..", "scripts", "run_mv08s_same_axis_baseline_entry.R"
  ), warn = FALSE), collapse = "\n")
  expect_match(baseline, 'source("R/toy_baseline.R")', fixed = TRUE)
  expect_match(runner, "full_ph_jobs_authorized = 0L")
  expect_match(runner, "landscape_groups_authorized = 0L")
  expect_match(runner, "label_jobs_authorized = 0L")
  expect_match(runner, "outcome_jobs_authorized = 0L")
  expect_false(grepl("run_landscape", runner, fixed = TRUE))
})

test_that("MV8-SA binds the output-free helper recovery", {
  root <- file.path("..", "..", "docs", "audits",
                    "mv08sa-baseline-helper-recovery-prefreeze-v1")
  read <- function(name) utils::read.csv(
    file.path(root, name), check.names = FALSE, stringsAsFactors = FALSE
  )
  failure <- read("mv08sa-failure-receipt.csv")
  decision <- read("mv08sa-decision.csv")
  checks <- read("mv08sa-validation.csv")
  implementation <- read("mv08sa-implementation-bindings.csv")
  manifest <- read("mv08sa-artifact-manifest.csv")
  expect_equal(nrow(failure), 1L)
  expect_identical(failure$disposition, "failed")
  expect_false(failure$output_published)
  expect_equal(failure$ph_records_computed, 0L)
  expect_false(decision$scientific_contract_changed)
  expect_false(decision$resource_contract_changed)
  expect_true(decision$replacement_roots_required)
  expect_equal(decision$within_replacement_retries, 0L)
  expect_true(all(checks$passed))
  observed_implementation <- vapply(implementation$file, function(path) {
    digest::digest(file = file.path("..", "..", path), algo = "sha256",
                   serialize = FALSE)
  }, character(1L))
  expect_identical(unname(observed_implementation), implementation$sha256)
  observed_manifest <- vapply(file.path(root, manifest$artifact), function(path) {
    digest::digest(file = path, algo = "sha256", serialize = FALSE)
  }, character(1L))
  expect_identical(unname(observed_manifest), manifest$sha256)
})
