test_that("MV8-ZT freezes fresh engine-v2 production from order zero", {
  root_v2 <- testthat::test_path(
    "..", "..", "docs", "audits",
    "mv08zt-engine-v2-full-landscape-prefreeze-v2"
  )
  root_v3 <- testthat::test_path(
    "..", "..", "docs", "audits",
    "mv08zt-engine-v2-full-landscape-prefreeze-v3"
  )
  root_v1 <- testthat::test_path(
    "..", "..", "docs", "audits",
    "mv08zt-engine-v2-full-landscape-prefreeze-v1"
  )
  root <- if (dir.exists(root_v3)) root_v3 else
    if (dir.exists(root_v2)) root_v2 else root_v1
  skip_if_not(dir.exists(root), "MV8-ZT prefreeze has not been produced")
  read <- function(name) utils::read.csv(
    file.path(root, name), check.names = FALSE, stringsAsFactors = FALSE
  )
  contract <- read("mv08zt-contract.csv")
  queue <- read("mv08zt-production-queue.csv")
  checks <- read("mv08zt-validation.csv")
  firewall <- read("mv08zt-old-root-firewall.csv")
  closure <- read("mv08zu-prospective-closure.csv")
  decision <- read("mv08zt-decision.csv")
  manifest <- read("mv08zt-artifact-manifest.csv")

  expect_equal(nrow(checks), 36L)
  expect_true(all(checks$passed))
  expect_equal(contract$scientific_engine_version, 2L)
  expect_equal(contract$production_start_order, 1L)
  expect_equal(contract$old_completed_chunks_reused, 0L)
  expect_equal(contract$old_completed_pairs_reused, 0L)
  expect_identical(contract$level_policy, "all_consecutive_active_levels")
  expect_identical(contract$integration_policy,
                   "exact_streamed_squared_l2_on_dimension_support")
  expect_identical(contract$dimension_policy,
                   "h0_h1_separate_primary_outputs")
  expect_identical(contract$grid_policy, "no_universal_fixed_grid")
  expect_identical(contract$level_cap_policy, "no_universal_level_cap")
  expect_equal(nrow(queue), 628L)
  expect_identical(as.integer(queue$production_order), 1:628)
  expect_equal(sum(queue$pair_count), 152744L)
  expect_true(all(queue$scientific_engine_version == 2L))
  expect_true(all(queue$workers == 1L))
  expect_true(all(queue$retries == 0L))
  expect_false(firewall$resume)
  expect_false(firewall$reuse)
  expect_false(firewall$overwrite)
  expect_false(firewall$delete)
  expect_equal(closure$required_chunks, 628L)
  expect_equal(closure$required_pairs, 152744L)
  expect_equal(closure$required_engine_version, 2L)
  expect_true(decision$fresh_production_authorized_after_commit)
  expect_false(decision$old_root_reuse)
  expect_equal(decision$comparison_jobs_authorized, 0L)
  expect_equal(decision$clustering_jobs_authorized, 0L)
  expect_equal(decision$fusion_jobs_authorized, 0L)
  expect_equal(decision$label_jobs_authorized, 0L)
  expect_equal(decision$outcome_jobs_authorized, 0L)
  expect_identical(decision$outcome_label_state, "closed")
  expect_false(decision$biological_outcomes_computed)

  observed <- vapply(file.path(root, manifest$artifact), function(path) {
    digest::digest(file = path, algo = "sha256", serialize = FALSE)
  }, character(1L))
  expect_identical(unname(observed), manifest$sha256)
})

test_that("MV8-ZT version-aware execution cannot mix or reuse old outputs", {
  worker <- testthat::test_path(
    "..", "..", "scripts", "run_mv08z_landscape_chunk.R"
  )
  runner <- testthat::test_path(
    "..", "..", "scripts", "run_mv08zt_full_landscape_production.R"
  )
  builder <- testthat::test_path(
    "..", "..", "scripts", "build_mv08zt_engine_v2_full_landscape_prefreeze.R"
  )
  expect_silent(parse(worker))
  expect_silent(parse(runner))
  expect_silent(parse(builder))
  worker_text <- paste(readLines(worker, warn = FALSE), collapse = "\n")
  runner_text <- paste(readLines(runner, warn = FALSE), collapse = "\n")
  builder_text <- paste(readLines(builder, warn = FALSE), collapse = "\n")
  expect_match(worker_text, "as.integer(value$engine_version) != expected_engine_version",
               fixed = TRUE)
  expect_match(worker_text, "scientific_engine_version = expected_engine_version",
               fixed = TRUE)
  expect_match(worker_text, "version_aware_chunk_worker", fixed = TRUE)
  expect_match(worker_text, "sentinel_R_oracle_runner", fixed = TRUE)
  expect_match(runner_text, "status$scientific_engine_version == 2L", fixed = TRUE)
  expect_match(runner_text, "MV08ZT_PREFREEZE = zt_root", fixed = TRUE)
  expect_false(grepl("recovery_roots|MV08Z_RECOVERY_CHAIN", runner_text))
  expect_match(builder_text, "old_completed_chunks_reused = 0L", fixed = TRUE)
  expect_false(grepl("system2|processx::process|process$new", builder_text))
})
