test_that("MV8-ZT interruption recovery freezes engine-v2 no-retry adoption", {
  root_v2 <- testthat::test_path(
    "..", "..", "docs", "audits",
    "mv08zt-engine-v2-interruption-recovery-prefreeze-v2"
  )
  root_v1 <- testthat::test_path(
    "..", "..", "docs", "audits",
    "mv08zt-engine-v2-interruption-recovery-prefreeze-v1"
  )
  root <- if (dir.exists(root_v2)) root_v2 else root_v1
  skip_if_not(dir.exists(root), "MV8-ZT recovery prefreeze has not been produced")
  read <- function(name) utils::read.csv(
    file.path(root, name), stringsAsFactors = FALSE, check.names = FALSE
  )
  validation <- read("mv08zt-r1-validation.csv")
  snapshot <- read("mv08zt-r1-stopped-snapshot.csv")
  orphan <- read("mv08zt-r1-orphan-binding.csv")
  policy <- read("mv08zt-r1-recovery-policy.csv")
  decision <- read("mv08zt-r1-decision.csv")
  ledger <- read("mv08zt-r1-prefix-ledger.csv")
  completed <- read("mv08zt-r1-prefix-completions.csv")

  expect_equal(nrow(validation), 28L)
  expect_true(all(validation$passed))
  expect_equal(snapshot$completed_prefix, 325L)
  expect_equal(snapshot$unreceipted_complete_order, 326L)
  expect_equal(snapshot$partial_files, 0L)
  expect_equal(snapshot$active_runner_processes, 0L)
  expect_equal(nrow(ledger), 325L)
  expect_equal(nrow(completed), 325L)
  expect_true(all(ledger$scientific_engine_version == 2L))
  expect_true(all(completed$scientific_engine_version == 2L))
  expect_equal(orphan$production_order, 326L)
  expect_equal(orphan$pair_count, 250L)
  expect_equal(orphan$scientific_engine_version, 2L)
  expect_true(orphan$upper_bound_not_measurement)
  expect_equal(orphan$parent_elapsed_upper_bound_seconds, 3600)
  expect_equal(orphan$parent_rss_upper_bound_bytes, 4 * 1024^3)
  expect_false(policy$landscape_recomputation)
  expect_false(policy$automatic_retry)
  expect_false(policy$deletion_allowed)
  expect_false(policy$overwrite_allowed)
  expect_false(policy$scientific_contract_changed)
  expect_true(decision$orphan_adoption_authorized)
  expect_equal(decision$resume_at_production_order, 327L)
  expect_equal(decision$automatic_retries, 0L)
})

test_that("MV8-ZT recovery implementation cannot recompute or delete", {
  root <- testthat::test_path("..", "..")
  builder <- file.path(root, "scripts", "build_mv08zt_interruption_recovery_prefreeze.R")
  recovery <- file.path(root, "scripts", "recover_mv08zt_landscape_interruption.R")
  expect_silent(parse(builder))
  expect_silent(parse(recovery))
  text <- paste(readLines(recovery, warn = FALSE), collapse = "\n")
  expect_match(text, "accepted-prefix private evidence drift", fixed = TRUE)
  expect_match(text, "orphan scientific validation failed", fixed = TRUE)
  expect_match(text, "conservative upper bounds", fixed = TRUE)
  expect_match(text, "resume_at=327", fixed = TRUE)
  expect_match(text, "public recovery receipts are non-monotone", fixed = TRUE)
  expect_match(text, "MV08ZT_RECOVERY_GIT_HEAD", fixed = TRUE)
  expect_false(grepl("system2\\(\"git\"", text, perl = TRUE))
  expect_false(grepl("run_mv08z_landscape_chunk", text, fixed = TRUE))
  expect_false(grepl("unlink\\(|file.remove\\(", text, perl = TRUE))
})
