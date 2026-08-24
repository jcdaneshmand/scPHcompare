test_that("MV8-ZH prospectively freezes no-retry orphan adoption", {
  root <- file.path("..", "..", "docs", "audits",
                    "mv08zh-landscape-launcher-recovery-prefreeze-v1")
  read <- function(name) utils::read.csv(file.path(root, name),
                                          stringsAsFactors = FALSE,
                                          check.names = FALSE)
  validation <- read("mv08zh-validation.csv")
  snapshot <- read("mv08zh-stopped-snapshot.csv")
  orphan <- read("mv08zh-orphan-binding.csv")
  policy <- read("mv08zh-recovery-policy.csv")
  decision <- read("mv08zh-decision.csv")

  expect_equal(nrow(validation), 25L)
  expect_true(all(validation$passed))
  expect_equal(snapshot$completed_prefix, 163L)
  expect_equal(snapshot$unreceipted_complete_order, 164L)
  expect_equal(snapshot$partial_files, 0L)
  expect_equal(orphan$production_order, 164L)
  expect_equal(orphan$pair_count, 250L)
  expect_true(orphan$upper_bound_not_measurement)
  expect_equal(orphan$parent_elapsed_upper_bound_seconds, 3600)
  expect_equal(orphan$parent_rss_upper_bound_bytes, 4 * 1024^3)
  expect_false(policy$landscape_recomputation)
  expect_false(policy$automatic_retry)
  expect_false(policy$deletion_allowed)
  expect_false(policy$overwrite_allowed)
  expect_true(decision$orphan_adoption_authorized)
  expect_equal(decision$resume_at_production_order, 165L)
  expect_false(decision$scientific_contract_changed)
  expect_equal(decision$automatic_retries, 0L)
})

test_that("MV8-ZH recovery implementation is bounded and auditable", {
  root <- testthat::test_path("..", "..")
  recovery <- file.path(root, "scripts",
                        "recover_mv08zh_landscape_launcher_interruption.R")
  builder <- file.path(root, "scripts",
                       "build_mv08zh_landscape_launcher_recovery_prefreeze.R")
  expect_silent(parse(recovery))
  expect_silent(parse(builder))
  text <- paste(readLines(recovery, warn = FALSE), collapse = "\n")
  expect_match(text, "accepted-prefix private evidence drift", fixed = TRUE)
  expect_match(text, "orphan scientific validation failed", fixed = TRUE)
  expect_match(text, "conservative upper bounds", fixed = TRUE)
  expect_match(text, "resume_at=165", fixed = TRUE)
  expect_false(grepl("run_mv08z_landscape_chunk", text, fixed = TRUE))
  expect_false(grepl("unlink\\(|file.remove\\(", text, perl = TRUE))
})
