test_that("MV8-P freezes bounded serial source production without opening topology", {
  audit <- file.path("..", "..", "docs", "audits", "mv08p-full-source-production-prefreeze-v2")
  report_path <- file.path(audit, "MV08P_FULL_SOURCE_PRODUCTION_PREFREEZE_2026-08-21.md")
  builder_path <- file.path("..", "..", "scripts", "build_mv08p_full_source_production_prefreeze.R")
  expect_true(all(file.exists(c(report_path, builder_path))))

  read_audit <- function(name) read.csv(file.path(audit, name), check.names = FALSE, stringsAsFactors = FALSE)
  contract <- read_audit("mv08p-contract.csv")
  queue <- read_audit("mv08p-remaining-source-queue.csv")
  validation <- read_audit("mv08p-validation.csv")
  manifest <- read_audit("mv08p-artifact-manifest.csv")
  expect_equal(nrow(contract), 1L)
  expect_equal(unname(as.integer(contract[1, c("all_source_fits", "completed_sentinel_primary_fits",
    "remaining_source_fits", "remaining_internal_fits", "remaining_external_fits",
    "maximum_sentinel_fit_cells", "maximum_remaining_fit_cells", "workers", "retries")])),
    c(132L, 3L, 129L, 122L, 7L, 11475L, 9071L, 1L, 0L))
  expect_true(contract$observed_memory_margin_bytes > 0)
  expect_identical(contract$source_execution_state, "authorized_after_commit")
  expect_identical(contract$topology_execution_state, "closed")

  expect_equal(nrow(queue), 129L)
  expect_identical(as.integer(queue$job_order), seq_len(129L))
  expect_true(all(diff(queue$fit_cells) >= 0L))
  expect_true(all(queue$authorization_state == "authorized_after_mv08p_commit"))
  expect_true(all(queue$workers == 1L & queue$retries == 0L))
  expect_false(any(queue$ph_authorized | queue$landscapes_authorized | queue$comparisons_authorized |
    queue$clustering_authorized | queue$fusion_authorized | queue$labels_authorized | queue$outcomes_authorized))
  expect_true(all(queue$outcome_label_state == "closed" & !queue$biological_outcomes_computed))
  expect_equal(nrow(validation), 8L)
  expect_true(all(validation$passed))
  expect_equal(nrow(manifest), 4L)
  live_hashes <- vapply(file.path(audit, manifest$artifact), function(path) {
    digest::digest(file = path, algo = "sha256", serialize = FALSE)
  }, character(1L))
  expect_identical(unname(live_hashes), manifest$sha256)
  public_text <- paste(unlist(lapply(list.files(audit, pattern = "\\.csv$", full.names = TRUE), readLines, warn = FALSE)), collapse = "\n")
  expect_false(grepl("(?:^|[,\"'])tmp[/\\\\]|/mnt/[a-z]/|[A-Za-z]:[/\\\\]", public_text, perl = TRUE))
})
