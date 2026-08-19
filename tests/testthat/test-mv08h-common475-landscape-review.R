test_that("MV8-H common475 landscapes retain all active levels separately", {
  root <- file.path("..", "..", "docs", "audits", "mv08h-common475-landscape-v1")
  execution_path <- file.path(root, "mv08h-common475-landscape-execution.csv")
  levels_path <- file.path(root, "mv08h-common475-landscape-level-summary.csv")
  ledger_path <- file.path(root, "mv08h-common475-landscape-repeat-ledger.csv")
  skip_if_not(file.exists(execution_path))
  skip_if_not(file.exists(levels_path))
  skip_if_not(file.exists(ledger_path))
  execution <- read.csv(execution_path, check.names = FALSE, stringsAsFactors = FALSE)
  levels <- read.csv(levels_path, check.names = FALSE, stringsAsFactors = FALSE)
  ledger <- read.csv(ledger_path, check.names = FALSE, stringsAsFactors = FALSE)

  expect_equal(nrow(execution), 4L)
  expect_setequal(execution$view_id, c("cell_topology_v1", "gene_topology_v1"))
  expect_setequal(execution$homology_dimension, c("H0", "H1"))
  expect_equal(execution$finite_positive_intervals, c(383L, 316L, 474L, 1073L))
  expect_equal(execution$active_level_count, c(383L, 81L, 474L, 287L))
  expect_true(all(execution$breakpoint_count >= 2L))
  expect_true(all(execution$total_integrated_lambda_squared > 0))
  expect_true(all(execution$integration_method == "exact_critical_breakpoint_v1"))
  expect_true(all(execution$level_policy == "all_active_consecutive_levels"))
  expect_true(all(execution$grid_policy == "none"))
  expect_true(all(!execution$labels_outcomes_opened))
  expect_true(all(!execution$fusion_computed))

  expect_equal(nrow(levels), sum(execution$active_level_count))
  expect_true(all(levels$level >= 1L))
  expect_true(all(levels$integrated_lambda_squared >= 0))
  for (i in seq_len(nrow(execution))) {
    rows <- levels[levels$view_id == execution$view_id[[i]] &
                     levels$homology_dimension == execution$homology_dimension[[i]], , drop = FALSE]
    expect_equal(rows$level, seq_len(execution$active_level_count[[i]]))
  }
  forbidden <- c("label", "outcome", "barcode", "expression", "gene_id")
  expect_false(any(tolower(names(levels)) %in% forbidden))
  expect_equal(nrow(ledger), 2L)
  expect_true(all(ledger$byte_identical))
  expect_equal(ledger$primary_sha256, ledger$repeat_sha256)
})
