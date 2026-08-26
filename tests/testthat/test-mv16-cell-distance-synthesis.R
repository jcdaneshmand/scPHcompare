test_that("MV16 summarization retains every prospectively fixed family", {
  root <- testthat::test_path("..", "..")
  source_root <- file.path(root, "tmp", "mv15-cell-distance-comparison-public-v1")
  skip_if_not(dir.exists(source_root), "private MV15 production is absent")
  global <- read.csv(file.path(source_root, "mv15-global-summary.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)
  neighbor <- read.csv(file.path(source_root, "mv15-neighbor-summary.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
  result <- mv16_build_descriptive_synthesis_v1(global, neighbor)
  expect_equal(nrow(result$complete_global), 36L)
  expect_equal(nrow(result$complete_neighbor), 42L)
  expect_equal(nrow(result$global_summary), 10L)
  expect_equal(nrow(result$neighbor_summary), 16L)
  expect_setequal(result$complete_global$contrast_family, c(
    "cell_seed_stability", "cell_panel_sensitivity",
    "cell_gene_view_agreement"))
  expect_setequal(result$complete_global$homology_dimension, c("H0", "H1"))
})

test_that("MV16 standalone scripts parse", {
  root <- testthat::test_path("..", "..")
  paths <- file.path(root, "scripts", c(
    "build_mv16_descriptive_synthesis_prefreeze.R",
    "run_mv16_descriptive_synthesis.R",
    "build_mv16_descriptive_synthesis_closure.R"
  ))
  expect_true(all(file.exists(paths)))
  expect_silent(lapply(paths, parse))
})
