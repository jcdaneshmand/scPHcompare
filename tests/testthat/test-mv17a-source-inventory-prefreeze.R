test_that("MV17-A committed audit is aggregate-only and fully gated", {
  root <- testthat::test_path("..", "..")
  audit <- file.path(root, "docs", "audits",
                     "mv17a-source-inventory-estimand-prefreeze-v4")
  skip_if_not(dir.exists(audit), "MV17-A audit has not been built yet")
  read_at <- function(name) read.csv(file.path(audit, name),
    stringsAsFactors = FALSE, check.names = FALSE)
  contract <- read_at("mv17a-contract.csv")
  inventory <- read_at("mv17a-source-inventory.csv")
  axes <- read_at("mv17a-axis-inventory.csv")
  validation <- read_at("mv17a-validation.csv")
  decision <- read_at("mv17a-decision.csv")
  manifest <- read_at("mv17a-artifact-manifest.csv")
  expect_equal(nrow(validation), 20L)
  expect_true(all(validation$passed))
  expect_equal(inventory$records[inventory$source_family == "gene_PH"], 2528L)
  expect_equal(inventory$records[inventory$source_family == "gene_landscape"],
               152688L)
  expect_equal(inventory$axes[inventory$source_family == "gene_landscape"],
               26L)
  expect_equal(inventory$records[inventory$source_family == "cell_PH"], 1272L)
  expect_equal(inventory$records[inventory$source_family == "cell_landscape"],
               76372L)
  expect_equal(axes$PCA_models[axes$source_family == "cell"], 7L)
  expect_true(all(axes$H0_H1_separate))
  expect_true(all(axes$essential_H0_excluded))
  expect_false(any(contract[c(
    "values_published", "null_computation_authorized",
    "localization_computation_authorized", "H2_computation_authorized",
    "labels_authorized", "outcomes_authorized", "clustering_authorized",
    "fusion_authorized", "biology_authorized", "manuscript_claims_authorized"
  )]))
  expect_true(decision$MV17A_closed_after_commit)
  expect_equal(decision$next_scientific_execution,
               "none_without_distinct_exact_head_prefreeze")
  expect_false(any(grepl("unit_id|sample_id|private_cache_path|output_file",
                         names(inventory))))
  files <- file.path(audit, manifest$artifact)
  expect_equal(as.numeric(file.info(files)$size), manifest$bytes)
  expect_equal(unname(vapply(files, function(path) digest::digest(
    file = path, algo = "sha256", serialize = FALSE
  ), character(1L))), manifest$sha256)
})
