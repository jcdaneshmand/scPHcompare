test_that("MV8-ZW correction machinery is additive and fail closed", {
  root <- test_path("..", "..")
  files <- file.path(root, c(
    "scripts/run_mv08zw_correction_worker.R",
    "scripts/run_mv08zw_landscape_correction.R",
    "scripts/build_mv08zw_correction_execution_prefreeze.R",
    "scripts/build_mv08zx_correction_closure.R"
  ))
  expect_true(all(file.exists(files)))
  text <- paste(vapply(files, function(path)
    paste(readLines(path, warn = FALSE), collapse = "\n"), character(1L)),
    collapse = "\n")
  expect_match(text, "gene_topology_v1", fixed = TRUE)
  expect_match(text, "all_consecutive_active_levels", fixed = TRUE)
  expect_match(text, "exact_streamed_squared_L2", fixed = TRUE)
  expect_match(text, "rust_scph_landscape_kernel_v2", fixed = TRUE)
  expect_match(text, "strict prefix", ignore.case = TRUE)
  expect_match(text, "MV8-ZU", fixed = TRUE)
  expect_match(text, "comparison_jobs = 0L", fixed = TRUE)
  expect_false(grepl("unlink\\(|file.remove\\(|remove-item|rm -", text,
                     ignore.case = TRUE))
})

test_that("MV8-ZW prefreeze binds exactly two correction groups when present", {
  audit <- test_path("..", "..", "docs", "audits",
                     "mv08zw-correction-execution-prefreeze-v1")
  skip_if_not(dir.exists(audit), "MV8-ZW prefreeze has not been generated")
  read_audit <- function(name) read.csv(file.path(audit, name),
                                        stringsAsFactors = FALSE,
                                        check.names = FALSE)
  contract <- read_audit("mv08zw-contract.csv")
  queue <- read_audit("mv08zw-correction-queue.csv")
  implementation <- read_audit("mv08zw-implementation-bindings.csv")
  validation <- read_audit("mv08zw-validation.csv")
  decision <- read_audit("mv08zw-decision.csv")
  manifest <- read_audit("mv08zw-artifact-manifest.csv")
  expect_equal(contract$groups, 2L)
  expect_equal(contract$total_pairs, 56L)
  expect_equal(contract$input_view_kind, "gene_topology_v1")
  expect_equal(contract$scientific_engine_version, 2L)
  expect_equal(contract$workers, 1L)
  expect_equal(contract$retries, 0L)
  expect_false(contract$mv08zu_mutation_authorized)
  expect_identical(queue$homology_dimension, c("H0", "H1"))
  expect_true(all(queue$unordered_pairs == 28L))
  expect_true(all(queue$input_view_kind == "gene_topology_v1"))
  expect_equal(nrow(implementation), 12L)
  expect_true(all(nchar(implementation$sha256) == 64L))
  expect_true(all(validation$passed))
  expect_equal(nrow(validation), 13L)
  expect_true(decision$correction_authorized_after_commit)
  expect_false(decision$comparisons_authorized)
  expect_false(decision$clustering_authorized)
  expect_false(decision$fusion_authorized)
  expect_false(decision$labels_authorized)
  expect_false(decision$outcomes_authorized)
  expect_false(decision$manuscript_claims_authorized)
  expect_equal(nrow(manifest), 7L)
})
