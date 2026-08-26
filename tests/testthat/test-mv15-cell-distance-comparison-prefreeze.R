test_that("MV15 prospectively freezes exactly 36 label-closed comparisons", {
  root <- testthat::test_path("..", "..")
  audit <- file.path(root, "docs", "audits",
                     "mv15-cell-distance-comparison-prefreeze-v1")
  read_at <- function(name) read.csv(
    file.path(audit, name), stringsAsFactors = FALSE, check.names = FALSE
  )
  contract <- read_at("mv15-contract.csv")
  bindings <- read_at("mv15-stack-bindings.csv")
  queue <- read_at("mv15-comparison-queue.csv")
  validation <- read_at("mv15-validation.csv")
  decision <- read_at("mv15-decision.csv")
  implementation <- read_at("mv15-implementation-bindings.csv")
  manifest <- read_at("mv15-artifact-manifest.csv")

  expect_equal(nrow(validation), 24L)
  expect_true(all(validation$passed))
  expect_equal(nrow(bindings), 28L)
  expect_equal(as.integer(table(bindings$source_kind)), c(14L, 14L))
  expect_true(all(bindings$source_rehashed))
  expect_true(all(nzchar(bindings$payload_set_sha256)))
  expect_true(all(nzchar(bindings$pair_axis_sha256)))
  expect_equal(nrow(queue), 36L)
  expect_equal(sum(queue$contrast_family == "cell_seed_stability"), 20L)
  expect_equal(sum(queue$contrast_family == "cell_panel_sensitivity"), 2L)
  expect_equal(sum(queue$contrast_family == "cell_gene_view_agreement"), 14L)
  expect_equal(as.integer(table(queue$homology_dimension)), c(18L, 18L))
  expect_true(all(queue$neighbor_k[queue$dataset_scope == "internal124"] ==
                    "10"))
  expect_true(all(queue$neighbor_k[queue$dataset_scope == "external8"] ==
                    "2;3"))
  expect_equal(contract$workers, 1L)
  expect_equal(contract$retries, 0L)
  expect_true(contract$strict_prefix_resume)
  expect_true(contract$independent_recomputation_required)
  expect_false(any(contract[c(
    "labels_authorized", "outcomes_authorized", "clustering_authorized",
    "fusion_authorized", "inference_authorized",
    "biological_claims_authorized", "manuscript_claims_authorized"
  )]))
  expect_true(decision$execution_authorized_after_commit)
  expect_false(decision$values_inspected_during_prefreeze)
  implementation_paths <- file.path(root, implementation$file)
  expect_true(all(file.exists(implementation_paths)))
  expect_equal(unname(vapply(implementation_paths, function(path) {
    digest::digest(file = path, algo = "sha256", serialize = FALSE)
  }, character(1L))), implementation$sha256)

  files <- file.path(audit, manifest$artifact)
  expect_equal(as.numeric(file.info(files)$size), manifest$bytes)
  expect_equal(unname(vapply(files, function(path) digest::digest(
    file = path, algo = "sha256", serialize = FALSE
  ), character(1L))), manifest$sha256)
})
