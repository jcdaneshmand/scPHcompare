test_that("MV10-H/I closes the label-closed clustering review synthesis", {
  root <- testthat::test_path("..", "..")
  synthesis <- file.path(root, "docs", "audits",
                         "mv10h-clustering-review-synthesis-v1")
  closure <- file.path(root, "docs", "audits",
                       "mv10i-clustering-review-closure-v1")
  verify_manifest <- function(path, name) {
    manifest <- read.csv(file.path(path, name), stringsAsFactors = FALSE,
                         check.names = FALSE)
    files <- file.path(path, manifest$artifact)
    expect_true(all(file.exists(files)))
    expect_equal(as.numeric(file.info(files)$size), as.numeric(manifest$bytes))
    expect_equal(unname(vapply(files, function(file) digest::digest(
      file = file, algo = "sha256", serialize = FALSE
    ), character(1L))), manifest$sha256)
  }
  verify_manifest(synthesis, "mv10h-artifact-manifest.csv")
  verify_manifest(closure, "mv10i-artifact-manifest.csv")
  read_synthesis <- function(name) read.csv(
    file.path(synthesis, name), stringsAsFactors = FALSE, check.names = FALSE
  )
  receipt <- read_synthesis("mv10h-terminal-receipt.csv")
  stability <- read_synthesis("mv10h-stability-grid.csv")
  quality <- read_synthesis("mv10h-quality-summary.csv")
  agreement <- read_synthesis("mv10h-method-agreement-summary.csv")
  selection <- read_synthesis("mv10h-primary-selection.csv")
  primary_stability <- read_synthesis("mv10h-primary-stability.csv")
  primary_quality <- read_synthesis("mv10h-primary-quality.csv")
  validation <- read.csv(file.path(closure, "mv10i-validation.csv"),
                         stringsAsFactors = FALSE, check.names = FALSE)
  decision <- read.csv(file.path(closure, "mv10i-decision.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
  rehash <- read.csv(file.path(closure, "mv10i-synthesis-rehash.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)

  expect_equal(receipt$completion_state, "complete")
  expect_equal(receipt$output_tables, 6L)
  expect_identical(
    c(nrow(stability), nrow(quality), nrow(agreement), nrow(selection),
      nrow(primary_stability), nrow(primary_quality)),
    c(270L, 270L, 540L, 2L, 18L, 90L)
  )
  expect_identical(sort(unique(stability$stack_id)), c(
    "allqc_data_exact500", "allqc_residual_exact500",
    "existing_selectedfit_data_exact500"
  ))
  expect_identical(sort(unique(stability$homology_dimension)), c("H0", "H1"))
  expect_equal(length(unique(stability$method_id)), 5L)
  expect_identical(sort(unique(stability$k)), 2:10)
  expect_identical(selection$selected_k, c(2L, 3L))
  expect_equal(nrow(rehash), 6L)
  expect_true(all(rehash$numeric_repeat))
  expect_equal(nrow(validation), 23L)
  expect_true(all(validation$passed))
  expect_true(decision$figure_render_authorized_after_commit)
  expect_equal(decision$interpretation_state,
               "closed_pending_exact_figures_and_visual_review")
  expect_equal(decision$biological_claims_state, "closed")
  expect_equal(decision$manuscript_claims_state, "closed")
  expect_false(receipt$labels_used)
  expect_false(receipt$outcomes_used)
  expect_false(receipt$inference_performed)
  expect_false(receipt$ranking_performed)
  expect_false(receipt$combined_score)
  expect_false(receipt$H0_H1_combined)
  expect_false(receipt$cell_gene_combined)
  expect_false(receipt$biological_claims)
  expect_false(receipt$manuscript_claims)
})
