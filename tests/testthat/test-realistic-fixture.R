test_that("realistic fixture inputs are deterministic and satisfy fixed QC gates", {
  first_dir <- tempfile("realistic-input-a-")
  second_dir <- tempfile("realistic-input-b-")
  set.seed(909L)
  rng_before <- .Random.seed
  first <- scPHcompare:::.generate_realistic_fixture_inputs(first_dir, 20260805L)
  second <- scPHcompare:::.generate_realistic_fixture_inputs(second_dir, 20260805L)

  expect_identical(.Random.seed, rng_before)
  expect_identical(first$matrices, second$matrices)
  expect_equal(length(first$matrices), 2L)
  expect_equal(unname(vapply(first$matrices, nrow, integer(1))), c(520L, 520L))
  expect_equal(unname(vapply(first$matrices, ncol, integer(1))), c(20L, 20L))
  expect_true(all(vapply(first$matrices, function(x) {
    all(Matrix::colSums(x > 0) >= 500L)
  }, logical(1))))
  expect_true(all(vapply(first$matrices, function(x) {
    all(Matrix::rowSums(x > 0) > 3L)
  }, logical(1))))
  expect_true(all(vapply(first$matrices, function(x) {
    ribo <- grep("^RP[SL]", rownames(x))
    all(Matrix::colSums(x[ribo, , drop = FALSE]) / Matrix::colSums(x) > 0.05)
  }, logical(1))))
  expect_true(all(vapply(first$matrices, function(x) {
    mito <- grep("^MT-", rownames(x))
    all(Matrix::colSums(x[mito, , drop = FALSE]) / Matrix::colSums(x) <= 0.20)
  }, logical(1))))
})

test_that("realistic fixture reference contract is packaged", {
  reference_path <- system.file(
    "extdata", "realistic_fixture_reference.csv", package = "scPHcompare"
  )
  expect_true(nzchar(reference_path))
  reference <- utils::read.csv(reference_path, check.names = FALSE)
  expect_equal(reference$ph_max_dimension, 0:1)
  expect_true(all(reference$samples_loaded == 2L))
  expect_true(all(reference$samples_eligible == 2L))
  expect_true(all(reference$loaded_features == "520;520"))
  expect_true(all(reference$post_qc_cells == "20;20"))
  expect_true(all(reference$ph_completed == 2L))
  expect_true(all(reference$finite_fraction_minimum <= 1))
  expect_true(all(reference$numeric_tolerance > 0))
  expect_equal(reference$h0_feature_counts, c("519;519", "519;519"))
  expect_equal(reference$h1_feature_counts, c("not_requested", "922;856"))
  expect_identical(eval(formals(run_realistic_fixture)$max_dimension), 0L)
  expect_equal(eval(formals(run_realistic_fixture)$ph_poll_interval), 0.25)
})

test_that("focused PH representation selection is explicit and validated", {
  expected_defaults <- c(
    "sct_individual", "raw", "sct_whole",
    "seurat_integration", "harmony_integration"
  )
  expect_identical(
    eval(formals(process_datasets_PH)$ph_representations),
    expected_defaults
  )
  expect_error(
    process_datasets_PH(data.frame(), ph_representations = "unknown"),
    "ph_representations must contain supported values"
  )
  expect_error(
    process_datasets_PH(data.frame(), ph_representations = character()),
    "ph_representations must contain supported values"
  )
})
