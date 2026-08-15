test_that("MV7-FP combines the exact 450 plus 170 cache axis", {
  primary <- data.frame(sample_id = sprintf("p%03d", rep(1:90, 5)),
    seed = rep(20260805:20260809, each = 90), selected_cell_sha256 = "a",
    normalization_cache_key = "b", payload_contract_id =
      "mv05d0_sct_data_matrix_v1", payload_sha256 = "c",
    private_cache_file = paste0("p", seq_len(450), ".rds"),
    private_cache_size_bytes = 1, private_cache_sha256 = "d",
    outcome_label_state = "closed", biological_outcomes_computed = FALSE)
  added <- data.frame(sample_id = sprintf("a%03d", rep(1:34, 5)),
    seed = rep(20260805:20260809, each = 34), selected_cell_sha256 = "e",
    normalization_cache_key = "f", payload_sha256 = "g",
    private_cache_file = paste0("a", seq_len(170), ".rds"),
    private_cache_bytes = 1, private_cache_sha256 = "h",
    finite_sct_payload = TRUE, outcome_label_state = "closed",
    biological_outcomes_computed = FALSE)
  result <- mv07fp_combined_cache_manifest_v1(primary, added)
  expect_equal(nrow(result), 620L)
  expect_equal(length(unique(result$sample_id)), 124L)
  expect_true(all(table(result$seed) == 124L))
  expect_equal(mv07fp_resource_caps_v1()$elapsed_cap_seconds, 2700)
})

test_that("MV7-FP variance calculation is sparse-safe and unchanged", {
  dense <- matrix(c(0, 1, 2, 3, 0, 4, 1, 1), nrow = 2L)
  sparse <- Matrix::Matrix(dense, sparse = TRUE)

  expect_equal(.mv03_row_variance(sparse), .mv03_row_variance(dense),
               tolerance = 1e-14)
  expect_equal(.mv03_row_variance(dense), apply(dense, 1L, stats::var),
               tolerance = 1e-14)
})
