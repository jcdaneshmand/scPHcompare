test_that("MV5-Q freezes the explicit ten-analysis alias registry", {
  registry <- mv05q_method_alias_registry_v1()
  expect_invisible(mv05q_validate_alias_registry_v1(registry))
  expect_equal(nrow(registry), 10L)
  expect_equal(sum(registry$training_source_policy ==
                     "shared_sct_training_pseudobulk"), 2L)
  expect_equal(registry$query_representation[registry$representation ==
    "sct_whole" & registry$distance_id == "cell_landscape_h0_v1"], "sct_fold")
  expect_equal(registry$query_representation[registry$representation ==
    "sct_whole" & registry$distance_id ==
      "pseudobulk_training_standardized_panel_v1"], "sct_fold")
  expect_error(mv05q_validate_alias_registry_v1(registry[-1, ]),
               "ten-analysis")
})

test_that("MV5-Q reconstructs complete matrices and raw H0-H1 composites", {
  ids <- c("C", "A", "B")
  pairs <- data.frame(first_sample_id = c("A", "A", "B"),
                      second_sample_id = c("B", "C", "C"),
                      distance = c(3, 4, 5),
                      outcome_label_state = "closed",
                      biological_outcomes_computed = FALSE)
  matrix <- mv05q_reconstruct_distance_matrix_v1(pairs, ids)
  expect_identical(rownames(matrix), c("A", "B", "C"))
  expect_equal(unname(matrix["B", "C"]), 5)
  combined <- mv05q_combine_h0_h1_v1(matrix, matrix)
  expect_equal(unname(combined["A", "B"]), sqrt(18))
  expect_error(mv05q_reconstruct_distance_matrix_v1(pairs[-1, ], ids),
               "complete")
})

test_that("MV5-Q rejects query slices that open labels or miss axes", {
  query <- expand.grid(query_sample_id = c("Q1", "Q2"),
                       training_sample_id = c("A", "B", "C"),
                       stringsAsFactors = FALSE)
  query$fold_id <- "fold"
  query$seed <- 20260805L
  query$representation <- "sct_fold"
  query$method_id <- "cell_landscape_h0_v1"
  query$distance <- seq_len(nrow(query))
  query$outcome_label_state <- "closed"
  query$biological_outcomes_computed <- FALSE
  expect_equal(nrow(mv05q_validate_query_slice_v1(query, c("A", "B", "C"))), 6L)
  query$tissue <- "forbidden"
  expect_error(mv05q_validate_query_slice_v1(query, c("A", "B", "C")),
               "prohibited outcome")
})

test_that("MV5-Q production wrapper preserves frozen selection and tie rules", {
  ids <- sprintf("S%02d", 1:12)
  base <- as.matrix(dist(cbind(c(rep(0, 6), rep(10, 6)), seq_along(ids) / 100)))
  dimnames(base) <- list(ids, ids)
  matrices <- stats::setNames(lapply(20260805:20260809, function(seed) base),
                              as.character(20260805:20260809))
  queries <- stats::setNames(lapply(20260805:20260809, function(seed) {
    data.frame(fold_id = "fold", seed = seed, representation = "sct_fold",
               method_id = "cell_landscape_h0_v1", query_sample_id = "Q",
               training_sample_id = ids, distance = rep(1, length(ids)),
               outcome_label_state = "closed",
               biological_outcomes_computed = FALSE,
               stringsAsFactors = FALSE)
  }), as.character(20260805:20260809))
  result <- mv05q_fit_analysis_group_v1(matrices, queries, "fold", "sct_whole",
                                        "cell_landscape_h0_v1")
  expect_invisible(mv05q_validate_group_result_v1(result, 12L, 1L))
  expect_equal(result$selected_k, 2L)
  pam <- result$heldout_assignments[result$heldout_assignments$algorithm_id ==
    "pam_stability_k_v1", ]
  expect_true(all(pam$assignment_reference == min(pam$assignment_reference)))
})
