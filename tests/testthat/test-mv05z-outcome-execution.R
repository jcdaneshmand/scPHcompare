test_that("MV5-Z ranks exact ties canonically without labels", {
  rows <- expand.grid(
    query_sample_id = c("q1", "q2"),
    training_sample_id = c("b", "a", "c"),
    method_id = c("cell_landscape_h0_v1", "cell_landscape_h1_v1",
                  "cell_landscape_h0_h1_raw_euclidean_v1",
                  "cell_distribution_energy_shared_pca_v1"),
    stringsAsFactors = FALSE)
  rows$distance <- c(2, 1, 1)[match(rows$training_sample_id, c("b", "a", "c"))]
  rows$robustness_group_id <- "g"
  rows$outcome_label_state <- "closed"; rows$outcomes_computed <- FALSE
  group <- data.frame(robustness_group_id = "g", fold_id = "large_loso_v1:S",
    seed = 20260805L, representation = "sct_whole",
    configuration_id = "cells384_pc20_euclidean_v1",
    method_rows = nrow(rows), stringsAsFactors = FALSE)
  ranked <- mv05z_rank_group_v1(rows, group)
  part <- ranked[ranked$query_sample_id == "q1" &
                   ranked$method_id == "cell_landscape_h0_v1", ]
  expect_equal(part$training_sample_id, c("a", "c", "b"))
  expect_equal(part$neighbor_rank, 1:3)
  expect_equal(part$distance_tie_size, c(2L, 2L, 1L))
  expect_false(any(ranked$labels_opened))
})

test_that("MV5-Z blocked inference freezes 24 intervals and four tests", {
  registry <- mv05y_estimand_registry_v1()
  study_sizes <- c(5L, 2L, 2L, 10L, 6L, 4L, 1L, 5L, 9L, 4L, 4L,
                   1L, 25L, 4L, 8L)
  studies <- paste0("S", seq_along(study_sizes))
  tissue <- rep(c("a", "a", "a", "b", "b", "c", "c", "d", "d", "d",
                  "d", "e", "e", "e", "e"), study_sizes)
  study <- rep(studies, study_sizes)
  samples <- paste0("sample", seq_along(study))
  value <- seq_along(study) / 1000
  summary <- do.call(rbind, lapply(seq_len(nrow(registry)), function(i) {
    data.frame(
      estimand_id = registry$estimand_id[[i]],
      estimand_order = registry$estimand_order[[i]],
      query_sample_id = samples, query_tissue = tissue,
      held_out_study = study, estimate = value * i,
      completed_seeds = 5L, stringsAsFactors = FALSE)
  }))
  result <- mv05z_block_inference_v1(
    summary, registry, bootstrap_replicates = 20L,
    randomization_replicates = 31L)
  expect_equal(nrow(result$intervals), 24L)
  expect_equal(nrow(result$primary), 4L)
  expect_equal(nrow(result$bootstrap_counts), 20L)
  expect_equal(nrow(result$sign_matrix), 31L)
  expect_true(all(result$primary$raw_p_value > 0))
})

test_that("MV5-Z boundary exceedances conservatively retain numerical ties", {
  observed <- 0.25
  null <- c(observed, -observed, observed - 2 * .Machine$double.eps, 0)
  expect_equal(.mv05z_boundary_exceedances(null, observed), 3L)
})
