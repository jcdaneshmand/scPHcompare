testthat::test_that("MV6-G training pair identities are order invariant", {
  first <- mv06g_training_pair_id_v1(
    "group", "sample_b", "sample_a", "cell_topology_v1", "H1"
  )
  second <- mv06g_training_pair_id_v1(
    "group", "sample_a", "sample_b", "cell_topology_v1", "H1"
  )
  testthat::expect_identical(first, second)
  testthat::expect_error(
    mv06g_training_pair_id_v1(
      "group", "sample_a", "sample_a", "cell_topology_v1", "H1"
    ), "distinct"
  )
})

testthat::test_that("MV6-G fits exactly four training-only median scales", {
  components <- c("cell_H0", "cell_H1", "gene_H0", "gene_H1")
  training <- do.call(rbind, lapply(seq_along(components), function(index) {
    data.frame(
      group_id = "group", view_id = if (index <= 2L) {
        "cell_topology_v1"
      } else "gene_topology_v1",
      homology_dimension = if (index %% 2L) "H0" else "H1",
      component_id = components[[index]],
      distance = seq_len(2080L) + index,
      outcome_label_state = "closed",
      biological_outcomes_computed = FALSE,
      stringsAsFactors = FALSE
    )
  }))
  scales <- mv06g_fit_component_scales_v1(training)
  testthat::expect_equal(nrow(scales), 4L)
  testthat::expect_identical(scales$component_id, components)
  testthat::expect_true(all(scales$query_rows_used == 0L))
  testthat::expect_true(all(scales$labels_used == 0L))
  testthat::expect_equal(scales$scale_value, 1040.5 + seq_along(components))
})

testthat::test_that("MV6-G produces nine deterministic label-closed rankings", {
  queries <- sprintf("query_%02d", 1:25)
  training <- sprintf("training_%02d", 1:65)
  components <- data.frame(
    view_id = rep(c("cell_topology_v1", "gene_topology_v1"), each = 2L),
    homology_dimension = rep(c("H0", "H1"), 2L),
    distance = c(2, 3, 5, 7), stringsAsFactors = FALSE
  )
  rows <- expand.grid(
    query_sample_id = queries, training_sample_id = training,
    component_order = 1:4, stringsAsFactors = FALSE
  )
  rows$group_id <- "group"
  rows$fold_id <- "fold"
  rows$seed <- 20260807L
  rows$view_id <- components$view_id[rows$component_order]
  rows$homology_dimension <-
    components$homology_dimension[rows$component_order]
  rows$distance <- components$distance[rows$component_order]
  rows$outcome_label_state <- "closed"
  rows$biological_outcomes_computed <- FALSE
  scales <- data.frame(
    component_id = c("cell_H0", "cell_H1", "gene_H0", "gene_H1"),
    scale_value = c(2, 3, 5, 7), stringsAsFactors = FALSE
  )
  result <- mv06g_build_rankings_v1(rows, scales)
  testthat::expect_equal(nrow(result), 14625L)
  testthat::expect_setequal(result$method_id, mv06g_method_panel_v1()$method_id)
  testthat::expect_true(all(result$normalized_distance == 1))
  block <- result[result$query_sample_id == queries[[1L]] &
                    result$method_id == "fusion_gene_weight_050", ]
  testthat::expect_identical(block$training_sample_id, training)
  testthat::expect_identical(as.integer(block$rank), 1:65)
  testthat::expect_true(all(result$outcome_label_state == "closed"))
  testthat::expect_false(any(result$biological_outcomes_computed))
})

testthat::test_that("MV6-G stage-one implementation inventory is exact", {
  paths <- mv06g_stage1_source_paths_v1()
  testthat::expect_length(paths, 10L)
  testthat::expect_identical(anyDuplicated(paths), 0L)
  testthat::expect_true(all(file.exists(testthat::test_path("..", "..", paths))))
})
