test_that("MV9-H proves k=7 degeneracy and freezes deterministic k2/k3", {
  units <- sprintf("unit_%02d", 1:8)
  axis <- utils::combn(units, 2L)
  left <- data.frame(
    first_unit_id = axis[1L, ], second_unit_id = axis[2L, ],
    distance = seq_len(ncol(axis)), stringsAsFactors = FALSE
  )
  right <- left
  right$distance <- rev(left$distance)
  result <- mv09h_external_neighbor_sensitivity_v1(
    left, right, "external_fixture_H0", c(2L, 3L)
  )
  repeat_result <- mv09h_external_neighbor_sensitivity_v1(
    left[nrow(left):1L, ], right[nrow(right):1L, ],
    "external_fixture_H0", c(2L, 3L)
  )

  expect_equal(nrow(result$summary), 2L)
  expect_equal(nrow(result$unit), 16L)
  expect_identical(sort(unique(result$summary$k)), c(2L, 3L))
  expect_true(all(result$unit$neighbor_jaccard >= 0 &
                    result$unit$neighbor_jaccard <= 1))
  expect_identical(result$summary, repeat_result$summary)
  expect_identical(result$unit, repeat_result$unit)

  proof <- result$degeneracy
  expect_true(proof$k_equals_all_other_units)
  expect_equal(proof$possible_neighbor_sets_per_unit, 1)
  expect_equal(proof$jaccard_for_any_two_complete_rankings, 1)
  expect_false(proof$informative_for_neighborhood_preservation)
  expect_match(proof$disposition, "structurally_noninformative")

  expect_error(mv09h_external_neighbor_sensitivity_v1(
    left, right, "external_fixture_H0", c(2L, 4L)
  ), "contract drift")
  expect_error(mv09h_external_neighbor_sensitivity_v1(
    left[left$first_unit_id != "unit_08" &
           left$second_unit_id != "unit_08", ],
    right[right$first_unit_id != "unit_08" &
            right$second_unit_id != "unit_08", ],
    "external_fixture_H0", c(2L, 3L)
  ), "contract drift")
})

test_that("MV9-H corrected review data separates global and local evidence", {
  mv09b <- testthat::test_path("..", "..", "docs", "audits",
                               "mv09b-robustness-synthesis-v1")
  base <- mv09d_prepare_review_figure_data_v1(mv09b)
  mapping <- unique(base$external[c("contrast_id", "homology_dimension")])
  neighbor <- do.call(rbind, lapply(seq_len(nrow(mapping)), function(i) {
    data.frame(
      execution_head = paste(rep("a", 40L), collapse = ""),
      sensitivity_order = i, comparison_order = 30L + i,
      comparison_id = paste0("external_fixture_", i),
      contrast_id = mapping$contrast_id[[i]],
      homology_dimension = as.character(mapping$homology_dimension[[i]]),
      k = c(2L, 3L),
      mean_neighbor_jaccard = c(0.5, 0.6),
      median_neighbor_jaccard = c(0.5, 0.6),
      p10_neighbor_jaccard = c(0.25, 1 / 3),
      stringsAsFactors = FALSE
    )
  }))
  root <- tempfile("mv09h-neighbor-")
  dir.create(root)
  utils::write.csv(neighbor, file.path(
    root, "mv09i-external-neighbor-summary.csv"
  ), row.names = FALSE)
  utils::write.csv(mv09h_neighbor_degeneracy_v1(8L, 7L), file.path(
    root, "mv09i-degeneracy-classification.csv"
  ), row.names = FALSE)

  result <- mv09h_prepare_corrected_review_data_v1(mv09b, root)
  expect_equal(nrow(result$internal), 120L)
  expect_equal(nrow(result$internal_summary), 24L)
  expect_equal(nrow(result$external_global), 30L)
  expect_equal(nrow(result$external_neighbor), 20L)
  expect_equal(nrow(result$global_delta), 60L)
  expect_identical(sort(unique(result$external_neighbor$k)), c(2L, 3L))
  expect_false(any(result$external_global$metric ==
                     "mean_neighbor_jaccard"))
  expect_false(any(result$global_delta$metric ==
                     "mean_neighbor_jaccard"))
  expect_identical(
    result$metric_contract$disposition[which(result$metric_contract$k == 7L)],
    "structurally_noninformative_exclude"
  )
})

test_that("MV9-H through MV9-L scripts parse before prospective execution", {
  scripts <- c(
    "build_mv09h_external_neighbor_prefreeze.R",
    "run_mv09i_external_neighbor_sensitivity.R",
    "build_mv09j_external_neighbor_closure.R",
    "render_mv09k_corrected_review_figures.R",
    "build_mv09l_corrected_review_figure_closure.R"
  )
  paths <- testthat::test_path("..", "..", "scripts", scripts)
  expect_true(all(file.exists(paths)))
  expect_silent(lapply(paths, parse))
})
