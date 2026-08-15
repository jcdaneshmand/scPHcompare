testthat::test_that("MV6-G dynamic scales support all production group sizes", {
  for (samples in c(65L, 80L, 89L)) {
    pairs <- samples * (samples - 1L) / 2L
    components <- c("cell_H0", "cell_H1", "gene_H0", "gene_H1")
    rows <- do.call(rbind, lapply(seq_along(components), function(index) {
      data.frame(group_id = "g", component_id = components[[index]],
                 distance = seq_len(pairs) + index,
                 outcome_label_state = "closed",
                 biological_outcomes_computed = FALSE)
    }))
    scales <- mv06g_fit_group_scales_v1(rows, pairs)
    testthat::expect_equal(nrow(scales), 4L)
    testthat::expect_true(all(scales$training_pairs == pairs))
  }
})

testthat::test_that("MV6-G production source inventory is exact", {
  paths <- mv06g_production_source_paths_v1()
  testthat::expect_length(paths, 8L)
  testthat::expect_identical(anyDuplicated(paths), 0L)
  testthat::expect_true(all(file.exists(testthat::test_path("..", "..", paths))))
})
