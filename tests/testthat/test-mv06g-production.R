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

testthat::test_that("MV6-G post-sort rank validation handles multiples of nine", {
  components <- data.frame(
    view_id = rep(c("cell_topology_v1", "gene_topology_v1"), each = 2L),
    homology_dimension = rep(c("H0", "H1"), 2L),
    stringsAsFactors = FALSE
  )
  scales <- data.frame(
    component_id = c("cell_H0", "cell_H1", "gene_H0", "gene_H1"),
    scale_value = c(2, 3, 5, 7)
  )
  for (training_samples in c(72L, 81L)) {
    training_ids <- sprintf("train_%03d", seq_len(training_samples))
    query <- do.call(rbind, lapply(seq_along(training_ids), function(index) {
      data.frame(
        group_id = "g", fold_id = "f", seed = 1L,
        query_sample_id = "query_1", training_sample_id = training_ids[[index]],
        view_id = components$view_id,
        homology_dimension = components$homology_dimension,
        distance = index + seq_len(4L) / 10,
        outcome_label_state = "closed", biological_outcomes_computed = FALSE,
        stringsAsFactors = FALSE
      )
    }))
    rankings <- mv06g_build_group_rankings_v1(
      query, scales, training_samples, training_samples
    )
    testthat::expect_equal(nrow(rankings), 9L * training_samples)
    testthat::expect_true(all(table(rankings$method_id) == training_samples))
    testthat::expect_true(all(vapply(split(rankings, rankings$method_id),
      function(block) identical(sort(as.integer(block$rank)),
                                 seq_len(training_samples)), logical(1L))))
  }
})
