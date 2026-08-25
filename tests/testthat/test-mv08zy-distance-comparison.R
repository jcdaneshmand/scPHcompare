test_that("MV8-ZY compares identical pair axes without labels or clustering", {
  units <- sprintf("u%d", 1:5)
  index <- combn(units, 2L)
  left <- data.frame(
    first_unit_id = index[1L, ], second_unit_id = index[2L, ],
    distance = seq_len(ncol(index)), stringsAsFactors = FALSE
  )
  right <- left
  right$distance <- 2 * left$distance
  result <- mv08zy_compare_distance_pairs_v1(left, right, "fixture")
  expect_equal(result$summary$units, 5L)
  expect_equal(result$summary$unordered_pairs, 10L)
  expect_equal(result$summary$pearson, 1)
  expect_equal(result$summary$spearman, 1)
  expect_equal(result$summary$nonnegative_scale, 2)
  expect_equal(result$summary$relative_stress, 0)
  expect_equal(result$summary$median_abs_median_scaled_change, 0)
  expect_equal(result$summary$p95_abs_median_scaled_change, 0)
  expect_equal(result$summary$mean_neighbor_jaccard, 1)
  expect_equal(result$summary$neighbor_k, 4L)
  expect_equal(nrow(result$neighbor), 5L)
  expect_true(all(result$neighbor$neighbor_jaccard == 1))
  expect_equal(result$summary$outcome_label_state, "closed")
  expect_false(result$summary$biological_outcomes_computed)
  expect_equal(result$summary$clustering_jobs, 0L)
  expect_equal(result$summary$fusion_jobs, 0L)
})

test_that("MV8-ZY neighbor ties are deterministic and pair axes fail closed", {
  units <- sprintf("u%d", 1:12)
  index <- combn(units, 2L)
  value <- data.frame(
    first_unit_id = index[2L, ], second_unit_id = index[1L, ],
    distance = rep(1:3, length.out = ncol(index)), stringsAsFactors = FALSE
  )
  result <- mv08zy_compare_distance_pairs_v1(value, value[sample(nrow(value)), ],
                                              "ties")
  expect_equal(result$summary$neighbor_k, 10L)
  expect_equal(result$summary$mean_neighbor_jaccard, 1)
  bad <- value[-1L, ]
  expect_error(mv08zy_compare_distance_pairs_v1(value, bad, "missing"),
               "incomplete")
  duplicated <- rbind(value, value[1L, ])
  expect_error(mv08zy_compare_distance_pairs_v1(duplicated, duplicated, "dup"),
               "duplicated")
  negative <- value
  negative$distance[[1L]] <- -1
  expect_error(mv08zy_compare_distance_pairs_v1(negative, value, "negative"),
               "invalid")
  constant <- value
  constant$distance <- 1
  expect_error(mv08zy_compare_distance_pairs_v1(constant, constant, "constant"),
               "degenerate")
})

test_that("MV8-ZY reads all three immutable stack layouts canonically", {
  root <- tempfile("mv08zy-layout-")
  dir.create(root)
  old <- file.path(root, "old")
  new <- file.path(root, "new")
  corrected <- file.path(root, "corrected")
  dir.create(file.path(old, "landscape", "old_group"), recursive = TRUE)
  dir.create(file.path(new, "production", "group_03", "chunk_001"),
             recursive = TRUE)
  dir.create(file.path(corrected, "groups", "h0"), recursive = TRUE)
  pairs <- data.frame(first_unit_id = c("a", "a", "b"),
                      second_unit_id = c("b", "c", "c"),
                      distance = c(1, 2, 3))
  old_pairs <- pairs
  names(old_pairs)[1:2] <- c("first_sample_id", "second_sample_id")
  write.csv(old_pairs, file.path(old, "landscape", "old_group",
                                 "distances.csv"), row.names = FALSE)
  write.csv(pairs, file.path(new, "production", "group_03", "chunk_001",
                             "distances.csv"), row.names = FALSE)
  write.csv(pairs, file.path(corrected, "groups", "h0", "distances.csv"),
            row.names = FALSE)
  binding <- data.frame(
    source_stage = "MV7-H", source_group_id = "old:group",
    source_group_order = 1L, homology_dimension = "H0",
    unordered_pairs = 3L
  )
  old_result <- mv08zy_read_distance_stack_v1(binding, old, new, corrected)
  binding$source_stage <- "MV8-ZU"
  binding$source_group_order <- 3L
  new_result <- mv08zy_read_distance_stack_v1(binding, old, new, corrected)
  binding$source_stage <- "MV8-ZV-correction"
  corrected_result <- mv08zy_read_distance_stack_v1(binding, old, new,
                                                      corrected)
  expect_identical(old_result$pairs, new_result$pairs)
  expect_identical(new_result$pairs, corrected_result$pairs)
  expect_identical(old_result$pair_axis_sha256,
                   corrected_result$pair_axis_sha256)
  expect_equal(nrow(old_result$file_manifest), 1L)
  unlink(root, recursive = TRUE)
})
