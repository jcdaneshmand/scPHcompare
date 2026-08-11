test_that("MV5-T admits only the three outcome-independent robustness families", {
  candidates <- mv05t_candidate_registry_v1()
  expect_setequal(candidates$candidate_id[candidates$selected_for_admission],
                  c("nested_cell_count_192_256", "pc20_truncation",
                    "cosine_chord_geometry"))
  expect_false(any(candidates$mv05s_values_used_for_candidate_choice))
  expect_false(any(candidates$new_outcomes_computed))
})

test_that("MV5-T configurations are one factor at a time", {
  configurations <- mv05t_configuration_registry_v1()
  expect_equal(nrow(configurations), 4L)
  changed <- (configurations$cells != 384L) +
    (configurations$coordinates != 30L) +
    (configurations$point_metric != "euclidean")
  expect_true(all(changed == 1L))
  expect_false(any(configurations$outcomes_authorized))
})

test_that("MV5-T nested cell selections are deterministic and representation-neutral", {
  ids <- sprintf("cell%03d", 1:384)
  first <- mv05t_nested_point_ids_v1("sample", 20260805L, ids, 192L)
  second <- mv05t_nested_point_ids_v1("sample", 20260805L, rev(ids), 192L)
  larger <- mv05t_nested_point_ids_v1("sample", 20260805L, ids, 256L)
  expect_identical(first, second)
  expect_true(all(first %in% larger))
  expect_equal(length(unique(larger)), 256L)
})

test_that("MV5-T coordinate pairing requires exact cell identities", {
  matrix_value <- matrix(seq_len(384L * 30L), 384L, 30L,
                         dimnames = list(sprintf("cell%03d", 1:384),
                                         sprintf("PC%d", 1:30)))
  view <- list(payload = matrix_value, point_ids = rownames(matrix_value),
               coordinate_ids = colnames(matrix_value))
  expect_invisible(mv05t_validate_coordinate_pair_v1(view, matrix_value))
  changed <- matrix_value
  rownames(changed)[[1L]] <- "different"
  expect_error(mv05t_validate_coordinate_pair_v1(view, changed),
               "identity or shape")
})

test_that("MV5-T admission queue is exact and unopened", {
  queue <- mv05t_build_admission_queue_v1(paste(rep("a", 64), collapse = ""))
  expect_equal(nrow(queue), 24L)
  expect_equal(sum(queue$representation == "sct_whole"), 12L)
  expect_equal(sum(queue$representation == "inductive_integrated"), 12L)
  expect_invisible(mv05t_validate_admission_queue_v1(queue))
  queue$admission_executed[[1L]] <- TRUE
  expect_error(mv05t_validate_admission_queue_v1(queue), "24-unit")
})
