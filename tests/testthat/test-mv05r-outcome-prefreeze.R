test_that("MV5-R freezes primary and sensitivity algorithm roles", {
  algorithms <- mv05r_algorithm_registry_v1()
  expect_equal(algorithms$role, c("primary", "sensitivity"))
  expect_false(any(algorithms$refit_authorized))
  expect_false(any(algorithms$oracle_k_authorized))
})

test_that("MV5-R endpoint registry separates alignment and generalization", {
  endpoints <- mv05r_endpoint_registry_v1()
  expect_equal(nrow(endpoints), 8L)
  expect_equal(sum(endpoints$evaluation_scope ==
                     "overlapping_training_partition_alignment"), 6L)
  expect_equal(sum(endpoints$evaluation_scope ==
                     "heldout_label_prediction_from_frozen_training_cluster"), 2L)
  expect_false(any(endpoints$p_value_authorized))
  expect_true(all(endpoints$multiplicity_family ==
                    "none_secondary_clustering_complete_reporting"))
  expect_equal(endpoints$uncertainty_policy[[7L]],
               "2000_tissue_stratified_study_block_bootstrap_95_percentile")
  expect_equal(endpoints$uncertainty_policy[[8L]],
               "2000_global_study_block_bootstrap_for_approach_95_percentile")
})

test_that("MV5-R plurality mapping uses a lexical exact tie", {
  partition <- data.frame(sample_id = c("A", "B", "C", "D"),
                          cluster = c(1, 1, 2, 2))
  labels <- data.frame(sample_id = c("A", "B", "C", "D"),
                       tissue = c("zeta", "alpha", "beta", "beta"))
  mapped <- mv05r_plurality_map_v1(partition, labels, "tissue")
  expect_equal(mapped$predicted_label, c("alpha", "beta"))
  expect_equal(mapped$plurality_tie_size, c(2L, 1L))
  expect_error(mv05r_plurality_map_v1(partition[-1, ], labels, "tissue"),
               "misaligned")
})

test_that("MV5-R builds exactly 2,400 unopened evaluation units", {
  analysis <- data.frame(
    analysis_group_id = sprintf("group%03d", 1:150),
    fold_id = rep(sprintf("fold%02d", 1:15), each = 10),
    representation = rep(rep(c("sct_whole", "inductive_integrated"), each = 5), 15),
    distance_id = rep(sprintf("distance%d", 1:5), 30),
    training_samples = 80L, stringsAsFactors = FALSE)
  queue <- mv05r_build_evaluation_queue_v1(analysis,
    source_freeze_sha256 = paste(rep("a", 64), collapse = ""))
  expect_equal(nrow(queue), 2400L)
  expect_invisible(mv05r_validate_evaluation_queue_v1(queue))
  queue$evaluation_executed[[1L]] <- TRUE
  expect_error(mv05r_validate_evaluation_queue_v1(queue), "2,400-unit")
})

test_that("MV5-R preoutcome firewall rejects result columns", {
  expect_error(.mv05r_assert_preoutcome(data.frame(ari = 0.5)),
               "outcome-result")
  expect_invisible(.mv05r_assert_preoutcome(data.frame(
    label_source_state = "frozen_external", outcomes_computed = FALSE)))
})
