make_mv06h_fixture <- function() {
  labels <- data.frame(
    sample_id = c("a", "b", "c", "d"), study = c("s1", "s2", "s1", "s2"),
    tissue = c("pbmc", "pbmc", "colon", "colon"), stringsAsFactors = FALSE)
  rows <- list(); index <- 0L
  for (study in c("s1", "s2")) {
    queries <- labels$sample_id[labels$study == study]
    training <- labels$sample_id[labels$study != study]
    for (seed in 1:2) for (method in .mv06h_methods) for (query in queries) {
      same <- labels$tissue[match(training, labels$sample_id)] ==
        labels$tissue[match(query, labels$sample_id)]
      distance <- ifelse(same, 1, 2)
      order_index <- order(distance, training, method = "radix")
      group <- paste(study, seed, sep = "_")
      for (rank in seq_along(order_index)) {
        j <- order_index[[rank]]; index <- index + 1L
        rows[[index]] <- data.frame(
          group_id = group, fold_id = paste0("large_loso_v1:", study),
          seed = seed, query_sample_id = query,
          training_sample_id = training[[j]], method_id = method,
          normalized_distance = distance[[j]],
          tie_break = "canonical_training_sample_id", rank = rank,
          outcome_label_state = "closed", biological_outcomes_computed = FALSE,
          fusion_evaluations = 0L, outcome_jobs = 0L,
          stringsAsFactors = FALSE)
      }
    }
  }
  list(rankings = do.call(rbind, rows), labels = labels)
}

test_that("MV6-H evaluates only immutable canonical ranks", {
  fixture <- make_mv06h_fixture()
  outcomes <- mv06h_evaluate_rankings_v1(
    fixture$rankings, fixture$labels,
    expected_rows = nrow(fixture$rankings), expected_observations = 72L,
    expected_samples = 4L, expected_seeds = 1:2)
  expect_equal(nrow(outcomes), 72L)
  expect_true(all(outcomes$reciprocal_rank == 1))
  expect_true(all(outcomes$one_nn_correct))
  expect_false(any(outcomes$upstream_refit))
  expect_false(any(outcomes$reranked_after_label_open))
  bad <- fixture$rankings
  bad$rank[c(1, 2)] <- bad$rank[c(2, 1)]
  expect_error(mv06h_evaluate_rankings_v1(
    bad, fixture$labels, expected_rows = nrow(bad),
    expected_observations = 72L, expected_samples = 4L,
    expected_seeds = 1:2), "not canonical")
})

test_that("MV6-H averages seeds within biological samples", {
  fixture <- make_mv06h_fixture()
  outcomes <- mv06h_evaluate_rankings_v1(
    fixture$rankings, fixture$labels,
    expected_rows = nrow(fixture$rankings), expected_observations = 72L,
    expected_samples = 4L, expected_seeds = 1:2)
  summary <- mv06h_summarize_outcomes_v1(outcomes, expected_seeds = 1:2)
  expect_equal(nrow(summary$sample), 36L)
  expect_true(all(summary$sample$seeds == 2L))
  expect_equal(nrow(summary$tissue), 18L)
  expect_equal(nrow(summary$method), 9L)
})

test_that("MV6-H blocked inference is deterministic and compares fusion twice", {
  samples <- data.frame(
    query_sample_id = c("a", "b", "c", "d"),
    query_tissue = c("pbmc", "pbmc", "colon", "colon"),
    held_out_study = c("s1", "s2", "s1", "s2"), stringsAsFactors = FALSE)
  fixture <- do.call(rbind, lapply(seq_along(.mv06h_methods), function(i) {
    x <- samples; x$method_id <- .mv06h_methods[[i]]
    x$mean_reciprocal_rank <- 0.2 + i / 100 + seq_len(4) / 1000
    x$one_nn_balanced_accuracy <- 0.1 + i / 100
    x$seeds <- 2L
    x
  }))
  # Make the fixed primary better than both required comparators.
  fixture$mean_reciprocal_rank[fixture$method_id == .mv06h_primary] <-
    fixture$mean_reciprocal_rank[fixture$method_id == .mv06h_primary] + 0.2
  set.seed(42); before <- .Random.seed
  first <- mv06h_block_inference_v1(
    fixture, bootstrap_replicates = 20L, bootstrap_seed = 7L,
    randomization_replicates = 99L, randomization_seed = 8L,
    expected_samples = 4L, expected_seeds = 1:2)
  expect_identical(.Random.seed, before)
  second <- mv06h_block_inference_v1(
    fixture, bootstrap_replicates = 20L, bootstrap_seed = 7L,
    randomization_replicates = 99L, randomization_seed = 8L,
    expected_samples = 4L, expected_seeds = 1:2)
  expect_identical(first, second)
  expect_equal(nrow(first$method_intervals), 18L)
  expect_equal(nrow(first$contrasts), 2L)
  expect_setequal(first$contrasts$comparator_id, .mv06h_comparators)
  expect_true(all(first$contrasts$estimate > 0))
  expect_equal(first$contrasts$holm_p_value,
               p.adjust(first$contrasts$raw_p_value, "holm"))
  expect_false(any(first$bootstrap_audit$seeds_treated_as_independent))
})

test_that("MV6-H group manifest and label boundary fail closed", {
  manifest <- data.frame(
    group_id = sprintf("group%02d", 1:75),
    fold_id = rep(sprintf("fold%02d", 1:15), each = 5),
    held_out_study = rep(sprintf("study%02d", 1:15), each = 5),
    seed = rep(.mv06h_seeds, 15), execution_order = 1:75,
    group_root_kind = c("accepted_stage1_sentinel",
                        rep("corrected_serial_completion", 74)),
    group_locator = paste0("ignored_tmp:group", 1:75),
    rankings_sha256 = strrep("a", 64), rankings_bytes = 1,
    ranking_rows = 4242L, scales_sha256 = strrep("b", 64),
    status_sha256 = strrep("c", 64),
    production_implementation_root_sha256 = strrep("d", 64),
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    fusion_evaluations = 0L, outcome_jobs = 0L, stringsAsFactors = FALSE)
  expect_invisible(mv06h_validate_group_manifest_v1(manifest))
  manifest$outcome_jobs[[1L]] <- 1L
  expect_error(mv06h_validate_group_manifest_v1(manifest), "outcome-open")
  expect_error(mv06h_open_frozen_labels_v1(
    "missing.csv", data.frame(label_open_authorized = FALSE,
                              prediction_lock_status = "pending")),
    "only after a passing prediction lock")
})
