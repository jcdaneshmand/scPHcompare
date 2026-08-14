test_that("MV5-S ARI and max-normalized NMI match independent libraries", {
  first <- c(1, 1, 1, 2, 2, 3, 3, 3)
  second <- c("a", "a", "b", "b", "b", "c", "c", "a")
  expect_equal(mv05s_adjusted_rand_index_v1(first, second),
               mclust::adjustedRandIndex(first, second), tolerance = 1e-14)
  expect_equal(mv05s_nmi_max_v1(first, second),
               aricode::NMI(first, second, variant = "max"), tolerance = 1e-14)
  expect_true(is.na(mv05s_nmi_max_v1(rep(1, 8), rep("a", 8))))
})

mv05s_training_fixture <- function() {
  labels <- data.frame(
    sample_id = LETTERS[1:6], study = c("s1", "s1", "s2", "s2", "s3", "s3"),
    tissue = c("x", "x", "x", "y", "y", "y"),
    approach = rep(c("a", "b"), 3), stringsAsFactors = FALSE)
  selected <- do.call(rbind, lapply(20260805:20260809, function(seed) {
    data.frame(analysis_group_id = "group", algorithm_id = "pam_stability_k_v1",
               seed = seed, sample_id = labels$sample_id,
               cluster = c(1, 1, 1, 2, 2, 2), selected_k = 2L,
               stringsAsFactors = FALSE)
  }))
  queue <- data.frame(
    evaluation_unit_id = "unit", analysis_group_id = "group", fold_id = "fold",
    representation = "rep", distance_id = "distance", training_samples = 6L,
    algorithm_id = "pam_stability_k_v1", algorithm_role = "primary",
    endpoint_id = "training_tissue_ari_v1",
    evaluation_scope = "overlapping_training_partition_alignment",
    label_axis = "tissue", metric_id = "adjusted_rand_index",
    stringsAsFactors = FALSE)
  list(labels = labels, selected = selected, queue = queue)
}

test_that("MV5-S training units retain five seeds and descriptive uncertainty", {
  fixture <- mv05s_training_fixture()
  result <- mv05s_evaluate_training_unit_v1(
    fixture$queue, fixture$selected, fixture$labels)
  expect_equal(nrow(result$seed), 5L)
  expect_true(all(result$seed$status == "completed"))
  expect_equal(result$summary$completed_seeds, 5L)
  expect_equal(result$summary$seed_mean, result$seed$estimate[[1L]])
  expect_equal(result$summary$seed_jackknife_se, 0)
  expect_false(result$summary$p_value_computed)
})

mv05s_heldout_fixture <- function() {
  labels <- data.frame(
    sample_id = LETTERS[1:6], study = c("s1", "s1", "s2", "s2", "s3", "s3"),
    tissue = c("x", "x", "y", "y", "x", "y"),
    approach = c("a", "a", "b", "b", "a", "b"),
    stringsAsFactors = FALSE)
  selected <- do.call(rbind, lapply(20260805:20260809, function(seed) {
    data.frame(analysis_group_id = "group", algorithm_id = "pam_stability_k_v1",
               seed = seed, sample_id = LETTERS[1:4], cluster = c(1, 1, 2, 2),
               selected_k = 2L, stringsAsFactors = FALSE)
  }))
  heldout <- do.call(rbind, lapply(20260805:20260809, function(seed) {
    data.frame(analysis_group_id = "group", algorithm_id = "pam_stability_k_v1",
               seed = seed, query_sample_id = LETTERS[5:6], cluster = c(1, 2),
               stringsAsFactors = FALSE)
  }))
  queue <- data.frame(
    evaluation_unit_id = "unit", analysis_group_id = "group", fold_id = "fold:s3",
    representation = "rep", distance_id = "distance", training_samples = 4L,
    algorithm_id = "pam_stability_k_v1", algorithm_role = "primary",
    endpoint_id = "heldout_tissue_plurality_balanced_accuracy_v1",
    evaluation_scope = "heldout_label_prediction_from_frozen_training_cluster",
    label_axis = "tissue", metric_id = "balanced_accuracy",
    stringsAsFactors = FALSE)
  list(labels = labels, selected = selected, heldout = heldout, queue = queue)
}

test_that("MV5-S held-out mapping is learned only from training labels", {
  fixture <- mv05s_heldout_fixture()
  first <- mv05s_evaluate_heldout_unit_v1(
    fixture$queue, fixture$selected, fixture$heldout, fixture$labels)
  altered <- fixture$labels
  altered$tissue[altered$sample_id %in% LETTERS[5:6]] <- c("q", "q")
  second <- mv05s_evaluate_heldout_unit_v1(
    fixture$queue, fixture$selected, fixture$heldout, altered)
  expect_equal(first$private$predicted_label, second$private$predicted_label)
  expect_false(identical(first$private$true_label, second$private$true_label))
  expect_true(all(first$private$correct))
  expect_equal(nrow(first$seed), 5L)
  expect_false(any(first$seed$p_value_computed))
})

mv05s_bootstrap_fixture <- function() {
  study <- rep(sprintf("s%02d", 1:15), each = 2)
  sample_id <- sprintf("sample%02d", seq_along(study))
  tissue_by_study <- rep(sprintf("t%d", 1:5), each = 3)
  tissue <- rep(tissue_by_study, each = 2)
  approach <- rep("a", length(study))
  for (mixed in c("s01", "s06", "s11")) {
    approach[study == mixed] <- c("a", "b")
  }
  approach[study %in% c("s02", "s07", "s12")] <- "b"
  data.frame(sample_id = sample_id, study = study, tissue = tissue,
             approach = approach,
             seed_mean_correct = rep(c(0, 1), 15), stringsAsFactors = FALSE)
}

test_that("MV5-S bootstrap preserves tissue strata and mixed study blocks", {
  fixture <- mv05s_bootstrap_fixture()
  tissue_first <- mv05s_bootstrap_counts_v1(fixture, "tissue", 50L, 7L)
  tissue_second <- mv05s_bootstrap_counts_v1(fixture, "tissue", 50L, 7L)
  approach <- mv05s_bootstrap_counts_v1(fixture, "approach", 50L, 7L)
  expect_identical(tissue_first, tissue_second)
  expect_true(all(rowSums(tissue_first) == 15L))
  expect_true(all(rowSums(approach) == 15L))
  study_tissue <- unique(fixture[c("study", "tissue")])
  for (tissue_id in unique(study_tissue$tissue)) {
    columns <- study_tissue$study[study_tissue$tissue == tissue_id]
    expect_true(all(rowSums(tissue_first[, columns, drop = FALSE]) == 3L))
  }
  result <- mv05s_apply_bootstrap_v1(fixture, "tissue", tissue_first)
  expect_equal(result$status, "completed")
  expect_equal(result$estimable_replicates, 50L)
  expect_equal(result$ci_lower, 0.5)
  expect_equal(result$ci_upper, 0.5)
})

test_that("MV5-S bootstrap reports insufficient class support", {
  fixture <- mv05s_bootstrap_fixture()
  counts <- matrix(0L, nrow = 10L, ncol = 15L,
                   dimnames = list(NULL, sprintf("s%02d", 1:15)))
  counts[, "s03"] <- 15L
  result <- mv05s_apply_bootstrap_v1(
    fixture, "approach", counts, minimum_estimable_fraction = 0.95)
  expect_equal(result$status, "bootstrap_support_insufficient")
  expect_true(is.na(result$ci_lower))
})
