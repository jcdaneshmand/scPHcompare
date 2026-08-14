make_mv05e_ranking_fixture <- function() {
  labels <- data.frame(
    sample_id = c("a", "b", "c", "d"),
    study = c("s1", "s2", "s1", "s2"),
    tissue = c("pbmc", "pbmc", "colon", "colon"),
    stringsAsFactors = FALSE
  )
  methods <- c("cell_landscape_h0_v1", "cell_landscape_h1_v1")
  seeds <- 1:2
  rows <- list()
  index <- 0L
  for (study in c("s1", "s2")) {
    queries <- labels$sample_id[labels$study == study]
    training <- labels$sample_id[labels$study != study]
    for (seed in seeds) {
      for (method in methods) {
        for (query in queries) {
          query_tissue <- labels$tissue[match(query, labels$sample_id)]
          training_tissue <- labels$tissue[match(training, labels$sample_id)]
          if (method == "cell_landscape_h0_v1") {
            distance <- ifelse(training_tissue == query_tissue, 1, 2)
          } else {
            distance <- ifelse(training_tissue == query_tissue, 2, 1)
          }
          order_index <- order(distance, training, method = "radix")
          for (rank in seq_along(order_index)) {
            pair <- order_index[[rank]]
            index <- index + 1L
            rows[[index]] <- data.frame(
              fold_id = paste0("large_loso_v1:", study), seed = seed,
              method_id = method, method_role = "fixture",
              query_sample_id = query,
              training_sample_id = training[[pair]],
              distance = distance[[pair]], pair_id = paste0("pair_", index),
              neighbor_rank = rank, distance_tie_size = 1L,
              distance_tied = FALSE,
              tie_break_policy = "exact_distance_then_canonical_sample_id_v1",
              outcome_label_state = "closed",
              biological_outcomes_computed = FALSE,
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }
  }
  list(rankings = do.call(rbind, rows), labels = labels,
       methods = methods, seeds = seeds)
}

test_that("MV5-E uses immutable ranks for MRR and fixed 1-NN", {
  fixture <- make_mv05e_ranking_fixture()
  observed <- mv05e_evaluate_retrieval_v1(
    fixture$rankings, fixture$labels,
    expected_rows = nrow(fixture$rankings),
    expected_methods = fixture$methods,
    expected_seeds = fixture$seeds,
    expected_observations = 16L
  )
  h0 <- observed[observed$method_id == "cell_landscape_h0_v1", ]
  h1 <- observed[observed$method_id == "cell_landscape_h1_v1", ]
  expect_equal(h0$reciprocal_rank, rep(1, 8))
  expect_true(all(h0$one_nn_correct))
  expect_equal(h1$reciprocal_rank, rep(0.5, 8))
  expect_false(any(h1$one_nn_correct))
  expect_true(all(observed$endpoint_status == "estimable"))
  expect_false(any(observed$upstream_refit))
  expect_false(any(observed$reranked_after_label_open))
})

test_that("MV5-E rejects noncanonical rank order and premature label opening", {
  fixture <- make_mv05e_ranking_fixture()
  bad <- fixture$rankings
  bad$neighbor_rank[c(1, 2)] <- bad$neighbor_rank[c(2, 1)]
  expect_error(
    mv05e_evaluate_retrieval_v1(
      bad, fixture$labels, expected_rows = nrow(bad),
      expected_methods = fixture$methods, expected_seeds = fixture$seeds,
      expected_observations = 16L
    ),
    "not canonical"
  )
  expect_error(
    mv05e_open_frozen_labels_v1("missing.csv", data.frame(
      label_open_authorized = FALSE
    )),
    "only after a passing prediction lock"
  )
})

test_that("MV5-E records prespecified non-estimability without scoring an error", {
  fixture <- make_mv05e_ranking_fixture()
  fixture$labels$tissue[fixture$labels$sample_id == "b"] <- "liver"
  observed <- mv05e_evaluate_retrieval_v1(
    fixture$rankings, fixture$labels,
    expected_rows = nrow(fixture$rankings),
    expected_methods = fixture$methods,
    expected_seeds = fixture$seeds,
    expected_observations = 16L,
    require_all_estimable = FALSE
  )
  nonestimable <- observed[observed$query_sample_id %in% c("a", "b"), ]
  expect_true(all(
    nonestimable$endpoint_status == "training_tissue_absent_not_estimable"
  ))
  expect_true(all(is.na(nonestimable$reciprocal_rank)))
  expect_true(all(is.na(nonestimable$one_nn_correct)))
})

test_that("MV5-E summaries average technical seeds before tissue macro use", {
  fixture <- make_mv05e_ranking_fixture()
  observed <- mv05e_evaluate_retrieval_v1(
    fixture$rankings, fixture$labels,
    expected_rows = nrow(fixture$rankings),
    expected_methods = fixture$methods,
    expected_seeds = fixture$seeds,
    expected_observations = 16L
  )
  summary <- mv05e_summarize_retrieval_v1(observed)
  expect_equal(summary$sample$seeds, rep(2L, 8))
  expect_equal(
    summary$method$macro_mean_reciprocal_rank,
    c(1, 0.5)
  )
  expect_equal(
    summary$method$macro_one_nn_balanced_accuracy,
    c(1, 0)
  )
})

make_mv05e_inference_fixture <- function() {
  tissue_counts <- c(
    "bone marrow" = 31L, colon = 13L, liver = 6L, pbmc = 12L, testis = 28L
  )
  tissue_studies <- list(
    "bone marrow" = paste0("bm", 1:3), colon = paste0("co", 1:2),
    liver = paste0("li", 1:2), pbmc = paste0("pb", 1:4),
    testis = paste0("te", 1:4)
  )
  samples <- do.call(rbind, lapply(names(tissue_counts), function(tissue) {
    count <- tissue_counts[[tissue]]
    data.frame(
      query_sample_id = paste0(gsub(" ", "_", tissue), "_", seq_len(count)),
      query_tissue = tissue,
      held_out_study = rep(tissue_studies[[tissue]], length.out = count),
      stringsAsFactors = FALSE
    )
  }))
  methods <- c(
    "cell_landscape_h0_v1", "cell_landscape_h1_v1",
    "cell_landscape_h0_h1_raw_euclidean_v1",
    "cell_distribution_energy_shared_pca_v1",
    "pseudobulk_shared_panel_euclidean_v1"
  )
  do.call(rbind, lapply(seq_along(methods), function(method_index) {
    values <- samples
    values$contract_id <- "fixture"
    values$method_id <- methods[[method_index]]
    values$method_role <- "fixture"
    values$mean_reciprocal_rank <-
      0.2 + method_index / 100 + (seq_len(nrow(values)) %% 11) / 100
    values$one_nn_balanced_accuracy <-
      0.1 + method_index / 100 + (seq_len(nrow(values)) %% 7) / 100
    values$observations <- 5L
    values$samples <- 1L
    values$studies <- 1L
    values$seeds <- 5L
    values$failures <- 0L
    values[, c(
      "contract_id", "method_id", "method_role", "query_sample_id",
      "query_tissue", "held_out_study", "mean_reciprocal_rank",
      "one_nn_balanced_accuracy", "observations", "samples", "studies",
      "seeds", "failures"
    )]
  }))
}

test_that("MV5-E inference is paired, blocked, deterministic, and Holm adjusted", {
  fixture <- make_mv05e_inference_fixture()
  set.seed(42)
  rng_before <- .Random.seed
  first <- mv05e_block_inference_v1(
    fixture, bootstrap_replicates = 50L, bootstrap_seed = 7L,
    randomization_replicates = 99L, randomization_seed = 8L
  )
  expect_identical(.Random.seed, rng_before)
  second <- mv05e_block_inference_v1(
    fixture, bootstrap_replicates = 50L, bootstrap_seed = 7L,
    randomization_replicates = 99L, randomization_seed = 8L
  )
  expect_identical(first, second)
  expect_equal(nrow(first$method_intervals), 10L)
  expect_equal(nrow(first$contrasts), 4L)
  primary <- first$contrasts$family_id == "F1_primary_retrieval"
  expect_true(all(first$contrasts$study_blocks == 15L))
  expect_true(all(first$contrasts$raw_p_value[primary] > 0))
  expect_equal(
    first$contrasts$holm_p_value[primary],
    stats::p.adjust(first$contrasts$raw_p_value[primary], method = "holm")
  )
  expect_false(any(first$bootstrap_audit$seeds_treated_as_independent))
})
