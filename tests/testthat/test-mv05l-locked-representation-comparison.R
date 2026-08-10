make_mv05l_endpoint_fixture <- function() {
  samples <- data.frame(
    query_sample_id = c("a", "b", "c", "d"),
    held_out_study = c("s1", "s2", "s1", "s2"),
    query_tissue = c("pbmc", "pbmc", "colon", "colon"),
    stringsAsFactors = FALSE
  )
  map <- .mv05l_method_map
  build <- function(representation) {
    rows <- list()
    counter <- 0L
    for (method_index in seq_len(nrow(map))) {
      for (seed in 1:2) {
        for (sample_index in seq_len(nrow(samples))) {
          counter <- counter + 1L
          family <- map$family_id[[method_index]]
          base_rank <- c(h0 = 2L, h1 = 3L, raw_composite = 2L,
                         energy = 2L, pseudobulk = 1L)[[family]]
          rank <- base_rank
          if (representation == "integrated" && family == "h0") rank <- 1L
          nearest_correct <- rank == 1L
          rows[[counter]] <- data.frame(
            contract_id = paste0("fixture_", representation),
            fold_id = paste0("large_loso_v1:", samples$held_out_study[[sample_index]]),
            held_out_study = samples$held_out_study[[sample_index]],
            seed = seed,
            method_id = if (representation == "sct")
              map$sct_method_id[[method_index]] else
                map$integrated_method_id[[method_index]],
            method_role = "fixture",
            query_sample_id = samples$query_sample_id[[sample_index]],
            query_tissue = samples$query_tissue[[sample_index]],
            training_samples = 2L, training_studies = 1L,
            first_same_tissue_sample_id = if (sample_index %% 2L) "b" else "a",
            first_same_tissue_rank = rank,
            reciprocal_rank = 1 / rank,
            nearest_sample_id = if (nearest_correct)
              if (sample_index %% 2L) "b" else "a" else
                if (sample_index <= 2L) "c" else "a",
            nearest_tissue = if (nearest_correct)
              samples$query_tissue[[sample_index]] else
                if (samples$query_tissue[[sample_index]] == "pbmc") "colon" else "pbmc",
            one_nn_correct = nearest_correct,
            nearest_distance_tied = FALSE,
            endpoint_status = "estimable",
            labels_opened_after_prediction_lock = TRUE,
            upstream_refit = FALSE,
            reranked_after_label_open = FALSE,
            stringsAsFactors = FALSE
          )
        }
      }
    }
    do.call(rbind, rows)
  }
  list(sct = build("sct"), integrated = build("integrated"))
}

test_that("MV5-L pairs only exact locked endpoint identities", {
  fixture <- make_mv05l_endpoint_fixture()
  paired <- mv05l_pair_locked_endpoints_v1(
    fixture$sct, fixture$integrated, expected_rows_per_family = 8L
  )
  expect_equal(nrow(paired), 40L)
  expect_identical(unique(paired$family_id), .mv05l_method_map$family_id)
  expect_true(isTRUE(attr(paired, "pseudobulk_identical")))
  expect_true(all(paired$endpoint_status == "estimable"))
  expect_equal(
    paired$direct_reciprocal_rank_difference[paired$family_id == "h0"],
    rep(0.5, 8L)
  )
})

test_that("MV5-L difference-in-differences subtracts matched energy change", {
  fixture <- make_mv05l_endpoint_fixture()
  paired <- mv05l_pair_locked_endpoints_v1(
    fixture$sct, fixture$integrated, expected_rows_per_family = 8L
  )
  summaries <- mv05l_summarize_estimands_v1(paired)
  h0_mrr <- summaries$macro$estimate[
    summaries$macro$endpoint_id == "cross_study_tissue_mrr_v1" &
      summaries$macro$estimand_id == "did_h0_topology_minus_energy"
  ]
  h1_mrr <- summaries$macro$estimate[
    summaries$macro$endpoint_id == "cross_study_tissue_mrr_v1" &
      summaries$macro$estimand_id == "did_h1_topology_minus_energy"
  ]
  expect_equal(h0_mrr, 0.5)
  expect_equal(h1_mrr, 0)
  expect_equal(nrow(summaries$sample), 56L)
  expect_equal(nrow(summaries$tissue), 28L)
  expect_equal(nrow(summaries$macro), 14L)
  pseudo <- summaries$macro$estimate[
    summaries$macro$estimand_id ==
      "direct_pseudobulk_integrated_minus_sct"
  ]
  expect_equal(pseudo, c(0, 0))
})

test_that("MV5-L rejects a shared-pseudobulk identity violation", {
  fixture <- make_mv05l_endpoint_fixture()
  row <- which(fixture$integrated$method_id ==
                 "pseudobulk_training_standardized_panel_v1")[[1L]]
  fixture$integrated$nearest_sample_id[[row]] <- "unexpected"
  expect_error(
    mv05l_pair_locked_endpoints_v1(
      fixture$sct, fixture$integrated, expected_rows_per_family = 8L
    ),
    "pseudobulk identity control failed"
  )
})

make_mv05l_inference_fixture <- function() {
  tissue_counts <- c("bone marrow" = 31L, colon = 13L, liver = 6L,
                     pbmc = 12L, testis = 28L)
  studies <- list(
    "bone marrow" = paste0("bm", 1:3), colon = paste0("co", 1:2),
    liver = paste0("li", 1:2), pbmc = paste0("pb", 1:4),
    testis = paste0("te", 1:4)
  )
  samples <- do.call(rbind, lapply(names(tissue_counts), function(tissue) {
    count <- tissue_counts[[tissue]]
    data.frame(
      query_sample_id = paste0(gsub(" ", "_", tissue), "_", seq_len(count)),
      query_tissue = tissue,
      held_out_study = rep(studies[[tissue]], length.out = count),
      stringsAsFactors = FALSE
    )
  }))
  rows <- list()
  counter <- 0L
  for (endpoint_index in seq_along(.mv05l_endpoints)) {
    for (estimand_index in seq_along(.mv05l_estimands)) {
      counter <- counter + 1L
      rows[[counter]] <- data.frame(
        contract_id = "fixture",
        endpoint_id = .mv05l_endpoints[[endpoint_index]],
        estimand_id = .mv05l_estimands[[estimand_index]],
        estimand_role = unname(.mv05l_estimand_roles[[estimand_index]]),
        samples,
        estimate = estimand_index / 100 +
          (seq_len(nrow(samples)) %% 9L) / 1000,
        completed_seeds = 5L,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

test_that("MV5-L inference is paired, blocked, deterministic, and Holm adjusted", {
  fixture <- make_mv05l_inference_fixture()
  set.seed(91)
  rng_before <- .Random.seed
  first <- mv05l_block_inference_v1(
    fixture, bootstrap_replicates = 30L, bootstrap_seed = 7L,
    randomization_replicates = 49L, randomization_seed = 8L
  )
  expect_identical(.Random.seed, rng_before)
  second <- mv05l_block_inference_v1(
    fixture, bootstrap_replicates = 30L, bootstrap_seed = 7L,
    randomization_replicates = 49L, randomization_seed = 8L
  )
  expect_identical(first, second)
  expect_equal(nrow(first$intervals), 14L)
  expect_equal(nrow(first$primary_contrasts), 4L)
  primary <- first$primary_contrasts$family_id ==
    "F1_representation_did_mrr"
  expect_equal(sum(primary), 2L)
  expect_true(all(first$primary_contrasts$raw_p_value[primary] > 0))
  expect_equal(
    first$primary_contrasts$holm_p_value[primary],
    stats::p.adjust(first$primary_contrasts$raw_p_value[primary],
                    method = "holm")
  )
  expect_true(all(
    first$bootstrap_audit$paired_across_representations_methods_and_estimands
  ))
  expect_false(any(first$bootstrap_audit$seeds_treated_as_independent))
})

test_that("MV5-L conservative boundary rule counts numerical ties", {
  observed <- 0.25
  null <- c(observed - 8 * .Machine$double.eps, 0.1, -observed)
  expect_equal(.mv05l_boundary_exceedances(null, observed), 2L)
})

test_that("MV5-L input lock binds commit and every hash", {
  files <- vapply(1:2, function(index) {
    path <- tempfile(fileext = ".txt")
    writeLines(paste0("locked-", index), path)
    path
  }, character(1L))
  names(files) <- c("a", "b")
  hashes <- vapply(files, .mv05l_file_sha256, character(1L))
  lock <- mv05l_verify_input_lock_v1(
    files, source_commit = "fixture_commit", expected_hashes = hashes,
    expected_commit = "fixture_commit"
  )
  expect_equal(nrow(lock), 2L)
  expect_true(all(lock$lock_passed_before_endpoint_read))
  expect_true(all(lock$marginal_aggregate_outcomes_known_at_specification))
  expect_false(any(lock$joint_sample_contrasts_known_at_specification))
  expect_error(
    mv05l_verify_input_lock_v1(
      files, source_commit = "wrong", expected_hashes = hashes,
      expected_commit = "fixture_commit"
    ),
    "frozen pre-join commit"
  )
})
