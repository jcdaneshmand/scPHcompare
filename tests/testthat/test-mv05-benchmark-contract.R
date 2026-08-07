test_that("MV-05 contract is immutable and sample-level", {
  first <- new_mv05_benchmark_contract_v1()
  second <- new_mv05_benchmark_contract_v1()
  expect_s3_class(first, "scph_mv05_benchmark_contract_v1")
  expect_identical(first$cache_key, second$cache_key)
  expect_true(all(!grepl("^cell$", first$endpoint_registry$analysis_unit)))
  expect_true(all(c(
    "cell_distribution_energy_shared_pca_v1",
    "pseudobulk_shared_panel_euclidean_v1",
    "gene_correlation_frobenius_v1"
  ) %in% first$method_registry$method_id))
  excluded <- first$method_registry$role == "excluded"
  expect_true(all(first$method_registry$execution_status[excluded] == "prohibited"))
})

test_that("MV-05 label firewall opens labels only at declared stages", {
  expect_invisible(mv05_validate_label_use_v1("outer_split_assignment", "study"))
  expect_error(
    mv05_validate_label_use_v1("outer_split_assignment", "tissue"),
    "Prohibited label use"
  )
  expect_error(
    mv05_validate_label_use_v1("fit_transform", "study"),
    "Prohibited label use"
  )
  expect_invisible(mv05_validate_label_use_v1(
    "evaluate_technical", c("study", "approach")
  ))
})

test_that("existing-data feasibility and LOSO folds respect study blocks", {
  metadata <- data.frame(
    cohort = c(rep("large", 5), rep("bone", 2)),
    sample_id = letters[1:7],
    study = c("s1", "s1", "s2", "s3", "s3", "b1", "b1"),
    tissue = c("t1", "t2", "t1", "t2", "t2", "marrow", "marrow"),
    approach = c("sc", "sc", "sn", "sc", "sn", "sc", "sn"),
    filtered_cells = rep(500L, 7), stringsAsFactors = FALSE
  )
  feasibility <- mv05_design_feasibility_v1(metadata)
  expect_true(feasibility$cross_study_tissue_eligible[
    feasibility$cohort == "large" & feasibility$tissue == "t1"
  ])
  expect_true(feasibility$cross_study_tissue_eligible[
    feasibility$cohort == "large" & feasibility$tissue == "t2"
  ])
  expect_false(feasibility$cross_study_tissue_eligible[
    feasibility$cohort == "bone"
  ])
  folds <- mv05_loso_fold_summary_v1(metadata)
  expect_equal(nrow(folds), 3L)
  expect_true(all(folds$no_study_overlap))
  expect_equal(sum(folds$test_samples), 5L)
})

test_that("Monte Carlo p-values cannot be zero", {
  expect_equal(mv05_monte_carlo_p_v1(0, 9999), 1 / 10000)
  expect_equal(mv05_monte_carlo_p_v1(9999, 9999), 1)
  expect_error(mv05_monte_carlo_p_v1(10, 9), "valid integer counts")
})

test_that("stable-k selection uses repeated inputs without labels", {
  samples <- letters[1:6]
  seeds <- as.character(20260805:20260809)
  rows <- list()
  index <- 0L
  for (seed_index in seq_along(seeds)) {
    index <- index + 1L
    rows[[index]] <- data.frame(
      seed = seeds[[seed_index]], k = 2L, sample_id = samples,
      cluster = c(1, 1, 1, 2, 2, 2)
    )
    index <- index + 1L
    rows[[index]] <- data.frame(
      seed = seeds[[seed_index]], k = 3L, sample_id = samples,
      cluster = c(1, 1, 2, 2, 3, if (seed_index %% 2L) 3 else 2)
    )
  }
  result <- mv05_select_stable_k_v1(do.call(rbind, rows))
  expect_identical(result$status, "selected")
  expect_equal(result$selected_k, 2L)
  expect_equal(result$summary$mean_stability[result$summary$k == 2L], 1)

  degenerate <- do.call(rbind, lapply(seeds, function(seed) {
    data.frame(seed = seed, k = 2L, sample_id = samples, cluster = 1L)
  }))
  expect_identical(mv05_select_stable_k_v1(degenerate)$status, "no_stable_k")

  incomplete <- do.call(rbind, rows[seq_len(length(rows) - 2L)])
  expect_identical(mv05_select_stable_k_v1(incomplete)$status, "no_stable_k")
})

test_that("descriptive pilot fits cannot masquerade as blocked confirmation", {
  distance <- matrix(
    c(0, 1, 2, 1, 0, 3, 2, 3, 0), nrow = 3,
    dimnames = list(c("a", "b", "c"), c("a", "b", "c"))
  )
  expect_error(
    mv05_validate_matrix_record_v1(
      distance, "confirmatory_blocked", "descriptive_all_pilot_samples"
    ),
    "prohibited"
  )
  expect_invisible(mv05_validate_matrix_record_v1(
    distance, "technical_smoke_test", "descriptive_all_pilot_samples"
  ))
  expect_invisible(mv05_validate_matrix_record_v1(
    distance, "confirmatory_blocked", "large_loso_v1:s1:training"
  ))
})
