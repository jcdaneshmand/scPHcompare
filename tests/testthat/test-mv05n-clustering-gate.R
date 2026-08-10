test_that("MV5-N builds canonical unordered training requests with H0 and H1", {
  ids <- c("q", "d", "a", "c", "b")
  ph <- data.frame(
    job_id = paste0("job_", ids), group_id = "source_group", group_order = 1L,
    fold_id = "large_loso_v1:fixture", seed = 20260805L, sample_id = ids,
    execution_role = c("held_out", rep("training", 4L)),
    representation = "sct_whole", view_id = "cell_topology_v1",
    stringsAsFactors = FALSE
  )
  metrics <- data.frame(
    job_id = ph$job_id, record_cache_key = paste0("record_", ids),
    diagram_sha256 = rep(strrep("a", 64L), length(ids)),
    result_file = paste0(ids, ".rds"),
    result_file_sha256 = rep(strrep("b", 64L), length(ids)),
    stringsAsFactors = FALSE
  )
  observed <- mv05n_build_group_pair_manifest_v1(ph, metrics, max_pairs = 4L)
  expect_equal(nrow(observed), choose(4L, 2L) * 2L)
  expect_true(all(observed$first_sample_id < observed$second_sample_id))
  expect_setequal(observed$homology_dimension, c("H0", "H1"))
  expect_equal(length(unique(observed$pair_request_id)), nrow(observed))
  expect_true(all(table(observed$chunk_id) <= 4L))
  inventory <- mv05n_group_inventory_v1(observed)
  expect_equal(inventory$training_samples, 4L)
  expect_equal(inventory$unordered_training_pairs, 6L)
  expect_equal(inventory$h0_h1_request_rows, 12L)
})

test_that("MV5-N request identities change with representation and reject labels", {
  make_group <- function(representation) {
    ids <- c("a", "b", "c")
    ph <- data.frame(
      job_id = paste0(representation, ids), group_id = paste0("g_", representation),
      group_order = 1L, fold_id = "fold", seed = 20260805L,
      sample_id = ids, execution_role = "training", representation = representation,
      view_id = "cell_topology_v1", stringsAsFactors = FALSE
    )
    metrics <- data.frame(
      job_id = ph$job_id, record_cache_key = paste0("r_", representation, ids),
      diagram_sha256 = rep(strrep("c", 64L), 3L), result_file = paste0(ids, ".rds"),
      result_file_sha256 = rep(strrep("d", 64L), 3L), stringsAsFactors = FALSE
    )
    mv05n_build_group_pair_manifest_v1(ph, metrics)
  }
  sct <- make_group("sct_whole")
  integrated <- make_group("inductive_integrated")
  expect_false(any(sct$pair_request_id %in% integrated$pair_request_id))
  labelled <- data.frame(x = 1, tissue = "forbidden")
  expect_error(.mv05n_assert_no_outcome_columns(labelled), "prohibited")
})

test_that("MV5-N profile selection is deterministic at minimum median and maximum", {
  inventory <- data.frame(
    fold_id = paste0("fold_", letters[1:5]), seed = 20260805L,
    representation = "sct_whole", training_samples = c(65, 80, 86, 86, 89),
    stringsAsFactors = FALSE
  )
  observed <- mv05n_select_admission_profiles_v1(inventory)
  expect_identical(observed$profile, c("minimum", "representative", "maximum"))
  expect_identical(observed$training_samples, c(65, 86, 89))
  expect_identical(observed$fold_id, c("fold_a", "fold_c", "fold_e"))
})

test_that("MV5-N canonical labels ignore arbitrary source cluster names", {
  ids <- c("z", "a", "m", "b")
  first <- mv05n_canonicalize_clusters_v1(ids, c("x", "y", "x", "y"))
  second <- mv05n_canonicalize_clusters_v1(ids, c("7", "2", "7", "2"))
  expect_identical(first, second)
  expect_identical(unname(first[c("a", "b", "m", "z")]), c(1L, 1L, 2L, 2L))
})

test_that("MV5-N fits the frozen five-seed grid and selects k label-free", {
  ids <- paste0("s", 1:6)
  base <- as.matrix(stats::dist(cbind(c(0, 0.1, 0.2, 10, 10.1, 10.2), 0)))
  dimnames(base) <- list(ids, ids)
  matrices <- stats::setNames(lapply(20260805:20260809, function(seed) base),
                              as.character(20260805:20260809))
  assignments <- mv05n_fit_five_seed_partitions_v1(matrices, "pam")
  expect_setequal(unique(assignments$k), 2:5)
  selected <- mv05n_select_k_v1(assignments)
  expect_identical(selected$status, "selected")
  expect_identical(selected$selected_k, 2L)
  expect_equal(selected$summary$mean_stability, rep(1, 4L))
})

test_that("MV5-N PAM held-out assignment uses canonical medoid tie-breaking", {
  training <- data.frame(
    sample_id = c("a", "b", "c", "d"), cluster = c(1L, 1L, 2L, 2L),
    is_medoid = c(TRUE, FALSE, TRUE, FALSE), stringsAsFactors = FALSE
  )
  distances <- expand.grid(
    query_sample_id = c("q2", "q1"), training_sample_id = training$sample_id,
    stringsAsFactors = FALSE
  )
  distances$distance <- c(1, 4, 1, 4, 1, 5, 1, 5)
  # q1 and q2 tie at medoids a/c; canonical medoid a wins.
  observed <- mv05n_assign_pam_heldout_v1(distances, training)
  expect_identical(observed$query_sample_id, c("q1", "q2"))
  expect_identical(observed$cluster, c(1L, 1L))
  expect_identical(observed$assignment_reference, c("a", "a"))
})

test_that("MV5-N average-linkage held-out assignment uses cluster means", {
  training <- data.frame(
    sample_id = c("a", "b", "c", "d"), cluster = c(1L, 1L, 2L, 2L),
    stringsAsFactors = FALSE
  )
  distances <- data.frame(
    query_sample_id = "q", training_sample_id = training$sample_id,
    distance = c(1, 3, 2, 2), stringsAsFactors = FALSE
  )
  # Both means are 2, so the cluster containing canonical member a wins.
  observed <- mv05n_assign_average_heldout_v1(distances, training)
  expect_identical(observed$cluster, 1L)
  expect_match(observed$assignment_reference, "^a")
})

test_that("MV5-N pseudobulk training pairs preserve Euclidean distance", {
  vectors <- list(a = c(g1 = 0, g2 = 0), b = c(g1 = 3, g2 = 4),
                  c = c(g1 = 0, g2 = 4))
  pairs <- data.frame(first_sample_id = c("a", "b"),
                      second_sample_id = c("b", "c"), stringsAsFactors = FALSE)
  observed <- mv05n_training_pseudobulk_pairs_v1(vectors, pairs)
  expect_equal(observed$distance, c(5, 3))
})

test_that("MV5-N training energy preserves the accepted V-statistic formula", {
  make_view <- function(sample_id, shift) {
    cells <- paste0(sample_id, "_c", 1:4)
    expression <- matrix(seq_len(12L), nrow = 3L,
                         dimnames = list(paste0("g", 1:3), cells))
    coordinates <- cbind(
      PC1 = c(1, 2, 3, 4) + shift,
      PC2 = c(4, 3, 2, 1) + shift / 2
    )
    rownames(coordinates) <- cells
    source <- new_cell_projection_source(
      expression, sample_id = sample_id, cohort = "fixture",
      representation = "fixture", fit_scope_id = "fixture:training",
      subsample_seed = 20260805L, standardization_id = "fixture_scale",
      contract_profile = "analytical_fixture", expected_genes = 3L,
      expected_cells = 4L, expected_pcs = 2L
    )
    construct_frozen_cell_topology_view(
      source, coordinates, "fixture_coordinates", "fixture_fit"
    )
  }
  views <- list(a = make_view("a", 0), b = make_view("b", 1),
                c = make_view("c", 3))
  pairs <- data.frame(first_sample_id = c("a", "b"),
                      second_sample_id = c("b", "c"), stringsAsFactors = FALSE)
  observed <- mv05n_training_energy_pairs_v1(views, pairs)
  expected <- c(
    .mv05d5_energy_distance(
      views$a$payload, views$b$payload,
      .mv05d5_within_mean_distance(views$a$payload),
      .mv05d5_within_mean_distance(views$b$payload)
    ),
    .mv05d5_energy_distance(
      views$b$payload, views$c$payload,
      .mv05d5_within_mean_distance(views$b$payload),
      .mv05d5_within_mean_distance(views$c$payload)
    )
  )
  expect_equal(observed$distance, expected, tolerance = 1e-14)
})

test_that("MV5-N rejects incomplete and invalid clustering inputs", {
  ids <- c("a", "b", "c")
  matrix <- matrix(c(0, 1, 2, 1, 0, 1, 2, 1, 0), 3L,
                   dimnames = list(ids, ids))
  expect_error(mv05n_pam_partition_v1(matrix, 3L), "2:min")
  matrices <- stats::setNames(rep(list(matrix), 4L),
                              as.character(20260805:20260808))
  expect_error(mv05n_fit_five_seed_partitions_v1(matrices), "exactly the five")
  bad <- data.frame(query_sample_id = "q", training_sample_id = "a", distance = 1)
  partition <- data.frame(sample_id = ids, cluster = c(1L, 1L, 2L),
                          is_medoid = c(TRUE, FALSE, TRUE))
  expect_error(mv05n_assign_pam_heldout_v1(bad, partition), "every training")
})
