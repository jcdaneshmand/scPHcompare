mv06d_fixture_axes <- function() {
  studies <- paste0("S", sprintf("%02d", seq_len(15L)))
  candidate <- do.call(rbind, lapply(seq_along(studies), function(index) {
    data.frame(
      sample_id = paste0(studies[[index]], "_", seq_len(6L)),
      study = studies[[index]], outcome_label_state = "closed",
      biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
    )
  }))
  folds <- do.call(rbind, lapply(20260805:20260809, function(seed) {
    data.frame(
      fold_id = paste0("fold:", studies, ":", seed),
      fit_scope_id = paste0("fit:", studies, ":", seed),
      held_out_study = studies, seed = seed, training_samples = 84L,
      held_out_samples = 6L, outcome_label_state = "closed",
      biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
    )
  }))
  resources <- do.call(rbind, lapply(20260805:20260809, function(seed) {
    data.frame(
      sample_id = candidate$sample_id, seed = seed,
      selected_cell_sha256 = paste(rep("a", 64L), collapse = ""),
      normalization_cache_key = paste0("cache:", candidate$sample_id, ":", seed),
      private_cache_file = paste0(candidate$sample_id, "__", seed, ".rds"),
      private_cache_size_bytes = seq_len(90L),
      private_cache_sha256 = paste(rep("b", 64L), collapse = ""),
      disposition = "built_atomic", outcome_label_state = "closed",
      biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
    )
  }))
  list(candidate = candidate, folds = folds, resources = resources)
}

test_that("MV6-D sentinel selection is deterministic and label closed", {
  fixture <- mv06d_fixture_axes()
  first <- mv06d_select_sentinels_v1(
    fixture$candidate, fixture$folds, fixture$resources
  )
  second <- mv06d_select_sentinels_v1(
    fixture$candidate[sample(nrow(fixture$candidate)), ],
    fixture$folds[sample(nrow(fixture$folds)), ],
    fixture$resources[sample(nrow(fixture$resources)), ]
  )
  expect_identical(first, second)
  expect_equal(nrow(first), 10L)
  expect_equal(sort(unique(first$study_size_rank)), c(1L, 4L, 8L, 12L, 15L))
  expect_equal(sort(unique(first$seed)), 20260805:20260809)
  expect_equal(sum(first$stage == 1L), 2L)
  expect_setequal(first$role[first$stage == 1L], c("held_out", "training"))
  expect_true(all(first$outcome_label_state == "closed"))
  expect_false(any(first$biological_outcomes_computed))
})

mv06d_analytical_sources <- function() {
  genes <- paste0("g", seq_len(4L))
  cells <- paste0("c", seq_len(6L))
  ids <- paste0("sample", seq_len(4L))
  matrices <- lapply(seq_along(ids), function(sample_index) {
    value <- outer(seq_len(4L), seq_len(6L), function(gene, cell) {
      sin(gene * cell / 5 + sample_index) + gene / 7 + cell / 11
    })
    dimnames(value) <- list(paste0("f", seq_len(4L)), cells)
    value
  })
  names(matrices) <- ids
  panel <- data.frame(feature_id = paste0("f", seq_len(4L)), gene = genes)
  keys <- stats::setNames(paste0("key:", ids), ids)
  prepared <- mv06d_prepare_matched_sources_v1(
    matrices, panel, ids[1:3], "fold", "fit", 1L, keys,
    contract_profile = "analytical_fixture", expected_genes = 4L,
    expected_cells = 6L, expected_pcs = 2L
  )
  pca <- fit_cell_topology_pca(prepared$sources[ids[1:3]],
                               n_components = 2L, pca_seed = 1L)
  source <- prepared$sources[[ids[[4L]]]]
  cell <- construct_frozen_cell_topology_view(
    source, t(source$matrix) %*% pca$rotation, "fixture_pca", pca$cache_key
  )
  gene <- construct_gene_topology_view(source)
  list(cell = cell, gene = gene)
}

test_that("MV6-D validates cell and gene PH against their actual metrics", {
  views <- mv06d_analytical_sources()
  for (view_id in names(views)) {
    result <- run_topology_view_ph(views[[view_id]])
    oracle <- mv06d_validate_ph_result_v1(result, views[[view_id]])
    expect_true(oracle$passed)
    expect_equal(oracle$finite_h0_intervals,
                 length(views[[view_id]]$point_ids) - 1L)
    record <- mv06d_new_ph_record_v1(
      "source:key", paste0("sentinel:", view_id), "held_out",
      views[[view_id]], result
    )
    expect_invisible(mv06d_validate_ph_record_v1(record))
    expect_true(all(unlist(record$downstream_execution) == 0L))
  }
})

test_that("MV6-D landscape records retain all-level components without runtime", {
  views <- mv06d_analytical_sources()
  cell_result <- run_topology_view_ph(views$cell)
  first <- mv06d_new_ph_record_v1(
    "source:key", "held", "held_out", views$cell, cell_result
  )
  second <- mv06d_new_ph_record_v1(
    "source:key", "train", "training", views$cell, cell_result
  )
  landscape <- persistence_landscape_distance(
    first$topology_result$diagram, second$topology_result$diagram,
    method = "auto"
  )
  record <- mv06d_new_landscape_record_v1(first, second, landscape)
  repeated <- mv06d_new_landscape_record_v1(first, second, landscape)
  expect_identical(record, repeated)
  expect_null(record$result$runtime)
  expect_identical(record$result$specification,
                   "full_l2_error_controlled_v1")
  expect_identical(record$result$provenance$level_policy,
                   "all consecutive active levels; zero-pad missing depth")
  expect_equal(record$result$distances[c("H0", "H1")], c(H0 = 0, H1 = 0))
})

test_that("MV6-D projection preserves separate source, PH, and landscape units", {
  source_metrics <- data.frame(elapsed_seconds = 1:5)
  ph_metrics <- data.frame(
    elapsed_seconds = rep(c(2, 4), each = 10),
    view_id = rep(c("cell_topology_v1", "gene_topology_v1"), each = 10)
  )
  landscape_metrics <- data.frame(
    elapsed_seconds = rep(c(3, 6), each = 5),
    view_id = rep(c("cell_topology_v1", "gene_topology_v1"), each = 5)
  )
  projection <- mv06d_project_workload_v1(
    source_metrics, ph_metrics, landscape_metrics
  )
  expect_equal(nrow(projection), 15L)
  expect_setequal(unique(projection$component), c(
    "fold_source_pca", "cell_ph", "gene_ph", "cell_landscape_pair",
    "gene_landscape_pair"
  ))
  expect_true(all(projection$outcome_label_state == "closed"))
  expect_false(any(projection$biological_outcomes_computed))
})
