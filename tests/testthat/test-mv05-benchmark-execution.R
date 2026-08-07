make_mv05_sources <- function() {
  genes <- paste0("g", 1:4)
  cells <- paste0("c", 1:5)
  first <- matrix(
    c(1, 2, 4, 7, 11, 2, 5, 3, 8, 13, 4, 1, 6, 9, 15, 3, 8, 2, 12, 5),
    nrow = 4, byrow = TRUE, dimnames = list(genes, cells)
  )
  second <- first + matrix(
    c(0, 1, 0, 1, 0, 1, 0, 2, 0, 1, 0, 2, 1, 0, 1, 2, 0, 1, 1, 0),
    nrow = 4, byrow = TRUE
  )
  third <- first + matrix(
    c(1, 0, 2, 0, 1, 0, 1, 0, 2, 0, 2, 0, 1, 0, 2, 0, 1, 2, 0, 1),
    nrow = 4, byrow = TRUE
  )
  matrices <- list(sample_b = second, sample_a = first, sample_c = third)
  lapply(names(matrices), function(id) new_dual_view_source(
    matrices[[id]], sample_id = id, cohort = "fixture",
    representation = "sct_fixture", fit_scope_id = "fixture_fold:training",
    subsample_seed = 20260805L, standardization_id = "fixture_scaling_v1",
    contract_profile = "analytical_fixture", expected_genes = 4L,
    expected_cells = 5L, expected_pcs = 2L
  ))
}

test_that("label-closed LOSO manifests are immutable and partition studies", {
  metadata <- data.frame(
    cohort = rep("large", 5), sample_id = letters[1:5],
    study = c("s1", "s1", "s2", "s3", "s3"),
    tissue = c("t1", "t2", "t1", "t2", "t2"),
    approach = c("sc", "sc", "sn", "sc", "sn"),
    filtered_cells = rep(500L, 5), stringsAsFactors = FALSE
  )
  manifest <- mv05_build_loso_manifest_v1(metadata)
  expect_s3_class(manifest, "scph_mv05_loso_manifest_v1")
  expect_equal(nrow(manifest$table), 15L)
  expect_equal(length(unique(manifest$table$fold_id)), 3L)
  expect_false(any(c("tissue", "approach") %in% names(manifest$table)))
  expect_true(all(manifest$table$outcome_label_state == "closed"))
  summary <- mv05_loso_manifest_summary_v1(manifest)
  expect_equal(sum(summary$held_out_samples), 5L)
  expect_true(all(summary$total_samples == 5L))

  tampered <- manifest
  tampered$table$sample_id[[1L]] <- "changed"
  expect_error(mv05_validate_loso_manifest_v1(tampered), "partition|stale")
})

test_that("matched baselines have analytical and provenance contracts", {
  sources <- make_mv05_sources()
  pca <- fit_cell_topology_pca(sources, n_components = 2L)
  cell_views <- lapply(sources, construct_cell_topology_view, pca_model = pca)

  energy <- mv05_cell_energy_baseline_v1(cell_views)
  pseudobulk <- mv05_pseudobulk_baseline_v1(sources)
  gene <- mv05_gene_correlation_baseline_v1(sources)
  for (bundle in list(energy, pseudobulk, gene)) {
    expect_s3_class(bundle, "scph_mv05_baseline_bundle_v1")
    expect_invisible(mv05_validate_baseline_bundle_v1(bundle))
    expect_identical(rownames(bundle$distance_matrix), sort(c(
      "sample_a", "sample_b", "sample_c"
    )))
    expect_equal(unname(diag(bundle$distance_matrix)), rep(0, 3))
  }

  means <- do.call(rbind, lapply(sources[order(vapply(
    sources, `[[`, character(1L), "sample_id"
  ))], function(x) rowMeans(x$matrix)))
  expect_equal(unname(pseudobulk$distance_matrix),
               unname(as.matrix(stats::dist(means))))

  identical_energy <- .mv05_empirical_energy_distance(
    cell_views[[1L]]$payload, cell_views[[1L]]$payload
  )
  expect_equal(identical_energy, 0, tolerance = 1e-12)
  expect_equal(
    .mv05_empirical_energy_distance(matrix(c(0, 2), ncol = 1),
                                    matrix(c(1, 3), ncol = 1)),
    1, tolerance = 1e-12
  )

  canonical_sources <- sources[order(vapply(
    sources, `[[`, character(1L), "sample_id"
  ))]
  first_correlation <- stats::cor(t(canonical_sources[[1L]]$matrix))
  second_correlation <- stats::cor(t(canonical_sources[[2L]]$matrix))
  expect_equal(
    gene$distance_matrix[1, 2],
    sqrt(mean((first_correlation - second_correlation) ^ 2))
  )

  tampered <- energy
  tampered$distance_matrix[1, 2] <- tampered$distance_matrix[1, 2] + 1
  tampered$distance_matrix[2, 1] <- tampered$distance_matrix[1, 2]
  expect_error(mv05_validate_baseline_bundle_v1(tampered), "digest")
})

test_that("baseline constructors reject mismatched fit scopes", {
  sources <- make_mv05_sources()
  original <- sources[[2L]]
  sources[[2L]] <- new_dual_view_source(
    original$matrix, sample_id = original$sample_id, cohort = original$cohort,
    representation = original$representation,
    fit_scope_id = "leaking_other_fold",
    subsample_seed = original$subsample_seed,
    standardization_id = original$standardization_id,
    contract_profile = "analytical_fixture", expected_genes = 4L,
    expected_cells = 5L, expected_pcs = 2L
  )
  expect_error(mv05_pseudobulk_baseline_v1(sources), "fit_scope_id")
})

test_that("inductive mapping result contract is label-closed and immutable", {
  embedding <- matrix(
    c(1, 2, 3, 4, 5, 6), nrow = 3,
    dimnames = list(paste0("q", 1:3), c("PC1", "PC2"))
  )
  result <- .mv05_new_mapping_result(
    fold_id = "fixture:q", fit_scope_id = "fixture:q:training",
    reference_sample_ids = c("r2", "r1"), held_out_sample_id = "q",
    features = c("g1", "g2"), dimensions = 1:2, seed = 20260805L,
    reference_identity_sha256 = paste(rep("a", 64), collapse = ""),
    anchor_count = 4L, query_embeddings = embedding
  )
  expect_s3_class(result, "scph_mv05_inductive_mapping_v1")
  expect_identical(result$reference_sample_ids, c("r1", "r2"))
  expect_false(result$biological_outcomes_computed)
  expect_invisible(mv05_validate_inductive_mapping_v1(result))

  tampered <- result
  tampered$query_embeddings[1, 1] <- 99
  expect_error(mv05_validate_inductive_mapping_v1(tampered), "stale")

  unavailable <- mv05_try_inductive_mapping_v1(reference = NULL, query = NULL)
  expect_identical(unavailable$status, "held_out_mapping_unavailable")
  expect_match(unavailable$reason, "Seurat objects")
})

test_that("frozen shared coordinates retain typed source identity", {
  genes <- paste0("g", seq_len(4L))
  cells <- paste0("c", seq_len(5L))
  value <- matrix(
    seq_len(20L), nrow = 4L, dimnames = list(genes, cells)
  )
  source <- new_dual_view_source(
    value, sample_id = "sample_a", cohort = "fixture",
    representation = "inductive_integrated",
    fit_scope_id = "fold_a:training", subsample_seed = 20260805L,
    standardization_id = "fixture_standardization",
    contract_profile = "analytical_fixture", expected_genes = 4L,
    expected_cells = 5L, expected_pcs = 2L
  )
  coordinates <- matrix(
    c(1, 2, 2, 1, 3, 1, 1, 3, 4, 2), nrow = 5L, byrow = TRUE,
    dimnames = list(cells, c("PC1", "PC2"))
  )
  view <- construct_frozen_cell_topology_view(
    source, coordinates, "fixture_mapping_v1", "fixture_mapping_cache_key"
  )
  expect_silent(validate_topology_view(view))
  expect_identical(view$sample_id, source$sample_id)
  expect_identical(view$payload, coordinates)
  expect_identical(
    view$transformations$coordinate_contract_id, "fixture_mapping_v1"
  )
  expect_error(
    construct_frozen_cell_topology_view(
      source, coordinates[rev(rownames(coordinates)), ],
      "fixture_mapping_v1", "fixture_mapping_cache_key"
    ),
    "exact ordered cell IDs"
  )
})

test_that("cell-projection sources permit held-out constants without gene use", {
  cells <- paste0("c", seq_len(5L))
  first <- matrix(
    c(1, 1, 1, 1, 1, 1:15), nrow = 4L, byrow = TRUE,
    dimnames = list(paste0("g", 1:4), cells)
  )
  second <- first + 1
  sources <- lapply(list(sample_a = first, sample_b = second), function(value) {
    new_cell_projection_source(
      value, sample_id = paste0("sample_", if (identical(value, first)) "a" else "b"),
      cohort = "fixture", representation = "sct_fold",
      fit_scope_id = "fold:training", subsample_seed = 20260805L,
      standardization_id = "training_scale", contract_profile = "analytical_fixture",
      expected_genes = 4L, expected_cells = 5L, expected_pcs = 2L
    )
  })
  expect_s3_class(sources[[1L]], "scph_cell_projection_source_v1")
  expect_error(construct_gene_topology_view(sources[[1L]]), "dual_view_source")
  pseudobulk <- mv05_pseudobulk_baseline_v1(sources)
  expect_silent(mv05_validate_baseline_bundle_v1(pseudobulk))
  coordinates <- matrix(
    seq_len(10L), nrow = 5L,
    dimnames = list(cells, c("PC1", "PC2"))
  )
  view <- construct_frozen_cell_topology_view(
    sources[[1L]], coordinates, "training_pca_v1", "training_pca_cache"
  )
  expect_silent(validate_topology_view(view))
})
