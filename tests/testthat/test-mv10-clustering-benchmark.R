test_that("MV10 registries freeze the authorized and deferred method scope", {
  methods <- mv10_method_registry_v1()
  expect_equal(nrow(methods), 9L)
  expect_equal(sum(methods$authorized_for_mv10b), 5L)
  expect_identical(methods$method_id[methods$role == "primary"],
                   "pam_dissimilarity_v1")
  expect_true(all(!methods$authorized_for_mv10b[methods$role %in%
                                                  c("deferred", "excluded")]))

  distances <- mv10_distance_registry_v1()
  expect_equal(nrow(distances), 6L)
  expect_identical(distances$distance_id[distances$authorized_for_mv10b], c(
    "exact_all_level_landscape_l2_H0_v2",
    "exact_all_level_landscape_l2_H1_v2"
  ))
  expect_true(all(distances$level_policy[distances$authorized_for_mv10b] ==
                    "all_consecutive_active_levels"))
})

test_that("MV10 validates exactly 30 label-closed internal stack bindings", {
  catalog <- expand.grid(
    stack_id = c("existing_selectedfit_data_exact500",
                 "allqc_data_exact500", "allqc_residual_exact500"),
    seed = 20260805:20260809,
    homology_dimension = c("H0", "H1"),
    stringsAsFactors = FALSE
  )
  catalog$catalog_order <- seq_len(nrow(catalog))
  catalog$dataset_scope <- "internal124"
  catalog$representation_id <- paste0("representation_", catalog$stack_id)
  catalog$panel_id <- "exact500"
  catalog$view_kind <- "gene_topology_v1"
  catalog$units <- 124L
  catalog$unordered_pairs <- choose(124L, 2L)
  catalog$source_stage <- "test"
  catalog$payload_set_sha256 <- strrep("a", 64L)
  catalog$pair_axis_sha256 <- strrep("b", 64L)
  expect_invisible(mv10_validate_stack_catalog_v1(catalog))
  analyses <- mv10_analysis_registry_v1(catalog)
  expect_equal(nrow(analyses), 6L)
  expect_true(all(!analyses$H0_H1_combined))
  expect_true(all(!analyses$cell_gene_combined))

  catalog$tissue <- "prohibited"
  expect_error(mv10_validate_stack_catalog_v1(catalog),
               "prohibited label/outcome")
})

test_that("MV10 executes a deterministic complete label-closed partition grid", {
  points <- rbind(
    cbind(seq(0, 0.5, length.out = 4L), 0),
    cbind(seq(5, 5.5, length.out = 4L), 5),
    cbind(seq(10, 10.5, length.out = 4L), 0)
  )
  ids <- sprintf("unit_%02d", seq_len(nrow(points)))
  matrix <- as.matrix(stats::dist(points))
  rownames(matrix) <- colnames(matrix) <- ids

  first <- mv10_partition_grid_v1(matrix)
  second <- mv10_partition_grid_v1(matrix)
  expect_identical(first, second)
  expect_equal(nrow(first$partitions), 5L * 9L * 12L)
  expect_equal(nrow(first$quality), 5L * 9L)
  expect_setequal(unique(first$partitions$method_id),
                  mv10_method_registry_v1()$method_id[
                    mv10_method_registry_v1()$authorized_for_mv10b])
  expect_identical(sort(unique(first$partitions$k)), 2:10)
  expect_true(all(is.finite(first$quality$mean_silhouette)))
  expect_true(all(first$partitions$outcome_label_state == "closed"))
})

test_that("MV10 computes seed stability, method agreement, and primary K", {
  ids <- sprintf("unit_%02d", 1:12)
  base <- expand.grid(
    stack_id = "allqc_residual_exact500",
    homology_dimension = c("H0", "H1"),
    seed = 20260805:20260809,
    method_id = mv10_method_registry_v1()$method_id[
      mv10_method_registry_v1()$authorized_for_mv10b],
    k = 2:10, sample_id = ids, stringsAsFactors = FALSE
  )
  base$cluster <- ((match(base$sample_id, ids) - 1L) %% base$k) + 1L
  stability <- mv10_seed_stability_v1(base)
  expect_equal(nrow(stability), 2L * 5L * 9L)
  expect_true(all(stability$mean_adjusted_rand == 1))

  agreement_input <- base[
    base$homology_dimension == "H0" & base$seed == 20260805 & base$k == 2L,
    , drop = FALSE
  ]
  agreement <- mv10_method_agreement_v1(agreement_input)
  expect_equal(nrow(agreement), choose(5L, 2L))
  expect_true(all(agreement$adjusted_rand == 1))
  expect_false("sample_id" %in% names(agreement))

  selection <- mv10_select_primary_k_v1(base)
  expect_equal(nrow(selection), 2L)
  expect_identical(selection$homology_dimension, c("H0", "H1"))
  expect_true(all(selection$selected_k == 2L))
  expect_true(all(selection$outcome_label_state == "closed"))
})

test_that("MV10 prefreeze builder parses and creates output only after validation", {
  path <- testthat::test_path(
    "..", "..", "scripts", "build_mv10a_clustering_benchmark_prefreeze.R"
  )
  expect_silent(parse(path))
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")
  validation_position <- regexpr(
    "if (!all(validation$passed))", text, fixed = TRUE
  )[[1L]]
  create_position <- regexpr(
    "if (!dir.exists(output)) dir.create", text, fixed = TRUE
  )[[1L]]
  expect_gt(validation_position, 0L)
  expect_gt(create_position, validation_position)
})
