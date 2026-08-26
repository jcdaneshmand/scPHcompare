make_mv10g_fixture <- function() {
  stacks <- names(.mv10g_stack_labels)
  methods <- names(.mv10g_method_labels)
  quality <- expand.grid(
    stack_id = stacks, seed = 20260805:20260809,
    homology_dimension = c("H0", "H1"), method_id = methods, k = 2:10,
    stringsAsFactors = FALSE
  )
  index <- seq_len(nrow(quality))
  quality$mean_silhouette <- -0.2 + (index %% 100L) / 100
  quality$median_silhouette <- quality$mean_silhouette
  quality$minimum_silhouette <- pmax(-1, quality$mean_silhouette - 0.2)
  quality$negative_silhouette_fraction <- (index %% 20L) / 20
  quality$minimum_cluster_size <- 1L + index %% 10L
  quality$maximum_cluster_size <- 20L + index %% 50L
  quality$singleton_clusters <- index %% 4L
  quality$outcome_label_state <- "closed"
  quality$biological_outcomes_computed <- FALSE

  stability <- expand.grid(
    stack_id = stacks, homology_dimension = c("H0", "H1"),
    method_id = methods, k = 2:10, stringsAsFactors = FALSE
  )
  stability$mean_adjusted_rand <- 0.2 +
    (seq_len(nrow(stability)) %% 60L) / 100
  stability$minimum_adjusted_rand <- stability$mean_adjusted_rand - 0.1
  stability$maximum_adjusted_rand <- stability$mean_adjusted_rand + 0.1
  stability$outcome_label_state <- "closed"
  stability$biological_outcomes_computed <- FALSE

  primary_k <- data.frame(
    stack_id = "allqc_residual_exact500",
    homology_dimension = c("H0", "H1"),
    method_id = "pam_dissimilarity_v1", selected_k = c(3L, 4L),
    threshold = c(0.55, 0.60), outcome_label_state = "closed",
    biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
  )
  method_pairs <- t(utils::combn(methods, 2L))
  agreement <- expand.grid(
    stack_id = stacks, homology_dimension = c("H0", "H1"),
    seed = 20260805:20260809, k = 2:10,
    method_pair = seq_len(nrow(method_pairs)), stringsAsFactors = FALSE
  )
  agreement$first_method_id <- method_pairs[agreement$method_pair, 1L]
  agreement$second_method_id <- method_pairs[agreement$method_pair, 2L]
  agreement$adjusted_rand <- -0.1 +
    (seq_len(nrow(agreement)) %% 100L) / 100
  agreement$method_pair <- NULL
  agreement$outcome_label_state <- "closed"
  agreement$biological_outcomes_computed <- FALSE
  list(quality = quality, stability = stability, primary_k = primary_k,
       agreement = agreement)
}

test_that("MV10-G builds a complete deterministic review synthesis", {
  fixture <- make_mv10g_fixture()
  first <- do.call(mv10g_build_review_data_v1, fixture)
  second <- do.call(mv10g_build_review_data_v1, fixture)
  expect_identical(first, second)
  expect_identical(unname(vapply(first, nrow, integer(1L))),
                   c(270L, 270L, 540L, 2L, 18L, 90L))
  expect_setequal(first$stability_grid$representation_label,
                  unname(.mv10g_stack_labels))
  expect_setequal(first$stability_grid$method_label,
                  unname(.mv10g_method_labels))
  expect_identical(sort(unique(first$quality_summary$k)), 2:10)
  expect_true(all(first$quality_summary$seeds == 5L))
  expect_true(all(first$agreement_summary$seeds == 5L))
  expect_equal(length(unique(first$agreement_summary$method_pair_id)), 10L)
  expect_setequal(first$primary_selection$homology_dimension, c("H0", "H1"))
  expect_false(any(vapply(first, function(value) {
    any(c("tissue", "approach", "label", "outcome") %in% tolower(names(value)))
  }, logical(1L))))
})

test_that("MV10-G rejects labels and incomplete sources", {
  fixture <- make_mv10g_fixture()
  fixture$quality$tissue <- "prohibited"
  expect_error(do.call(mv10g_build_review_data_v1, fixture),
               "prohibited label/outcome")
  fixture <- make_mv10g_fixture()
  fixture$stability <- fixture$stability[-1L, ]
  expect_error(do.call(mv10g_build_review_data_v1, fixture),
               "closed source schema")
})

test_that("MV10-G through MV10-K scripts parse and freeze review ordering", {
  root <- testthat::test_path("..", "..")
  scripts <- file.path(root, "scripts", c(
    "build_mv10g_clustering_review_prefreeze.R",
    "run_mv10h_clustering_review_synthesis.R",
    "build_mv10i_clustering_review_closure.R",
    "render_mv10j_clustering_review_figures.R",
    "build_mv10k_clustering_review_figure_closure.R"
  ))
  expect_true(all(file.exists(scripts)))
  for (path in scripts) expect_silent(parse(path))
  text <- setNames(lapply(scripts, function(path) paste(
    readLines(path, warn = FALSE), collapse = "\n"
  )), basename(scripts))
  prefreeze <- text[["build_mv10g_clustering_review_prefreeze.R"]]
  expect_match(prefreeze, "nrows = 0L", fixed = TRUE)
  expect_match(prefreeze,
               "source_values_opened_before_metric_and_figure_freeze = FALSE",
               fixed = TRUE)
  expect_match(prefreeze, "synthesis_authorized_after_commit = TRUE",
               fixed = TRUE)
  expect_match(prefreeze, "figure_render_authorized = FALSE", fixed = TRUE)
  expect_match(prefreeze, "figures = 7L", fixed = TRUE)
  expect_match(text[["run_mv10h_clustering_review_synthesis.R"]],
               "mv10g_build_review_data_v1", fixed = TRUE)
  expect_match(text[["build_mv10i_clustering_review_closure.R"]],
               "independent_numeric_repeat", fixed = TRUE)
  expect_match(text[["render_mv10j_clustering_review_figures.R"]],
               "owner_review_state = \"pending\"", fixed = TRUE)
  figure_closure <- text[["build_mv10k_clustering_review_figure_closure.R"]]
  expect_match(figure_closure, "byte_identical_repeat", fixed = TRUE)
  expect_match(figure_closure, "exact_pixel_repeat", fixed = TRUE)
  expect_match(figure_closure,
               "closed_pending_original_resolution_visual_review", fixed = TRUE)
})
