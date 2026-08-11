mv05u_test_view <- function(payload = NULL) {
  if (is.null(payload)) {
    payload <- matrix(seq_len(384L * 30L) / 1000, nrow = 384L, ncol = 30L)
  }
  rownames(payload) <- sprintf("cell_%03d", seq_len(384L))
  colnames(payload) <- paste0("PC", seq_len(30L))
  source <- list(
    contract = list(profile = "scientific", scientific_eligible = TRUE),
    cache_key = "test_source_v1:abc", sample_id = "sample_a",
    cohort = "large_loso_v1:test", representation = "sct_whole",
    fit_scope_id = "training_only:test", subsample_seed = 20260805L
  )
  .new_topology_view(
    view_id = "cell_topology_v1", source = source,
    point_metric = "euclidean_frozen_shared_coordinates_v1",
    payload = payload, point_ids = rownames(payload),
    coordinate_ids = colnames(payload),
    transformations = list(source = "test"),
    payload_sha256 = .scientific_digest(payload)
  )
}

test_that("MV5-U transforms only the four admitted coordinate factors", {
  view <- mv05u_test_view()
  configurations <- mv05t_configuration_registry_v1()
  transformed <- lapply(configurations$configuration_id, function(id) {
    mv05u_transform_view_v1(view, id)
  })
  expect_equal(vapply(transformed, function(item) nrow(item$payload), integer(1L)),
               c(192L, 256L, 384L, 384L))
  expect_equal(vapply(transformed, function(item) ncol(item$payload), integer(1L)),
               c(30L, 30L, 20L, 30L))
  expect_true(all(transformed[[1L]]$point_ids %in% transformed[[2L]]$point_ids))
  expect_identical(transformed[[3L]]$payload, view$payload[, seq_len(20L)])
  expect_equal(unname(rowSums(transformed[[4L]]$payload^2)), rep(1, 384L),
               tolerance = 1e-12)
  expect_true(all(vapply(transformed, function(item) {
    item$transformations$mv05u_robustness$outcomes_authorized == FALSE
  }, logical(1L))))
})

test_that("MV5-U rejects unadmitted configurations and zero-norm cosine rows", {
  view <- mv05u_test_view()
  expect_error(mv05u_transform_view_v1(view, "factorial_search"),
               "configuration")
  payload <- view$payload
  payload[1L, ] <- 0
  zero <- mv05u_test_view(payload)
  expect_error(mv05u_transform_view_v1(
    zero, "cells384_pc30_cosine_chord_v1"
  ), "zero-norm")
})

test_that("MV5-U freezes balanced deterministic pair coverage", {
  training <- sprintf("train_%02d", 1:20)
  query <- sprintf("query_%02d", 1:4)
  first <- mv05u_pair_coverage_from_ids_v1(
    training, query, "large_loso_v1:test", 20260805L
  )
  second <- mv05u_pair_coverage_from_ids_v1(
    rev(training), rev(query), "large_loso_v1:test", 20260805L
  )
  expect_identical(first, second)
  expect_equal(nrow(first), 32L)
  expect_equal(as.integer(table(first$pair_scope)), c(16L, 16L))
  expect_equal(length(unique(first$pair_request_id)), 32L)
  expect_true(all(first$outcome_label_state == "closed"))
  expect_false(any(first$biological_outcomes_computed))
})

test_that("MV5-U variable-size PH validation uses the exact MST death multiset", {
  view <- mv05u_transform_view_v1(
    mv05u_test_view(), "nested_cells_192_pc30_euclidean_v1"
  )
  deaths <- mv05d2_h0_mst_deaths_v1(view$payload)
  diagram <- rbind(
    cbind(dimension = 0, birth = 0, death = deaths),
    c(dimension = 0, birth = 0, death = Inf)
  )
  result <- structure(list(
    diagram = diagram,
    provenance = list(
      view_cache_key = view$cache_key, point_count = 192L, max_dim = 1L,
      threshold = -1, field = 2L, invalid_interval_count = 0L,
      zero_persistence_count = 0L, essential_h0_count = 1L
    )
  ), class = c("scph_topology_result_v1", "scph_topology_result"))
  oracle <- mv05u_validate_ph_result_v1(result, view)
  expect_true(oracle$passed)
  expect_equal(oracle$finite_h0_intervals, 191L)
  result$diagram[1L, "death"] <- result$diagram[1L, "death"] + 1
  expect_error(mv05u_validate_ph_result_v1(result, view), "MST oracle")
})

test_that("MV5-U bound queue remains exact and unopened", {
  queue <- mv05t_build_admission_queue_v1(paste(rep("a", 64), collapse = ""))
  queue$mv05t_source_freeze_sha256 <- queue$source_freeze_sha256
  queue$mv05t_queue_sha256 <- paste(rep("b", 64), collapse = "")
  queue$implementation_sha256 <- paste(rep("c", 64), collapse = "")
  queue$prospective_head <- paste(rep("d", 40), collapse = "")
  queue$python_executable_sha256 <- paste(rep("e", 64), collapse = "")
  queue$python_version <- "3.10.20"
  queue$persim_version <- "0.3.8"
  queue$numpy_version <- "2.2.6"
  queue$scipy_version <- "1.15.3"
  queue$pair_coverage_per_scope <- 16L
  queue$view_count <- 90L
  queue$landscape_dimensions <- "H0|H1"
  queue$contract_id <- "mv05u_execution_queue_v1"
  expect_invisible(mv05u_validate_execution_queue_v1(queue))
  queue$outcomes_computed[[1L]] <- TRUE
  expect_error(mv05u_validate_execution_queue_v1(queue), "24-unit")
})

test_that("MV5-U resource decision never authorizes full robustness", {
  metrics <- data.frame(
    admission_unit_id = paste0("unit_", 1:24),
    disposition = "completed", elapsed_seconds = 10,
    peak_process_tree_rss_bytes = 1e8, completed_views = 90L,
    landscape_pair_rows = 64L, energy_pair_rows = 32L,
    labels_opened = FALSE, outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
  decision <- mv05u_resource_decision_v1(metrics, 1e6)
  expect_true(decision$admission_complete)
  expect_false(decision$full_robustness_authorized)
  metrics$elapsed_seconds[[1L]] <- 601
  failed <- mv05u_resource_decision_v1(metrics, 1e6)
  expect_false(failed$admission_complete)
  expect_false(failed$full_robustness_authorized)
})

test_that("MV5-U validator parses cross-language boolean fields strictly", {
  expect_identical(
    mv05u_parse_strict_boolean_v1(
      c(TRUE, FALSE, "True", "False", "TRUE", "false"), "flag"
    ),
    c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE)
  )
  expect_error(
    mv05u_parse_strict_boolean_v1(c("true", "yes"), "flag"),
    "must contain only true/false"
  )
  expect_error(
    mv05u_parse_strict_boolean_v1(c(TRUE, NA), "flag"),
    "must contain only true/false"
  )
})
