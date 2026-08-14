source(testthat::test_path("..", "..", "R", "mv06f_production.R"))

mv06f_fixture <- function() {
  study_sizes <- c(1L, 1L, 2L, 2L, 4L, 4L, 4L, 4L, 5L, 5L, 6L, 8L,
                   9L, 10L, 25L)
  studies <- paste0("study_", sprintf("%02d", seq_along(study_sizes)))
  candidate <- do.call(rbind, lapply(seq_along(studies), function(index) {
    data.frame(
      sample_id = paste0(studies[[index]], "_sample_",
                         seq_len(study_sizes[[index]])),
      study = studies[[index]], outcome_label_state = "closed",
      biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
    )
  }))
  folds <- do.call(rbind, lapply(studies, function(study) {
    held_out <- sum(candidate$study == study)
    data.frame(
      fold_id = paste0("fold:", study),
      fit_scope_id = paste0("fold:", study, ":training"),
      held_out_study = study, seed = 20260805:20260809,
      training_samples = 90L - held_out, held_out_samples = held_out,
      outcome_label_state = "closed", biological_outcomes_computed = FALSE,
      stringsAsFactors = FALSE
    )
  }))
  resources <- do.call(rbind, lapply(20260805:20260809, function(seed) {
    data.frame(
      sample_id = candidate$sample_id, seed = seed,
      selected_cell_sha256 = strrep("a", 64L),
      normalization_cache_key = paste0("cache:", candidate$sample_id, ":", seed),
      private_cache_file = paste0(candidate$sample_id, "__", seed, ".rds"),
      private_cache_size_bytes = 1000,
      private_cache_sha256 = strrep("b", 64L), disposition = "built_atomic",
      outcome_label_state = "closed", biological_outcomes_computed = FALSE,
      stringsAsFactors = FALSE
    )
  }))
  panel <- data.frame(
    panel_sha256 =
      "7be22cdb9056fed427c78d58be2b19e258c7c6807e6f7ac3900bd1bfa1d19eb8",
    panel_order = seq_len(500L), feature_id = paste0("feature_", seq_len(500L)),
    gene = paste0("gene_", seq_len(500L)), selected_all_cache_core = TRUE,
    stringsAsFactors = FALSE
  )
  list(candidate = candidate, folds = folds, resources = resources,
       panel = panel)
}

testthat::test_that("MV6-F reconstructs the frozen complete workload", {
  fixture <- mv06f_fixture()
  hashes <- stats::setNames(rep(strrep("c", 64L), 4L),
                            c("candidate", "folds", "resources", "panel"))
  queue <- mv06f_build_group_queue_v1(
    fixture$candidate, fixture$folds, fixture$resources, fixture$panel, hashes
  )
  testthat::expect_silent(mv06f_validate_queue_v1(queue))
  testthat::expect_equal(nrow(queue), 75L)
  testthat::expect_equal(sum(queue$cell_ph_jobs + queue$gene_ph_jobs), 13500L)
  testthat::expect_equal(sum(queue$diagram_dimension_records), 27000L)
  testthat::expect_equal(sum(queue$biological_pairs), 35350L)
  testthat::expect_equal(sum(queue$landscape_component_rows), 141400L)
  stage <- queue[queue$stage == "stage_1_maximum", ]
  testthat::expect_equal(stage$held_out_samples, 25L)
  testthat::expect_equal(stage$seed, 20260807L)
  testthat::expect_equal(stage$execution_order, 1L)
  testthat::expect_match(mv06f_queue_root_v1(queue), "^[0-9a-f]{64}$")
})

testthat::test_that("MV6-F fails closed on axis and label drift", {
  fixture <- mv06f_fixture()
  testthat::expect_silent(mv06f_validate_prefreeze_inputs_v1(
    fixture$candidate, fixture$folds, fixture$resources, fixture$panel
  ))
  opened <- fixture$candidate
  opened$outcome_label_state[[1L]] <- "open"
  testthat::expect_error(mv06f_validate_prefreeze_inputs_v1(
    opened, fixture$folds, fixture$resources, fixture$panel
  ), "closed-label")
  broken <- fixture$folds
  broken$training_samples[[1L]] <- broken$training_samples[[1L]] - 1L
  testthat::expect_error(mv06f_validate_prefreeze_inputs_v1(
    fixture$candidate, broken, fixture$resources, fixture$panel
  ), "identities")
})

testthat::test_that("MV6-F pair identities bind direction, view, and dimension", {
  base <- mv06f_pair_id_v1(
    "group", "query", "training", "cell_topology_v1", "H0"
  )
  testthat::expect_identical(base, mv06f_pair_id_v1(
    "group", "query", "training", "cell_topology_v1", "H0"
  ))
  testthat::expect_false(identical(base, mv06f_pair_id_v1(
    "group", "training", "query", "cell_topology_v1", "H0"
  )))
  testthat::expect_false(identical(base, mv06f_pair_id_v1(
    "group", "query", "training", "gene_topology_v1", "H0"
  )))
  testthat::expect_false(identical(base, mv06f_pair_id_v1(
    "group", "query", "training", "cell_topology_v1", "H1"
  )))
})

testthat::test_that("MV6-F completed group validation binds companion hashes", {
  root <- tempfile("mv06f-complete-")
  dir.create(root)
  records <- lapply(seq_len(180L), function(index) {
    identity <- list(
      contract_id = "mv06f_matched_ph_identity_v1",
      sample_id = paste0("sample_", index), outcome_label_state = "closed",
      biological_outcomes_computed = FALSE
    )
    topology_result <- list(diagram = matrix(numeric(), ncol = 3L))
    oracle <- list(index = index)
    payload_sha <- digest::digest(list(
      identity = identity, topology_result = topology_result,
      h0_mst_oracle = oracle
    ), algo = "sha256", serialize = TRUE)
    structure(list(
      contract_id = "mv06f_matched_ph_record_v1", identity = identity,
      topology_result = topology_result, h0_mst_oracle = oracle,
      payload_sha256 = payload_sha,
      cache_key = paste0("mv06f_matched_ph_record_v1:", payload_sha),
      downstream_execution = list(
        landscape_pairs = 0L, fusion_jobs = 0L, clustering_jobs = 0L,
        outcome_jobs = 0L
      )
    ), class = c("scph_mv06f_ph_record_v1", "list"))
  })
  manifest <- data.frame(
    ph_cache_key = vapply(records, `[[`, character(1L), "cache_key"),
    outcome_label_state = "closed", biological_outcomes_computed = FALSE
  )
  distances <- data.frame(
    pair_id = "pair", squared_distance = 0, exact = TRUE,
    all_active_levels = TRUE, level_cap_applied = FALSE, rust_status = 0L,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE
  )
  metrics <- data.frame(contract_id = "metrics", value = 1)
  saveRDS(records, file.path(root, "diagrams.rds"), compress = FALSE,
          version = 3)
  utils::write.csv(manifest, file.path(root, "diagram-manifest.csv"),
                   row.names = FALSE)
  utils::write.csv(distances, file.path(root, "distances.csv"),
                   row.names = FALSE)
  utils::write.csv(metrics, file.path(root, "metrics.csv"), row.names = FALSE)
  artifacts <- file.path(root, c(
    "diagrams.rds", "diagram-manifest.csv", "distances.csv", "metrics.csv"
  ))
  hashes <- vapply(artifacts, .mv06f_sha256, character(1L))
  bytes <- unname(file.info(artifacts)$size)
  queue_row <- data.frame(
    group_id = "group", biological_pairs = 1L,
    landscape_component_rows = 1L
  )
  status <- data.frame(
    contract_id = "mv06f_group_status_v1", group_id = "group",
    queue_root_sha256 = strrep("a", 64L),
    implementation_root_sha256 = strrep("b", 64L),
    rust_library_sha256 = strrep("c", 64L), completion_state = "complete",
    diagrams_sha256 = hashes[[1L]], diagrams_bytes = bytes[[1L]],
    diagram_manifest_sha256 = hashes[[2L]],
    diagram_manifest_bytes = bytes[[2L]],
    distances_sha256 = hashes[[3L]], distances_bytes = bytes[[3L]],
    metrics_sha256 = hashes[[4L]], metrics_bytes = bytes[[4L]],
    ph_jobs = 180L, diagram_dimension_records = 360L, biological_pairs = 1L,
    landscape_component_rows = 1L, outcome_label_state = "closed",
    biological_outcomes_computed = FALSE, fusion_jobs = 0L,
    clustering_jobs = 0L, outcome_jobs = 0L
  )
  utils::write.csv(status, file.path(root, "status.csv"), row.names = FALSE)
  testthat::expect_silent(mv06f_validate_group_directory_v1(
    root, queue_row, strrep("a", 64L), strrep("b", 64L), strrep("c", 64L)
  ))
  utils::write.csv(data.frame(contract_id = "changed", value = 2),
                   file.path(root, "metrics.csv"), row.names = FALSE)
  testthat::expect_error(mv06f_validate_group_directory_v1(
    root, queue_row, strrep("a", 64L), strrep("b", 64L), strrep("c", 64L)
  ), "hashes")
})
