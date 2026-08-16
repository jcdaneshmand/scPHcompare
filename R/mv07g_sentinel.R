# MV7-G full-124 global descriptive typed-view/PH sentinel helpers.

mv07g_resource_caps_v1 <- function() {
  data.frame(
    contract_id = "mv07g_resource_caps_v1",
    stage = c("global_fit_views", "cell_ph", "gene_ph", "aggregate"),
    elapsed_cap_seconds = c(1800, 600, 1800, 14400),
    rss_cap_bytes = c(8, 4, 8, 8) * 1024^3,
    storage_cap_bytes = c(NA, NA, NA, 2 * 1024^3),
    workers = 1L, retries = 0L, atomic_write = TRUE,
    stringsAsFactors = FALSE
  )
}

mv07g_sentinel_axis_v1 <- function(sentinels, cache_manifest) {
  required <- c("sample_id", "selection_boundary", "selected_cells")
  if (!is.data.frame(sentinels) || !all(required %in% names(sentinels)) ||
      nrow(sentinels) != 6L || anyDuplicated(sentinels$sample_id) ||
      !all(sentinels$selected_cells == 384L)) {
    stop("MV7-G requires the exact six frozen sentinels.", call. = FALSE)
  }
  if (!is.data.frame(cache_manifest) || nrow(cache_manifest) != 620L ||
      length(unique(cache_manifest$sample_id)) != 124L ||
      !identical(sort(unique(cache_manifest$seed)), 20260805:20260809) ||
      any(table(cache_manifest$seed) != 124L) ||
      !all(sentinels$sample_id %in% cache_manifest$sample_id)) {
    stop("MV7-G cache manifest is not the frozen 124 by five axis.",
         call. = FALSE)
  }
  axis <- merge(
    sentinels[c("sample_id", "selection_boundary", "selected_cells")],
    unique(cache_manifest[c("sample_id", "seed")]),
    by = "sample_id", sort = FALSE
  )
  axis <- axis[order(axis$seed, axis$sample_id, method = "radix"),,
               drop = FALSE]
  rownames(axis) <- NULL
  axis$contract_id <- "mv07g_sentinel_axis_v1"
  axis$sentinel_basis <- "frozen_mv07d_tissue_depth_extremes_identity_only"
  axis$outcome_label_state <- "closed"
  axis$biological_outcomes_computed <- FALSE
  axis[c("contract_id", "sample_id", "selection_boundary", "seed",
         "selected_cells", "sentinel_basis", "outcome_label_state",
         "biological_outcomes_computed")]
}

mv07g_queue_v1 <- function(axis) {
  if (!is.data.frame(axis) || nrow(axis) != 30L ||
      any(table(axis$seed) != 6L)) {
    stop("MV7-G sentinel axis must contain six samples by five seeds.",
         call. = FALSE)
  }
  caps <- mv07g_resource_caps_v1()
  source <- data.frame(
    stage = "global_fit_views", seed = sort(unique(axis$seed)),
    sample_id = NA_character_, view_id = NA_character_,
    stringsAsFactors = FALSE
  )
  ph <- merge(
    axis[c("sample_id", "seed")],
    data.frame(view_id = c("cell_topology_v1", "gene_topology_v1"),
               stringsAsFactors = FALSE),
    all = TRUE, sort = FALSE
  )
  ph$stage <- ifelse(ph$view_id == "cell_topology_v1", "cell_ph", "gene_ph")
  ph <- ph[c("stage", "seed", "sample_id", "view_id")]
  queue <- rbind(source, ph)
  queue <- queue[order(match(queue$stage,
                             c("global_fit_views", "cell_ph", "gene_ph")),
                       queue$seed, queue$sample_id, method = "radix"),,
                 drop = FALSE]
  queue$job_order <- seq_len(nrow(queue))
  queue$job_id <- ifelse(
    queue$stage == "global_fit_views",
    paste0("fit__", queue$seed),
    paste(queue$stage, queue$seed, queue$sample_id, sep = "__")
  )
  queue$output_file <- ifelse(
    queue$stage == "global_fit_views",
    paste0("source/mv07g__", queue$seed, "__source.rds"),
    paste0("ph/mv07g__", queue$seed, "__", queue$sample_id, "__",
           queue$view_id, "__ph.rds")
  )
  matched <- match(queue$stage, caps$stage)
  queue$elapsed_cap_seconds <- caps$elapsed_cap_seconds[matched]
  queue$rss_cap_bytes <- caps$rss_cap_bytes[matched]
  queue$workers <- 1L
  queue$retries <- 0L
  queue$outcome_label_state <- "closed"
  queue$biological_outcomes_computed <- FALSE
  queue[c("job_order", "job_id", "stage", "seed", "sample_id", "view_id",
          "output_file", "elapsed_cap_seconds", "rss_cap_bytes", "workers",
          "retries", "outcome_label_state", "biological_outcomes_computed")]
}

.mv07g_source_payload <- function(record) {
  list(identity = record$identity, panel = record$panel,
       center = record$center, scale = record$scale,
       pca_model = record$pca_model, views = record$views)
}

mv07g_new_source_record_v1 <- function(identity, panel, prepared, pca_model,
                                       views) {
  record <- list(
    contract_id = "mv07g_global_transform_sentinel_views_v1",
    identity = identity, panel = panel, center = prepared$center,
    scale = prepared$scale, pca_model = pca_model, views = views,
    payload_sha256 = NULL, cache_key = NULL,
    downstream_execution = list(ph_jobs = 0L, landscape_jobs = 0L,
      distance_jobs = 0L, clustering_jobs = 0L, label_jobs = 0L,
      biological_outcome_jobs = 0L)
  )
  record$payload_sha256 <- digest::digest(
    .mv07g_source_payload(record), algo = "sha256", serialize = TRUE)
  record$cache_key <- paste0(
    "mv07g_global_transform_sentinel_views_v1:", record$payload_sha256)
  class(record) <- c("scph_mv07g_source_record_v1", "list")
  mv07g_validate_source_record_v1(record)
  record
}

mv07g_validate_source_record_v1 <- function(record) {
  required <- c("contract_id", "identity", "panel", "center", "scale",
                "pca_model", "views", "payload_sha256", "cache_key",
                "downstream_execution")
  if (!is.list(record) || !all(required %in% names(record)) ||
      !identical(record$contract_id,
                 "mv07g_global_transform_sentinel_views_v1") ||
      !is.list(record$identity) || record$identity$fit_samples != 124L ||
      record$identity$fit_cells != 47616L ||
      record$identity$panel_size != 500L ||
      !(record$identity$seed %in% 20260805:20260809) ||
      length(record$identity$input_cache_keys) != 124L ||
      length(record$identity$sentinel_ids) != 6L ||
      !identical(names(record$views), record$identity$sentinel_ids) ||
      nrow(record$panel) != 500L || length(record$center) != 500L ||
      length(record$scale) != 500L || any(!is.finite(record$center)) ||
      any(!is.finite(record$scale)) ||
      any(record$scale <= sqrt(.Machine$double.eps))) {
    stop("MV7-G source record is incomplete or violates the frozen axis.",
         call. = FALSE)
  }
  .validate_cell_pca_model(record$pca_model)
  if (record$pca_model$n_components != 30L ||
      length(record$pca_model$fit_sample_ids) != 124L ||
      !identical(record$pca_model$cache_key,
                 record$identity$pca_model_cache_key) ||
      !identical(record$pca_model$standardization_id,
                 record$identity$standardization_id)) {
    stop("MV7-G PCA identity is stale.", call. = FALSE)
  }
  for (sample_id in names(record$views)) {
    pair <- record$views[[sample_id]]
    if (!identical(sort(names(pair)),
                   sort(c("cell_topology_v1", "gene_topology_v1")))) {
      stop("MV7-G source record lacks a typed view pair.", call. = FALSE)
    }
    lapply(pair, validate_topology_view)
    if (any(vapply(pair, `[[`, character(1L), "sample_id") != sample_id)) {
      stop("MV7-G typed-view sample identity is stale.", call. = FALSE)
    }
  }
  expected <- digest::digest(.mv07g_source_payload(record), algo = "sha256",
                             serialize = TRUE)
  zeros <- unlist(record$downstream_execution, use.names = FALSE)
  if (!identical(record$payload_sha256, expected) ||
      !identical(record$cache_key, paste0(
        "mv07g_global_transform_sentinel_views_v1:", expected)) ||
      any(zeros != 0)) {
    stop("MV7-G source payload or stop boundary is stale.", call. = FALSE)
  }
  invisible(TRUE)
}

mv07g_h0_mst_deaths_v1 <- function(view) {
  validate_topology_view(view)
  distances <- if (view$view_id == "cell_topology_v1") {
    as.matrix(stats::dist(view$payload, method = "euclidean"))
  } else {
    as.matrix(view$payload)
  }
  n <- nrow(distances)
  selected <- rep(FALSE, n)
  selected[[1L]] <- TRUE
  nearest <- distances[1L, ]
  nearest[[1L]] <- Inf
  edges <- numeric(n - 1L)
  for (index in seq_len(n - 1L)) {
    candidates <- which(!selected)
    chosen <- candidates[[which.min(nearest[candidates])]]
    edges[[index]] <- nearest[[chosen]]
    selected[[chosen]] <- TRUE
    nearest <- pmin(nearest, distances[chosen, ])
    nearest[selected] <- Inf
  }
  sort(edges, method = "radix")
}

mv07g_validate_ph_against_view_v1 <- function(result, view) {
  validate_topology_view(view)
  diagram <- result$diagram
  if (!inherits(result, "scph_topology_result_v1") || !is.matrix(diagram) ||
      !identical(result$provenance$view_cache_key, view$cache_key) ||
      result$provenance$invalid_interval_count != 0L ||
      result$provenance$zero_persistence_count != 0L ||
      result$provenance$essential_h0_count != 1L) {
    stop("MV7-G PH result violates the corrected typed contract.",
         call. = FALSE)
  }
  h0 <- diagram[diagram[, "dimension"] == 0, , drop = FALSE]
  h1 <- diagram[diagram[, "dimension"] == 1, , drop = FALSE]
  finite_h0 <- sort(h0[is.finite(h0[, "death"]), "death"], method = "radix")
  oracle <- mv07g_h0_mst_deaths_v1(view)
  tolerance <- max(1e-7, max(oracle) * 1e-7)
  error <- if (length(finite_h0) == length(oracle)) {
    max(abs(finite_h0 - oracle))
  } else Inf
  if (!is.finite(error) || error > tolerance ||
      length(finite_h0) != length(view$point_ids) - 1L ||
      any(!is.finite(h1[, "death"])) ||
      any(h1[, "death"] <= h1[, "birth"])) {
    stop("MV7-G PH intervals fail the MST/finite-positive gate.",
         call. = FALSE)
  }
  list(
    contract_id = "mv07g_h0_mst_oracle_v1",
    point_count = length(view$point_ids),
    finite_h0_intervals = length(finite_h0),
    finite_h1_intervals = nrow(h1),
    maximum_absolute_error = error, tolerance = tolerance, passed = TRUE
  )
}

.mv07g_ph_payload <- function(record) {
  list(identity = record$identity, topology_result = record$topology_result,
       h0_mst_oracle = record$h0_mst_oracle)
}

mv07g_new_ph_record_v1 <- function(source_record, sample_id, view_id, result) {
  mv07g_validate_source_record_v1(source_record)
  view <- source_record$views[[sample_id]][[view_id]]
  oracle <- mv07g_validate_ph_against_view_v1(result, view)
  identity <- list(
    contract_id = "mv07g_ph_identity_v1",
    source_cache_key = source_record$cache_key,
    sample_id = sample_id, seed = source_record$identity$seed,
    fit_scope_id = source_record$identity$fit_scope_id,
    panel_sha256 = source_record$identity$panel_sha256,
    view_id = view_id, view_cache_key = view$cache_key,
    point_axis_role = view$point_axis_role,
    coordinate_axis_role = view$coordinate_axis_role,
    point_metric = view$point_metric, point_count = length(view$point_ids),
    max_dim = 1L, threshold = -1, field = 2L,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE
  )
  record <- list(
    contract_id = "mv07g_ph_record_v1", identity = identity,
    topology_result = result, h0_mst_oracle = oracle,
    payload_sha256 = NULL, cache_key = NULL,
    downstream_execution = list(landscape_jobs = 0L, distance_jobs = 0L,
      clustering_jobs = 0L, label_jobs = 0L, biological_outcome_jobs = 0L)
  )
  record$payload_sha256 <- digest::digest(
    .mv07g_ph_payload(record), algo = "sha256", serialize = TRUE)
  record$cache_key <- paste0("mv07g_ph_record_v1:", record$payload_sha256)
  class(record) <- c("scph_mv07g_ph_record_v1", "list")
  mv07g_validate_ph_record_v1(record, view)
  record
}

mv07g_validate_ph_record_v1 <- function(record, view = NULL) {
  required <- c("contract_id", "identity", "topology_result",
                "h0_mst_oracle", "payload_sha256", "cache_key",
                "downstream_execution")
  if (!is.list(record) || !all(required %in% names(record)) ||
      !identical(record$contract_id, "mv07g_ph_record_v1") ||
      !(record$identity$view_id %in%
        c("cell_topology_v1", "gene_topology_v1")) ||
      !(record$identity$seed %in% 20260805:20260809) ||
      !identical(record$identity$outcome_label_state, "closed") ||
      isTRUE(record$identity$biological_outcomes_computed) ||
      !isTRUE(record$h0_mst_oracle$passed)) {
    stop("MV7-G PH record is incomplete.", call. = FALSE)
  }
  expected <- digest::digest(.mv07g_ph_payload(record), algo = "sha256",
                             serialize = TRUE)
  if (!identical(record$payload_sha256, expected) ||
      !identical(record$cache_key,
                 paste0("mv07g_ph_record_v1:", expected)) ||
      any(unlist(record$downstream_execution, use.names = FALSE) != 0)) {
    stop("MV7-G PH record payload or stop boundary is stale.", call. = FALSE)
  }
  if (!is.null(view)) {
    if (!identical(record$identity$view_cache_key, view$cache_key)) {
      stop("MV7-G PH record belongs to a different typed view.", call. = FALSE)
    }
    mv07g_validate_ph_against_view_v1(record$topology_result, view)
  }
  invisible(TRUE)
}
