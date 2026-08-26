.mv13_seeds <- 20260805:20260809
.mv13_residual_id <- "sct_pearson_residual_all_qc_fit_selected384"

mv13_cell_topology_queue_v1 <- function(
    internal_units, external_units,
    seeds = .mv13_seeds) {
  internal_units <- sort(unique(as.character(internal_units)), method = "radix")
  external_units <- sort(unique(as.character(external_units)), method = "radix")
  seeds <- as.integer(seeds)
  if (length(internal_units) != 124L || length(external_units) != 8L ||
      length(seeds) != 5L || anyDuplicated(seeds) ||
      anyNA(c(internal_units, external_units, seeds)) ||
      any(!nzchar(c(internal_units, external_units)))) {
    stop("MV13 queue axes must be internal124, external8, and five seeds.",
         call. = FALSE)
  }
  internal <- expand.grid(
    dataset_scope = "internal124", unit_id = internal_units, seed = seeds,
    panel_id = "exact500", stringsAsFactors = FALSE
  )
  external <- expand.grid(
    dataset_scope = "external8", unit_id = external_units,
    seed = seeds[[1L]], panel_id = c("common475", "exact500"),
    stringsAsFactors = FALSE
  )
  views <- rbind(internal, external)
  views <- views[order(views$dataset_scope, views$panel_id, views$seed,
                       views$unit_id, method = "radix"), , drop = FALSE]
  rownames(views) <- NULL
  views$view_order <- seq_len(nrow(views))
  views$representation_id <- .mv13_residual_id
  views$view_id <- "cell_topology_v1"
  views$selected_cells <- 384L
  views$pca_components <- 30L
  views$point_metric <- "euclidean_shared_pca_v1"
  views$outcome_label_state <- "closed"
  views$biological_outcomes_computed <- FALSE
  views <- views[, c(
    "view_order", "dataset_scope", "unit_id", "seed", "panel_id",
    "representation_id", "view_id", "selected_cells", "pca_components",
    "point_metric", "outcome_label_state", "biological_outcomes_computed"
  )]
  ph <- views[rep(seq_len(nrow(views)), each = 2L), , drop = FALSE]
  ph$homology_dimension <- rep(c("H0", "H1"), times = nrow(views))
  ph$ph_order <- seq_len(nrow(ph))
  ph <- ph[, c("ph_order", names(views), "homology_dimension")]
  pca <- unique(views[, c("dataset_scope", "seed", "panel_id",
                          "representation_id", "pca_components")])
  pca <- pca[order(pca$dataset_scope, pca$panel_id, pca$seed,
                   method = "radix"), , drop = FALSE]
  rownames(pca) <- NULL
  pca$pca_order <- seq_len(nrow(pca))
  pca <- pca[, c("pca_order", setdiff(names(pca), "pca_order"))]
  list(pca = pca, views = views, ph = ph)
}

mv13_validate_residual_cache_v1 <- function(cache, expected_unit = NULL,
                                             expected_scope = NULL) {
  required <- c("contract_id", "identity", "panels", "representations",
                "payload_sha256", "cache_key", "downstream_execution")
  if (!is.list(cache) || !all(required %in% names(cache)) ||
      !identical(cache$contract_id, "mv08o_residual_source_cache_v1") ||
      !is.data.frame(cache$panels) || nrow(cache$panels) != 500L ||
      !identical(as.integer(cache$panels$panel_order), seq_len(500L)) ||
      anyDuplicated(cache$panels$feature_id) ||
      !identical(cache$identity$outcome_label_state, "closed") ||
      isTRUE(cache$identity$biological_outcomes_computed) ||
      any(unlist(cache$downstream_execution, use.names = FALSE) != 0)) {
    stop("MV13 source cache violates the closed MV8-O/Q contract.",
         call. = FALSE)
  }
  if (!is.null(expected_unit) &&
      !identical(cache$identity$unit_id, as.character(expected_unit))) {
    stop("MV13 source cache unit identity drift.", call. = FALSE)
  }
  if (!is.null(expected_scope) &&
      !identical(cache$identity$dataset_scope, as.character(expected_scope))) {
    stop("MV13 source cache scope drift.", call. = FALSE)
  }
  expected_seeds <- if (cache$identity$dataset_scope == "internal124") {
    as.character(.mv13_seeds)
  } else if (cache$identity$dataset_scope == "external8") {
    as.character(.mv13_seeds[[1L]])
  } else stop("MV13 cache has an unsupported dataset scope.", call. = FALSE)
  if (!identical(sort(names(cache$representations)), sort(expected_seeds))) {
    stop("MV13 cache seed axis drift.", call. = FALSE)
  }
  for (seed in expected_seeds) {
    matrix <- cache$representations[[seed]][[.mv13_residual_id]]
    if (!is.matrix(matrix) || !is.numeric(matrix) ||
        !identical(dim(matrix), c(500L, 384L)) ||
        any(!is.finite(matrix)) || anyDuplicated(rownames(matrix)) ||
        anyDuplicated(colnames(matrix)) ||
        !identical(rownames(matrix), cache$panels$feature_id)) {
      stop("MV13 residual matrix is malformed at seed ", seed, ".",
           call. = FALSE)
    }
  }
  invisible(TRUE)
}

mv13_source_from_cache_v1 <- function(cache, seed, panel,
                                       common475_features = NULL) {
  mv13_validate_residual_cache_v1(cache)
  seed <- as.character(as.integer(seed))
  if (!(seed %in% names(cache$representations))) {
    stop("MV13 requested seed is absent from the cache.", call. = FALSE)
  }
  panel <- match.arg(panel, c("exact500", "common475"))
  matrix <- cache$representations[[seed]][[.mv13_residual_id]]
  profile <- "scientific"
  if (panel == "common475") {
    common475_features <- as.character(common475_features)
    if (length(common475_features) != 475L || anyNA(common475_features) ||
        anyDuplicated(common475_features) ||
        !all(common475_features %in% rownames(matrix))) {
      stop("MV13 common475 feature axis is invalid.", call. = FALSE)
    }
    matrix <- matrix[common475_features, , drop = FALSE]
    profile <- "scientific_common475"
  }
  new_dual_view_source(
    matrix, sample_id = cache$identity$unit_id,
    cohort = cache$identity$dataset_scope,
    representation = .mv13_residual_id,
    fit_scope_id = "all_qc_per_unit_sct_fit_selected384",
    subsample_seed = as.integer(seed),
    standardization_id = "pearson_residual_all_qc_fit_v1",
    contract_profile = profile
  )
}

mv13_fit_shared_cell_pca_v1 <- function(caches, seed, panel,
                                         common475_features = NULL) {
  if (!is.list(caches) || length(caches) < 2L) {
    stop("MV13 PCA requires a multi-unit cache list.", call. = FALSE)
  }
  sources <- lapply(caches, mv13_source_from_cache_v1, seed = seed,
                    panel = panel, common475_features = common475_features)
  fit_cell_topology_pca(sources, n_components = 30L,
                        pca_seed = as.integer(seed))
}

mv13_load_group_sources_v1 <- function(locator, dataset_scope, seed, panel,
                                        common475_features = NULL,
                                        verify_hash = NULL) {
  required <- c("dataset_scope", "unit_id", "private_cache_path",
                "cache_sha256")
  if (!is.data.frame(locator) || !all(required %in% names(locator))) {
    stop("MV13 private locator schema drift.", call. = FALSE)
  }
  dataset_scope <- match.arg(dataset_scope, c("internal124", "external8"))
  panel <- match.arg(panel, c("exact500", "common475"))
  rows <- locator[locator$dataset_scope == dataset_scope, , drop = FALSE]
  rows <- rows[order(rows$unit_id, method = "radix"), , drop = FALSE]
  expected <- if (dataset_scope == "internal124") 124L else 8L
  if (nrow(rows) != expected || anyDuplicated(rows$unit_id)) {
    stop("MV13 group locator population drift.", call. = FALSE)
  }
  lapply(seq_len(nrow(rows)), function(i) {
    path <- rows$private_cache_path[[i]]
    if (!is.null(verify_hash) &&
        !identical(tolower(verify_hash(path)),
                   tolower(rows$cache_sha256[[i]]))) {
      stop("MV13 group cache hash drift at unit order ", i, ".",
           call. = FALSE)
    }
    cache <- readRDS(path)
    mv13_validate_residual_cache_v1(cache, rows$unit_id[[i]], dataset_scope)
    mv13_source_from_cache_v1(
      cache, seed = seed, panel = panel,
      common475_features = common475_features
    )
  })
}

mv13_compute_cell_group_v1 <- function(
    sources, dataset_scope, seed, panel, model = NULL,
    adopted_unit = NULL, adopted_view = NULL, adopted_result = NULL) {
  dataset_scope <- match.arg(dataset_scope, c("internal124", "external8"))
  panel <- match.arg(panel, c("exact500", "common475"))
  expected <- if (dataset_scope == "internal124") 124L else 8L
  if (!is.list(sources) || length(sources) != expected) {
    stop("MV13 group source population drift.", call. = FALSE)
  }
  sample_ids <- vapply(sources, `[[`, character(1L), "sample_id")
  if (anyDuplicated(sample_ids)) stop("MV13 group sample IDs are duplicated.")
  if (is.null(model)) {
    model <- fit_cell_topology_pca(
      sources, n_components = 30L, pca_seed = as.integer(seed)
    )
    model_origin <- "new_fit"
  } else {
    .validate_cell_pca_model(model)
    if (!identical(sort(model$fit_sample_ids), sort(sample_ids)) ||
        model$pca_seed != as.integer(seed) || model$n_components != 30L) {
      stop("MV13 adopted PCA model axis drift.", call. = FALSE)
    }
    model_origin <- "independently_closed_adoption"
  }
  adopting <- !is.null(adopted_unit)
  if (adopting && (is.null(adopted_view) || is.null(adopted_result) ||
                   !(adopted_unit %in% sample_ids))) {
    stop("MV13 adopted sentinel evidence is incomplete.", call. = FALSE)
  }
  records <- vector("list", length(sources)); names(records) <- sample_ids
  for (i in seq_along(sources)) {
    unit <- sample_ids[[i]]
    if (adopting && identical(unit, adopted_unit)) {
      validate_topology_view(adopted_view)
      if (adopted_view$sample_id != unit ||
          adopted_view$transformations$pca_model_cache_key != model$cache_key ||
          adopted_result$provenance$view_cache_key != adopted_view$cache_key) {
        stop("MV13 adopted sentinel identity drift.", call. = FALSE)
      }
      view <- adopted_view; result <- adopted_result
      record_origin <- "independently_closed_adoption"
    } else {
      view <- construct_cell_topology_view(
        sources[[i]], model, n_components = 30L
      )
      result <- run_topology_view_ph(
        view, max_dim = 1L, threshold = -1, field = 2L
      )
      record_origin <- "new_computation"
    }
    oracle <- mv07g_validate_ph_against_view_v1(result, view)
    if (!isTRUE(oracle$passed)) stop("MV13 H0 MST oracle failed.")
    records[[i]] <- list(
      unit_id = unit, origin = record_origin, view = view,
      result = result, oracle = oracle
    )
  }
  group_id <- paste(dataset_scope, panel, as.integer(seed), sep = "__")
  artifact <- list(
    contract_id = "mv13d_allqc_cell_group_v1", group_id = group_id,
    dataset_scope = dataset_scope, panel_id = panel,
    seed = as.integer(seed), representation_id = .mv13_residual_id,
    model_origin = model_origin, model = model, records = records,
    labels_used = FALSE, outcomes_used = FALSE, downstream_jobs = 0L
  )
  class(artifact) <- c("scph_mv13d_cell_group_v1", "list")
  mv13_validate_cell_group_v1(artifact)
  artifact
}

mv13_validate_cell_group_v1 <- function(artifact) {
  required <- c("contract_id", "group_id", "dataset_scope", "panel_id",
                "seed", "representation_id", "model_origin", "model",
                "records", "labels_used", "outcomes_used", "downstream_jobs")
  expected <- if (identical(artifact$dataset_scope, "internal124")) 124L else 8L
  if (!is.list(artifact) || !all(required %in% names(artifact)) ||
      !identical(artifact$contract_id, "mv13d_allqc_cell_group_v1") ||
      !(artifact$dataset_scope %in% c("internal124", "external8")) ||
      !(artifact$panel_id %in% c("exact500", "common475")) ||
      artifact$representation_id != .mv13_residual_id ||
      length(artifact$records) != expected || artifact$labels_used ||
      artifact$outcomes_used || artifact$downstream_jobs != 0L) {
    stop("MV13 cell group violates its closed contract.", call. = FALSE)
  }
  .validate_cell_pca_model(artifact$model)
  units <- vapply(artifact$records, `[[`, character(1L), "unit_id")
  if (anyDuplicated(units) ||
      !identical(sort(unname(units), method = "radix"),
                 sort(unname(artifact$model$fit_sample_ids), method = "radix"))) {
    stop("MV13 cell group unit/model axis drift.", call. = FALSE)
  }
  for (record in artifact$records) {
    validate_topology_view(record$view)
    if (!isTRUE(record$oracle$passed) ||
        record$result$provenance$view_cache_key != record$view$cache_key ||
        record$view$sample_id != record$unit_id) {
      stop("MV13 cell group record drift.", call. = FALSE)
    }
  }
  invisible(TRUE)
}
