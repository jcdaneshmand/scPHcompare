.dual_view_contract_version <- "1.0.0"
.dual_view_scientific_shape <- c(genes = 500L, cells = 384L, pcs = 30L)

.one_nonempty_string <- function(value, label) {
  if (length(value) != 1L || is.na(value) || !is.character(value) ||
      !nzchar(trimws(value))) {
    stop(label, " must be one non-empty string.", call. = FALSE)
  }
  trimws(value)
}

.one_integer <- function(value, label, minimum = 1L) {
  if (length(value) != 1L || is.na(value) || !is.numeric(value) ||
      !is.finite(value) || value != as.integer(value) || value < minimum) {
    stop(label, " must be one integer of at least ", minimum, ".", call. = FALSE)
  }
  as.integer(value)
}

.validate_named_numeric_matrix <- function(x, label) {
  if (!is.matrix(x) || !is.numeric(x) || length(dim(x)) != 2L) {
    stop(label, " must be a numeric matrix.", call. = FALSE)
  }
  if (any(!is.finite(x))) {
    stop(label, " must contain only finite values.", call. = FALSE)
  }
  axes <- dimnames(x)
  if (is.null(axes) || is.null(axes[[1L]]) || is.null(axes[[2L]])) {
    stop(label, " must have explicit row and column identifiers.", call. = FALSE)
  }
  for (index in seq_len(2L)) {
    ids <- as.character(axes[[index]])
    axis_label <- if (index == 1L) "row" else "column"
    if (length(ids) != dim(x)[[index]] || anyNA(ids) || any(!nzchar(ids)) ||
        anyDuplicated(ids)) {
      stop(label, " has missing, empty, or duplicated ", axis_label,
           " identifiers.", call. = FALSE)
    }
  }
  x
}

.resolve_dual_view_shape <- function(contract_profile, expected_genes,
                                     expected_cells, expected_pcs) {
  contract_profile <- match.arg(
    contract_profile, c("scientific", "analytical_fixture")
  )
  defaults <- .dual_view_scientific_shape
  supplied <- c(
    genes = if (is.null(expected_genes)) defaults[["genes"]] else expected_genes,
    cells = if (is.null(expected_cells)) defaults[["cells"]] else expected_cells,
    pcs = if (is.null(expected_pcs)) defaults[["pcs"]] else expected_pcs
  )
  supplied <- c(
    genes = .one_integer(supplied[["genes"]], "expected_genes", 2L),
    cells = .one_integer(supplied[["cells"]], "expected_cells", 3L),
    pcs = .one_integer(supplied[["pcs"]], "expected_pcs", 1L)
  )
  if (supplied[["pcs"]] >= min(supplied[["genes"]], supplied[["cells"]])) {
    stop("expected_pcs must be smaller than both point-axis dimensions.",
         call. = FALSE)
  }
  if (contract_profile == "scientific" && !identical(supplied, defaults)) {
    stop(
      "The scientific contract requires exactly 500 genes, 384 cells, and 30 PCs; ",
      "use contract_profile = 'analytical_fixture' for reduced tests.",
      call. = FALSE
    )
  }
  list(
    profile = contract_profile,
    version = .dual_view_contract_version,
    expected_genes = unname(supplied[["genes"]]),
    expected_cells = unname(supplied[["cells"]]),
    expected_pcs = unname(supplied[["pcs"]]),
    scientific_eligible = identical(contract_profile, "scientific")
  )
}

.scientific_digest <- function(value) {
  digest::digest(value, algo = "sha256", serialize = TRUE)
}

.dual_view_source_identity <- function(source) {
  list(
    object_type = "dual_view_source",
    contract_version = .dual_view_contract_version,
    contract_profile = source$contract$profile,
    axis_role = "genes_by_cells",
    genes = source$gene_ids,
    cells = source$cell_ids,
    subsample_seed = source$subsample_seed,
    matrix_sha256 = source$input_sha256,
    sample_id = source$sample_id,
    cohort = source$cohort,
    representation = source$representation,
    fit_scope_id = source$fit_scope_id,
    standardization_id = source$standardization_id
  )
}

.cell_projection_source_identity <- function(source) {
  list(
    object_type = "cell_projection_source",
    contract_version = .dual_view_contract_version,
    contract_profile = source$contract$profile,
    axis_role = "genes_by_cells",
    genes = source$gene_ids,
    cells = source$cell_ids,
    subsample_seed = source$subsample_seed,
    matrix_sha256 = source$input_sha256,
    sample_id = source$sample_id,
    cohort = source$cohort,
    representation = source$representation,
    fit_scope_id = source$fit_scope_id,
    standardization_id = source$standardization_id
  )
}

#' Select a deterministic matched-cell subset
#'
#' Sorts and validates cell identifiers before sampling without replacement,
#' restores the caller's random-number state, and returns the selected IDs in
#' canonical sorted order. The same result can therefore be reused across
#' representations and topology views.
#'
#' @param cell_ids Unique non-empty cell identifiers.
#' @param n Number of cells to select.
#' @param seed One finite integer-compatible seed.
#'
#' @return A sorted character vector with attributes recording the seed, source
#'   count, and SHA-256 digest.
#' @export
select_matched_cells <- function(cell_ids, n = 384L, seed = 20260805L) {
  cell_ids <- as.character(cell_ids)
  if (anyNA(cell_ids) || any(!nzchar(cell_ids)) || anyDuplicated(cell_ids)) {
    stop("cell_ids must be unique, non-missing, and non-empty.", call. = FALSE)
  }
  n <- .one_integer(n, "n", 1L)
  seed <- .one_integer(seed, "seed", 0L)
  if (length(cell_ids) < n) {
    stop("cell_ids contains fewer than n eligible cells.", call. = FALSE)
  }
  eligible <- sort(cell_ids, method = "radix")
  selected <- .with_preserved_seed(seed, {
    sort(sample(eligible, size = n, replace = FALSE), method = "radix")
  })
  attr(selected, "subsample_seed") <- seed
  attr(selected, "eligible_cell_count") <- length(eligible)
  attr(selected, "selected_cell_sha256") <- .scientific_digest(selected)
  selected
}

#' Construct a typed dual-view source matrix
#'
#' Wraps an already standardized genes-by-cells matrix in the frozen MV-01
#' identity and provenance contract. This constructor deliberately does not
#' infer axis roles from dimensions and does not perform representation-specific
#' preprocessing. Those fitted transformations must be completed and recorded
#' before this boundary.
#'
#' @param x Numeric matrix with genes in rows and selected cells in columns.
#' @param sample_id,cohort,representation,fit_scope_id Non-empty provenance IDs.
#' @param subsample_seed Integer-compatible cell-subsampling seed.
#' @param standardization_id Stable ID or digest for the fitted scaling operation.
#' @param contract_profile `"scientific"` or `"analytical_fixture"`.
#' @param expected_genes,expected_cells,expected_pcs Expected shape. Overrides
#'   are allowed only for analytical fixtures.
#'
#' @return A versioned `scph_dual_view_source_v1` object.
#' @export
new_dual_view_source <- function(x, sample_id, cohort, representation,
                                 fit_scope_id, subsample_seed,
                                 standardization_id,
                                 contract_profile = "scientific",
                                 expected_genes = NULL,
                                 expected_cells = NULL,
                                 expected_pcs = NULL) {
  x <- .validate_named_numeric_matrix(x, "x")
  shape <- .resolve_dual_view_shape(
    contract_profile, expected_genes, expected_cells, expected_pcs
  )
  if (nrow(x) != shape$expected_genes || ncol(x) != shape$expected_cells) {
    stop(
      "x must be genes_by_cells with shape ", shape$expected_genes, " x ",
      shape$expected_cells, "; observed ", nrow(x), " x ", ncol(x), ".",
      call. = FALSE
    )
  }
  gene_sd <- apply(x, 1L, stats::sd)
  if (any(!is.finite(gene_sd)) ||
      any(gene_sd <= sqrt(.Machine$double.eps))) {
    stop("x contains nonfinite or effectively zero-variance genes.",
         call. = FALSE)
  }
  identifiers <- list(
    sample_id = .one_nonempty_string(sample_id, "sample_id"),
    cohort = .one_nonempty_string(cohort, "cohort"),
    representation = .one_nonempty_string(representation, "representation"),
    fit_scope_id = .one_nonempty_string(fit_scope_id, "fit_scope_id"),
    standardization_id = .one_nonempty_string(
      standardization_id, "standardization_id"
    )
  )
  subsample_seed <- .one_integer(subsample_seed, "subsample_seed", 0L)
  result <- structure(
    list(
      matrix = x,
      contract = shape,
      axis_roles = c(rows = "genes", columns = "cells"),
      gene_ids = rownames(x),
      cell_ids = colnames(x),
      sample_id = identifiers$sample_id,
      cohort = identifiers$cohort,
      representation = identifiers$representation,
      fit_scope_id = identifiers$fit_scope_id,
      subsample_seed = subsample_seed,
      standardization_id = identifiers$standardization_id,
      input_sha256 = .scientific_digest(x),
      cache_key = NULL
    ),
    class = c("scph_dual_view_source_v1", "scph_dual_view_source")
  )
  result$cache_key <- paste0(
    "dual_view_source_v1:",
    .scientific_digest(.dual_view_source_identity(result))
  )
  result
}

#' Construct a typed source for held-out cell projection
#'
#' This source has the same scientific shape and immutable provenance as a
#' dual-view source but permits genes that are constant within a held-out
#' sample. Such genes remain valid inputs to a training-fitted PCA projection;
#' the source is intentionally ineligible for gene topology.
#'
#' @inheritParams new_dual_view_source
#'
#' @return A versioned `scph_cell_projection_source_v1` object.
#' @export
new_cell_projection_source <- function(
    x, sample_id, cohort, representation, fit_scope_id, subsample_seed,
    standardization_id, contract_profile = "scientific",
    expected_genes = NULL, expected_cells = NULL, expected_pcs = NULL) {
  x <- .validate_named_numeric_matrix(x, "x")
  shape <- .resolve_dual_view_shape(
    contract_profile, expected_genes, expected_cells, expected_pcs
  )
  if (nrow(x) != shape$expected_genes || ncol(x) != shape$expected_cells) {
    stop(
      "x must be genes_by_cells with shape ", shape$expected_genes, " x ",
      shape$expected_cells, "; observed ", nrow(x), " x ", ncol(x), ".",
      call. = FALSE
    )
  }
  result <- structure(
    list(
      matrix = x, contract = shape,
      axis_roles = c(rows = "genes", columns = "cells"),
      gene_ids = rownames(x), cell_ids = colnames(x),
      sample_id = .one_nonempty_string(sample_id, "sample_id"),
      cohort = .one_nonempty_string(cohort, "cohort"),
      representation = .one_nonempty_string(representation, "representation"),
      fit_scope_id = .one_nonempty_string(fit_scope_id, "fit_scope_id"),
      subsample_seed = .one_integer(subsample_seed, "subsample_seed", 0L),
      standardization_id = .one_nonempty_string(
        standardization_id, "standardization_id"
      ),
      input_sha256 = .scientific_digest(x), cache_key = NULL
    ),
    class = c("scph_cell_projection_source_v1", "scph_cell_projection_source")
  )
  result$cache_key <- paste0(
    "cell_projection_source_v1:",
    .scientific_digest(.cell_projection_source_identity(result))
  )
  result
}

.validate_cell_projection_source <- function(source) {
  if (!inherits(source, "scph_cell_projection_source_v1") ||
      !is.list(source)) {
    stop("source must be a scph_cell_projection_source_v1 object.",
         call. = FALSE)
  }
  x <- .validate_named_numeric_matrix(source$matrix, "source$matrix")
  if (!identical(source$axis_roles, c(rows = "genes", columns = "cells")) ||
      !identical(rownames(x), source$gene_ids) ||
      !identical(colnames(x), source$cell_ids) ||
      nrow(x) != source$contract$expected_genes ||
      ncol(x) != source$contract$expected_cells ||
      !identical(.scientific_digest(x), source$input_sha256)) {
    stop("Cell-projection source payload or axes are stale.", call. = FALSE)
  }
  expected <- paste0(
    "cell_projection_source_v1:",
    .scientific_digest(.cell_projection_source_identity(source))
  )
  if (!identical(expected, source$cache_key)) {
    stop("Cell-projection source provenance is stale.", call. = FALSE)
  }
  source
}

.validate_dual_view_source <- function(source) {
  if (!inherits(source, "scph_dual_view_source_v1") || !is.list(source)) {
    stop("source must be a scph_dual_view_source_v1 object.", call. = FALSE)
  }
  x <- .validate_named_numeric_matrix(source$matrix, "source$matrix")
  if (!identical(source$axis_roles, c(rows = "genes", columns = "cells")) ||
      !identical(rownames(x), source$gene_ids) ||
      !identical(colnames(x), source$cell_ids)) {
    stop("source axis roles or identifiers are inconsistent.", call. = FALSE)
  }
  if (nrow(x) != source$contract$expected_genes ||
      ncol(x) != source$contract$expected_cells) {
    stop("source dimensions no longer match its contract.", call. = FALSE)
  }
  expected_input <- .scientific_digest(x)
  if (!identical(expected_input, source$input_sha256)) {
    stop("source matrix digest is stale or has been modified.", call. = FALSE)
  }
  expected_cache <- paste0(
    "dual_view_source_v1:",
    .scientific_digest(.dual_view_source_identity(source))
  )
  if (!identical(expected_cache, source$cache_key)) {
    stop("source provenance or cache identity is stale.", call. = FALSE)
  }
  source
}

#' Fit the shared PCA model for cell topology
#'
#' Fits deterministic `stats::prcomp` scores in the conventional cells-by-genes
#' orientation over equally sized, typed source objects. The model is shared by
#' every sample in one cohort, representation, fit scope, and repetition.
#'
#' @param sources Named list of compatible `scph_dual_view_source_v1` objects.
#' @param n_components Number of PCs; defaults to the source contract.
#' @param pca_seed Recorded integer-compatible seed. The selected SVD algorithm
#'   is deterministic, but the seed remains part of provenance/cache identity.
#'
#' @return A versioned `scph_cell_pca_model_v1` object.
#' @export
fit_cell_topology_pca <- function(sources, n_components = NULL,
                                  pca_seed = 20260805L) {
  if (!is.list(sources) || length(sources) == 0L) {
    stop("sources must be a non-empty list.", call. = FALSE)
  }
  sources <- lapply(sources, .validate_dual_view_source)
  sample_ids <- vapply(sources, `[[`, character(1L), "sample_id")
  if (anyDuplicated(sample_ids)) {
    stop("sources must have unique sample IDs.", call. = FALSE)
  }
  canonical_order <- order(sample_ids, method = "radix")
  sources <- sources[canonical_order]
  sample_ids <- sample_ids[canonical_order]
  reference <- sources[[1L]]
  comparable_fields <- c("cohort", "representation", "fit_scope_id", "subsample_seed")
  for (field in comparable_fields) {
    values <- vapply(sources, function(source) as.character(source[[field]]), character(1L))
    if (length(unique(values)) != 1L) {
      stop("sources must share ", field, ".", call. = FALSE)
    }
  }
  if (any(!vapply(sources, function(source) {
    identical(source$gene_ids, reference$gene_ids) &&
      identical(source$contract, reference$contract)
  }, logical(1L)))) {
    stop("sources must share the exact ordered genes and contract.", call. = FALSE)
  }
  n_components <- if (is.null(n_components)) {
    reference$contract$expected_pcs
  } else {
    .one_integer(n_components, "n_components", 1L)
  }
  if (n_components > reference$contract$expected_pcs) {
    stop("n_components exceeds the source contract.", call. = FALSE)
  }
  pca_seed <- .one_integer(pca_seed, "pca_seed", 0L)
  pooled <- do.call(rbind, lapply(sources, function(source) {
    values <- t(source$matrix)
    rownames(values) <- paste(source$sample_id, rownames(values), sep = "::")
    values
  }))
  fit <- .with_preserved_seed(pca_seed, stats::prcomp(
    pooled, center = FALSE, scale. = FALSE, rank. = n_components
  ))
  if (ncol(fit$rotation) < n_components || length(fit$sdev) < n_components ||
      any(!is.finite(fit$sdev[seq_len(n_components)])) ||
      any(fit$sdev[seq_len(n_components)] <= sqrt(.Machine$double.eps))) {
    stop("Fewer than the requested number of usable PCs were fitted.", call. = FALSE)
  }
  rotation <- fit$rotation[, seq_len(n_components), drop = FALSE]
  colnames(rotation) <- paste0("PC", seq_len(n_components))
  identity <- list(
    object_type = "cell_pca_model",
    contract_version = .dual_view_contract_version,
    contract_profile = reference$contract$profile,
    source_cache_keys = vapply(sources, `[[`, character(1L), "cache_key"),
    sample_ids = sample_ids,
    genes = rownames(rotation),
    coordinates = colnames(rotation),
    n_components = n_components,
    pca_seed = pca_seed,
    algorithm = "stats_prcomp_svd_v1",
    rotation = unname(rotation),
    singular_values = fit$sdev[seq_len(n_components)]
  )
  structure(
    list(
      view_id = "cell_topology_v1",
      contract = reference$contract,
      cohort = reference$cohort,
      representation = reference$representation,
      fit_scope_id = reference$fit_scope_id,
      subsample_seed = reference$subsample_seed,
      fit_sample_ids = sample_ids,
      source_cache_keys = vapply(sources, `[[`, character(1L), "cache_key"),
      gene_ids = rownames(rotation),
      rotation = rotation,
      singular_values = fit$sdev[seq_len(n_components)],
      n_components = n_components,
      pca_seed = pca_seed,
      algorithm = "stats_prcomp_svd_v1",
      cache_key = paste0("cell_pca_model_v1:", .scientific_digest(identity))
    ),
    class = c("scph_cell_pca_model_v1", "scph_topology_fit")
  )
}

.validate_cell_pca_model <- function(pca_model) {
  if (!inherits(pca_model, "scph_cell_pca_model_v1") ||
      !is.matrix(pca_model$rotation) ||
      !identical(rownames(pca_model$rotation), pca_model$gene_ids)) {
    stop("pca_model must be a valid scph_cell_pca_model_v1 object.",
         call. = FALSE)
  }
  identity <- list(
    object_type = "cell_pca_model",
    contract_version = .dual_view_contract_version,
    contract_profile = pca_model$contract$profile,
    source_cache_keys = pca_model$source_cache_keys,
    sample_ids = pca_model$fit_sample_ids,
    genes = pca_model$gene_ids,
    coordinates = colnames(pca_model$rotation),
    n_components = pca_model$n_components,
    pca_seed = pca_model$pca_seed,
    algorithm = pca_model$algorithm,
    rotation = unname(pca_model$rotation),
    singular_values = pca_model$singular_values
  )
  expected <- paste0("cell_pca_model_v1:", .scientific_digest(identity))
  if (!identical(expected, pca_model$cache_key)) {
    stop("pca_model payload, provenance, or cache identity is stale.",
         call. = FALSE)
  }
  pca_model
}

.topology_view_identity <- function(view) {
  list(
    object_type = "topology_view",
    view_id = view$view_id,
    contract_version = view$contract_version,
    contract_profile = view$contract_profile,
    source_cache_key = view$source_cache_key,
    point_metric = view$point_metric,
    point_ids = view$point_ids,
    coordinate_ids = view$coordinate_ids,
    transformations = view$transformations,
    payload_sha256 = view$payload_sha256
  )
}

.new_topology_view <- function(view_id, source, point_metric, payload,
                               point_ids, coordinate_ids, transformations,
                               payload_sha256, diagnostics = list()) {
  identity <- list(
    object_type = "topology_view",
    view_id = view_id,
    contract_version = .dual_view_contract_version,
    contract_profile = source$contract$profile,
    source_cache_key = source$cache_key,
    point_metric = point_metric,
    point_ids = point_ids,
    coordinate_ids = coordinate_ids,
    transformations = transformations,
    payload_sha256 = payload_sha256
  )
  class_name <- if (view_id == "cell_topology_v1") {
    "scph_cell_topology_view_v1"
  } else {
    "scph_gene_topology_view_v1"
  }
  structure(
    list(
      view_id = view_id,
      contract_version = .dual_view_contract_version,
      contract_profile = source$contract$profile,
      scientific_eligible = source$contract$scientific_eligible,
      source_cache_key = source$cache_key,
      sample_id = source$sample_id,
      cohort = source$cohort,
      representation = source$representation,
      fit_scope_id = source$fit_scope_id,
      subsample_seed = source$subsample_seed,
      input_axis_role = "genes_by_cells",
      point_axis_role = if (view_id == "cell_topology_v1") "cells" else "genes",
      coordinate_axis_role = if (view_id == "cell_topology_v1") "shared_pcs" else "matched_cells",
      point_metric = point_metric,
      point_ids = point_ids,
      coordinate_ids = coordinate_ids,
      transformations = transformations,
      payload = payload,
      payload_sha256 = payload_sha256,
      diagnostics = diagnostics,
      cache_key = paste0(view_id, ":", .scientific_digest(identity))
    ),
    class = c(class_name, "scph_topology_view_v1", "scph_topology_view")
  )
}

#' Construct the corrected cell-topology view
#'
#' Projects cells into an explicitly fitted shared PCA model and returns a typed
#' cells-by-PC point cloud. Bare assay matrices are not accepted.
#'
#' @param source A typed dual-view source.
#' @param pca_model A compatible model from `fit_cell_topology_pca()`.
#' @param n_components Number of leading fitted PCs to use.
#'
#' @return A versioned `scph_cell_topology_view_v1` object.
#' @export
construct_cell_topology_view <- function(source, pca_model,
                                         n_components = NULL) {
  source <- .validate_dual_view_source(source)
  pca_model <- .validate_cell_pca_model(pca_model)
  if (!identical(source$gene_ids, pca_model$gene_ids) ||
      !identical(source$cohort, pca_model$cohort) ||
      !identical(source$representation, pca_model$representation) ||
      !identical(source$fit_scope_id, pca_model$fit_scope_id) ||
      !identical(source$subsample_seed, pca_model$subsample_seed)) {
    stop("source is incompatible with the shared PCA model.", call. = FALSE)
  }
  n_components <- if (is.null(n_components)) pca_model$n_components else {
    .one_integer(n_components, "n_components", 1L)
  }
  if (n_components > pca_model$n_components) {
    stop("n_components exceeds the fitted PCA model.", call. = FALSE)
  }
  rotation <- pca_model$rotation[, seq_len(n_components), drop = FALSE]
  scores <- t(source$matrix) %*% rotation
  rownames(scores) <- source$cell_ids
  colnames(scores) <- colnames(rotation)
  if (any(!is.finite(scores))) {
    stop("Cell PCA scores contain nonfinite values.", call. = FALSE)
  }
  zero_norm <- sqrt(rowSums(scores ^ 2)) <= sqrt(.Machine$double.eps)
  if (any(zero_norm)) {
    stop("Cell PCA scores contain one or more zero-norm cells.", call. = FALSE)
  }
  .new_topology_view(
    view_id = "cell_topology_v1",
    source = source,
    point_metric = "euclidean_shared_pca_v1",
    payload = scores,
    point_ids = rownames(scores),
    coordinate_ids = colnames(scores),
    transformations = list(
      standardization_id = source$standardization_id,
      pca_model_cache_key = pca_model$cache_key,
      pca_algorithm = pca_model$algorithm,
      pca_seed = pca_model$pca_seed,
      n_components = n_components
    ),
    payload_sha256 = .scientific_digest(scores),
    diagnostics = list(duplicated_point_rows = sum(duplicated(scores)))
  )
}

#' Construct a cell-topology view from frozen shared coordinates
#'
#' Wraps cell-by-coordinate embeddings produced by a separately validated
#' training-reference transform (for example, inductive query mapping). The
#' ordered cell IDs must exactly match the typed source, and the fit/cache
#' identities are included in the immutable topology-view provenance.
#'
#' @param source A typed dual-view source defining the sample and cell identity.
#' @param coordinates Named cells-by-coordinates numeric matrix.
#' @param coordinate_contract_id Non-empty identifier for the coordinate method.
#' @param coordinate_fit_cache_key Non-empty immutable fit or mapping cache key.
#'
#' @return A versioned `scph_cell_topology_view_v1` object.
#' @export
construct_frozen_cell_topology_view <- function(
    source, coordinates, coordinate_contract_id, coordinate_fit_cache_key) {
  source <- if (inherits(source, "scph_cell_projection_source_v1")) {
    .validate_cell_projection_source(source)
  } else {
    .validate_dual_view_source(source)
  }
  coordinates <- .validate_named_numeric_matrix(coordinates, "coordinates")
  coordinate_contract_id <- .one_nonempty_string(
    coordinate_contract_id, "coordinate_contract_id"
  )
  coordinate_fit_cache_key <- .one_nonempty_string(
    coordinate_fit_cache_key, "coordinate_fit_cache_key"
  )
  if (!identical(rownames(coordinates), source$cell_ids)) {
    stop("coordinates must have the source's exact ordered cell IDs.",
         call. = FALSE)
  }
  zero_norm <- sqrt(rowSums(coordinates ^ 2)) <= sqrt(.Machine$double.eps)
  if (any(zero_norm)) {
    stop("Frozen cell coordinates contain one or more zero-norm cells.",
         call. = FALSE)
  }
  .new_topology_view(
    view_id = "cell_topology_v1",
    source = source,
    point_metric = "euclidean_frozen_shared_coordinates_v1",
    payload = coordinates,
    point_ids = rownames(coordinates),
    coordinate_ids = colnames(coordinates),
    transformations = list(
      standardization_id = source$standardization_id,
      coordinate_contract_id = coordinate_contract_id,
      coordinate_fit_cache_key = coordinate_fit_cache_key,
      n_components = ncol(coordinates)
    ),
    payload_sha256 = .scientific_digest(coordinates),
    diagnostics = list(duplicated_point_rows = sum(duplicated(coordinates)))
  )
}

#' Construct the intentional gene-topology view
#'
#' Centers and unit-normalizes each gene over matched cells, then constructs an
#' explicit named Euclidean chord distance object equal to
#' `sqrt(2 * (1 - Pearson correlation))`.
#'
#' @param source A typed dual-view source.
#'
#' @return A versioned `scph_gene_topology_view_v1` object.
#' @export
construct_gene_topology_view <- function(source) {
  source <- .validate_dual_view_source(source)
  centered <- source$matrix - rowMeans(source$matrix)
  norms <- sqrt(rowSums(centered ^ 2))
  if (any(!is.finite(norms)) || any(norms <= sqrt(.Machine$double.eps))) {
    failed <- source$gene_ids[!is.finite(norms) |
                                norms <= sqrt(.Machine$double.eps)]
    stop(
      "gene_topology_v1 requires nonconstant finite genes; failed: ",
      paste(utils::head(failed, 10L), collapse = ", "),
      call. = FALSE
    )
  }
  unit_vectors <- centered / norms
  rownames(unit_vectors) <- source$gene_ids
  colnames(unit_vectors) <- source$cell_ids
  distances <- stats::dist(unit_vectors, method = "euclidean")
  attr(distances, "Labels") <- source$gene_ids
  .new_topology_view(
    view_id = "gene_topology_v1",
    source = source,
    point_metric = "pearson_correlation_chord_v1",
    payload = distances,
    point_ids = source$gene_ids,
    coordinate_ids = source$cell_ids,
    transformations = list(
      standardization_id = source$standardization_id,
      row_centering = "within_sample_gene_mean",
      row_normalization = "l2_unit_norm",
      distance_representation = "explicit_dist"
    ),
    payload_sha256 = .scientific_digest(distances),
    diagnostics = list(
      minimum_gene_norm = min(norms),
      maximum_gene_norm = max(norms),
      duplicated_point_rows = sum(duplicated(unit_vectors))
    )
  )
}

#' Validate a corrected topology-view object
#'
#' @param view Object returned by a corrected view constructor.
#'
#' @return The validated object, invisibly suitable for internal composition.
#' @export
validate_topology_view <- function(view) {
  if (!inherits(view, "scph_topology_view_v1") || !is.list(view)) {
    stop(
      "view must be a typed cell_topology_v1 or gene_topology_v1 object; ",
      "bare matrices and dist objects are prohibited.",
      call. = FALSE
    )
  }
  allowed <- c("cell_topology_v1", "gene_topology_v1")
  if (!(view$view_id %in% allowed) ||
      !identical(view$contract_version, .dual_view_contract_version)) {
    stop("Unknown or incompatible topology-view contract.", call. = FALSE)
  }
  if (anyNA(view$point_ids) || any(!nzchar(view$point_ids)) ||
      anyDuplicated(view$point_ids)) {
    stop("view point identifiers are invalid.", call. = FALSE)
  }
  if (view$view_id == "cell_topology_v1") {
    payload <- .validate_named_numeric_matrix(view$payload, "view$payload")
    if (!identical(rownames(payload), view$point_ids) ||
        !identical(colnames(payload), view$coordinate_ids) ||
        !identical(view$point_axis_role, "cells") ||
        !identical(view$coordinate_axis_role, "shared_pcs")) {
      stop("cell topology payload axes are inconsistent.", call. = FALSE)
    }
  } else {
    if (!inherits(view$payload, "dist") ||
        !identical(attr(view$payload, "Labels"), view$point_ids) ||
        any(!is.finite(unclass(view$payload))) ||
        !identical(view$point_axis_role, "genes") ||
        !identical(view$coordinate_axis_role, "matched_cells")) {
      stop("gene topology distance payload is inconsistent.", call. = FALSE)
    }
  }
  if (!identical(.scientific_digest(view$payload), view$payload_sha256)) {
    stop("topology-view payload digest is stale or has been modified.", call. = FALSE)
  }
  expected_cache <- paste0(view$view_id, ":", .scientific_digest(
    .topology_view_identity(view)
  ))
  if (!identical(expected_cache, view$cache_key)) {
    stop("topology-view provenance or cache identity is stale.", call. = FALSE)
  }
  invisible(view)
}

.as_diagram_matrix <- function(value) {
  result <- as.matrix(value)
  if (ncol(result) != 3L) {
    stop("ripserr returned an unexpected persistence-diagram schema.",
         call. = FALSE)
  }
  storage.mode(result) <- "double"
  colnames(result) <- c("dimension", "birth", "death")
  result
}

#' Run PH through a corrected typed topology view
#'
#' This is the only MV-02 corrected PH entry point. It refuses bare matrices and
#' bare distance objects, validates payload identity, and dispatches cells as a
#' named point cloud or genes as an explicit `dist` object.
#'
#' @param view Typed topology view.
#' @param max_dim Maximum homology dimension. The v1 contract requires `1`.
#' @param threshold Filtration threshold. The v1 primary contract requires `-1`.
#' @param field Prime coefficient field passed to `ripserr` as `p`.
#'
#' @return A versioned result containing the diagram and immutable provenance.
#' @export
run_topology_view_ph <- function(view, max_dim = 1L, threshold = -1,
                                 field = 2L) {
  validate_topology_view(view)
  max_dim <- .one_integer(max_dim, "max_dim", 0L)
  field <- .one_integer(field, "field", 2L)
  if (!identical(max_dim, 1L)) {
    stop("The v1 corrected PH contract requires max_dim = 1.", call. = FALSE)
  }
  if (length(threshold) != 1L || is.na(threshold) || !is.numeric(threshold) ||
      !is.finite(threshold) || !identical(as.numeric(threshold), -1)) {
    stop("The v1 corrected PH contract requires threshold = -1.", call. = FALSE)
  }
  raw <- ripserr::vietoris_rips(
    dataset = view$payload,
    max_dim = max_dim,
    threshold = threshold,
    p = field,
    return_format = "mat"
  )
  diagram <- .as_diagram_matrix(raw)
  has_essential_h0 <- any(
    diagram[, "dimension"] == 0 & is.infinite(diagram[, "death"])
  )
  essential_h0_added <- !has_essential_h0
  if (essential_h0_added) {
    diagram <- rbind(
      diagram,
      c(dimension = 0, birth = 0, death = Inf)
    )
  }
  invalid_intervals <- !is.finite(diagram[, "dimension"]) |
    !is.finite(diagram[, "birth"]) |
    is.na(diagram[, "death"]) |
    diagram[, "death"] < diagram[, "birth"]
  zero_persistence <- is.finite(diagram[, "death"]) &
    diagram[, "death"] == diagram[, "birth"]
  provenance <- list(
    result_contract_id = "corrected_topology_result_v1",
    view_id = view$view_id,
    view_cache_key = view$cache_key,
    contract_version = view$contract_version,
    contract_profile = view$contract_profile,
    scientific_eligible = view$scientific_eligible,
    sample_id = view$sample_id,
    point_axis_role = view$point_axis_role,
    coordinate_axis_role = view$coordinate_axis_role,
    point_metric = view$point_metric,
    point_count = length(view$point_ids),
    max_dim = max_dim,
    threshold = as.numeric(threshold),
    field = field,
    ph_engine = "ripserr",
    ph_engine_version = as.character(utils::packageVersion("ripserr")),
    essential_h0_added = essential_h0_added,
    essential_h0_count = sum(
      diagram[, "dimension"] == 0 & is.infinite(diagram[, "death"])
    ),
    finite_interval_count = sum(is.finite(diagram[, "death"])),
    infinite_interval_count = sum(is.infinite(diagram[, "death"])),
    zero_persistence_count = sum(zero_persistence),
    invalid_interval_count = sum(invalid_intervals),
    diagram_sha256 = .scientific_digest(diagram)
  )
  structure(
    list(
      diagram = diagram,
      provenance = provenance,
      cache_key = paste0(
        "corrected_topology_result_v1:", .scientific_digest(provenance)
      )
    ),
    class = c("scph_topology_result_v1", "scph_topology_result")
  )
}

#' Run the explicitly acknowledged historical matrix route
#'
#' This compatibility helper preserves the historical rows-as-points behavior.
#' It is not scientifically eligible for the corrected analysis and its result
#' contract cannot collide with corrected topology-view cache keys.
#'
#' @param x Numeric rows-by-coordinates matrix.
#' @param acknowledge_legacy Must be `TRUE` to confirm historical-only use.
#' @param max_dim,threshold,field Parameters passed to `ripserr`.
#'
#' @return A legacy-stamped persistence result.
#' @export
run_legacy_matrix_ph <- function(x, acknowledge_legacy = FALSE,
                                 max_dim = 1L, threshold = -1,
                                 field = 2L) {
  warning(
    "legacy_gene_view_v0 treats matrix rows as points and is ineligible for ",
    "corrected cell/gene topology analyses.",
    call. = FALSE
  )
  if (!isTRUE(acknowledge_legacy)) {
    stop("Set acknowledge_legacy = TRUE for historical reproduction only.",
         call. = FALSE)
  }
  x <- .validate_named_numeric_matrix(x, "x")
  raw <- ripserr::vietoris_rips(
    dataset = x, max_dim = max_dim, threshold = threshold, p = field,
    return_format = "mat"
  )
  diagram <- .as_diagram_matrix(raw)
  provenance <- list(
    result_contract_id = "legacy_topology_result_v0",
    view_id = "legacy_gene_view_v0",
    scientific_eligible = FALSE,
    axis_role = "rows_as_points_historical",
    input_sha256 = .scientific_digest(x),
    diagram_sha256 = .scientific_digest(diagram)
  )
  structure(
    list(
      diagram = diagram,
      provenance = provenance,
      cache_key = paste0(
        "legacy_topology_result_v0:", .scientific_digest(provenance)
      )
    ),
    class = c("scph_legacy_topology_result_v0", "scph_topology_result")
  )
}
