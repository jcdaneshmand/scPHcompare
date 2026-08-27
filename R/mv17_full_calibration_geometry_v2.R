# View-aware MV17-G calibration geometry.
#
# Version 1 of the full-calibration worker passed raw matrices directly to
# Ripserr.  That is correct for the frozen cell PCA coordinates but bypasses
# the frozen gene correlation-chord construction.  Version 2 keeps null
# generation in its qualified raw input space, then applies the accepted
# view-specific topology transform before PH.

mv17g_topology_coordinates_v2 <- function(x, view) {
  view <- match.arg(view, c("cell", "gene"))
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  if (nrow(x) < 4L || ncol(x) < 2L || anyNA(x) || any(!is.finite(x))) {
    stop("invalid MV17-G topology-coordinate input", call. = FALSE)
  }
  if (view == "cell") {
    return(x)
  }

  centered <- x - rowMeans(x)
  norms <- sqrt(rowSums(centered ^ 2))
  if (any(!is.finite(norms)) ||
      any(norms <= sqrt(.Machine$double.eps))) {
    stop(
      "MV17-G gene correlation-chord transform has a constant row",
      call. = FALSE
    )
  }
  output <- centered / norms
  dimnames(output) <- dimnames(x)
  if (max(abs(rowMeans(output))) > 1e-12 ||
      max(abs(rowSums(output ^ 2) - 1)) > 1e-12) {
    stop("MV17-G gene correlation-chord invariant failed", call. = FALSE)
  }
  output
}

mv17g_geometry_id_v2 <- function(view) {
  view <- match.arg(view, c("cell", "gene"))
  if (view == "cell") {
    "euclidean_shared_pca_v1"
  } else {
    "euclidean_correlation_chord_v1"
  }
}

mv17g_group_metrics_v2 <- function(x, view, null_family, seeds) {
  view <- match.arg(view, c("cell", "gene"))
  x <- as.matrix(x)
  seeds <- as.integer(seeds)
  if (anyNA(x) || any(!is.finite(x)) || !length(seeds) || anyNA(seeds)) {
    stop("finite MV17-G matrix and seeds required", call. = FALSE)
  }

  registry <- mv17b_null_registry_v1()
  if (identical(null_family, "observed")) {
    if (!identical(seeds, 0L)) {
      stop("observed MV17-G group requires seed zero", call. = FALSE)
    }
    generator <- NULL
  } else {
    hit <- registry[registry$null_family == null_family, , drop = FALSE]
    compatible <- nrow(hit) == 1L &&
      (hit$compatible_view == "cell_and_gene" ||
       (view == "gene" && hit$compatible_view == "gene_only"))
    if (!isTRUE(compatible)) {
      stop("incompatible MV17-G null family/view", call. = FALSE)
    }
    generator <- get(hit$function_name, mode = "function")
  }

  rows <- lapply(seq_along(seeds), function(i) {
    seed <- seeds[[i]]
    raw_null <- if (is.null(generator)) x else generator(x, seed)
    topology_coordinates <- mv17g_topology_coordinates_v2(raw_null, view)
    diagram <- as.data.frame(ripserr::vietoris_rips(
      topology_coordinates, max_dim = 1L, threshold = Inf
    ))
    metrics <- do.call(rbind, lapply(0:1, function(dimension) {
      mv17c_diagram_metrics_v1(diagram, dimension)
    }))
    h0 <- diagram[
      diagram[, 1L] == 0L & is.finite(diagram[, 3L]) &
        diagram[, 3L] > diagram[, 2L], , drop = FALSE
    ]
    mst <- mv17d_h0_mergers_v1(topology_coordinates)
    h0_error <- if (nrow(h0) == nrow(mst)) {
      max(abs(sort(h0[, 3L]) - sort(mst$death)))
    } else {
      Inf
    }
    if (!is.finite(h0_error) || h0_error > 1e-8) {
      stop("MV17-G H0 MST oracle failed", call. = FALSE)
    }
    metrics$replicate <- if (identical(null_family, "observed")) 0L else i
    metrics$seed <- seed
    metrics$h0_mst_maximum_absolute_error <- h0_error
    metrics$geometry <- mv17g_geometry_id_v2(view)
    metrics[c(
      "replicate", "seed", "homology_dimension", "summary_id", "value",
      "h0_mst_maximum_absolute_error", "geometry"
    )]
  })

  output <- do.call(rbind, rows)
  rownames(output) <- NULL
  if (any(!is.finite(output$value))) {
    stop("non-finite MV17-G metrics", call. = FALSE)
  }
  output
}
