# Exact, resource-bounded H1 profiling for corrected MV17-G gene geometry.

mv17gr_materialize_gene_case_v1 <- function(x, null_family, seed) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  seed <- as.integer(seed)
  if (anyNA(x) || any(!is.finite(x)) || length(seed) != 1L || is.na(seed)) {
    stop("invalid MV17-GR source matrix/seed", call. = FALSE)
  }
  if (identical(null_family, "observed")) {
    if (seed != 0L) stop("observed MV17-GR case requires seed zero", call. = FALSE)
    raw <- x
  } else {
    registry <- mv17b_null_registry_v1()
    hit <- registry[registry$null_family == null_family, , drop = FALSE]
    if (nrow(hit) != 1L ||
        !(hit$compatible_view %in% c("cell_and_gene", "gene_only"))) {
      stop("invalid MV17-GR gene null family", call. = FALSE)
    }
    raw <- get(hit$function_name, mode = "function")(x, seed)
  }
  mv17g_topology_coordinates_v2(raw, "gene")
}

mv17gr_distance_contract_v1 <- function(coordinates) {
  coordinates <- as.matrix(coordinates)
  distances <- as.matrix(stats::dist(coordinates))
  eccentricity <- apply(distances, 1L, max)
  anchor <- which.min(eccentricity)[[1L]]
  raw_threshold <- eccentricity[[anchor]]
  numerical_guard <- 64 * .Machine$double.eps * max(1, raw_threshold)
  threshold <- raw_threshold + numerical_guard
  upper <- distances[upper.tri(distances)]
  if (!is.finite(threshold) || threshold <= 0 ||
      max(distances[anchor, ]) > threshold) {
    stop("MV17-GR cone threshold invariant failed", call. = FALSE)
  }
  list(
    distances = distances,
    cone_anchor_index = anchor,
    cone_threshold_raw = raw_threshold,
    cone_threshold = threshold,
    diameter = max(upper),
    admitted_edges = sum(upper <= threshold),
    possible_edges = length(upper),
    edge_fraction = mean(upper <= threshold)
  )
}

mv17gr_canonical_h1_v1 <- function(diagram) {
  diagram <- as.data.frame(diagram)
  if (ncol(diagram) < 3L) stop("invalid MV17-GR diagram", call. = FALSE)
  names(diagram)[1:3] <- c("dimension", "birth", "death")
  h1 <- diagram[
    diagram$dimension == 1L & is.finite(diagram$birth) &
      is.finite(diagram$death) & diagram$death > diagram$birth,
    c("dimension", "birth", "death"), drop = FALSE
  ]
  h1 <- h1[order(h1$birth, h1$death, method = "radix"), , drop = FALSE]
  rownames(h1) <- NULL
  h1
}

mv17gr_run_exact_h1_v1 <- function(coordinates, engine) {
  engine <- match.arg(
    engine, c("ripserr_infinite_v1", "ripserr_cone_exact_v1",
              "gudhi_cone_exact_v1")
  )
  contract <- mv17gr_distance_contract_v1(coordinates)
  if (engine == "ripserr_infinite_v1") {
    raw <- as.data.frame(ripserr::vietoris_rips(
      stats::as.dist(contract$distances), max_dim = 1L, threshold = Inf, p = 2L
    ))
    threshold <- Inf
    engine_version <- as.character(utils::packageVersion("ripserr"))
  } else if (engine == "ripserr_cone_exact_v1") {
    raw <- as.data.frame(ripserr::vietoris_rips(
      stats::as.dist(contract$distances), max_dim = 1L,
      threshold = contract$cone_threshold, p = 2L
    ))
    threshold <- contract$cone_threshold
    engine_version <- as.character(utils::packageVersion("ripserr"))
  } else {
    if (!requireNamespace("TDA", quietly = TRUE)) {
      stop("TDA/GUDHI is unavailable for MV17-GR", call. = FALSE)
    }
    gudhi <- TDA::ripsDiag(
      X = contract$distances, maxdimension = 1L,
      maxscale = contract$cone_threshold, dist = "arbitrary",
      library = "GUDHI", location = FALSE, printProgress = FALSE
    )[["diagram"]]
    raw <- data.frame(
      dimension = as.numeric(gudhi[, 1L]),
      birth = as.numeric(gudhi[, 2L]),
      death = as.numeric(gudhi[, 3L])
    )
    threshold <- contract$cone_threshold
    engine_version <- as.character(utils::packageVersion("TDA"))
  }
  h1 <- mv17gr_canonical_h1_v1(raw)
  if (nrow(h1) && max(h1$death) > contract$cone_threshold + 1e-10) {
    stop("MV17-GR H1 death exceeds the exact cone threshold", call. = FALSE)
  }
  list(
    contract_id = "mv17gr_exact_h1_result_v1",
    engine = engine,
    engine_version = engine_version,
    filtration = "complete_vietoris_rips_field_2",
    exact_cone_argument = TRUE,
    threshold = threshold,
    cone_threshold = contract$cone_threshold,
    cone_threshold_raw = contract$cone_threshold_raw,
    diameter = contract$diameter,
    admitted_edges = contract$admitted_edges,
    possible_edges = contract$possible_edges,
    edge_fraction = contract$edge_fraction,
    h1 = h1,
    h1_metrics = mv17c_diagram_metrics_v1(h1, 1L),
    finite = all(is.finite(unlist(mv17c_diagram_metrics_v1(h1, 1L)["value"])))
  )
}

mv17gr_maximum_h1_difference_v1 <- function(a, b) {
  a <- mv17gr_canonical_h1_v1(a)
  b <- mv17gr_canonical_h1_v1(b)
  if (nrow(a) != nrow(b)) return(Inf)
  if (!nrow(a)) return(0)
  max(abs(as.matrix(a[c("birth", "death")]) -
          as.matrix(b[c("birth", "death")])))
}
