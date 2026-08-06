# Internal, versioned persistence-landscape contract helpers.

landscape_specification_registry <- function() {
  data.frame(
    specification = c(
      "legacy_k1_unit_grid_v0",
      "paper_k1_common_grid_v1",
      "full_l2_common_grid_v1",
      "full_l2_error_controlled_v1"
    ),
    levels = c(
      "1", "1", "explicit positive integer sequence",
      "all consecutive active levels"
    ),
    grid = c(
      "100 points on [0,1]",
      "dimension-specific common absolute finite filtration range",
      "dimension-specific common absolute finite filtration range",
      "dimension-specific shared adaptive grid or exact integration"
    ),
    level_summary = c(
      "not applicable",
      "first level",
      "pointwise L2 across levels",
      "pointwise L2 across levels"
    ),
    distance = c(
      "unweighted discrete Euclidean",
      "trapezoidal L2",
      "trapezoidal L2 over levels and filtration",
      "exact or error-controlled L2 over levels and filtration"
    ),
    status = c(
      "historical compatibility",
      "paper sensitivity",
      "superseded fixed-grid proposal",
      "owner-approved target; not activated"
    ),
    stringsAsFactors = FALSE
  )
}

validate_landscape_grid <- function(grid) {
  if (!is.numeric(grid) || length(grid) < 2L || any(!is.finite(grid)) ||
      is.unsorted(grid, strictly = TRUE)) {
    stop("grid must contain at least two strictly increasing finite values.",
         call. = FALSE)
  }
  as.numeric(grid)
}

validate_landscape_levels <- function(levels) {
  if (!is.numeric(levels) || length(levels) < 1L || any(!is.finite(levels)) ||
      any(levels < 1) || any(levels != as.integer(levels)) ||
      anyDuplicated(levels)) {
    stop("levels must be distinct positive integers.", call. = FALSE)
  }
  as.integer(levels)
}

canonicalize_landscape_matrix <- function(values, grid, levels) {
  grid <- validate_landscape_grid(grid)
  levels <- validate_landscape_levels(levels)
  if (!is.numeric(values)) {
    stop("landscape values must be numeric.", call. = FALSE)
  }
  if (is.null(dim(values))) {
    if (length(values) != length(grid) ||
        (any(values != 0) && length(levels) != 1L)) {
      stop("Vector landscape output is incompatible with grid/levels.",
           call. = FALSE)
    }
    values <- if (length(levels) == 1L) {
      matrix(values, ncol = 1L)
    } else {
      matrix(0, nrow = length(grid), ncol = length(levels))
    }
  } else {
    values <- as.matrix(values)
    time_by_level <- identical(dim(values), c(length(grid), length(levels)))
    level_by_time <- identical(dim(values), c(length(levels), length(grid)))
    if (time_by_level && level_by_time && length(grid) > 1L) {
      stop("Square landscape matrix orientation is ambiguous.", call. = FALSE)
    }
    if (level_by_time) {
      values <- t(values)
    } else if (!time_by_level) {
      stop(
        "Landscape matrix dimensions do not match the declared grid and levels.",
        call. = FALSE
      )
    }
  }
  storage.mode(values) <- "double"
  rownames(values) <- NULL
  colnames(values) <- paste0("lambda_", levels)
  values
}

finite_landscape_diagram <- function(pd, dimension) {
  if (!is.matrix(pd) && !is.data.frame(pd)) {
    stop("pd must be a matrix or data frame.", call. = FALSE)
  }
  pd <- as.matrix(pd)
  if (!is.numeric(pd) || ncol(pd) < 3L) {
    stop("pd must have at least three numeric columns: dimension, birth, death.",
         call. = FALSE)
  }
  keep <- pd[, 1] == dimension & is.finite(pd[, 2]) & is.finite(pd[, 3]) &
    pd[, 2] < pd[, 3]
  pd[keep, seq_len(3L), drop = FALSE]
}

compute_landscape_values <- function(pd, dimension, grid, levels = 1L) {
  grid <- validate_landscape_grid(grid)
  levels <- validate_landscape_levels(levels)
  finite_pd <- finite_landscape_diagram(pd, dimension)
  if (nrow(finite_pd) == 0L) {
    return(canonicalize_landscape_matrix(
      numeric(length(grid)), grid = grid, levels = levels
    ))
  }
  values <- TDA::landscape(
    Diag = finite_pd,
    dimension = dimension,
    KK = levels,
    tseq = grid
  )
  canonicalize_landscape_matrix(values, grid = grid, levels = levels)
}

derive_common_landscape_grid <- function(pd_list, grid_points = 500L,
                                         dimension = NULL) {
  if (!is.list(pd_list) || length(pd_list) < 1L) {
    stop("pd_list must be a non-empty list.", call. = FALSE)
  }
  if (length(grid_points) != 1L || is.na(grid_points) || grid_points < 2L ||
      grid_points != as.integer(grid_points)) {
    stop("grid_points must be one integer of at least two.", call. = FALSE)
  }
  if (!is.null(dimension) && (length(dimension) != 1L || is.na(dimension) ||
      dimension < 0 || dimension != as.integer(dimension))) {
    stop("dimension must be NULL or one non-negative integer.", call. = FALSE)
  }
  ranges <- lapply(pd_list, function(pd) {
    pd <- as.matrix(pd)
    if (!is.numeric(pd) || ncol(pd) < 3L) {
      stop("Every persistence diagram must have three numeric columns.",
           call. = FALSE)
    }
    keep <- is.finite(pd[, 2]) & is.finite(pd[, 3]) & pd[, 2] < pd[, 3]
    if (!is.null(dimension)) keep <- keep & pd[, 1] == dimension
    if (!any(keep)) return(NULL)
    range(c(pd[keep, 2], pd[keep, 3]))
  })
  ranges <- Filter(Negate(is.null), ranges)
  if (length(ranges) == 0L) {
    stop("No finite positive-persistence intervals are available.", call. = FALSE)
  }
  limits <- range(unlist(ranges, use.names = FALSE))
  seq(limits[[1]], limits[[2]], length.out = as.integer(grid_points))
}

normalize_landscape_grids <- function(grids) {
  if (!is.list(grids)) grids <- list(dim0 = grids, dim1 = grids)
  if (!identical(names(grids), c("dim0", "dim1"))) {
    stop("grids must be numeric or a dim0/dim1 list.", call. = FALSE)
  }
  lapply(grids, validate_landscape_grid)
}

normalize_landscape_level_sets <- function(levels) {
  if (!is.list(levels)) levels <- list(dim0 = levels, dim1 = levels)
  if (!identical(names(levels), c("dim0", "dim1"))) {
    stop("levels must be numeric or a dim0/dim1 list.", call. = FALSE)
  }
  lapply(levels, validate_landscape_levels)
}

new_landscape_contract <- function(pd, grids, levels = 1L,
                                   specification = "paper_k1_common_grid_v1") {
  registry <- landscape_specification_registry()
  if (length(specification) != 1L ||
      !specification %in% registry$specification) {
    stop("Unknown landscape specification: ", specification, call. = FALSE)
  }
  grids <- normalize_landscape_grids(grids)
  levels <- normalize_landscape_level_sets(levels)
  if (specification %in% c(
      "legacy_k1_unit_grid_v0", "paper_k1_common_grid_v1"
    ) && (!identical(levels$dim0, 1L) || !identical(levels$dim1, 1L))) {
    stop(specification, " requires landscape level 1 only.", call. = FALSE)
  }
  if (specification %in% c(
      "full_l2_common_grid_v1", "full_l2_error_controlled_v1"
    )) {
    consecutive <- vapply(levels, function(x) {
      identical(x, seq_len(max(x)))
    }, logical(1))
    if (!all(consecutive)) {
      stop(
        specification, " requires consecutive levels from 1 through K.",
        call. = FALSE
      )
    }
  }
  list(
    dim0 = compute_landscape_values(pd, 0L, grids$dim0, levels$dim0),
    dim1 = compute_landscape_values(pd, 1L, grids$dim1, levels$dim1),
    grids = grids,
    levels = levels,
    specification = specification
  )
}

landscape_level_summary <- function(values, grid, levels,
                                    method = c("l2", "first", "mean")) {
  method <- match.arg(method)
  values <- canonicalize_landscape_matrix(values, grid = grid, levels = levels)
  switch(
    method,
    l2 = sqrt(rowSums(values ^ 2)),
    first = values[, 1L],
    mean = rowMeans(values)
  )
}

aggregate_landscape_contract <- function(landscape,
                                         level_method = c("l2", "first", "mean")) {
  level_method <- match.arg(level_method)
  if (!identical(landscape$grids$dim0, landscape$grids$dim1)) {
    stop("Pointwise H0/H1 aggregation requires one shared filtration grid.",
         call. = FALSE)
  }
  curve0 <- landscape_level_summary(
    landscape$dim0, landscape$grids$dim0, landscape$levels$dim0, level_method
  )
  curve1 <- landscape_level_summary(
    landscape$dim1, landscape$grids$dim1, landscape$levels$dim1, level_method
  )
  sqrt(curve0 ^ 2 + curve1 ^ 2)
}

trapezoidal_l2 <- function(pointwise_squared_difference, grid) {
  grid <- validate_landscape_grid(grid)
  if (!is.numeric(pointwise_squared_difference) ||
      length(pointwise_squared_difference) != length(grid) ||
      any(!is.finite(pointwise_squared_difference)) ||
      any(pointwise_squared_difference < 0)) {
    stop("Squared differences must be finite, non-negative, and match grid.",
         call. = FALSE)
  }
  delta <- diff(grid)
  sqrt(sum(delta * (
    pointwise_squared_difference[-length(grid)] +
      pointwise_squared_difference[-1L]
  ) / 2))
}

landscape_distance_components <- function(first, second,
                                          dimension_weights = c(dim0 = 1, dim1 = 1)) {
  if (!identical(first$grids, second$grids) ||
      !identical(first$levels, second$levels)) {
    stop("Landscapes must use identical dimension-specific grids and level sets.",
         call. = FALSE)
  }
  if (!is.numeric(dimension_weights) ||
      !identical(names(dimension_weights), c("dim0", "dim1")) ||
      any(!is.finite(dimension_weights)) || any(dimension_weights < 0)) {
    stop("dimension_weights must be finite non-negative dim0/dim1 weights.",
         call. = FALSE)
  }
  per_dimension <- vapply(c("dim0", "dim1"), function(dimension) {
    a <- canonicalize_landscape_matrix(
      first[[dimension]], first$grids[[dimension]], first$levels[[dimension]]
    )
    b <- canonicalize_landscape_matrix(
      second[[dimension]], second$grids[[dimension]], second$levels[[dimension]]
    )
    trapezoidal_l2(rowSums((a - b) ^ 2), first$grids[[dimension]])
  }, numeric(1))
  c(
    per_dimension,
    combined = sqrt(sum(dimension_weights * per_dimension ^ 2))
  )
}
