# Internal exact and error-controlled persistence-landscape reference engine.

landscape_reference_intervals <- function(pd, dimension) {
  finite_pd <- finite_landscape_diagram(pd, dimension)
  if (!nrow(finite_pd)) {
    return(matrix(numeric(), nrow = 0L, ncol = 2L,
                  dimnames = list(NULL, c("birth", "death"))))
  }
  intervals <- finite_pd[, 2:3, drop = FALSE]
  colnames(intervals) <- c("birth", "death")
  intervals
}

landscape_reference_values <- function(intervals, location) {
  if (!nrow(intervals)) return(numeric())
  values <- pmin(location - intervals[, "birth"],
                 intervals[, "death"] - location)
  sort(values[values > 0], decreasing = TRUE)
}

landscape_reference_difference_squared <- function(first, second, location) {
  a <- landscape_reference_values(first, location)
  b <- landscape_reference_values(second, location)
  depth <- max(length(a), length(b))
  if (!depth) return(0)
  length(a) <- depth
  length(b) <- depth
  a[is.na(a)] <- 0
  b[is.na(b)] <- 0
  sum((a - b) ^ 2)
}

landscape_reference_breakpoints <- function(intervals) {
  if (!nrow(intervals)) return(numeric())
  births <- intervals[, "birth"]
  deaths <- intervals[, "death"]
  midpoints <- (births + deaths) / 2
  segments <- rbind(
    data.frame(left = births, right = midpoints, slope = 1, intercept = -births),
    data.frame(left = midpoints, right = deaths, slope = -1, intercept = deaths)
  )
  crossings <- numeric()
  if (nrow(segments) > 1L) {
    pairs <- utils::combn(seq_len(nrow(segments)), 2L)
    for (column in seq_len(ncol(pairs))) {
      i <- pairs[1L, column]
      j <- pairs[2L, column]
      if (segments$slope[[i]] == segments$slope[[j]]) next
      left <- max(segments$left[[i]], segments$left[[j]])
      right <- min(segments$right[[i]], segments$right[[j]])
      if (left >= right) next
      location <- (segments$intercept[[j]] - segments$intercept[[i]]) /
        (segments$slope[[i]] - segments$slope[[j]])
      if (location > left && location < right) {
        crossings <- c(crossings, location)
      }
    }
  }
  sort(unique(c(births, midpoints, deaths, crossings)))
}

landscape_reference_exact_dimension <- function(first, second, dimension,
                                                exact_max_intervals = 500L) {
  a <- landscape_reference_intervals(first, dimension)
  b <- landscape_reference_intervals(second, dimension)
  largest <- max(nrow(a), nrow(b))
  if (largest > exact_max_intervals) {
    stop(
      "Exact landscape reference guard exceeded: ", largest,
      " intervals versus exact_max_intervals=", exact_max_intervals, ".",
      call. = FALSE
    )
  }
  nodes <- sort(unique(c(
    landscape_reference_breakpoints(a),
    landscape_reference_breakpoints(b)
  )))
  if (length(nodes) < 2L) {
    squared <- 0
  } else {
    squared <- 0
    previous <- NULL
    for (index in seq_along(nodes)) {
      location <- nodes[[index]]
      current_a <- landscape_reference_values(a, location)
      current_b <- landscape_reference_values(b, location)
      depth <- max(length(current_a), length(current_b),
                   if (is.null(previous)) 0L else length(previous))
      current <- numeric(depth)
      if (length(current_a)) current[seq_along(current_a)] <- current_a
      if (length(current_b)) current[seq_along(current_b)] <-
        current[seq_along(current_b)] - current_b
      if (!is.null(previous)) {
        length(previous) <- depth
        previous[is.na(previous)] <- 0
        width <- location - nodes[[index - 1L]]
        squared <- squared + width * sum(
          previous ^ 2 + previous * current + current ^ 2
        ) / 3
      }
      previous <- current
    }
  }
  list(
    distance = sqrt(max(0, squared)),
    squared_distance = max(0, squared),
    method = "exact_breakpoint_stream_v1",
    exact = TRUE,
    requested_absolute_tolerance = 0,
    requested_relative_tolerance = 0,
    achieved_absolute_error_estimate = 0,
    refinement_delta = 0,
    within_requested_tolerance = TRUE,
    integration_nodes = length(nodes),
    first_finite_intervals = nrow(a),
    second_finite_intervals = nrow(b)
  )
}

landscape_reference_integrate_partition <- function(
    integrand, lower, upper, subdivisions, rel_tol, abs_tol,
    integrate_fn = stats::integrate, split_depth = 0L,
    max_split_depth = 20L) {
  result <- tryCatch(
    integrate_fn(
      integrand, lower = lower, upper = upper,
      subdivisions = subdivisions, rel.tol = rel_tol,
      abs.tol = max(.Machine$double.eps, abs_tol),
      stop.on.error = TRUE
    ),
    error = identity
  )
  if (!inherits(result, "error")) {
    result$fallback_splits <- 0L
    return(result)
  }
  recoverable <- any(vapply(
    c("extremely bad integrand behaviour", "roundoff error is detected"),
    grepl, logical(1), x = conditionMessage(result), fixed = TRUE
  ))
  if (!recoverable || split_depth >= max_split_depth) {
    stop(result)
  }
  midpoint <- lower + (upper - lower) / 2
  if (!is.finite(midpoint) || midpoint <= lower || midpoint >= upper) {
    stop(result)
  }
  left <- landscape_reference_integrate_partition(
    integrand, lower, midpoint, subdivisions, rel_tol, abs_tol / 2,
    integrate_fn, split_depth + 1L, max_split_depth
  )
  right <- landscape_reference_integrate_partition(
    integrand, midpoint, upper, subdivisions, rel_tol, abs_tol / 2,
    integrate_fn, split_depth + 1L, max_split_depth
  )
  list(
    value = left$value + right$value,
    abs.error = left$abs.error + right$abs.error,
    subdivisions = left$subdivisions + right$subdivisions,
    fallback_splits = 1L + left$fallback_splits + right$fallback_splits
  )
}

landscape_reference_adaptive_pass <- function(first, second, breaks,
                                              abs_tol, rel_tol,
                                              subdivisions) {
  if (length(breaks) < 2L) {
    return(list(
      value = 0, absolute.error = 0, evaluations = 0L,
      pilot_value = 0, global_error_budget = abs_tol,
      local_absolute_tolerance = abs_tol, fallback_splits = 0L
    ))
  }
  widths <- diff(breaks)
  partition_count <- length(widths)
  integrand <- function(location) {
    vapply(location, function(point) {
      landscape_reference_difference_squared(first, second, point)
    }, numeric(1))
  }
  # The midpoint pilot allocates one global error budget; it is not used as
  # the answer or as certification. Equal allocation avoids demanding
  # near-machine precision from short partitions containing unresolved
  # landscape-order crossings. The summed QUADPACK errors and an independent
  # quarter-tolerance pass still certify the returned global integral.
  midpoints <- (breaks[-length(breaks)] + breaks[-1L]) / 2
  pilot_value <- sum(widths * integrand(midpoints))
  global_error_budget <- max(abs_tol, rel_tol * abs(pilot_value))
  local_abs_tol <- global_error_budget / partition_count
  value <- 0
  absolute.error <- 0
  evaluations <- 0L
  fallback_splits <- 0L
  for (index in seq_along(widths)) {
    result <- landscape_reference_integrate_partition(
      integrand, breaks[[index]], breaks[[index + 1L]],
      subdivisions, rel_tol, local_abs_tol
    )
    value <- value + result$value
    absolute.error <- absolute.error + result$abs.error
    evaluations <- evaluations + result$subdivisions
    fallback_splits <- fallback_splits + result$fallback_splits
  }
  list(
    value = value, absolute.error = absolute.error,
    evaluations = evaluations, pilot_value = pilot_value,
    global_error_budget = global_error_budget,
    local_absolute_tolerance = local_abs_tol,
    fallback_splits = fallback_splits
  )
}

landscape_reference_adaptive_dimension <- function(
    first, second, dimension, abs_tol = 1e-8, rel_tol = 1e-8,
    subdivisions = 200L) {
  a <- landscape_reference_intervals(first, dimension)
  b <- landscape_reference_intervals(second, dimension)
  intervals <- rbind(a, b)
  breaks <- if (nrow(intervals)) {
    sort(unique(c(
      intervals[, "birth"],
      rowMeans(intervals),
      intervals[, "death"]
    )))
  } else numeric()
  coarse <- landscape_reference_adaptive_pass(
    a, b, breaks, abs_tol, rel_tol, subdivisions
  )
  fine <- landscape_reference_adaptive_pass(
    a, b, breaks, abs_tol / 4, rel_tol / 4, subdivisions
  )
  refinement_delta <- abs(fine$value - coarse$value)
  # Conservatively combine the fine-pass quadrature estimate with the
  # independent coarse/fine refinement delta; taking only their maximum would
  # understate the joint numerical uncertainty.
  achieved <- fine$absolute.error + refinement_delta
  threshold <- max(abs_tol, rel_tol * abs(fine$value))
  list(
    distance = sqrt(max(0, fine$value)),
    squared_distance = max(0, fine$value),
    method = "adaptive_quadpack_partitioned_v3",
    exact = FALSE,
    requested_absolute_tolerance = abs_tol,
    requested_relative_tolerance = rel_tol,
    achieved_absolute_error_estimate = achieved,
    error_estimate_policy = "fine_quadrature_error_plus_refinement_delta_v2",
    refinement_delta = refinement_delta,
    within_requested_tolerance = achieved <= threshold,
    integration_nodes = length(breaks),
    integration_subdivisions = fine$evaluations,
    partition_failure_policy =
      "recursive_bisection_bad_integrand_or_roundoff_v3",
    coarse_fallback_splits = coarse$fallback_splits,
    fine_fallback_splits = fine$fallback_splits,
    tolerance_allocation = "global_midpoint_pilot_equal_partition_v2",
    coarse_pilot_squared_distance = coarse$pilot_value,
    fine_pilot_squared_distance = fine$pilot_value,
    coarse_global_error_budget = coarse$global_error_budget,
    fine_global_error_budget = fine$global_error_budget,
    coarse_summed_quadrature_error = coarse$absolute.error,
    fine_summed_quadrature_error = fine$absolute.error,
    first_finite_intervals = nrow(a),
    second_finite_intervals = nrow(b)
  )
}

landscape_reference_provenance <- function(first, second, method,
                                           exact_max_intervals,
                                           abs_tol, rel_tol) {
  list(
    specification = "full_l2_error_controlled_v1",
    engine_version = "landscape_reference_v3",
    method_requested = method,
    exact_max_intervals = as.integer(exact_max_intervals),
    absolute_tolerance = abs_tol,
    relative_tolerance = rel_tol,
    level_policy = "all consecutive active levels",
    dimension_policy = "H0 and H1 separate; Euclidean combination secondary",
    infinite_interval_policy = "exclude before landscape construction",
    first_diagram_sha256 = digest::digest(first, algo = "sha256"),
    second_diagram_sha256 = digest::digest(second, algo = "sha256"),
    activated_as_scientific_default = FALSE,
    r_version = as.character(getRversion())
  )
}

#' Compute a non-default all-level persistence-landscape reference distance
#'
#' This reference interface implements the owner-approved, dissertation-aligned
#' target contract without changing the package's public scientific defaults.
#' H0 and H1 are integrated separately and are only then combined by an
#' unweighted Euclidean norm. Exact breakpoint integration is used for tractable
#' diagrams; larger diagrams use partitioned adaptive quadrature with an
#' independent refinement check and recorded numerical error estimate.
#'
#' @param first,second Persistence diagrams with dimension, birth, and death
#'   in their first three columns.
#' @param method One of `"auto"`, `"exact"`, or `"adaptive"`.
#' @param exact_max_intervals Maximum finite intervals in either diagram and
#'   dimension for exact breakpoint integration.
#' @param abs_tol,rel_tol Requested tolerances for adaptive squared-distance
#'   integration.
#' @param subdivisions Maximum QUADPACK subdivisions per feature partition.
#' @return A versioned reference result containing separate H0/H1 distances,
#'   secondary combined distance, H1 squared-distance contribution, per-
#'   dimension diagnostics, and provenance.
#' @keywords internal
landscape_reference_distance <- function(
    first, second, method = c("auto", "exact", "adaptive"),
    exact_max_intervals = 500L, abs_tol = 1e-8, rel_tol = 1e-8,
    subdivisions = 200L) {
  method <- match.arg(method)
  if (length(exact_max_intervals) != 1L || is.na(exact_max_intervals) ||
      exact_max_intervals < 1L || exact_max_intervals != as.integer(exact_max_intervals)) {
    stop("exact_max_intervals must be one positive integer.", call. = FALSE)
  }
  if (any(!is.finite(c(abs_tol, rel_tol))) || abs_tol <= 0 || rel_tol <= 0) {
    stop("abs_tol and rel_tol must be finite and positive.", call. = FALSE)
  }
  if (length(subdivisions) != 1L || is.na(subdivisions) ||
      subdivisions < 1L || subdivisions != as.integer(subdivisions)) {
    stop("subdivisions must be one positive integer.", call. = FALSE)
  }
  dimensions <- lapply(0:1, function(dimension) {
    counts <- c(
      nrow(landscape_reference_intervals(first, dimension)),
      nrow(landscape_reference_intervals(second, dimension))
    )
    selected <- if (method == "auto") {
      if (max(counts) <= exact_max_intervals) "exact" else "adaptive"
    } else method
    if (selected == "exact") {
      landscape_reference_exact_dimension(
        first, second, dimension, as.integer(exact_max_intervals)
      )
    } else {
      landscape_reference_adaptive_dimension(
        first, second, dimension, abs_tol, rel_tol,
        as.integer(subdivisions)
      )
    }
  })
  names(dimensions) <- c("H0", "H1")
  squared <- vapply(dimensions, `[[`, numeric(1), "squared_distance")
  combined_squared <- sum(squared)
  distances <- c(
    H0 = sqrt(squared[["H0"]]),
    H1 = sqrt(squared[["H1"]]),
    combined = sqrt(combined_squared)
  )
  structure(list(
    distances = distances,
    h1_squared_distance_fraction = if (combined_squared > 0) {
      squared[["H1"]] / combined_squared
    } else 0,
    dimensions = dimensions,
    provenance = landscape_reference_provenance(
      first, second, method, exact_max_intervals, abs_tol, rel_tol
    )
  ), class = "scph_landscape_reference")
}
