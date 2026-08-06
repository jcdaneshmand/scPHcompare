# Internal diagnostics for persistence-landscape convergence profiling.

landscape_interval_depth <- function(pd, dimension) {
  finite_pd <- finite_landscape_diagram(pd, dimension)
  if (nrow(finite_pd) == 0L) return(0L)
  events <- data.frame(
    location = c(finite_pd[, 2], finite_pd[, 3]),
    change = c(rep.int(1L, nrow(finite_pd)), rep.int(-1L, nrow(finite_pd)))
  )
  # At tied coordinates, deaths precede births because landscape tents are zero
  # at both endpoints. The active count immediately after the coordinate then
  # contains new intervals but not intervals that have just died.
  events <- events[order(events$location, events$change), , drop = FALSE]
  as.integer(max(c(0L, cumsum(events$change))))
}

landscape_interval_inventory <- function(pd_list, stratum = NA_character_) {
  if (!is.list(pd_list) || length(pd_list) < 1L) {
    stop("pd_list must be a non-empty list.", call. = FALSE)
  }
  sample_ids <- names(pd_list)
  if (is.null(sample_ids)) sample_ids <- as.character(seq_along(pd_list))
  rows <- lapply(seq_along(pd_list), function(i) {
    do.call(rbind, lapply(0:1, function(dimension) {
      finite_pd <- finite_landscape_diagram(pd_list[[i]], dimension)
      data.frame(
        stratum = stratum,
        sample_id = sample_ids[[i]],
        dimension = paste0("H", dimension),
        finite_intervals = nrow(finite_pd),
        maximum_active_levels = landscape_interval_depth(
          pd_list[[i]], dimension
        ),
        minimum_birth = if (nrow(finite_pd)) min(finite_pd[, 2]) else NA_real_,
        maximum_death = if (nrow(finite_pd)) max(finite_pd[, 3]) else NA_real_,
        stringsAsFactors = FALSE
      )
    }))
  })
  do.call(rbind, rows)
}

landscape_exact_total_energy <- function(pd, dimension) {
  finite_pd <- finite_landscape_diagram(pd, dimension)
  if (!nrow(finite_pd)) return(0)
  persistence <- finite_pd[, 3] - finite_pd[, 2]
  sum(persistence ^ 3 / 12)
}

landscape_quadrature_weights <- function(grid) {
  grid <- validate_landscape_grid(grid)
  delta <- diff(grid)
  c(delta[[1]] / 2, (delta[-1L] + delta[-length(delta)]) / 2,
    delta[[length(delta)]] / 2)
}

landscape_energy_capture <- function(values, total_squared, grid, levels) {
  grid <- validate_landscape_grid(grid)
  values <- as.matrix(values)
  if (!is.numeric(values) || nrow(values) != length(grid) ||
      !is.numeric(total_squared) || length(total_squared) != length(grid) ||
      any(!is.finite(values)) || any(!is.finite(total_squared)) ||
      any(values < 0) || any(total_squared < 0)) {
    stop("Landscape values and total squared energy must match the grid.",
         call. = FALSE)
  }
  levels <- validate_landscape_levels(levels)
  if (max(levels) > ncol(values)) {
    stop("Requested levels exceed the evaluated landscape columns.",
         call. = FALSE)
  }
  weights <- landscape_quadrature_weights(grid)
  total_energy <- sum(weights * total_squared)
  captured <- vapply(levels, function(level) {
    sum(weights * rowSums(values[, seq_len(level), drop = FALSE] ^ 2))
  }, numeric(1))
  tail <- pmax(0, total_energy - captured)
  data.frame(
    levels = levels,
    captured_energy = captured,
    total_energy = total_energy,
    tail_energy = tail,
    tail_energy_fraction = if (total_energy > 0) tail / total_energy else 0,
    stringsAsFactors = FALSE
  )
}

landscape_distance_matrix_from_values <- function(values_list, grid, levels) {
  grid <- validate_landscape_grid(grid)
  if (!is.list(values_list) || length(values_list) < 2L) {
    stop("values_list must contain at least two landscape matrices.",
         call. = FALSE)
  }
  if (length(levels) != 1L || is.na(levels) || levels < 1L ||
      levels != as.integer(levels)) {
    stop("levels must be one positive integer.", call. = FALSE)
  }
  weights <- sqrt(landscape_quadrature_weights(grid))
  features <- t(vapply(values_list, function(values) {
    values <- as.matrix(values)
    if (!is.numeric(values) || nrow(values) != length(grid) ||
        ncol(values) < levels || any(!is.finite(values))) {
      stop("Every landscape matrix must match grid and requested levels.",
           call. = FALSE)
    }
    as.vector(values[, seq_len(levels), drop = FALSE] * weights)
  }, numeric(length(grid) * levels)))
  sample_ids <- names(values_list)
  if (is.null(sample_ids)) sample_ids <- as.character(seq_along(values_list))
  result <- as.matrix(stats::dist(features, method = "euclidean"))
  dimnames(result) <- list(sample_ids, sample_ids)
  result
}

landscape_distance_matrix_chunked <- function(values_list, grid, levels,
                                              level_chunk = 250L) {
  grid <- validate_landscape_grid(grid)
  if (!is.list(values_list) || length(values_list) < 2L) {
    stop("values_list must contain at least two landscape matrices.",
         call. = FALSE)
  }
  if (length(levels) != 1L || is.na(levels) || levels < 1L ||
      levels != as.integer(levels) || length(level_chunk) != 1L ||
      is.na(level_chunk) || level_chunk < 1L ||
      level_chunk != as.integer(level_chunk)) {
    stop("levels and level_chunk must be positive integers.", call. = FALSE)
  }
  sample_ids <- names(values_list)
  if (is.null(sample_ids)) sample_ids <- as.character(seq_along(values_list))
  for (values in values_list) {
    if (!is.numeric(values) || nrow(values) != length(grid) ||
        ncol(values) < levels || any(!is.finite(values))) {
      stop("Every landscape matrix must match grid and requested levels.",
           call. = FALSE)
    }
  }
  weights <- sqrt(landscape_quadrature_weights(grid))
  squared <- matrix(0, length(values_list), length(values_list))
  starts <- seq.int(1L, levels, by = as.integer(level_chunk))
  for (start in starts) {
    end <- min(levels, start + level_chunk - 1L)
    width <- end - start + 1L
    features <- t(vapply(values_list, function(values) {
      as.vector(values[, start:end, drop = FALSE] * weights)
    }, numeric(length(grid) * width)))
    squared <- squared + as.matrix(stats::dist(features)) ^ 2
  }
  result <- sqrt(squared)
  dimnames(result) <- list(sample_ids, sample_ids)
  result
}

combine_landscape_distance_matrices <- function(distance_h0, distance_h1) {
  distance_h0 <- as.matrix(distance_h0)
  distance_h1 <- as.matrix(distance_h1)
  if (!identical(dim(distance_h0), dim(distance_h1)) ||
      !identical(dimnames(distance_h0), dimnames(distance_h1)) ||
      any(!is.finite(distance_h0)) || any(!is.finite(distance_h1))) {
    stop("H0 and H1 distance matrices must be finite and aligned.",
         call. = FALSE)
  }
  sqrt(distance_h0 ^ 2 + distance_h1 ^ 2)
}

landscape_distance_stability <- function(candidate, reference) {
  candidate <- as.matrix(candidate)
  reference <- as.matrix(reference)
  if (!identical(dim(candidate), dim(reference)) ||
      !identical(dimnames(candidate), dimnames(reference)) ||
      nrow(candidate) < 2L || any(!is.finite(candidate)) ||
      any(!is.finite(reference))) {
    stop("Candidate and reference distances must be finite and aligned.",
         call. = FALSE)
  }
  lower <- lower.tri(reference)
  candidate_values <- candidate[lower]
  reference_values <- reference[lower]
  denominator <- sqrt(sum(reference_values ^ 2))
  data.frame(
    spearman = if (stats::sd(candidate_values) > 0 &&
        stats::sd(reference_values) > 0) {
      suppressWarnings(stats::cor(
        candidate_values, reference_values, method = "spearman"
      ))
    } else if (isTRUE(all.equal(candidate_values, reference_values))) {
      1
    } else {
      NA_real_
    },
    relative_frobenius_error = if (denominator > 0) {
      sqrt(sum((candidate_values - reference_values) ^ 2)) / denominator
    } else if (all(candidate_values == 0)) 0 else Inf,
    maximum_absolute_error = max(abs(candidate_values - reference_values)),
    stringsAsFactors = FALSE
  )
}

landscape_clustering_stability <- function(candidate, reference,
                                           cluster_counts = c(2L, 5L, 8L),
                                           methods = c("average", "ward.D2")) {
  candidate <- as.matrix(candidate)
  reference <- as.matrix(reference)
  if (!identical(dim(candidate), dim(reference)) ||
      !identical(dimnames(candidate), dimnames(reference)) ||
      nrow(candidate) < 3L) {
    stop("Candidate and reference distances must be aligned and at least 3x3.",
         call. = FALSE)
  }
  cluster_counts <- unique(as.integer(cluster_counts))
  cluster_counts <- cluster_counts[cluster_counts >= 2L &
    cluster_counts < nrow(candidate)]
  if (!length(cluster_counts)) {
    stop("No cluster count is valid for the supplied matrices.", call. = FALSE)
  }
  rows <- lapply(methods, function(method) {
    candidate_tree <- stats::hclust(stats::as.dist(candidate), method = method)
    reference_tree <- stats::hclust(stats::as.dist(reference), method = method)
    cophenetic_correlation <- suppressWarnings(stats::cor(
      stats::cophenetic(candidate_tree), stats::cophenetic(reference_tree),
      method = "pearson"
    ))
    do.call(rbind, lapply(cluster_counts, function(k) {
      data.frame(
        method = method,
        clusters = k,
        adjusted_rand_index = mclust::adjustedRandIndex(
          stats::cutree(candidate_tree, k), stats::cutree(reference_tree, k)
        ),
        cophenetic_correlation = cophenetic_correlation,
        stringsAsFactors = FALSE
      )
    }))
  })
  do.call(rbind, rows)
}

landscape_h1_energy_summary <- function(distance_h0, distance_h1) {
  distance_h0 <- as.matrix(distance_h0)
  distance_h1 <- as.matrix(distance_h1)
  if (!identical(dim(distance_h0), dim(distance_h1)) ||
      !identical(dimnames(distance_h0), dimnames(distance_h1))) {
    stop("H0 and H1 distance matrices must be aligned.", call. = FALSE)
  }
  lower <- lower.tri(distance_h0)
  h0 <- distance_h0[lower]
  h1 <- distance_h1[lower]
  total <- h0 ^ 2 + h1 ^ 2
  fraction <- ifelse(total > 0, h1 ^ 2 / total, NA_real_)
  quantiles <- stats::quantile(
    fraction, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE,
    names = FALSE
  )
  data.frame(
    pairs = sum(is.finite(fraction)),
    h1_energy_min = quantiles[[1]],
    h1_energy_q25 = quantiles[[2]],
    h1_energy_median = quantiles[[3]],
    h1_energy_q75 = quantiles[[4]],
    h1_energy_max = quantiles[[5]],
    stringsAsFactors = FALSE
  )
}
