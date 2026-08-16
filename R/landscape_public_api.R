# Public, versioned persistence-landscape distance contract.

.scph_landscape_contract_id_v1 <- "full_l2_error_controlled_v1"
.scph_landscape_pair_class_v1 <- "scph_landscape_distance_v1"
.scph_landscape_matrix_class_v1 <- "scph_landscape_distance_matrix_v1"

.canonical_public_diagram_v1 <- function(pd, label = "diagram") {
  if (!is.matrix(pd) && !is.data.frame(pd)) {
    stop(label, " must be a matrix or data frame.", call. = FALSE)
  }
  pd <- as.matrix(pd)
  if (!is.numeric(pd) || ncol(pd) < 3L) {
    stop(label, " must have at least three numeric columns: dimension, birth, death.",
         call. = FALSE)
  }
  pd <- pd[, seq_len(3L), drop = FALSE]
  storage.mode(pd) <- "double"
  colnames(pd) <- c("dimension", "birth", "death")
  rownames(pd) <- NULL
  if (!nrow(pd)) return(pd)
  dimension <- pd[, "dimension"]
  birth <- pd[, "birth"]
  death <- pd[, "death"]
  if (any(!is.finite(dimension)) ||
      any(dimension != as.integer(dimension)) ||
      any(!dimension %in% 0:1)) {
    stop(label, " dimensions must be finite integers 0 (H0) or 1 (H1).",
         call. = FALSE)
  }
  if (any(!is.finite(birth)) || anyNA(death) || any(death == -Inf)) {
    stop(label, " births must be finite and deaths must be finite or +Inf.",
         call. = FALSE)
  }
  finite <- is.finite(death)
  if (any(death[finite] <= birth[finite])) {
    stop(label, " finite intervals must have strictly positive persistence.",
         call. = FALSE)
  }
  order_index <- order(dimension, birth, death, method = "radix")
  pd <- pd[order_index, , drop = FALSE]
  rownames(pd) <- NULL
  pd
}

.public_diagram_summary_v1 <- function(pd, source_id) {
  finite <- is.finite(pd[, "death"])
  data.frame(
    source_id = source_id,
    diagram_sha256 = digest::digest(pd, algo = "sha256", serialize = TRUE),
    total_intervals = nrow(pd),
    finite_h0_intervals = sum(finite & pd[, "dimension"] == 0),
    finite_h1_intervals = sum(finite & pd[, "dimension"] == 1),
    excluded_essential_h0 = sum(!finite & pd[, "dimension"] == 0),
    excluded_essential_h1 = sum(!finite & pd[, "dimension"] == 1),
    stringsAsFactors = FALSE
  )
}

.validate_public_source_id_v1 <- function(value, fallback, label) {
  if (is.null(value)) return(fallback)
  value <- as.character(value)
  if (length(value) != 1L || is.na(value) || !nzchar(value)) {
    stop(label, " must be NULL or one non-empty string.", call. = FALSE)
  }
  value
}

.legacy_k1_unit_grid_distance_v0 <- function(first, second) {
  grid <- seq(0, 1, length.out = 100L)
  dimensions <- lapply(0:1, function(dimension) {
    first_values <- compute_landscape_values(first, dimension, grid, 1L)
    second_values <- compute_landscape_values(second, dimension, grid, 1L)
    squared <- sum((first_values - second_values) ^ 2)
    list(
      distance = sqrt(squared),
      squared_distance = squared,
      method = "legacy_k1_unit_grid_discrete_euclidean_v0",
      exact = FALSE,
      requested_absolute_tolerance = NA_real_,
      requested_relative_tolerance = NA_real_,
      achieved_absolute_error_estimate = NA_real_,
      refinement_delta = NA_real_,
      within_requested_tolerance = NA,
      integration_nodes = 100L,
      first_finite_intervals = nrow(landscape_reference_intervals(
        first, dimension
      )),
      second_finite_intervals = nrow(landscape_reference_intervals(
        second, dimension
      ))
    )
  })
  names(dimensions) <- c("H0", "H1")
  dimensions
}

.landscape_result_identity_v1 <- function(result) {
  list(
    contract_id = result$contract_id,
    specification = result$specification,
    mode = result$mode,
    first_source_id = result$provenance$first_source_id,
    second_source_id = result$provenance$second_source_id,
    first_diagram_sha256 = result$provenance$first_diagram_sha256,
    second_diagram_sha256 = result$provenance$second_diagram_sha256,
    distances = result$distances,
    dimension_methods = vapply(
      result$dimensions, `[[`, character(1L), "method"
    )
  )
}

#' Versioned persistence-landscape distance between two diagrams
#'
#' Computes H0 and H1 persistence-landscape distances separately using every
#' finite interval and every consecutive active landscape level. The default
#' scientific mode uses exact critical-breakpoint integration and never uses a
#' fixed grid or level cap. The combined value is a secondary, unweighted
#' Euclidean summary of the two dimension-specific distances.
#'
#' `mode = "legacy_k1_unit_grid_v0"` is an explicit reproduction mode for the
#' historical first-level, 100-point `[0, 1]` discrete calculation. It is not an
#' approximation to, or silent fallback for, the scientific contract.
#'
#' @param first,second Numeric persistence diagrams whose first three columns
#'   are dimension, birth, and death. Only H0 and H1 are accepted. Finite
#'   intervals must have positive persistence; `+Inf` deaths are retained in
#'   provenance and excluded before landscape construction.
#' @param method Scientific integration method: `"auto"` (default),
#'   `"exact"`, or `"adaptive"`. Auto uses exact integration through the
#'   computational interval guard and certified adaptive integration above it.
#'   Ignored in explicit legacy mode.
#' @param exact_max_intervals Positive exact-engine safety guard per diagram and
#'   dimension. The default of 500 covers the observed H0 range in the accepted
#'   MV-04 corpus. This routes numerical engines; it never removes intervals or
#'   caps landscape levels.
#' @param abs_tol,rel_tol Positive tolerances for squared-distance integration
#'   in adaptive mode.
#' @param subdivisions Positive QUADPACK subdivision limit per partition.
#' @param mode `"scientific"` or the explicit historical reproduction mode
#'   `"legacy_k1_unit_grid_v0"`.
#' @param first_id,second_id Optional non-empty source identifiers. Canonical
#'   diagram hashes are used when omitted.
#' @return A `scph_landscape_distance_v1` object containing H0, H1, and
#'   secondary combined distances; per-dimension diagnostics; immutable
#'   provenance; a deterministic cache key; and runtime metadata.
#' @examples
#' empty <- matrix(numeric(), nrow = 0, ncol = 3)
#' one_loop <- matrix(c(1, 0, 2), ncol = 3)
#' persistence_landscape_distance(one_loop, empty)$distances
#' @export
persistence_landscape_distance <- function(
    first, second, method = c("auto", "exact", "adaptive"),
    exact_max_intervals = 500L, abs_tol = 1e-8, rel_tol = 1e-8,
    subdivisions = 200L,
    mode = c("scientific", "legacy_k1_unit_grid_v0"),
    first_id = NULL, second_id = NULL) {
  started <- proc.time()[["elapsed"]]
  method <- match.arg(method)
  mode <- match.arg(mode)
  first <- .canonical_public_diagram_v1(first, "first")
  second <- .canonical_public_diagram_v1(second, "second")
  first_hash <- digest::digest(first, algo = "sha256", serialize = TRUE)
  second_hash <- digest::digest(second, algo = "sha256", serialize = TRUE)
  first_id <- .validate_public_source_id_v1(first_id, first_hash, "first_id")
  second_id <- .validate_public_source_id_v1(second_id, second_hash, "second_id")

  if (mode == "scientific") {
    reference <- landscape_reference_distance(
      first, second, method = method,
      exact_max_intervals = exact_max_intervals,
      abs_tol = abs_tol, rel_tol = rel_tol, subdivisions = subdivisions
    )
    dimensions <- reference$dimensions
    tolerance_ok <- vapply(
      dimensions, `[[`, logical(1L), "within_requested_tolerance"
    )
    if (!all(tolerance_ok)) {
      stop("Adaptive landscape integration failed its requested error bound.",
           call. = FALSE)
    }
    specification <- .scph_landscape_contract_id_v1
  } else {
    dimensions <- .legacy_k1_unit_grid_distance_v0(first, second)
    specification <- "legacy_k1_unit_grid_v0"
  }
  squared <- vapply(dimensions, `[[`, numeric(1L), "squared_distance")
  combined_squared <- sum(squared)
  distances <- c(
    H0 = sqrt(squared[["H0"]]),
    H1 = sqrt(squared[["H1"]]),
    combined = sqrt(combined_squared)
  )
  summaries <- rbind(
    .public_diagram_summary_v1(first, first_id),
    .public_diagram_summary_v1(second, second_id)
  )
  result <- structure(list(
    contract_id = "scph_public_landscape_pair_v1",
    result_version = 1L,
    specification = specification,
    mode = mode,
    distances = distances,
    h1_squared_distance_fraction = if (combined_squared > 0) {
      squared[["H1"]] / combined_squared
    } else 0,
    dimensions = dimensions,
    diagram_provenance = summaries,
    provenance = list(
      api = "persistence_landscape_distance",
      api_version = 1L,
      scientific_contract = .scph_landscape_contract_id_v1,
      specification = specification,
      method_requested = if (mode == "scientific") method else NA_character_,
      exact_max_intervals = if (mode == "scientific") {
        as.integer(exact_max_intervals)
      } else NA_integer_,
      absolute_tolerance = if (mode == "scientific") abs_tol else NA_real_,
      relative_tolerance = if (mode == "scientific") rel_tol else NA_real_,
      subdivisions = if (mode == "scientific") {
        as.integer(subdivisions)
      } else NA_integer_,
      level_policy = if (mode == "scientific") {
        "all consecutive active levels; zero-pad missing depth"
      } else "level 1 only",
      grid_policy = if (mode == "scientific") {
        "none: exact critical breakpoints or adaptive feature partitions"
      } else "100 equally spaced points on [0,1]",
      dimension_policy = "H0 and H1 separate; Euclidean combination secondary",
      infinite_interval_policy = "exclude +Inf deaths; report by dimension",
      first_source_id = first_id,
      second_source_id = second_id,
      first_diagram_sha256 = first_hash,
      second_diagram_sha256 = second_hash,
      legacy_reproduction = mode != "scientific",
      existing_workflow_default_changed = FALSE,
      engine_version = if (mode == "scientific") {
        "landscape_reference_v4"
      } else "legacy_k1_unit_grid_v0"
    ),
    runtime = list(elapsed_seconds = unname(
      proc.time()[["elapsed"]] - started
    ))
  ), class = c(.scph_landscape_pair_class_v1, "list"))
  result$cache_key <- paste0(
    "scph_landscape_distance_v1:",
    digest::digest(.landscape_result_identity_v1(result),
                   algo = "sha256", serialize = TRUE)
  )
  result
}

.validate_public_diagram_list_v1 <- function(diagrams) {
  if (!is.list(diagrams) || length(diagrams) < 1L) {
    stop("diagrams must be a non-empty named list.", call. = FALSE)
  }
  ids <- names(diagrams)
  if (is.null(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop("diagrams must be a named list with unique non-empty names.",
         call. = FALSE)
  }
  order_index <- order(ids, method = "radix")
  diagrams <- diagrams[order_index]
  lapply(seq_along(diagrams), function(index) {
    .canonical_public_diagram_v1(diagrams[[index]], names(diagrams)[[index]])
  }) |>
    stats::setNames(names(diagrams))
}

#' Complete versioned persistence-landscape distance matrices
#'
#' Applies [persistence_landscape_distance()] to every unordered pair in a
#' named diagram list. Sample order is canonical lexical order. H0 and H1
#' matrices remain primary and separate; `combined` is a secondary descriptive
#' matrix. Historical workflows and artifacts are not read, redirected, or
#' overwritten.
#'
#' @inheritParams persistence_landscape_distance
#' @param diagrams Non-empty named list of persistence diagrams with unique,
#'   non-empty names.
#' @return A `scph_landscape_distance_matrix_v1` object with separate named H0,
#'   H1, and combined matrices, pair diagnostics, diagram provenance, complete
#'   method provenance, deterministic cache key, and runtime metadata.
#' @examples
#' diagrams <- list(
#'   sample_a = matrix(c(0, 0, 1), ncol = 3),
#'   sample_b = matrix(c(0, 0, 2, 1, 0, 1), ncol = 3, byrow = TRUE)
#' )
#' persistence_landscape_distance_matrix(diagrams)$matrices$H0
#' @export
persistence_landscape_distance_matrix <- function(
    diagrams, method = c("auto", "exact", "adaptive"),
    exact_max_intervals = 500L, abs_tol = 1e-8, rel_tol = 1e-8,
    subdivisions = 200L,
    mode = c("scientific", "legacy_k1_unit_grid_v0")) {
  started <- proc.time()[["elapsed"]]
  method <- match.arg(method)
  mode <- match.arg(mode)
  diagrams <- .validate_public_diagram_list_v1(diagrams)
  ids <- names(diagrams)
  matrices <- lapply(c("H0", "H1", "combined"), function(unused) {
    matrix(0, nrow = length(ids), ncol = length(ids),
           dimnames = list(ids, ids))
  })
  names(matrices) <- c("H0", "H1", "combined")
  pair_rows <- list()
  pair_cache_keys <- character()
  if (length(ids) > 1L) {
    pairs <- utils::combn(seq_along(ids), 2L)
    for (column in seq_len(ncol(pairs))) {
      left <- pairs[1L, column]
      right <- pairs[2L, column]
      pair <- persistence_landscape_distance(
        diagrams[[left]], diagrams[[right]], method = method,
        exact_max_intervals = exact_max_intervals,
        abs_tol = abs_tol, rel_tol = rel_tol, subdivisions = subdivisions,
        mode = mode, first_id = ids[[left]], second_id = ids[[right]]
      )
      for (dimension in names(matrices)) {
        matrices[[dimension]][left, right] <- pair$distances[[dimension]]
        matrices[[dimension]][right, left] <- pair$distances[[dimension]]
      }
      pair_rows[[column]] <- data.frame(
        first_source_id = ids[[left]],
        second_source_id = ids[[right]],
        h0_distance = unname(pair$distances[["H0"]]),
        h1_distance = unname(pair$distances[["H1"]]),
        combined_distance = unname(pair$distances[["combined"]]),
        h1_squared_distance_fraction = pair$h1_squared_distance_fraction,
        h0_method = pair$dimensions$H0$method,
        h1_method = pair$dimensions$H1$method,
        h0_error_estimate = pair$dimensions$H0$achieved_absolute_error_estimate,
        h1_error_estimate = pair$dimensions$H1$achieved_absolute_error_estimate,
        cache_key = pair$cache_key,
        stringsAsFactors = FALSE
      )
      pair_cache_keys[[column]] <- pair$cache_key
    }
  }
  pair_diagnostics <- if (length(pair_rows)) do.call(rbind, pair_rows) else {
    data.frame(
      first_source_id = character(), second_source_id = character(),
      h0_distance = numeric(), h1_distance = numeric(),
      combined_distance = numeric(),
      h1_squared_distance_fraction = numeric(),
      h0_method = character(), h1_method = character(),
      h0_error_estimate = numeric(), h1_error_estimate = numeric(),
      cache_key = character(), stringsAsFactors = FALSE
    )
  }
  diagram_provenance <- do.call(rbind, lapply(ids, function(id) {
    .public_diagram_summary_v1(diagrams[[id]], id)
  }))
  rownames(diagram_provenance) <- NULL
  identity <- list(
    contract_id = "scph_public_landscape_matrix_v1",
    specification = if (mode == "scientific") {
      .scph_landscape_contract_id_v1
    } else "legacy_k1_unit_grid_v0",
    mode = mode,
    sample_ids = ids,
    diagram_sha256 = diagram_provenance$diagram_sha256,
    pair_cache_keys = pair_cache_keys,
    matrices = matrices
  )
  structure(list(
    contract_id = identity$contract_id,
    result_version = 1L,
    specification = identity$specification,
    mode = mode,
    sample_ids = ids,
    matrices = matrices,
    pair_diagnostics = pair_diagnostics,
    diagram_provenance = diagram_provenance,
    provenance = list(
      api = "persistence_landscape_distance_matrix",
      api_version = 1L,
      scientific_contract = .scph_landscape_contract_id_v1,
      method_requested = if (mode == "scientific") method else NA_character_,
      exact_max_intervals = if (mode == "scientific") {
        as.integer(exact_max_intervals)
      } else NA_integer_,
      absolute_tolerance = if (mode == "scientific") abs_tol else NA_real_,
      relative_tolerance = if (mode == "scientific") rel_tol else NA_real_,
      subdivisions = if (mode == "scientific") {
        as.integer(subdivisions)
      } else NA_integer_,
      sample_order_policy = "lexical radix order",
      pair_count = nrow(pair_diagnostics),
      dimension_policy = "H0 and H1 separate; Euclidean combination secondary",
      legacy_reproduction = mode != "scientific",
      existing_workflow_default_changed = FALSE
    ),
    runtime = list(elapsed_seconds = unname(
      proc.time()[["elapsed"]] - started
    )),
    cache_key = paste0(
      "scph_landscape_distance_matrix_v1:",
      digest::digest(identity, algo = "sha256", serialize = TRUE)
    )
  ), class = c(.scph_landscape_matrix_class_v1, "list"))
}

.detect_landscape_artifact_schema_v1 <- function(value) {
  classification <- "unrecognized"
  confidence <- "none"
  safe_action <- "reject"
  if (inherits(value, .scph_landscape_pair_class_v1)) {
    classification <- "scph_landscape_distance_v1"
    confidence <- "exact"
    safe_action <- "read_versioned"
  } else if (inherits(value, .scph_landscape_matrix_class_v1)) {
    classification <- "scph_landscape_distance_matrix_v1"
    confidence <- "exact"
    safe_action <- "read_versioned"
  } else if (is.list(value) && length(value) > 0L && all(vapply(
    value, function(item) {
      is.list(item) && identical(sort(names(item)), c("dim0", "dim1")) &&
        all(vapply(item, function(component) {
          is.numeric(component) && length(component) == 100L
        }, logical(1L)))
    }, logical(1L)))) {
    classification <- "legacy_k1_unit_grid_landscape_list_candidate"
    confidence <- "shape_only"
    safe_action <- "read_only_require_explicit_legacy_confirmation"
  } else if (is.matrix(value) && is.numeric(value) &&
             nrow(value) == ncol(value) && nrow(value) > 0L &&
             all(is.finite(value)) && isTRUE(all.equal(value, t(value)))) {
    classification <- "legacy_unversioned_combined_matrix_candidate"
    confidence <- "shape_only"
    safe_action <- "read_only_never_convert_silently"
  }
  data.frame(
    detector_version = "landscape_artifact_schema_detector_v1",
    classification = classification,
    confidence = confidence,
    safe_action = safe_action,
    silent_conversion_allowed = FALSE,
    stringsAsFactors = FALSE
  )
}
