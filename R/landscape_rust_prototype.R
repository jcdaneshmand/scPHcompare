# Internal optional MV5-BC Rust prototype. The accepted R engine remains
# canonical; these helpers are deliberately not exported.

.scph_rust_prototype_state <- new.env(parent = emptyenv())

.canonical_rust_intervals <- function(intervals) {
  if (is.null(intervals)) {
    return(matrix(numeric(), nrow = 0L, ncol = 2L,
                  dimnames = list(NULL, c("birth", "death"))))
  }
  if (!is.matrix(intervals) && !is.data.frame(intervals)) {
    stop("Rust prototype intervals must be a two-column matrix or data frame.",
         call. = FALSE)
  }
  intervals <- as.matrix(intervals)
  storage.mode(intervals) <- "double"
  if (ncol(intervals) != 2L || any(!is.finite(intervals)) ||
      any(intervals[, 2L] <= intervals[, 1L])) {
    stop("Rust prototype intervals require finite birth < death rows.",
         call. = FALSE)
  }
  if (nrow(intervals)) {
    intervals <- intervals[order(intervals[, 1L], -intervals[, 2L]), ,
                           drop = FALSE]
  }
  colnames(intervals) <- c("birth", "death")
  intervals
}

.load_scph_rust_prototype <- function(library) {
  if (!is.character(library) || length(library) != 1L || !nzchar(library) ||
      !file.exists(library)) {
    return(FALSE)
  }
  normalized <- normalizePath(library, winslash = "/", mustWork = TRUE)
  if (identical(.scph_rust_prototype_state$library, normalized) &&
      inherits(.scph_rust_prototype_state$symbol, "NativeSymbolInfo")) {
    return(TRUE)
  }
  loaded <- tryCatch(dyn.load(normalized), error = identity)
  if (inherits(loaded, "error")) return(FALSE)
  symbol <- tryCatch(
    getNativeSymbolInfo("scph_landscape_l2_r_v1", PACKAGE = loaded),
    error = identity
  )
  if (inherits(symbol, "error")) return(FALSE)
  .scph_rust_prototype_state$library <- normalized
  .scph_rust_prototype_state$dll <- loaded
  .scph_rust_prototype_state$symbol <- symbol
  TRUE
}

landscape_rust_prototype_dimension <- function(
    first, second, dimension,
    library = Sys.getenv("SCPH_RUST_LANDSCAPE_LIB", unset = "")) {
  first <- .canonical_rust_intervals(first)
  second <- .canonical_rust_intervals(second)
  if (!is.numeric(dimension) || length(dimension) != 1L ||
      is.na(dimension) || !dimension %in% 0:1) {
    stop("Rust prototype dimension must be 0 or 1.", call. = FALSE)
  }
  if (!.load_scph_rust_prototype(library)) {
    return(list(
      squared_distance = NA_real_, active_levels = 0, event_segments = 0,
      first_finite_intervals = nrow(first),
      second_finite_intervals = nrow(second), engine_version = 1L,
      status = 9001L, rust_used = FALSE
    ))
  }
  result <- tryCatch(
    base::.C(
      .scph_rust_prototype_state$symbol,
      first_births = as.double(first[, 1L]),
      first_deaths = as.double(first[, 2L]),
      first_len = as.integer(nrow(first)),
      second_births = as.double(second[, 1L]),
      second_deaths = as.double(second[, 2L]),
      second_len = as.integer(nrow(second)),
      dimension = as.integer(dimension),
      squared_distance = double(1L), active_levels = double(1L),
      event_segments = double(1L), first_finite_intervals = double(1L),
      second_finite_intervals = double(1L), engine_version = integer(1L),
      status = integer(1L), NAOK = FALSE
    ),
    error = identity
  )
  if (inherits(result, "error")) {
    return(list(
      squared_distance = NA_real_, active_levels = 0, event_segments = 0,
      first_finite_intervals = nrow(first),
      second_finite_intervals = nrow(second), engine_version = 1L,
      status = 9002L, rust_used = FALSE,
      error_message = conditionMessage(result)
    ))
  }
  list(
    squared_distance = result$squared_distance,
    active_levels = result$active_levels,
    event_segments = result$event_segments,
    first_finite_intervals = result$first_finite_intervals,
    second_finite_intervals = result$second_finite_intervals,
    engine_version = result$engine_version, status = result$status,
    rust_used = identical(result$status, 0L) &&
      is.finite(result$squared_distance) && result$squared_distance >= 0
  )
}

landscape_rust_prototype_with_fallback <- function(
    first, second, dimension, reference_squared,
    library = Sys.getenv("SCPH_RUST_LANDSCAPE_LIB", unset = "")) {
  if (!is.function(reference_squared)) {
    stop("reference_squared must be a zero-argument R fallback function.",
         call. = FALSE)
  }
  candidate <- landscape_rust_prototype_dimension(
    first, second, dimension, library = library
  )
  if (isTRUE(candidate$rust_used)) {
    candidate$engine <- "rust_prototype_v1"
    candidate$fallback_used <- FALSE
    return(candidate)
  }
  reference <- reference_squared()
  if (!is.numeric(reference) || length(reference) != 1L ||
      !is.finite(reference) || reference < 0) {
    stop("Canonical R fallback returned an invalid squared distance.",
         call. = FALSE)
  }
  candidate$squared_distance <- as.double(reference)
  candidate$engine <- "r_reference_fallback"
  candidate$fallback_used <- TRUE
  candidate
}
