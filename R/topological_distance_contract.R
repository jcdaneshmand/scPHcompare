.validate_distance_matrix_v1 <- function(value, label = "distance matrix") {
  if (!is.matrix(value) || !is.numeric(value) || nrow(value) != ncol(value) ||
      is.null(rownames(value)) || is.null(colnames(value)) ||
      !identical(rownames(value), colnames(value)) ||
      anyDuplicated(rownames(value)) || any(!is.finite(value)) ||
      any(value < 0)) {
    stop(label, " must be a finite non-negative named square matrix.",
         call. = FALSE)
  }
  tolerance <- 100 * .Machine$double.eps * max(1, max(value))
  if (max(abs(value - t(value))) > tolerance ||
      any(abs(diag(value)) > tolerance)) {
    stop(label, " must be symmetric with an exact numerical zero diagonal.",
         call. = FALSE)
  }
  value <- (value + t(value)) / 2
  diag(value) <- 0
  value
}

distance_pairs_to_matrix_v1 <- function(pairs, sample_ids,
                                        value_column = "distance") {
  required <- c("first_sample_id", "second_sample_id", value_column)
  if (!is.data.frame(pairs) || !all(required %in% names(pairs))) {
    stop("pairs lacks the required sample/value columns.", call. = FALSE)
  }
  sample_ids <- sort(as.character(sample_ids), method = "radix")
  if (length(sample_ids) < 2L || anyNA(sample_ids) ||
      any(!nzchar(sample_ids)) || anyDuplicated(sample_ids)) {
    stop("sample_ids must contain at least two unique non-empty IDs.",
         call. = FALSE)
  }
  expected_pairs <- length(sample_ids) * (length(sample_ids) - 1L) / 2L
  if (nrow(pairs) != expected_pairs) {
    stop("pairs does not contain exactly one row per unordered sample pair.",
         call. = FALSE)
  }
  left <- as.character(pairs$first_sample_id)
  right <- as.character(pairs$second_sample_id)
  values <- as.numeric(pairs[[value_column]])
  if (anyNA(values) || any(!is.finite(values)) || any(values < 0) ||
      any(!(left %in% sample_ids)) || any(!(right %in% sample_ids)) ||
      any(left == right)) {
    stop("pairs contains invalid IDs or distance values.", call. = FALSE)
  }
  keys <- vapply(seq_along(left), function(index) {
    paste(sort(c(left[[index]], right[[index]]), method = "radix"),
          collapse = "::")
  }, character(1L))
  if (anyDuplicated(keys)) {
    stop("pairs contains a duplicated unordered sample pair.", call. = FALSE)
  }
  expected_keys <- utils::combn(sample_ids, 2L, function(value) {
    paste(value, collapse = "::")
  })
  if (!setequal(keys, expected_keys)) {
    stop("pairs is missing one or more required sample combinations.",
         call. = FALSE)
  }
  result <- matrix(
    0, nrow = length(sample_ids), ncol = length(sample_ids),
    dimnames = list(sample_ids, sample_ids)
  )
  for (index in seq_along(values)) {
    result[left[[index]], right[[index]]] <- values[[index]]
    result[right[[index]], left[[index]]] <- values[[index]]
  }
  .validate_distance_matrix_v1(result)
}

fit_distance_scale_v1 <- function(distance_matrix, fit_sample_ids = NULL,
                                  fit_scope_id) {
  distance_matrix <- .validate_distance_matrix_v1(distance_matrix)
  fit_scope_id <- as.character(fit_scope_id)
  if (length(fit_scope_id) != 1L || is.na(fit_scope_id) ||
      !nzchar(fit_scope_id)) {
    stop("fit_scope_id must be one non-empty string.", call. = FALSE)
  }
  if (is.null(fit_sample_ids)) fit_sample_ids <- rownames(distance_matrix)
  fit_sample_ids <- sort(as.character(fit_sample_ids), method = "radix")
  if (length(fit_sample_ids) < 2L || anyDuplicated(fit_sample_ids) ||
      any(!(fit_sample_ids %in% rownames(distance_matrix)))) {
    stop("fit_sample_ids must identify at least two matrix samples.",
         call. = FALSE)
  }
  fit <- distance_matrix[fit_sample_ids, fit_sample_ids, drop = FALSE]
  values <- fit[lower.tri(fit)]
  scale <- stats::median(values)
  if (!is.finite(scale) || scale <= sqrt(.Machine$double.eps)) {
    stop("Distance component is degenerate in the declared fit scope.",
         call. = FALSE)
  }
  identity <- list(
    contract_id = "median_offdiagonal_distance_scale_v1",
    fit_scope_id = fit_scope_id,
    fit_sample_ids = fit_sample_ids,
    scale = unname(scale),
    matrix_sha256 = digest::digest(
      fit, algo = "sha256", serialize = TRUE
    )
  )
  structure(
    c(identity, list(
      cache_key = paste0(
        "distance_scale_v1:",
        digest::digest(identity, algo = "sha256", serialize = TRUE)
      )
    )),
    class = "scph_distance_scale_v1"
  )
}

apply_distance_scale_v1 <- function(distance_matrix, scale_fit) {
  distance_matrix <- .validate_distance_matrix_v1(distance_matrix)
  if (!inherits(scale_fit, "scph_distance_scale_v1") ||
      !identical(scale_fit$contract_id,
                 "median_offdiagonal_distance_scale_v1") ||
      !is.numeric(scale_fit$scale) || length(scale_fit$scale) != 1L ||
      !is.finite(scale_fit$scale) ||
      scale_fit$scale <= sqrt(.Machine$double.eps)) {
    stop("scale_fit is not a valid scph_distance_scale_v1 object.",
         call. = FALSE)
  }
  result <- distance_matrix / scale_fit$scale
  diag(result) <- 0
  .validate_distance_matrix_v1(result)
}

new_mv04_distance_bundle <- function(h0, h1, stratum_id, method_id,
                                     input_diagram_ids,
                                     fit_scope_id =
                                       "descriptive_all_pilot_samples") {
  h0 <- .validate_distance_matrix_v1(h0, "h0")
  h1 <- .validate_distance_matrix_v1(h1, "h1")
  if (!identical(dimnames(h0), dimnames(h1))) {
    stop("h0 and h1 matrices must have identical ordered sample axes.",
         call. = FALSE)
  }
  stratum_id <- as.character(stratum_id)
  method_id <- as.character(method_id)
  if (length(stratum_id) != 1L || !nzchar(stratum_id) ||
      length(method_id) != 1L || !nzchar(method_id)) {
    stop("stratum_id and method_id must be non-empty strings.",
         call. = FALSE)
  }
  input_diagram_ids <- as.character(input_diagram_ids)
  if (length(input_diagram_ids) != nrow(h0) || anyNA(input_diagram_ids) ||
      any(!nzchar(input_diagram_ids)) || anyDuplicated(input_diagram_ids)) {
    stop("input_diagram_ids must bind one unique diagram to every sample.",
         call. = FALSE)
  }
  combined <- sqrt(h0 ^ 2 + h1 ^ 2)
  diag(combined) <- 0
  combined <- .validate_distance_matrix_v1(combined, "combined")
  h0_scale <- fit_distance_scale_v1(h0, fit_scope_id = fit_scope_id)
  h1_scale <- fit_distance_scale_v1(h1, fit_scope_id = fit_scope_id)
  identity <- list(
    contract_id = "mv04_topological_distance_bundle_v1",
    stratum_id = stratum_id,
    method_id = method_id,
    sample_ids = rownames(h0),
    input_diagram_ids = input_diagram_ids,
    h0_sha256 = digest::digest(h0, algo = "sha256", serialize = TRUE),
    h1_sha256 = digest::digest(h1, algo = "sha256", serialize = TRUE),
    combined_sha256 = digest::digest(
      combined, algo = "sha256", serialize = TRUE
    ),
    h0_scale_key = h0_scale$cache_key,
    h1_scale_key = h1_scale$cache_key
  )
  structure(
    list(
      contract_id = identity$contract_id,
      stratum_id = stratum_id,
      method_id = method_id,
      fit_scope_id = fit_scope_id,
      sample_ids = rownames(h0),
      input_diagram_ids = input_diagram_ids,
      matrices = list(H0 = h0, H1 = h1, combined = combined),
      scales = list(H0 = h0_scale, H1 = h1_scale),
      normalized = list(
        H0 = apply_distance_scale_v1(h0, h0_scale),
        H1 = apply_distance_scale_v1(h1, h1_scale)
      ),
      cache_key = paste0(
        "mv04_distance_bundle_v1:",
        digest::digest(identity, algo = "sha256", serialize = TRUE)
      )
    ),
    class = "scph_mv04_distance_bundle_v1"
  )
}

mv04_distance_contributions <- function(bundle) {
  if (!inherits(bundle, "scph_mv04_distance_bundle_v1")) {
    stop("bundle must be a scph_mv04_distance_bundle_v1 object.",
         call. = FALSE)
  }
  h0 <- bundle$matrices$H0
  h1 <- bundle$matrices$H1
  combined_squared <- h0 ^ 2 + h1 ^ 2
  indices <- which(lower.tri(h0), arr.ind = TRUE)
  data.frame(
    stratum_id = bundle$stratum_id,
    method_id = bundle$method_id,
    first_sample_id = rownames(h0)[indices[, "row"]],
    second_sample_id = colnames(h0)[indices[, "col"]],
    h0_distance = h0[indices],
    h1_distance = h1[indices],
    combined_distance = sqrt(combined_squared[indices]),
    h1_squared_fraction = ifelse(
      combined_squared[indices] > 0,
      h1[indices] ^ 2 / combined_squared[indices], 0
    ),
    stringsAsFactors = FALSE
  )
}
