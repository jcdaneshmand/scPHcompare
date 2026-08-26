# Prospectively qualified null generators for MV17-B.

.mv17b_matrix <- function(x) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  if (nrow(x) < 4L || ncol(x) < 2L || anyNA(x) || any(!is.finite(x))) {
    stop("MV17-B requires a finite matrix with >=4 points and >=2 coordinates",
         call. = FALSE)
  }
  x
}

mv17b_coordinate_permutation_v1 <- function(x, seed) {
  x <- .mv17b_matrix(x); set.seed(as.integer(seed))
  out <- vapply(seq_len(ncol(x)), function(j) sample(x[, j], replace = FALSE),
                numeric(nrow(x)))
  dimnames(out) <- dimnames(x); out
}

mv17b_covariance_gaussian_v1 <- function(x, seed, tolerance = 1e-12) {
  x <- .mv17b_matrix(x); n <- nrow(x); d <- ncol(x)
  if (n <= d) stop("MV17-B exact covariance Gaussian requires points > coordinates")
  target <- stats::cov(x); center <- colMeans(x)
  set.seed(as.integer(seed)); z <- matrix(stats::rnorm(n * d), n, d)
  z <- sweep(z, 2L, colMeans(z), "-")
  ze <- eigen(crossprod(z) / (n - 1), symmetric = TRUE)
  te <- eigen(target, symmetric = TRUE)
  if (min(ze$values) <= tolerance || min(te$values) < -tolerance) {
    stop("MV17-B covariance is singular beyond the frozen tolerance")
  }
  zw <- ze$vectors %*% diag(1 / sqrt(ze$values), d) %*% t(ze$vectors)
  ts <- te$vectors %*% diag(sqrt(pmax(te$values, 0)), d) %*% t(te$vectors)
  out <- z %*% zw %*% ts
  out <- sweep(out, 2L, center, "+"); dimnames(out) <- dimnames(x); out
}

mv17b_radial_density_cloud_v1 <- function(x, seed) {
  x <- .mv17b_matrix(x); n <- nrow(x); d <- ncol(x); center <- colMeans(x)
  radius <- sqrt(rowSums(sweep(x, 2L, center, "-")^2))
  set.seed(as.integer(seed)); direction <- matrix(stats::rnorm(n * d), n, d)
  norm <- sqrt(rowSums(direction^2))
  if (any(norm == 0)) stop("MV17-B generated a zero direction")
  direction <- direction / norm
  out <- direction * sample(radius, replace = FALSE)
  out <- sweep(out, 2L, center, "+"); dimnames(out) <- dimnames(x); out
}

mv17b_within_row_axis_shuffle_v1 <- function(x, seed) {
  x <- .mv17b_matrix(x); set.seed(as.integer(seed))
  out <- t(apply(x, 1L, sample, replace = FALSE))
  dimnames(out) <- dimnames(x); out
}

mv17b_null_registry_v1 <- function() {
  data.frame(
    contract_id = "mv17b_null_family_v1",
    null_family = c("coordinate_permutation", "covariance_gaussian",
                    "radial_density_cloud", "within_row_axis_shuffle"),
    function_name = c("mv17b_coordinate_permutation_v1",
                      "mv17b_covariance_gaussian_v1",
                      "mv17b_radial_density_cloud_v1",
                      "mv17b_within_row_axis_shuffle_v1"),
    compatible_view = c("cell_and_gene", "cell_and_gene", "cell_and_gene",
                        "gene_only"),
    preserves = c("point_count;coordinate_count;coordinatewise_marginals",
                  "point_count;coordinate_count;mean;sample_covariance",
                  "point_count;coordinate_count;center;radial_multiset",
                  "point_count;coordinate_count;within_point_value_multiset"),
    destroys = c("cross_coordinate_dependence;local_neighborhoods",
                 "higher_moments;non_gaussian_geometry;local_topology",
                 "angular_structure;coordinate_dependence;local_topology",
                 "coordinate_identity;between_point_feature_alignment"),
    primary_invariant = c("sorted_coordinate_values", "mean_and_covariance",
                          "sorted_centered_radii", "sorted_row_values"),
    real_corpus_authorized = FALSE, stringsAsFactors = FALSE
  )
}

mv17b_fixture_registry_v1 <- function() {
  expand.grid(
    fixture = c("correlated_gaussian", "circle_in_3d", "independent_gaussian"),
    seed = c(17001L, 17002L, 17003L), stringsAsFactors = FALSE
  ) |>
    transform(contract_id = "mv17b_fixture_v1", points = 128L,
              coordinates = 3L, values_inspected = FALSE,
              execution_authorized = FALSE)
}
