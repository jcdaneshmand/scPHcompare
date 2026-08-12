.corrected_consumer_bundle_class_v1 <-
  "scph_corrected_landscape_consumer_bundle_v1"
.corrected_average_trees_class_v1 <-
  "scph_corrected_landscape_average_trees_v1"

.one_nonempty_string_v1 <- function(value, label) {
  value <- as.character(value)
  if (length(value) != 1L || is.na(value) || !nzchar(value)) {
    stop(label, " must be one non-empty string.", call. = FALSE)
  }
  value
}

#' Read a verified corrected persistence-landscape sidecar
#'
#' Loads a completed corrected-landscape artifact tree without mutating it.
#' The caller must explicitly bind the iteration and topology view; view identity
#' is never guessed from paths or matrix dimensions.
#'
#' @param sidecar A sidecar list returned by the corrected-landscape workflow.
#' @param iteration_id One non-empty iteration identifier.
#' @param view_id Exactly `"cell_topology_v1"` or `"gene_topology_v1"`.
#'
#' @return A versioned read-only bundle with separate H0 and H1 matrices and a
#'   descriptive combined matrix.
#' @export
read_corrected_landscape_bundle <- function(sidecar, iteration_id, view_id) {
  required <- c("artifact_dir", "input_set_sha256", "matrix_cache_key")
  if (!is.list(sidecar) || !all(required %in% names(sidecar))) {
    stop("sidecar lacks required corrected-landscape identities.", call. = FALSE)
  }
  iteration_id <- .one_nonempty_string_v1(iteration_id, "iteration_id")
  view_id <- .one_nonempty_string_v1(view_id, "view_id")
  allowed_views <- c("cell_topology_v1", "gene_topology_v1")
  if (!view_id %in% allowed_views) {
    stop("view_id must be cell_topology_v1 or gene_topology_v1.", call. = FALSE)
  }
  artifact_dir <- .one_nonempty_string_v1(sidecar$artifact_dir, "sidecar$artifact_dir")
  if (!dir.exists(artifact_dir)) stop("sidecar artifact directory does not exist.", call. = FALSE)
  completion_path <- file.path(artifact_dir, "completion-v1.csv")
  if (!file.exists(completion_path)) stop("corrected sidecar is incomplete.", call. = FALSE)
  completion <- utils::read.csv(completion_path, stringsAsFactors = FALSE,
                                check.names = FALSE)
  .verify_completion_v1(artifact_dir, completion)
  input_manifest <- utils::read.csv(file.path(artifact_dir, "input-manifest-v1.csv"),
                                    stringsAsFactors = FALSE, check.names = FALSE)
  matrix_value <- readRDS(file.path(artifact_dir, "distance-matrix-v1.rds"))
  if (!inherits(matrix_value, .scph_landscape_matrix_class_v1) ||
      !identical(matrix_value$contract_id, "scph_public_landscape_matrix_v1") ||
      !identical(matrix_value$mode, "scientific") ||
      !identical(matrix_value$provenance$workflow_artifact_contract,
                 "scph_corrected_landscape_artifact_v1") ||
      isTRUE(matrix_value$provenance$legacy_reproduction)) {
    stop("corrected sidecar matrix contract is invalid.", call. = FALSE)
  }
  input_key <- .one_nonempty_string_v1(sidecar$input_set_sha256,
                                       "sidecar$input_set_sha256")
  matrix_key <- .one_nonempty_string_v1(sidecar$matrix_cache_key,
                                        "sidecar$matrix_cache_key")
  if (!all(input_manifest$input_set_sha256 == input_key)) {
    stop("sidecar input-set key does not match its manifest.", call. = FALSE)
  }
  if (!identical(matrix_value$cache_key, matrix_key)) {
    stop("sidecar matrix cache key does not match its matrix.", call. = FALSE)
  }
  matrices <- lapply(matrix_value$matrices[c("H0", "H1", "combined")],
                     .validate_distance_matrix_v1)
  if (length(matrices) != 3L || !identical(dimnames(matrices$H0), dimnames(matrices$H1)) ||
      !identical(dimnames(matrices$H0), dimnames(matrices$combined)) ||
      !identical(rownames(matrices$H0), as.character(input_manifest$source_id))) {
    stop("corrected matrix axes do not match the input manifest.", call. = FALSE)
  }
  expected_combined <- sqrt(matrices$H0 ^ 2 + matrices$H1 ^ 2)
  diag(expected_combined) <- 0
  tolerance <- 100 * .Machine$double.eps * max(1, max(expected_combined))
  if (max(abs(expected_combined - matrices$combined)) > tolerance) {
    stop("descriptive combined matrix is inconsistent with H0 and H1.", call. = FALSE)
  }
  identity <- list(contract_id = .corrected_consumer_bundle_class_v1,
    iteration_id = iteration_id, view_id = view_id, input_set_sha256 = input_key,
    source_matrix_cache_key = matrix_key, sample_ids = rownames(matrices$H0),
    h0_sha256 = digest::digest(matrices$H0, algo = "sha256", serialize = TRUE),
    h1_sha256 = digest::digest(matrices$H1, algo = "sha256", serialize = TRUE))
  structure(list(contract_id = identity$contract_id, result_version = 1L,
    iteration_id = iteration_id, view_id = view_id,
    scientific_contract = matrix_value$specification,
    sample_ids = identity$sample_ids, matrices = matrices,
    dimension_policy = "H0_H1_separate_combined_descriptive_not_consumed",
    source = list(artifact_dir = normalizePath(artifact_dir, winslash = "/"),
      input_set_sha256 = input_key, matrix_cache_key = matrix_key),
    cache_key = paste0("scph_corrected_consumer_bundle_v1:", digest::digest(
      identity, algo = "sha256", serialize = TRUE))),
    class = c(.corrected_consumer_bundle_class_v1, "list"))
}

#' Build separate average-linkage trees from corrected landscapes
#'
#' @param bundle A bundle returned by [read_corrected_landscape_bundle()].
#'
#' @return A versioned result containing independent H0 and H1 `hclust` trees.
#'   No partition or combined-distance tree is produced.
#' @export
corrected_landscape_average_trees <- function(bundle) {
  if (!inherits(bundle, .corrected_consumer_bundle_class_v1) ||
      !identical(bundle$contract_id, .corrected_consumer_bundle_class_v1) ||
      !identical(names(bundle$matrices), c("H0", "H1", "combined"))) {
    stop("bundle is not a valid corrected landscape consumer bundle.", call. = FALSE)
  }
  matrices <- lapply(bundle$matrices, .validate_distance_matrix_v1)
  trees <- lapply(matrices[c("H0", "H1")], function(value) {
    stats::hclust(stats::as.dist(value), method = "average")
  })
  identity <- list(contract_id = .corrected_average_trees_class_v1,
    source_bundle_cache_key = bundle$cache_key, iteration_id = bundle$iteration_id,
    view_id = bundle$view_id, dimensions = c("H0", "H1"), method = "average",
    h0_cophenetic_sha256 = digest::digest(stats::cophenetic(trees$H0),
      algo = "sha256", serialize = TRUE),
    h1_cophenetic_sha256 = digest::digest(stats::cophenetic(trees$H1),
      algo = "sha256", serialize = TRUE))
  structure(list(contract_id = identity$contract_id, result_version = 1L,
    iteration_id = bundle$iteration_id, view_id = bundle$view_id,
    method = "average", dimensions = c("H0", "H1"), trees = trees,
    combined_tree = NULL, partitions = NULL, selected_k = NULL,
    provenance = list(label_free = TRUE, outcome_free = TRUE,
      combined_consumed = FALSE, source_bundle_cache_key = bundle$cache_key),
    cache_key = paste0("scph_corrected_average_trees_v1:", digest::digest(
      identity, algo = "sha256", serialize = TRUE))),
    class = c(.corrected_average_trees_class_v1, "list"))
}
