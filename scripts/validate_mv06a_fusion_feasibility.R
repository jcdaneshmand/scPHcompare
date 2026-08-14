args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop(
    "Usage: validate_mv06a_fusion_feasibility.R <input-dir> <output-dir>",
    call. = FALSE
  )
}

input_dir <- normalizePath(args[[1L]], mustWork = TRUE)
output_dir <- normalizePath(args[[2L]], mustWork = TRUE)
contract_id <- "mv06a_label_closed_fusion_feasibility_v1"
fit_scope_id <- "mv06a_descriptive_all_pilot_samples"
weights_frozen <- c(0, 0.25, 0.5, 0.75, 1)

read_output <- function(name) {
  utils::read.csv(
    file.path(output_dir, name),
    stringsAsFactors = FALSE, check.names = FALSE
  )
}

assert_matrix <- function(value) {
  stopifnot(
    is.matrix(value), is.numeric(value), nrow(value) == ncol(value),
    identical(rownames(value), colnames(value)),
    all(is.finite(value)), all(value >= 0),
    max(abs(value - t(value))) <=
      100 * .Machine$double.eps * max(1, max(value)),
    all(abs(diag(value)) <=
      100 * .Machine$double.eps * max(1, max(value)))
  )
  value <- (value + t(value)) / 2
  diag(value) <- 0
  value
}

pair_index <- function(ids) {
  value <- utils::combn(ids, 2L)
  data.frame(
    first_sample_id = value[1L, ],
    second_sample_id = value[2L, ],
    stringsAsFactors = FALSE
  )
}

pair_values <- function(matrix, pairs) {
  matrix[cbind(
    match(pairs$first_sample_id, rownames(matrix)),
    match(pairs$second_sample_id, colnames(matrix))
  )]
}

weight_id <- function(weight) {
  paste0("gene_weight_", sprintf("%03d", as.integer(round(100 * weight))))
}

ordered_neighbors <- function(matrix, sample_id, k) {
  candidates <- setdiff(rownames(matrix), sample_id)
  candidates[order(matrix[sample_id, candidates], candidates, method = "radix")][
    seq_len(k)
  ]
}

compare_frame <- function(observed, expected, keys, tolerance = 1e-12) {
  stopifnot(setequal(names(observed), names(expected)))
  observed <- observed[names(expected)]
  order_args_observed <- lapply(keys, function(key) observed[[key]])
  order_args_expected <- lapply(keys, function(key) expected[[key]])
  observed <- observed[do.call(order, c(order_args_observed, list(method = "radix"))), , drop = FALSE]
  expected <- expected[do.call(order, c(order_args_expected, list(method = "radix"))), , drop = FALSE]
  rownames(observed) <- NULL
  rownames(expected) <- NULL
  stopifnot(nrow(observed) == nrow(expected))
  for (name in names(expected)) {
    if (is.numeric(expected[[name]])) {
      stopifnot(
        is.numeric(observed[[name]]),
        all(is.finite(observed[[name]])),
        max(abs(observed[[name]] - expected[[name]])) <= tolerance
      )
    } else if (is.logical(expected[[name]])) {
      stopifnot(identical(as.logical(observed[[name]]), expected[[name]]))
    } else {
      stopifnot(identical(as.character(observed[[name]]),
                          as.character(expected[[name]])))
    }
  }
  invisible(TRUE)
}

sources <- read_output("mv06a-source-inventory.csv")
stopifnot(
  nrow(sources) == 8L,
  all(sources$contract_id == contract_id),
  !anyDuplicated(sources$file)
)
for (index in seq_len(nrow(sources))) {
  path <- file.path(input_dir, sources$file[[index]])
  stopifnot(
    file.exists(path),
    as.numeric(file.info(path)$size) == sources$bytes[[index]],
    identical(
      digest::digest(path, algo = "sha256", file = TRUE, serialize = FALSE),
      sources$sha256[[index]]
    )
  )
}

strata <- c(
  "bone__integrated", "bone__sct_whole",
  "large__sct_whole", "large__seurat_integration"
)
expected_scales <- list()
expected_pairs <- list()
expected_weight_rows <- list()
expected_correlations <- list()
expected_neighbors <- list()
expected_hashes <- list()

for (stratum_index in seq_along(strata)) {
  stratum_id <- strata[[stratum_index]]
  cell <- readRDS(file.path(
    input_dir, paste0(stratum_id, "__cell_topology_v1.rds")
  ))
  gene <- readRDS(file.path(
    input_dir, paste0(stratum_id, "__gene_topology_v1.rds")
  ))
  stopifnot(
    identical(cell$method_id, "full_l2_exact_critical_pairs_v1"),
    identical(gene$method_id, "full_l2_exact_critical_pairs_v1"),
    identical(cell$stratum_id, paste0(stratum_id, "__cell_topology_v1")),
    identical(gene$stratum_id, paste0(stratum_id, "__gene_topology_v1")),
    identical(cell$sample_ids, gene$sample_ids),
    identical(sort(cell$sample_ids, method = "radix"), cell$sample_ids)
  )
  raw <- list(
    cell_H0 = assert_matrix(cell$matrices$H0),
    cell_H1 = assert_matrix(cell$matrices$H1),
    gene_H0 = assert_matrix(gene$matrices$H0),
    gene_H1 = assert_matrix(gene$matrices$H1)
  )
  scale_values <- vapply(raw, function(matrix) {
    value <- stats::median(matrix[lower.tri(matrix)])
    stopifnot(is.finite(value), value > sqrt(.Machine$double.eps))
    value
  }, numeric(1L))
  normalized <- Map(`/`, raw, scale_values)
  normalized <- lapply(normalized, assert_matrix)
  composites <- list(
    cell = assert_matrix(0.5 * normalized$cell_H0 +
                           0.5 * normalized$cell_H1),
    gene = assert_matrix(0.5 * normalized$gene_H0 +
                           0.5 * normalized$gene_H1)
  )
  fusion <- lapply(weights_frozen, function(weight) {
    assert_matrix((1 - weight) * composites$cell + weight * composites$gene)
  })
  names(fusion) <- vapply(weights_frozen, weight_id, character(1L))

  expected_scales[[stratum_index]] <- data.frame(
    contract_id = contract_id,
    stratum_id = stratum_id,
    component_id = names(scale_values),
    fit_scope_id = fit_scope_id,
    fit_sample_count = length(cell$sample_ids),
    scale = unname(scale_values),
    stringsAsFactors = FALSE
  )
  pairs <- pair_index(cell$sample_ids)
  pair_frame <- data.frame(
    contract_id = contract_id,
    stratum_id = stratum_id,
    pairs,
    cell_H0 = pair_values(normalized$cell_H0, pairs),
    cell_H1 = pair_values(normalized$cell_H1, pairs),
    gene_H0 = pair_values(normalized$gene_H0, pairs),
    gene_H1 = pair_values(normalized$gene_H1, pairs),
    cell_composite = pair_values(composites$cell, pairs),
    gene_composite = pair_values(composites$gene, pairs),
    equal_weight_fusion = pair_values(fusion$gene_weight_050, pairs),
    stringsAsFactors = FALSE
  )
  for (component_id in c("cell_H0", "cell_H1", "gene_H0", "gene_H1")) {
    pair_frame[[paste0(component_id, "_contribution")]] <-
      0.25 * pair_frame[[component_id]]
  }
  expected_pairs[[stratum_index]] <- pair_frame

  expected_weight_rows[[stratum_index]] <- do.call(
    rbind, lapply(seq_along(weights_frozen), function(weight_index) {
      weight <- weights_frozen[[weight_index]]
      data.frame(
        contract_id = contract_id,
        stratum_id = stratum_id,
        weight_id = weight_id(weight),
        gene_weight = weight,
        cell_weight = 1 - weight,
        pairs,
        distance = pair_values(fusion[[weight_index]], pairs),
        stringsAsFactors = FALSE
      )
    })
  )
  axes <- list(
    H0 = list(cell = normalized$cell_H0, gene = normalized$gene_H0),
    H1 = list(cell = normalized$cell_H1, gene = normalized$gene_H1),
    composite = composites
  )
  expected_correlations[[stratum_index]] <- do.call(
    rbind, lapply(names(axes), function(axis_id) {
      cell_values <- pair_values(axes[[axis_id]]$cell, pairs)
      gene_values <- pair_values(axes[[axis_id]]$gene, pairs)
      data.frame(
        contract_id = contract_id,
        stratum_id = stratum_id,
        axis_id = axis_id,
        pair_count = length(cell_values),
        pearson = stats::cor(cell_values, gene_values),
        spearman = stats::cor(cell_values, gene_values, method = "spearman"),
        mean_absolute_difference = mean(abs(cell_values - gene_values)),
        stringsAsFactors = FALSE
      )
    })
  )
  k <- min(3L, length(cell$sample_ids) - 1L)
  neighbor_rows <- list()
  neighbor_index <- 1L
  for (axis_id in names(axes)) {
    for (sample_id in cell$sample_ids) {
      cell_neighbors <- ordered_neighbors(
        axes[[axis_id]]$cell, sample_id, k
      )
      gene_neighbors <- ordered_neighbors(
        axes[[axis_id]]$gene, sample_id, k
      )
      overlap <- length(intersect(cell_neighbors, gene_neighbors))
      neighbor_rows[[neighbor_index]] <- data.frame(
        contract_id = contract_id,
        stratum_id = stratum_id,
        axis_id = axis_id,
        sample_id = sample_id,
        k = k,
        cell_neighbors = paste(cell_neighbors, collapse = ";"),
        gene_neighbors = paste(gene_neighbors, collapse = ";"),
        overlap_count = overlap,
        jaccard = overlap / length(union(cell_neighbors, gene_neighbors)),
        exact_neighbor_set = identical(
          sort(cell_neighbors), sort(gene_neighbors)
        ),
        stringsAsFactors = FALSE
      )
      neighbor_index <- neighbor_index + 1L
    }
  }
  expected_neighbors[[stratum_index]] <- do.call(rbind, neighbor_rows)
  matrices <- c(normalized, composites, fusion)
  matrix_ids <- c(
    paste0("normalized_", names(normalized)),
    paste0("composite_", names(composites)),
    paste0("fusion_", names(fusion))
  )
  expected_hashes[[stratum_index]] <- data.frame(
    contract_id = contract_id,
    stratum_id = stratum_id,
    matrix_id = matrix_ids,
    sample_count = length(cell$sample_ids),
    minimum = vapply(matrices, min, numeric(1L)),
    maximum = vapply(matrices, max, numeric(1L)),
    sha256 = vapply(
      matrices, digest::digest, character(1L),
      algo = "sha256", serialize = TRUE
    ),
    stringsAsFactors = FALSE
  )
}

observed_scales <- read_output("mv06a-scales.csv")
expected_scales <- do.call(rbind, expected_scales)
stopifnot(all(grepl("^distance_scale_v1:", observed_scales$scale_cache_key)))
observed_scales$scale_cache_key <- NULL
compare_frame(
  observed_scales, expected_scales,
  c("stratum_id", "component_id")
)
compare_frame(
  read_output("mv06a-pairs.csv"), do.call(rbind, expected_pairs),
  c("stratum_id", "first_sample_id", "second_sample_id")
)
compare_frame(
  read_output("mv06a-weight-grid.csv"),
  do.call(rbind, expected_weight_rows),
  c("stratum_id", "weight_id", "first_sample_id", "second_sample_id")
)
compare_frame(
  read_output("mv06a-correlations.csv"),
  do.call(rbind, expected_correlations),
  c("stratum_id", "axis_id")
)
compare_frame(
  read_output("mv06a-neighbors.csv"), do.call(rbind, expected_neighbors),
  c("stratum_id", "axis_id", "sample_id")
)
compare_frame(
  read_output("mv06a-matrix-hashes.csv"), do.call(rbind, expected_hashes),
  c("stratum_id", "matrix_id")
)

contract <- read_output("mv06a-contract.csv")
stopifnot(
  nrow(contract) == 1L,
  identical(contract$contract_id, contract_id),
  contract$stratum_count == 4L,
  contract$source_count == 8L,
  contract$scale_count == 16L,
  contract$unordered_pair_count == 102L,
  contract$weight_row_count == 510L,
  contract$correlation_count == 12L,
  contract$neighbor_row_count == 84L,
  contract$matrix_hash_count == 44L,
  contract$biological_values_read == 0L,
  contract$endpoint_rows_computed == 0L,
  contract$clusterings_computed == 0L,
  contract$advanced_fusion_computed == 0L
)

manifest <- read_output("mv06a-artifact-manifest.csv")
stopifnot(nrow(manifest) == 10L, !anyDuplicated(manifest$file))
for (index in seq_len(nrow(manifest))) {
  path <- file.path(output_dir, manifest$file[[index]])
  stopifnot(
    file.exists(path),
    as.numeric(file.info(path)$size) == manifest$bytes[[index]],
    identical(
      digest::digest(path, algo = "sha256", file = TRUE, serialize = FALSE),
      manifest$sha256[[index]]
    )
  )
}

forbidden_names <- c("tissue", "approach", "study", "outcome", "ari", "nmi")
all_outputs <- list.files(output_dir, pattern = "^mv06a-.*[.]csv$", full.names = TRUE)
for (path in all_outputs) {
  headers <- names(read_output(basename(path)))
  stopifnot(!any(tolower(headers) %in% forbidden_names))
}

cat("MV6-A independent validation PASS: 12 reconstructed categories\n")
