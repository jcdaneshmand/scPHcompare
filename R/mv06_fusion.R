.mv06_weight_id <- function(weight) {
  if (!is.numeric(weight) || length(weight) != 1L || !is.finite(weight) ||
      weight < 0 || weight > 1) {
    stop("weight must be one finite value in [0, 1].", call. = FALSE)
  }
  paste0("gene_weight_", sprintf("%03d", as.integer(round(100 * weight))))
}

.mv06_validate_bundle <- function(bundle, expected_method_id) {
  if (!inherits(bundle, "scph_mv04_distance_bundle_v1") ||
      !identical(bundle$contract_id, "mv04_topological_distance_bundle_v1") ||
      !identical(bundle$method_id, expected_method_id) ||
      !is.character(bundle$sample_ids) || length(bundle$sample_ids) < 3L ||
      anyNA(bundle$sample_ids) || any(!nzchar(bundle$sample_ids)) ||
      anyDuplicated(bundle$sample_ids) ||
      !identical(sort(bundle$sample_ids, method = "radix"),
                 bundle$sample_ids) ||
      !all(c("H0", "H1") %in% names(bundle$matrices))) {
    stop("bundle does not satisfy the frozen MV-04 identity contract.",
         call. = FALSE)
  }
  h0 <- .validate_distance_matrix_v1(bundle$matrices$H0, "bundle H0")
  h1 <- .validate_distance_matrix_v1(bundle$matrices$H1, "bundle H1")
  if (!identical(dimnames(h0), dimnames(h1)) ||
      !identical(rownames(h0), bundle$sample_ids)) {
    stop("bundle H0/H1 sample axes do not match its declared sample IDs.",
         call. = FALSE)
  }
  list(H0 = h0, H1 = h1)
}

mv06_fit_fusion_components_v1 <- function(
    cell_bundle, gene_bundle, fit_sample_ids = NULL,
    fit_scope_id = "mv06a_descriptive_all_pilot_samples",
    gene_weights = c(0, 0.25, 0.5, 0.75, 1)) {
  cell <- .mv06_validate_bundle(cell_bundle, "cell_topology_v1")
  gene <- .mv06_validate_bundle(gene_bundle, "gene_topology_v1")
  if (!identical(cell_bundle$stratum_id, gene_bundle$stratum_id) ||
      !identical(cell_bundle$sample_ids, gene_bundle$sample_ids) ||
      !identical(dimnames(cell$H0), dimnames(gene$H0))) {
    stop("cell and gene bundles must have one identical stratum/sample axis.",
         call. = FALSE)
  }
  gene_weights <- as.numeric(gene_weights)
  frozen_weights <- c(0, 0.25, 0.5, 0.75, 1)
  if (!identical(gene_weights, frozen_weights)) {
    stop("gene_weights must equal the complete frozen MV6-A grid.",
         call. = FALSE)
  }
  if (is.null(fit_sample_ids)) fit_sample_ids <- cell_bundle$sample_ids
  fit_sample_ids <- sort(as.character(fit_sample_ids), method = "radix")

  raw <- list(
    cell_H0 = cell$H0,
    cell_H1 = cell$H1,
    gene_H0 = gene$H0,
    gene_H1 = gene$H1
  )
  scales <- lapply(names(raw), function(component_id) {
    fit_distance_scale_v1(
      raw[[component_id]],
      fit_sample_ids = fit_sample_ids,
      fit_scope_id = fit_scope_id
    )
  })
  names(scales) <- names(raw)
  normalized <- Map(apply_distance_scale_v1, raw, scales)
  cell_composite <- 0.5 * normalized$cell_H0 +
    0.5 * normalized$cell_H1
  gene_composite <- 0.5 * normalized$gene_H0 +
    0.5 * normalized$gene_H1
  cell_composite <- .validate_distance_matrix_v1(
    cell_composite, "cell composite"
  )
  gene_composite <- .validate_distance_matrix_v1(
    gene_composite, "gene composite"
  )
  fusion <- lapply(gene_weights, function(weight) {
    value <- (1 - weight) * cell_composite + weight * gene_composite
    .validate_distance_matrix_v1(value, "fusion matrix")
  })
  names(fusion) <- vapply(gene_weights, .mv06_weight_id, character(1L))

  identity <- list(
    contract_id = "mv06a_label_closed_fusion_feasibility_v1",
    stratum_id = cell_bundle$stratum_id,
    sample_ids = cell_bundle$sample_ids,
    fit_scope_id = fit_scope_id,
    fit_sample_ids = fit_sample_ids,
    cell_bundle_key = cell_bundle$cache_key,
    gene_bundle_key = gene_bundle$cache_key,
    scale_keys = vapply(scales, `[[`, character(1L), "cache_key"),
    gene_weights = gene_weights,
    formula_id = "four_component_convex_arithmetic_v1"
  )
  structure(
    list(
      contract_id = identity$contract_id,
      stratum_id = identity$stratum_id,
      sample_ids = identity$sample_ids,
      fit_scope_id = fit_scope_id,
      fit_sample_ids = fit_sample_ids,
      raw = raw,
      scales = scales,
      normalized = normalized,
      composites = list(cell = cell_composite, gene = gene_composite),
      gene_weights = gene_weights,
      fusion = fusion,
      cache_key = paste0(
        "mv06a_fusion_v1:",
        digest::digest(identity, algo = "sha256", serialize = TRUE)
      )
    ),
    class = "scph_mv06a_fusion_v1"
  )
}

.mv06_pair_index <- function(sample_ids) {
  pairs <- utils::combn(sample_ids, 2L)
  data.frame(
    first_sample_id = pairs[1L, ],
    second_sample_id = pairs[2L, ],
    stringsAsFactors = FALSE
  )
}

.mv06_pair_values <- function(matrix, pairs) {
  matrix[cbind(
    match(pairs$first_sample_id, rownames(matrix)),
    match(pairs$second_sample_id, colnames(matrix))
  )]
}

mv06_pair_diagnostics_v1 <- function(fit) {
  if (!inherits(fit, "scph_mv06a_fusion_v1")) {
    stop("fit must be a scph_mv06a_fusion_v1 object.", call. = FALSE)
  }
  pairs <- .mv06_pair_index(fit$sample_ids)
  result <- data.frame(
    contract_id = fit$contract_id,
    stratum_id = fit$stratum_id,
    pairs,
    cell_H0 = .mv06_pair_values(fit$normalized$cell_H0, pairs),
    cell_H1 = .mv06_pair_values(fit$normalized$cell_H1, pairs),
    gene_H0 = .mv06_pair_values(fit$normalized$gene_H0, pairs),
    gene_H1 = .mv06_pair_values(fit$normalized$gene_H1, pairs),
    cell_composite = .mv06_pair_values(fit$composites$cell, pairs),
    gene_composite = .mv06_pair_values(fit$composites$gene, pairs),
    equal_weight_fusion = .mv06_pair_values(
      fit$fusion$gene_weight_050, pairs
    ),
    stringsAsFactors = FALSE
  )
  for (component_id in c("cell_H0", "cell_H1", "gene_H0", "gene_H1")) {
    result[[paste0(component_id, "_contribution")]] <-
      0.25 * result[[component_id]]
  }
  result
}

mv06_weight_diagnostics_v1 <- function(fit) {
  if (!inherits(fit, "scph_mv06a_fusion_v1")) {
    stop("fit must be a scph_mv06a_fusion_v1 object.", call. = FALSE)
  }
  pairs <- .mv06_pair_index(fit$sample_ids)
  rows <- lapply(seq_along(fit$gene_weights), function(index) {
    weight <- fit$gene_weights[[index]]
    matrix <- fit$fusion[[index]]
    data.frame(
      contract_id = fit$contract_id,
      stratum_id = fit$stratum_id,
      weight_id = .mv06_weight_id(weight),
      gene_weight = weight,
      cell_weight = 1 - weight,
      pairs,
      distance = .mv06_pair_values(matrix, pairs),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

mv06_correlation_diagnostics_v1 <- function(fit) {
  if (!inherits(fit, "scph_mv06a_fusion_v1")) {
    stop("fit must be a scph_mv06a_fusion_v1 object.", call. = FALSE)
  }
  pairs <- .mv06_pair_index(fit$sample_ids)
  axes <- list(
    H0 = list(cell = fit$normalized$cell_H0,
              gene = fit$normalized$gene_H0),
    H1 = list(cell = fit$normalized$cell_H1,
              gene = fit$normalized$gene_H1),
    composite = fit$composites
  )
  do.call(rbind, lapply(names(axes), function(axis_id) {
    cell <- .mv06_pair_values(axes[[axis_id]]$cell, pairs)
    gene <- .mv06_pair_values(axes[[axis_id]]$gene, pairs)
    data.frame(
      contract_id = fit$contract_id,
      stratum_id = fit$stratum_id,
      axis_id = axis_id,
      pair_count = length(cell),
      pearson = stats::cor(cell, gene, method = "pearson"),
      spearman = stats::cor(cell, gene, method = "spearman"),
      mean_absolute_difference = mean(abs(cell - gene)),
      stringsAsFactors = FALSE
    )
  }))
}

.mv06_ordered_neighbors <- function(matrix, sample_id, k) {
  candidates <- setdiff(rownames(matrix), sample_id)
  order_index <- order(
    matrix[sample_id, candidates], candidates,
    method = "radix"
  )
  candidates[order_index][seq_len(k)]
}

mv06_neighbor_diagnostics_v1 <- function(fit) {
  if (!inherits(fit, "scph_mv06a_fusion_v1")) {
    stop("fit must be a scph_mv06a_fusion_v1 object.", call. = FALSE)
  }
  axes <- list(
    H0 = list(cell = fit$normalized$cell_H0,
              gene = fit$normalized$gene_H0),
    H1 = list(cell = fit$normalized$cell_H1,
              gene = fit$normalized$gene_H1),
    composite = fit$composites
  )
  k <- min(3L, length(fit$sample_ids) - 1L)
  rows <- list()
  index <- 1L
  for (axis_id in names(axes)) {
    for (sample_id in fit$sample_ids) {
      cell <- .mv06_ordered_neighbors(
        axes[[axis_id]]$cell, sample_id, k
      )
      gene <- .mv06_ordered_neighbors(
        axes[[axis_id]]$gene, sample_id, k
      )
      overlap <- length(intersect(cell, gene))
      union <- length(union(cell, gene))
      rows[[index]] <- data.frame(
        contract_id = fit$contract_id,
        stratum_id = fit$stratum_id,
        axis_id = axis_id,
        sample_id = sample_id,
        k = k,
        cell_neighbors = paste(cell, collapse = ";"),
        gene_neighbors = paste(gene, collapse = ";"),
        overlap_count = overlap,
        jaccard = overlap / union,
        exact_neighbor_set = identical(sort(cell), sort(gene)),
        stringsAsFactors = FALSE
      )
      index <- index + 1L
    }
  }
  do.call(rbind, rows)
}

mv06_scale_diagnostics_v1 <- function(fit) {
  if (!inherits(fit, "scph_mv06a_fusion_v1")) {
    stop("fit must be a scph_mv06a_fusion_v1 object.", call. = FALSE)
  }
  data.frame(
    contract_id = fit$contract_id,
    stratum_id = fit$stratum_id,
    component_id = names(fit$scales),
    fit_scope_id = fit$fit_scope_id,
    fit_sample_count = length(fit$fit_sample_ids),
    scale = vapply(fit$scales, `[[`, numeric(1L), "scale"),
    scale_cache_key = vapply(
      fit$scales, `[[`, character(1L), "cache_key"
    ),
    stringsAsFactors = FALSE
  )
}

mv06_matrix_hashes_v1 <- function(fit) {
  if (!inherits(fit, "scph_mv06a_fusion_v1")) {
    stop("fit must be a scph_mv06a_fusion_v1 object.", call. = FALSE)
  }
  matrices <- c(fit$normalized, fit$composites, fit$fusion)
  ids <- c(
    paste0("normalized_", names(fit$normalized)),
    paste0("composite_", names(fit$composites)),
    paste0("fusion_", names(fit$fusion))
  )
  data.frame(
    contract_id = fit$contract_id,
    stratum_id = fit$stratum_id,
    matrix_id = ids,
    sample_count = length(fit$sample_ids),
    minimum = vapply(matrices, min, numeric(1L)),
    maximum = vapply(matrices, max, numeric(1L)),
    sha256 = vapply(
      matrices, digest::digest, character(1L),
      algo = "sha256", serialize = TRUE
    ),
    stringsAsFactors = FALSE
  )
}
