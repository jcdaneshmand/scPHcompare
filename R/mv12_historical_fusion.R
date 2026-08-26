# Label-closed historical cell/gene distance-fusion feasibility for MV12.

.mv12_weights <- c(0, 0.25, 0.5, 0.75, 1)
.mv12_k <- c(2L, 3L)
.mv12_primary_method <- "pam_dissimilarity_v1"

.mv12_weight_id <- function(weight) {
  paste0("fusion_gene_weight_", sprintf("%03d", as.integer(100 * weight)))
}

mv12_normalize_distance_v1 <- function(matrix) {
  matrix <- .mv05n_validate_distance_matrix(matrix)
  values <- matrix[upper.tri(matrix)]
  positive <- values[is.finite(values) & values > 0]
  scale <- stats::median(positive)
  if (!length(positive) || !is.finite(scale) || scale <= 0) {
    stop("MV12 distance scale is not finite and positive", call. = FALSE)
  }
  list(matrix = .mv05n_validate_distance_matrix(matrix / scale), scale = scale,
       positive_pairs = length(positive), total_pairs = length(values),
       scale_rule = "median_positive_offdiagonal_full124_label_closed_v1")
}

mv12_build_fusion_matrices_v1 <- function(bundle) {
  mv11_validate_matrix_bundle_v1(bundle)
  matrices <- list(); scale_rows <- list(); catalog_rows <- list()
  matrix_cursor <- 0L; scale_cursor <- 0L
  for (dimension in c("H0", "H1")) {
    component <- paste0("cell_", dimension)
    gene_component <- paste0("gene_", dimension)
    if (!gene_component %in% names(bundle$seed_specific)) {
      stop("MV12 historical gene component is absent", call. = FALSE)
    }
    for (seed in .mv11_required_seeds) {
      cell <- bundle$seed_specific[[component]][[as.character(seed)]]
      cell <- .mv05n_validate_distance_matrix(cell)
      gene <- bundle$seed_specific[[gene_component]][[as.character(seed)]]
      gene <- .mv05n_validate_distance_matrix(gene)
      if (!identical(dimnames(cell), dimnames(gene))) {
        stop("MV12 cell/gene matrix axes differ", call. = FALSE)
      }
      normalized <- list(
        cell = mv12_normalize_distance_v1(cell),
        gene = mv12_normalize_distance_v1(gene)
      )
      for (view in names(normalized)) {
        scale_cursor <- scale_cursor + 1L
        scale_rows[[scale_cursor]] <- data.frame(
          contract_id = "mv12_distance_scale_v1",
          homology_dimension = dimension, seed = seed, view = view,
          scale = normalized[[view]]$scale,
          positive_pairs = normalized[[view]]$positive_pairs,
          total_pairs = normalized[[view]]$total_pairs,
          scale_rule = normalized[[view]]$scale_rule,
          labels_used = FALSE, outcomes_used = FALSE,
          stringsAsFactors = FALSE
        )
      }
      for (weight in .mv12_weights) {
        matrix_cursor <- matrix_cursor + 1L
        stack_id <- .mv12_weight_id(weight)
        value <- (1 - weight) * normalized$cell$matrix +
          weight * normalized$gene$matrix
        value <- .mv05n_validate_distance_matrix(value)
        key <- paste(dimension, seed, stack_id, sep = "__")
        matrices[[key]] <- value
        catalog_rows[[matrix_cursor]] <- data.frame(
          contract_id = "mv12_fusion_matrix_catalog_v1",
          matrix_order = matrix_cursor, matrix_key = key,
          stack_id = stack_id, homology_dimension = dimension, seed = seed,
          gene_weight = weight, cell_weight = 1 - weight,
          role = if (weight == 0.5) "fusion_primary" else if (
            weight %in% c(0.25, 0.75)) "fusion_sensitivity" else
              "component_comparator",
          samples = nrow(value), unordered_pairs = choose(nrow(value), 2L),
          H0_H1_combined = FALSE, labels_used = FALSE, outcomes_used = FALSE,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  scales <- do.call(rbind, scale_rows); rownames(scales) <- NULL
  catalog <- do.call(rbind, catalog_rows); rownames(catalog) <- NULL
  if (nrow(scales) != 20L || nrow(catalog) != 50L ||
      length(matrices) != 50L || any(catalog$samples != 124L) ||
      any(catalog$unordered_pairs != choose(124L, 2L))) {
    stop("MV12 fusion matrix cardinality drift", call. = FALSE)
  }
  list(matrices = matrices, scales = scales, catalog = catalog)
}

mv12_fit_fusion_grid_v1 <- function(bundle) {
  built <- mv12_build_fusion_matrices_v1(bundle)
  registry <- mv10_method_registry_v1()
  methods <- registry$method_id[registry$authorized_for_mv10b]
  assignments <- list(); quality <- list(); cursor <- 0L
  for (index in seq_len(nrow(built$catalog))) {
    metadata <- built$catalog[index, , drop = FALSE]
    matrix <- built$matrices[[metadata$matrix_key]]
    for (method in methods) for (k in .mv12_k) {
      cursor <- cursor + 1L
      fit <- mv10_fit_partition_v1(matrix, k, method)
      score <- mv10_partition_quality_v1(matrix, fit)
      common <- metadata[c(
        "matrix_order", "matrix_key", "stack_id", "homology_dimension",
        "seed", "gene_weight", "cell_weight", "role"
      )]
      assignments[[cursor]] <- cbind(
        common[rep(1L, nrow(fit)), , drop = FALSE], fit
      )
      quality[[cursor]] <- cbind(common, score)
    }
  }
  assignments <- do.call(rbind, assignments); rownames(assignments) <- NULL
  quality <- do.call(rbind, quality); rownames(quality) <- NULL
  if (nrow(assignments) != 62000L || nrow(quality) != 500L ||
      anyDuplicated(assignments[c(
        "stack_id", "homology_dimension", "seed", "method_id", "k",
        "sample_id"
      )])) {
    stop("MV12 fusion clustering grid drift", call. = FALSE)
  }
  list(assignments = assignments, quality = quality,
       scales = built$scales, catalog = built$catalog)
}

mv12_consensus_diagnostics_v1 <- function(assignments) {
  required <- c("stack_id", "homology_dimension", "seed", "method_id", "k",
                "sample_id", "cluster", "gene_weight")
  if (!is.data.frame(assignments) || nrow(assignments) != 62000L ||
      !all(required %in% names(assignments))) {
    stop("MV12 consensus assignments are malformed", call. = FALSE)
  }
  rows <- list(); cursor <- 0L
  methods <- mv10_method_registry_v1()
  methods <- methods$method_id[methods$authorized_for_mv10b]
  for (dimension in c("H0", "H1")) for (seed in .mv11_required_seeds) {
    for (method in methods) for (k in .mv12_k) {
      select <- function(weight) {
        value <- assignments[
          assignments$homology_dimension == dimension &
            assignments$seed == seed & assignments$method_id == method &
            assignments$k == k & abs(assignments$gene_weight - weight) < 1e-12,
          c("sample_id", "cluster"), drop = FALSE
        ]
        value[order(value$sample_id, method = "radix"), , drop = FALSE]
      }
      cell <- select(0); gene <- select(1)
      if (nrow(cell) != 124L || !identical(cell$sample_id, gene$sample_id)) {
        stop("MV12 component partition axis drift", call. = FALSE)
      }
      baseline <- mclust::adjustedRandIndex(cell$cluster, gene$cluster)
      for (weight in c(0.25, 0.5, 0.75)) {
        fusion <- select(weight)
        if (!identical(cell$sample_id, fusion$sample_id)) {
          stop("MV12 fusion partition axis drift", call. = FALSE)
        }
        cell_ari <- mclust::adjustedRandIndex(fusion$cluster, cell$cluster)
        gene_ari <- mclust::adjustedRandIndex(fusion$cluster, gene$cluster)
        cursor <- cursor + 1L
        rows[[cursor]] <- data.frame(
          contract_id = "mv12_consensus_diagnostic_v1",
          homology_dimension = dimension, seed = seed, method_id = method,
          k = k, gene_weight = weight,
          cell_gene_adjusted_rand = baseline,
          fusion_cell_adjusted_rand = cell_ari,
          fusion_gene_adjusted_rand = gene_ari,
          minimum_component_adjusted_rand = min(cell_ari, gene_ari),
          mean_component_adjusted_rand = mean(c(cell_ari, gene_ari)),
          balanced_gain_over_cell_gene = min(cell_ari, gene_ari) - baseline,
          labels_used = FALSE, outcomes_used = FALSE,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  result <- do.call(rbind, rows); rownames(result) <- NULL
  if (nrow(result) != 300L || any(!is.finite(unlist(result[c(
    "cell_gene_adjusted_rand", "fusion_cell_adjusted_rand",
    "fusion_gene_adjusted_rand", "minimum_component_adjusted_rand",
    "mean_component_adjusted_rand", "balanced_gain_over_cell_gene"
  )])))) stop("MV12 consensus diagnostic drift", call. = FALSE)
  result
}

mv12_fusion_decision_v1 <- function(stability, consensus) {
  required_stability <- c("stack_id", "homology_dimension", "method_id", "k",
                          "mean_adjusted_rand")
  if (!is.data.frame(stability) || nrow(stability) != 100L ||
      !all(required_stability %in% names(stability)) ||
      !is.data.frame(consensus) || nrow(consensus) != 300L) {
    stop("MV12 decision inputs are malformed", call. = FALSE)
  }
  pam_stability <- stability[
    stability$homology_dimension == "H1" &
      stability$method_id == .mv12_primary_method &
      stability$k %in% .mv12_k, , drop = FALSE
  ]
  detail <- lapply(.mv12_k, function(k) {
    value <- pam_stability[pam_stability$k == k, , drop = FALSE]
    lookup <- function(weight) value$mean_adjusted_rand[
      value$stack_id == .mv12_weight_id(weight)
    ]
    cell_stability <- lookup(0); gene_stability <- lookup(1)
    equal_stability <- lookup(0.5)
    equal_consensus <- consensus[
      consensus$homology_dimension == "H1" &
        consensus$method_id == .mv12_primary_method & consensus$k == k &
        abs(consensus$gene_weight - 0.5) < 1e-12, , drop = FALSE
    ]
    if (length(cell_stability) != 1L || length(gene_stability) != 1L ||
        length(equal_stability) != 1L || nrow(equal_consensus) != 5L) {
      stop("MV12 primary decision axis drift", call. = FALSE)
    }
    balanced_positive <- sum(equal_consensus$balanced_gain_over_cell_gene > 0)
    data.frame(
      contract_id = "mv12_primary_fusion_decision_detail_v1", k = k,
      cell_mean_seed_stability = cell_stability,
      gene_mean_seed_stability = gene_stability,
      equal_fusion_mean_seed_stability = equal_stability,
      stability_gain_over_better_component = equal_stability -
        max(cell_stability, gene_stability),
      balanced_gain_positive_seeds = balanced_positive,
      required_positive_seeds = 3L,
      primary_k_pass = equal_stability > max(cell_stability, gene_stability) &&
        balanced_positive >= 3L,
      stringsAsFactors = FALSE
    )
  })
  detail <- do.call(rbind, detail); rownames(detail) <- NULL
  sensitivity_signal <- FALSE
  for (k in .mv12_k) for (weight in c(0.25, 0.75)) {
    value <- pam_stability[pam_stability$k == k, , drop = FALSE]
    fusion_stability <- value$mean_adjusted_rand[
      value$stack_id == .mv12_weight_id(weight)
    ]
    components <- value$mean_adjusted_rand[value$stack_id %in% c(
      .mv12_weight_id(0), .mv12_weight_id(1)
    )]
    diagnostic <- consensus[
      consensus$homology_dimension == "H1" &
        consensus$method_id == .mv12_primary_method & consensus$k == k &
        abs(consensus$gene_weight - weight) < 1e-12, , drop = FALSE
    ]
    sensitivity_signal <- sensitivity_signal ||
      (length(fusion_stability) == 1L && length(components) == 2L &&
         fusion_stability > max(components) &&
         sum(diagnostic$balanced_gain_over_cell_gene > 0) >= 3L)
  }
  passes <- sum(detail$primary_k_pass)
  disposition <- if (passes == 2L) {
    "credible_equal_weight_signal_both_H1_common_K"
  } else if (passes == 1L) {
    "targeted_equal_weight_signal_one_H1_common_K"
  } else if (sensitivity_signal) {
    "weight_sensitive_ambiguous_historical_signal"
  } else {
    "clear_historical_fusion_negative"
  }
  decision <- data.frame(
    contract_id = "mv12_fusion_disposition_v1", disposition = disposition,
    primary_K_passes = passes,
    sensitivity_signal = sensitivity_signal,
    option2_new_allqc_cell_topology_required =
      disposition != "clear_historical_fusion_negative",
    labels_used = FALSE, outcomes_used = FALSE,
    method_or_weight_selected = FALSE, biological_claims = FALSE,
    manuscript_claims = FALSE, stringsAsFactors = FALSE
  )
  list(detail = detail, decision = decision)
}
