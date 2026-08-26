# Label-closed matched historical cell/gene partition agreement for MV11-G.

mv11g_select_common_partitions_v1 <- function(gene_partitions,
                                               cell_partitions) {
  required <- c(
    "stack_id", "representation_id", "panel_id", "seed",
    "homology_dimension", "method_id", "k", "sample_id", "cluster",
    "outcome_label_state", "biological_outcomes_computed"
  )
  for (value in list(gene_partitions, cell_partitions)) {
    if (!is.data.frame(value) || !all(required %in% names(value)) ||
        any(value$outcome_label_state != "closed") ||
        any(as.logical(value$biological_outcomes_computed))) {
      stop("MV11-G partition source is malformed or not label closed",
           call. = FALSE)
    }
    .mv10_assert_label_closed(value, "MV11-G partition source")
  }
  gene <- gene_partitions[
    gene_partitions$stack_id == "existing_selectedfit_data_exact500" &
      gene_partitions$k %in% c(2L, 3L), , drop = FALSE
  ]
  cell <- cell_partitions[
    cell_partitions$stack_id == "historical_selectedfit_cell_exact500" &
      cell_partitions$k %in% c(2L, 3L), , drop = FALSE
  ]
  methods <- mv10_method_registry_v1()
  methods <- methods$method_id[methods$authorized_for_mv10b]
  validate <- function(value, label) {
    key <- c("homology_dimension", "seed", "method_id", "k", "sample_id")
    if (nrow(value) != 12400L || anyDuplicated(value[key]) ||
        !identical(sort(unique(as.integer(value$seed))),
                   .mv11_required_seeds) ||
        !setequal(value$homology_dimension, c("H0", "H1")) ||
        !setequal(value$method_id, methods) ||
        !identical(sort(unique(as.integer(value$k))), c(2L, 3L)) ||
        any(table(value$homology_dimension, value$seed, value$method_id,
                  value$k) != 124L)) {
      stop("MV11-G ", label, " common-K axis drift", call. = FALSE)
    }
  }
  validate(gene, "gene"); validate(cell, "cell")
  if (!setequal(gene$sample_id, cell$sample_id)) {
    stop("MV11-G cell/gene sample axes differ", call. = FALSE)
  }
  list(gene = gene, cell = cell)
}

mv11g_cross_view_agreement_v1 <- function(gene_partitions,
                                           cell_partitions) {
  selected <- mv11g_select_common_partitions_v1(
    gene_partitions, cell_partitions
  )
  registry <- mv10_method_registry_v1()
  registry <- registry[registry$authorized_for_mv10b,
                       c("method_id", "role"), drop = FALSE]
  queue <- expand.grid(
    homology_dimension = c("H0", "H1"),
    method_id = registry$method_id,
    k = c(2L, 3L), seed = .mv11_required_seeds,
    stringsAsFactors = FALSE
  )
  queue <- queue[order(queue$homology_dimension,
                       match(queue$method_id, registry$method_id),
                       queue$k, queue$seed, method = "radix"), , drop = FALSE]
  seed_rows <- lapply(seq_len(nrow(queue)), function(index) {
    unit <- queue[index, , drop = FALSE]
    choose <- function(value) {
      result <- value[
        value$homology_dimension == unit$homology_dimension &
          value$method_id == unit$method_id & value$k == unit$k &
          value$seed == unit$seed, c("sample_id", "cluster"), drop = FALSE
      ]
      result[order(result$sample_id, method = "radix"), , drop = FALSE]
    }
    gene <- choose(selected$gene); cell <- choose(selected$cell)
    if (nrow(gene) != 124L || nrow(cell) != 124L ||
        !identical(gene$sample_id, cell$sample_id) ||
        anyDuplicated(gene$sample_id) || anyDuplicated(cell$sample_id)) {
      stop("MV11-G comparison unit sample axis drift", call. = FALSE)
    }
    adjusted_rand <- mclust::adjustedRandIndex(cell$cluster, gene$cluster)
    data.frame(
      contract_id = "mv11g_cross_view_seed_agreement_v1",
      homology_dimension = unit$homology_dimension,
      method_id = unit$method_id,
      method_role = registry$role[match(unit$method_id, registry$method_id)],
      k = as.integer(unit$k), seed = as.integer(unit$seed), samples = 124L,
      adjusted_rand = adjusted_rand,
      exact_partition_agreement = isTRUE(all.equal(adjusted_rand, 1,
                                                    tolerance = 1e-12)),
      label_state = "closed", outcome_state = "closed",
      comparison_semantics = "symmetric_agreement_not_view_ranking",
      stringsAsFactors = FALSE
    )
  })
  seed_agreement <- do.call(rbind, seed_rows); rownames(seed_agreement) <- NULL
  group <- interaction(
    seed_agreement$homology_dimension, seed_agreement$method_id,
    seed_agreement$k, drop = TRUE, lex.order = TRUE
  )
  summary_rows <- lapply(split(seed_agreement, group), function(value) {
    data.frame(
      contract_id = "mv11g_cross_view_summary_v1",
      homology_dimension = value$homology_dimension[[1L]],
      method_id = value$method_id[[1L]], method_role = value$method_role[[1L]],
      k = value$k[[1L]], seeds = nrow(value), samples_per_seed = 124L,
      mean_adjusted_rand = mean(value$adjusted_rand),
      median_adjusted_rand = stats::median(value$adjusted_rand),
      minimum_adjusted_rand = min(value$adjusted_rand),
      maximum_adjusted_rand = max(value$adjusted_rand),
      exact_agreement_seeds = sum(value$exact_partition_agreement),
      label_state = "closed", outcome_state = "closed",
      comparison_semantics = "symmetric_agreement_not_view_ranking",
      stringsAsFactors = FALSE
    )
  })
  summary <- do.call(rbind, summary_rows); rownames(summary) <- NULL
  summary <- summary[order(
    summary$homology_dimension,
    match(summary$method_id, registry$method_id), summary$k,
    method = "radix"
  ), , drop = FALSE]
  if (nrow(seed_agreement) != 100L || nrow(summary) != 20L ||
      any(!is.finite(seed_agreement$adjusted_rand)) ||
      any(seed_agreement$adjusted_rand < -1 |
            seed_agreement$adjusted_rand > 1) ||
      any(summary$seeds != 5L) || "sample_id" %in% names(seed_agreement) ||
      "sample_id" %in% names(summary)) {
    stop("MV11-G public agreement contract failed", call. = FALSE)
  }
  list(seed_agreement = seed_agreement, summary = summary)
}
