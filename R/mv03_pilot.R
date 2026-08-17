.mv03_variance_floor <- .Machine$double.eps

canonical_mv03_gene_ids <- function(ids) {
  ids <- as.character(ids)
  if (anyNA(ids) || any(!nzchar(ids))) {
    stop("Gene identifiers must be non-missing and non-empty.", call. = FALSE)
  }
  sub("-ENSG[0-9]+(?:\\.[0-9]+)?$", "", ids, perl = TRUE)
}

mv03_feature_category <- function(gene_ids) {
  genes <- canonical_mv03_gene_ids(gene_ids)
  result <- rep("retained_candidate", length(genes))
  result[grepl("^MT-", genes)] <- "mitochondrial_^MT-"
  result[grepl("^RP[SL]", genes)] <- "ribosomal_protein_^RP[SL]"
  result[grepl("^HB(?!P)", genes, perl = TRUE)] <-
    "hemoglobin_^HB(?!P)"
  result
}

.mv03_canonical_matrix <- function(value, label) {
  if ((!is.matrix(value) && !inherits(value, "Matrix")) ||
      is.null(rownames(value)) || is.null(colnames(value))) {
    stop(label, " must be a named matrix.", call. = FALSE)
  }
  genes <- canonical_mv03_gene_ids(rownames(value))
  if (anyDuplicated(genes)) {
    stop(label, " has duplicated canonical gene symbols.", call. = FALSE)
  }
  rownames(value) <- genes
  value
}

.mv03_row_variance <- function(value) {
  count <- ncol(value)
  if (count < 2L) {
    return(rep(NA_real_, nrow(value)))
  }
  if (inherits(value, "Matrix")) {
    means <- Matrix::rowMeans(value)
    centered_ss <- Matrix::rowSums(value * value) - count * means * means
  } else {
    means <- rowMeans(value)
    centered_ss <- rowSums(value * value) - count * means * means
  }
  pmax(as.numeric(centered_ss) / (count - 1), 0)
}

fit_mv03_gene_panel <- function(residual_samples, integrated_samples,
                                cohort, fit_scope_id,
                                panel_size = 500L) {
  if (!is.list(residual_samples) || !is.list(integrated_samples) ||
      is.null(names(residual_samples)) ||
      !identical(sort(names(residual_samples)),
                 sort(names(integrated_samples)))) {
    stop("Residual and integrated samples must be paired named lists.",
         call. = FALSE)
  }
  panel_size <- as.integer(panel_size)
  sample_ids <- sort(names(residual_samples), method = "radix")
  residual_samples <- lapply(sample_ids, function(id) {
    .mv03_canonical_matrix(residual_samples[[id]], paste0(id, " residual"))
  })
  names(residual_samples) <- sample_ids
  integrated_samples <- lapply(sample_ids, function(id) {
    .mv03_canonical_matrix(integrated_samples[[id]], paste0(id, " integrated"))
  })
  names(integrated_samples) <- sample_ids

  common <- Reduce(intersect, c(
    lapply(residual_samples, rownames), lapply(integrated_samples, rownames)
  ))
  common <- sort(common, method = "radix")
  if (length(common) < panel_size) {
    stop("Fewer than panel_size genes are common to every pilot matrix.",
         call. = FALSE)
  }

  rank_columns <- list()
  eligibility <- data.frame(
    gene = common,
    feature_category = mv03_feature_category(common),
    finite_all_samples_representations = TRUE,
    variance_above_floor_all_samples_representations = TRUE,
    stringsAsFactors = FALSE
  )
  for (sample_id in sample_ids) {
    for (representation in c("sct_whole", "integrated")) {
      value <- if (representation == "sct_whole") {
        residual_samples[[sample_id]][common, , drop = FALSE]
      } else {
        integrated_samples[[sample_id]][common, , drop = FALSE]
      }
      variance <- .mv03_row_variance(value)
      finite <- is.finite(variance)
      eligible_variance <- finite & variance > .mv03_variance_floor
      eligibility$finite_all_samples_representations <-
        eligibility$finite_all_samples_representations & finite
      eligibility$variance_above_floor_all_samples_representations <-
        eligibility$variance_above_floor_all_samples_representations &
        eligible_variance
      if (representation == "sct_whole") {
        rank_columns[[sample_id]] <- rank(
          -variance, ties.method = "min", na.last = "keep"
        )
      }
    }
  }
  rank_matrix <- do.call(cbind, rank_columns)
  rownames(rank_matrix) <- common
  median_rank <- apply(rank_matrix, 1L, stats::median, na.rm = TRUE)
  eligibility$median_sct_variance_rank <- median_rank
  eligibility$eligible <-
    eligibility$feature_category == "retained_candidate" &
    eligibility$finite_all_samples_representations &
    eligibility$variance_above_floor_all_samples_representations
  eligible <- eligibility[eligibility$eligible, , drop = FALSE]
  eligible <- eligible[order(
    eligible$median_sct_variance_rank, eligible$gene, method = "radix"
  ), , drop = FALSE]
  if (nrow(eligible) < panel_size) {
    stop("Fewer than panel_size genes pass the frozen eligibility rules.",
         call. = FALSE)
  }
  panel <- eligible$gene[seq_len(panel_size)]
  rank_table <- data.frame(
    cohort = cohort,
    fit_scope_id = fit_scope_id,
    gene = common,
    rank_matrix,
    median_sct_variance_rank = median_rank,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  result <- list(
    contract_id = "mv03_gene_panel_v1",
    cohort = cohort,
    fit_scope_id = fit_scope_id,
    panel_size = panel_size,
    panel = panel,
    panel_sha256 = digest::digest(panel, algo = "sha256", serialize = TRUE),
    variance_floor = .mv03_variance_floor,
    feature_rule_id = "technical_feature_rules_v1",
    eligibility = eligibility,
    ranks = rank_table,
    residual_samples = residual_samples,
    integrated_samples = integrated_samples
  )
  class(result) <- "scph_mv03_gene_panel_v1"
  result
}

prepare_mv03_sources <- function(samples, panel, cohort, representation,
                                 fit_scope_id, seed,
                                 selected_cells = NULL,
                                 contract_profile = "scientific",
                                 expected_genes = NULL,
                                 expected_cells = NULL,
                                 expected_pcs = NULL) {
  sample_ids <- sort(names(samples), method = "radix")
  if (is.null(selected_cells)) {
    selected_cells <- lapply(sample_ids, function(sample_id) {
      select_matched_cells(colnames(samples[[sample_id]]), n = 384L,
                           seed = seed)
    })
    names(selected_cells) <- sample_ids
  }
  selected <- lapply(sample_ids, function(sample_id) {
    value <- samples[[sample_id]]
    if (!all(panel %in% rownames(value)) ||
        !all(selected_cells[[sample_id]] %in% colnames(value))) {
      stop("Panel or matched cells are unavailable for ", sample_id,
           call. = FALSE)
    }
    as.matrix(value[panel, selected_cells[[sample_id]], drop = FALSE])
  })
  names(selected) <- sample_ids
  pooled <- do.call(cbind, selected)
  center <- rowMeans(pooled)
  scale <- apply(pooled, 1L, stats::sd)
  if (any(!is.finite(center)) || any(!is.finite(scale)) ||
      any(scale <= sqrt(.Machine$double.eps))) {
    stop("Fit-scope standardization produced invalid gene parameters.",
         call. = FALSE)
  }
  standardization_identity <- list(
    contract_id = "mv03_fit_scope_standardization_v1",
    cohort = cohort,
    representation = representation,
    fit_scope_id = fit_scope_id,
    seed = as.integer(seed),
    sample_ids = sample_ids,
    panel = panel,
    selected_cells = selected_cells,
    center = center,
    scale = scale
  )
  standardization_id <- paste0(
    "mv03_standardization_v1:",
    digest::digest(standardization_identity, algo = "sha256", serialize = TRUE)
  )
  standardized <- lapply(selected, function(value) {
    sweep(sweep(value, 1L, center, "-"), 1L, scale, "/")
  })
  sources <- lapply(sample_ids, function(sample_id) {
    new_dual_view_source(
      standardized[[sample_id]], sample_id = sample_id, cohort = cohort,
      representation = representation, fit_scope_id = fit_scope_id,
      subsample_seed = seed, standardization_id = standardization_id,
      contract_profile = contract_profile, expected_genes = expected_genes,
      expected_cells = expected_cells, expected_pcs = expected_pcs
    )
  })
  names(sources) <- sample_ids
  list(
    contract_id = "mv03_prepared_sources_v1",
    sources = sources,
    selected_cells = selected_cells,
    center = center,
    scale = scale,
    standardization_id = standardization_id
  )
}
