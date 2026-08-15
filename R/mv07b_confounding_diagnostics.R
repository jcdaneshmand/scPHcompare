# MV7-B no-new-PH confounding diagnostic helpers.

mv07b_methods_v1 <- function() data.frame(
  method_order = 1:6,
  method_id = c("cell_H0", "cell_H1", "gene_H0", "gene_H1",
    "cell_composite", "gene_composite"),
  view = c("cell", "cell", "gene", "gene", "cell", "gene"),
  dimension = c("H0", "H1", "H0", "H1", "H0_H1_descriptive", "H0_H1_descriptive"),
  stringsAsFactors = FALSE)

mv07b_endpoints_v1 <- function() data.frame(
  endpoint_order = 1:2,
  endpoint_id = c("mean_reciprocal_rank", "one_nn_balanced_accuracy"),
  role = c("primary_diagnostic", "supportive_diagnostic"),
  bounded_0_1 = TRUE, stringsAsFactors = FALSE)

mv07b_contrasts_v1 <- function() data.frame(
  contrast_order = 1:3,
  contrast_id = c("cell_H0_minus_gene_H0", "cell_H1_minus_gene_H1",
    "cell_composite_minus_gene_composite"),
  left_method = c("cell_H0", "cell_H1", "cell_composite"),
  right_method = c("gene_H0", "gene_H1", "gene_composite"),
  selection_role = "fixed_diagnostic_not_method_selection",
  stringsAsFactors = FALSE)

mv07b_macro <- function(x, endpoint) {
  means <- tapply(x[[endpoint]], x$query_tissue, mean)
  if (length(means) != 5L || any(!is.finite(means))) return(NA_real_)
  mean(means)
}

mv07b_within_study_rank_correlation <- function(x, endpoint) {
  parts <- split(x, x$held_out_study)
  centered <- lapply(parts, function(z) {
    if (nrow(z) < 2L || length(unique(z$retained_cells)) < 2L) return(NULL)
    a <- rank(log2(z$retained_cells), ties.method = "average")
    b <- rank(z[[endpoint]], ties.method = "average")
    data.frame(a = a - mean(a), b = b - mean(b))
  })
  centered <- do.call(rbind, centered)
  if (is.null(centered) || nrow(centered) < 3L || sd(centered$a) == 0 ||
      sd(centered$b) == 0) return(NA_real_)
  cor(centered$a, centered$b)
}

mv07b_percentile <- function(x) {
  if (sum(is.finite(x)) < 0.95 * length(x)) return(c(NA_real_, NA_real_))
  as.numeric(quantile(x[is.finite(x)], c(0.025, 0.975), names = FALSE, type = 7))
}
