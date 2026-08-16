#!/usr/bin/env Rscript

options(warn = 2)
for (package in c("digest", "Matrix")) {
  if (!requireNamespace(package, quietly = TRUE)) stop(package, " required.")
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop(paste(
    "usage: validate_mv07fp_seed_stability.R PREFREEZE PRIMARY_CACHE",
    "ADDED_CACHE EVIDENCE DETAIL_OUTPUT CHECK_OUTPUT"
  ))
}
prefreeze <- args[[1]]
primary_cache <- args[[2]]
added_cache <- args[[3]]
evidence <- args[[4]]
detail_output <- args[[5]]
check_output <- args[[6]]
source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/mv03_pilot.R")
source("R/mv05_resource_safe_execution.R")
sha <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}
manifest <- read.csv(
  file.path(prefreeze, "mv07fp-cache-manifest.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
paths <- ifelse(
  manifest$source_tier == "primary90",
  file.path(primary_cache, manifest$private_cache_file),
  file.path(added_cache, manifest$private_cache_file)
)
base <- NULL
presence <- NULL
variances <- NULL
cache_ok <- logical(nrow(manifest))
for (i in seq_len(nrow(manifest))) {
  record <- readRDS(paths[[i]])
  mv05d0_validate_normalization_cache_record_v2(record)
  value <- mv05d0_sct_matrix_from_cache_v1(record)
  cache_ok[[i]] <-
    sha(paths[[i]]) == manifest$private_cache_sha256[[i]] &&
    record$identity$sample_id == manifest$sample_id[[i]] &&
    record$identity$seed == manifest$seed[[i]] &&
    record$payload_sha256 == manifest$payload_sha256[[i]] &&
    ncol(value) == 384L && all(is.finite(value@x))
  if (!cache_ok[[i]]) stop("MV7-FP seed cache failure at row ", i)
  features <- rownames(value)
  if (is.null(base)) {
    base <- features
    presence <- integer(length(base))
    variances <- matrix(NA_real_, length(base), nrow(manifest))
  }
  matched <- match(base, features)
  present <- which(!is.na(matched))
  presence[present] <- presence[present] + 1L
  variances[present, i] <- .mv03_row_variance(
    value[matched[present], , drop = FALSE]
  )
  rm(record, value)
  if (i %% 25L == 0L) invisible(gc())
}
common <- sort(base[presence == nrow(manifest)], method = "radix")
variances <- variances[match(common, base), , drop = FALSE]
positive <- is.finite(variances) & variances > .Machine$double.eps
ranks <- apply(variances, 2L, function(value) {
  rank(-value, ties.method = "min", na.last = "keep")
})
global_rank <- apply(ranks, 1L, stats::median, na.rm = TRUE)
canonical <- canonical_mv03_gene_ids(common)
category <- mv03_feature_category(common)
global <- data.frame(
  feature_id = common, gene = canonical, category = category,
  positive = rowSums(positive) == ncol(positive),
  median_rank = global_rank, stringsAsFactors = FALSE
)
global <- global[global$category == "retained_candidate" & global$positive,
                 , drop = FALSE]
global <- global[order(global$median_rank, global$gene, global$feature_id,
                       method = "radix"), , drop = FALSE]
global <- global[!duplicated(global$gene), , drop = FALSE]
panel <- read.csv(file.path(evidence, "mv07fp-panel.csv"),
                  stringsAsFactors = FALSE, check.names = FALSE)
if (!identical(panel$feature_id, head(global$feature_id, 500L))) {
  stop("MV7-FP supplemental global panel mismatch.")
}
expected <- do.call(rbind, lapply(20260805:20260809, function(seed) {
  columns <- manifest$seed == seed
  seed_rank <- apply(ranks[, columns, drop = FALSE], 1L,
                     stats::median, na.rm = TRUE)
  candidate <- data.frame(
    feature_id = common, gene = canonical, category = category,
    positive = rowSums(positive[, columns, drop = FALSE]) == sum(columns),
    median_rank = seed_rank, stringsAsFactors = FALSE
  )
  candidate <- candidate[
    candidate$category == "retained_candidate" & candidate$positive,
    , drop = FALSE
  ]
  candidate <- candidate[order(candidate$median_rank, candidate$gene,
                               candidate$feature_id, method = "radix"),
                         , drop = FALSE]
  candidate <- candidate[!duplicated(candidate$gene), , drop = FALSE]
  seed_top <- head(candidate$feature_id, 500L)
  overlap <- length(intersect(panel$feature_id, seed_top))
  union_n <- length(union(panel$feature_id, seed_top))
  shared <- match(global$feature_id, candidate$feature_id)
  correlation <- suppressWarnings(stats::cor(
    global$median_rank, candidate$median_rank[shared],
    method = "spearman", use = "complete.obs"
  ))
  data.frame(
    contract_id = "mv07fp_seed_stability_independent_v1",
    seed = seed, cache_count = sum(columns),
    eligible_unique_genes = nrow(candidate),
    seed_top_panel_size = length(seed_top),
    overlap_with_global_panel = overlap,
    jaccard_with_global_panel = overlap / union_n,
    spearman_global_vs_seed_candidate_ranks = correlation,
    stringsAsFactors = FALSE
  )
}))
observed <- read.csv(file.path(evidence, "mv07fp-seed-stability.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)
fields <- setdiff(names(expected), "contract_id")
exact_fields <- fields[!fields %in% c(
  "jaccard_with_global_panel", "spearman_global_vs_seed_candidate_ranks"
)]
exact_ok <- all(vapply(exact_fields, function(field) {
  identical(expected[[field]], observed[[field]])
}, logical(1)))
numeric_ok <- isTRUE(all.equal(
  expected$jaccard_with_global_panel, observed$jaccard_with_global_panel,
  tolerance = 1e-14, check.attributes = FALSE
)) && isTRUE(all.equal(
  expected$spearman_global_vs_seed_candidate_ranks,
  observed$spearman_global_vs_seed_candidate_ranks,
  tolerance = 1e-14, check.attributes = FALSE
))
checks <- data.frame(
  contract_id = "mv07fp_seed_stability_validation_v1",
  category = c("cache_content", "global_panel", "seed_axes",
               "seed_exact_metrics", "seed_numeric_metrics"),
  passed = c(all(cache_ok), TRUE, identical(expected$seed, 20260805:20260809),
             exact_ok, numeric_ok),
  detail = c("620 payloads", "exact ordered 500", "five by 124",
             "counts and overlaps", "Jaccard and Spearman at 1e-14"),
  stringsAsFactors = FALSE
)
if (!all(checks$passed)) {
  stop("MV7-FP supplemental seed-stability validation failed: ",
       paste(checks$category[!checks$passed], collapse = ", "))
}
write.csv(expected, detail_output, row.names = FALSE, na = "")
write.csv(checks, check_output, row.names = FALSE, na = "")
message("MV7-FP supplemental seed-stability validation: 5/5 pass")
