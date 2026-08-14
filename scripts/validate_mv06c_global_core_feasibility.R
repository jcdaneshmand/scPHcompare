#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop(paste(
    "Usage: validate_mv06c_global_core_feasibility.R",
    "<repo-root> <d0-ledger.csv> <d0-cache-dir>",
    "<evidence-dir> <validation-output.csv>"
  ), call. = FALSE)
}

repo <- normalizePath(args[[1L]], mustWork = TRUE)
ledger_path <- normalizePath(args[[2L]], mustWork = TRUE)
cache_dir <- normalizePath(args[[3L]], mustWork = TRUE)
evidence_dir <- normalizePath(args[[4L]], mustWork = TRUE)
validation_path <- args[[5L]]
devtools::load_all(repo, quiet = TRUE)
file_sha <- function(path) digest::digest(
  file = path, algo = "sha256", serialize = FALSE
)
read_evidence <- function(name) utils::read.csv(
  file.path(evidence_dir, name), stringsAsFactors = FALSE, check.names = FALSE
)
same_numeric <- function(first, second, tolerance = 1e-12) {
  length(first) == length(second) &&
    all(abs(as.numeric(first) - as.numeric(second)) <= tolerance)
}

expected_ledger_sha <-
  "73f757a91c202a8e38dfa746fc8816f8e272caf912c33b4b959ef55102c68308"
ledger <- utils::read.csv(
  ledger_path, stringsAsFactors = FALSE, check.names = FALSE
)
ledger <- ledger[order(ledger$sample_id, ledger$seed, method = "radix"),,
                 drop = FALSE]
contract <- read_evidence("mv06c-contract.csv")
source <- read_evidence("mv06c-source-summary.csv")
eligibility <- read_evidence("mv06c-eligibility-summary.csv")
panel <- read_evidence("mv06c-panel.csv")
seed_observed <- read_evidence("mv06c-seed-stability.csv")
workload <- read_evidence("mv06c-workload.csv")
decision <- read_evidence("mv06c-decision.csv")
resource <- read_evidence("mv06c-resource.csv")
manifest <- read_evidence("mv06c-artifact-manifest.csv")

checks <- list()
record_check <- function(category, passed, detail) {
  checks[[length(checks) + 1L]] <<- data.frame(
    contract_id = "mv06c_independent_validation_v1",
    category = category, passed = isTRUE(passed), detail = detail,
    stringsAsFactors = FALSE
  )
}

deterministic_files <- c(
  "mv06c-contract.csv", "mv06c-source-summary.csv",
  "mv06c-eligibility-summary.csv", "mv06c-panel.csv",
  "mv06c-seed-stability.csv", "mv06c-workload.csv", "mv06c-decision.csv"
)
manifest <- manifest[match(deterministic_files, manifest$file), , drop = FALSE]
paths <- file.path(evidence_dir, deterministic_files)
manifest_ok <- nrow(manifest) == length(deterministic_files) &&
  !anyNA(manifest$file) && all(manifest$file == deterministic_files) &&
  all(as.numeric(manifest$bytes) == as.numeric(file.info(paths)$size)) &&
  all(manifest$sha256 == unname(vapply(paths, file_sha, character(1L))))
record_check("artifact_manifest", manifest_ok,
             "seven immutable outputs independently rehashed")

ledger_ok <- identical(file_sha(ledger_path), expected_ledger_sha) &&
  nrow(ledger) == 450L && length(unique(ledger$sample_id)) == 90L &&
  identical(sort(unique(as.integer(ledger$seed))), 20260805:20260809) &&
  all(table(ledger$seed) == 90L) &&
  anyDuplicated(paste(ledger$sample_id, ledger$seed, sep = "\r")) == 0L &&
  all(ledger$outcome_label_state == "closed") &&
  !any(as.logical(ledger$biological_outcomes_computed))
record_check("accepted_ledger", ledger_ok,
             "frozen ledger hash and 90-sample/five-seed axis verified")

first_features <- NULL
feature_union <- character()
variance_matrix <- NULL
presence_count <- NULL
for (index in seq_len(nrow(ledger))) {
  row <- ledger[index, , drop = FALSE]
  path <- file.path(cache_dir, row$private_cache_file)
  if (!file.exists(path) || !identical(file_sha(path), row$private_cache_sha256)) {
    stop("MV6-C independent validator found a missing or stale cache.",
         call. = FALSE)
  }
  record <- readRDS(path)
  mv05d0_validate_normalization_cache_record_v2(record)
  value <- mv05d0_sct_matrix_from_cache_v1(record)
  if (!identical(record$cache_key, row$normalization_cache_key) ||
      ncol(value) != 384L) {
    stop("MV6-C independent cache identity or cell axis mismatch.",
         call. = FALSE)
  }
  features <- rownames(value)
  feature_union <- union(feature_union, features)
  if (is.null(first_features)) {
    first_features <- features
    variance_matrix <- matrix(
      NA_real_, nrow = length(first_features), ncol = nrow(ledger),
      dimnames = list(first_features, NULL)
    )
    presence_count <- integer(length(first_features))
  }
  matched <- match(first_features, features)
  present <- which(!is.na(matched))
  presence_count[present] <- presence_count[present] + 1L
  variance_matrix[present, index] <- .mv03_row_variance(
    value[matched[present], , drop = FALSE]
  )
  rm(record, value)
}
common <- sort(first_features[presence_count == nrow(ledger)], method = "radix")
variance_matrix <- variance_matrix[
  match(common, first_features), , drop = FALSE
]
source_ok <- nrow(source) == 1L && source$records_verified == 450L &&
  source$biological_samples == 90L && source$seeds == 5L &&
  source$cells_per_record == 384L &&
  source$union_features == length(feature_union) &&
  source$common_features == length(common) &&
  source$source_set_sha256 == digest::digest(
    ledger[c("sample_id", "seed", "private_cache_sha256")],
    algo = "sha256", serialize = TRUE
  )
record_check("source_reconstruction", source_ok,
             "450 identities plus union/intersection axes independently rebuilt")

finite <- is.finite(variance_matrix)
positive <- finite & variance_matrix > .Machine$double.eps
finite_all <- rowSums(finite) == ncol(variance_matrix)
positive_all <- rowSums(positive) == ncol(variance_matrix)
ranks <- apply(variance_matrix, 2L, function(value) {
  rank(-value, ties.method = "min", na.last = "keep")
})
median_rank <- apply(ranks, 1L, stats::median, na.rm = TRUE)
minimum_variance <- apply(variance_matrix, 1L, function(value) {
  if (all(is.finite(value))) min(value) else NA_real_
})
canonical <- canonical_mv03_gene_ids(common)
category <- mv03_feature_category(common)
candidate <- data.frame(
  feature_id = common, gene = canonical, category = category,
  finite_cache_count = rowSums(finite),
  positive_cache_count = rowSums(positive),
  finite_all_caches = finite_all, positive_all_caches = positive_all,
  median_variance_rank = median_rank, minimum_variance = minimum_variance,
  stringsAsFactors = FALSE
)
retained <- candidate[
  candidate$category == "retained_candidate" &
    candidate$finite_all_caches & candidate$positive_all_caches,
  , drop = FALSE
]
retained <- retained[order(
  retained$median_variance_rank, retained$gene, retained$feature_id,
  method = "radix"
), , drop = FALSE]
before_dedup <- nrow(retained)
retained <- retained[!duplicated(retained$gene), , drop = FALSE]
expected_eligibility <- c(
  common_features = length(common),
  retained_category_features = sum(category == "retained_candidate"),
  technical_category_features = sum(category != "retained_candidate"),
  nonfinite_any_cache_features = sum(!finite_all),
  nonpositive_any_cache_features = sum(finite_all & !positive_all),
  eligible_features_before_canonical_deduplication = before_dedup,
  duplicate_canonical_features_removed = before_dedup - nrow(retained),
  eligible_unique_canonical_genes = nrow(retained), requested_panel_size = 500L,
  selected_panel_size = min(500L, nrow(retained)),
  eligibility_margin = nrow(retained) - 500L
)
eligibility_ok <- nrow(eligibility) == 1L && all(vapply(
  names(expected_eligibility), function(name) {
    as.numeric(eligibility[[name]]) == as.numeric(expected_eligibility[[name]])
  }, logical(1L)
))
record_check("eligibility_reconstruction", eligibility_ok,
             "all common/category/variance/deduplication counts rebuilt")

expected_panel <- head(retained, 500L)
expected_panel$panel_order <- seq_len(nrow(expected_panel))
expected_panel_sha <- digest::digest(
  expected_panel[c("panel_order", "feature_id", "gene")],
  algo = "sha256", serialize = TRUE
)
panel_ok <- nrow(panel) == 500L &&
  identical(panel$feature_id, expected_panel$feature_id) &&
  identical(panel$gene, expected_panel$gene) &&
  identical(as.integer(panel$panel_order), expected_panel$panel_order) &&
  identical(unique(panel$panel_sha256), expected_panel_sha) &&
  same_numeric(panel$median_variance_rank,
               expected_panel$median_variance_rank) &&
  same_numeric(panel$minimum_variance, expected_panel$minimum_variance) &&
  all(panel$finite_cache_count == 450L) &&
  all(panel$positive_cache_count == 450L)
record_check("ordered_panel_reconstruction", panel_ok,
             "exact ordered 500-gene panel and aggregate values rebuilt")

seed_expected <- lapply(20260805:20260809, function(seed) {
  columns <- as.integer(ledger$seed) == seed
  seed_positive <- rowSums(positive[, columns, drop = FALSE]) == sum(columns)
  seed_rank <- apply(ranks[, columns, drop = FALSE], 1L,
                     stats::median, na.rm = TRUE)
  value <- data.frame(
    feature_id = common, gene = canonical, category = category,
    positive = seed_positive, median_rank = seed_rank,
    stringsAsFactors = FALSE
  )
  value <- value[value$category == "retained_candidate" & value$positive,,
                 drop = FALSE]
  value <- value[order(value$median_rank, value$gene, value$feature_id,
                       method = "radix"), , drop = FALSE]
  value <- value[!duplicated(value$gene), , drop = FALSE]
  top <- head(value$feature_id, 500L)
  overlap <- length(intersect(expected_panel$feature_id, top))
  union_count <- length(union(expected_panel$feature_id, top))
  matched <- match(retained$feature_id, value$feature_id)
  data.frame(
    seed = seed, cache_count = sum(columns), eligible_unique_genes = nrow(value),
    seed_top_panel_size = length(top), overlap_with_global_panel = overlap,
    jaccard_with_global_panel = overlap / union_count,
    spearman_global_vs_seed_candidate_ranks = suppressWarnings(stats::cor(
      retained$median_variance_rank, value$median_rank[matched],
      method = "spearman", use = "complete.obs"
    )), stringsAsFactors = FALSE
  )
})
seed_expected <- do.call(rbind, seed_expected)
seed_ok <- identical(as.integer(seed_observed$seed), seed_expected$seed) &&
  identical(as.integer(seed_observed$eligible_unique_genes),
            seed_expected$eligible_unique_genes) &&
  identical(as.integer(seed_observed$overlap_with_global_panel),
            seed_expected$overlap_with_global_panel) &&
  same_numeric(seed_observed$jaccard_with_global_panel,
               seed_expected$jaccard_with_global_panel) &&
  same_numeric(seed_observed$spearman_global_vs_seed_candidate_ranks,
               seed_expected$spearman_global_vs_seed_candidate_ranks)
record_check("seed_stability_reconstruction", seed_ok,
             "all five seed panels, overlaps, Jaccards, and correlations rebuilt")

workload_ok <- nrow(workload) == 1L && workload$matched_cell_views == 6750L &&
  workload$matched_gene_views == 6750L &&
  workload$h0_h1_diagram_components == 27000L &&
  workload$directed_query_training_pairs == 35350L &&
  workload$four_component_landscape_distances == 141400L &&
  workload$training_fitted_component_scales == 300L &&
  workload$five_weight_fusion_pair_rows == 176750L &&
  !workload$execution_authorized
record_check("workload_inventory", workload_ok,
             "future matched-SCT workload is complete and unauthorized")

decision_ok <- nrow(decision) == 1L &&
  decision$decision == "go_bounded_matched_sct_profile" &&
  decision$panel_sha256 == expected_panel_sha &&
  decision$selected_panel_size == 500L &&
  decision$eligible_unique_canonical_genes == nrow(retained) &&
  decision$eligibility_margin == nrow(retained) - 500L &&
  decision$all_selected_present_finite_nonconstant &&
  decision$resource_caps_pass
record_check("decision_rule", decision_ok,
             "prefrozen go rule follows from panel and resource evidence")

resource_ok <- nrow(resource) == 1L && resource$elapsed_seconds <= 1800 &&
  resource$peak_process_rss_bytes <= 8 * 1024^3 &&
  resource$elapsed_cap_pass && resource$rss_cap_pass
record_check("resource_gate", resource_ok,
             "inventory elapsed and RSS remain within frozen caps")

zero_columns <- c(
  "cell_source_jobs_executed", "gene_source_jobs_executed",
  "pca_jobs_executed", "ph_jobs_executed", "landscape_jobs_executed",
  "fusion_jobs_executed"
)
zero_ok <- all(as.numeric(decision[zero_columns]) == 0) &&
  !decision$biological_outcomes_computed &&
  decision$outcome_label_state == "closed" &&
  !contract$biological_outcomes_computed &&
  contract$outcome_label_state == "closed"
record_check("stop_boundary", zero_ok,
             "source, PCA, PH, landscape, fusion, and outcome jobs remain zero")

public_frames <- list(
  contract, source, eligibility, panel, seed_observed, workload, decision,
  resource, manifest
)
forbidden <- c("tissue", "approach", "endpoint", "cell_id", "expression")
privacy_ok <- !any(vapply(public_frames, function(value) {
  any(vapply(tolower(names(value)), function(name) {
    any(startsWith(name, forbidden))
  }, logical(1L)))
}, logical(1L)))
record_check("public_safe_schema", privacy_ok,
             "no tissue, approach, endpoint, cell, or expression-value columns")

finite_ok <- all(vapply(public_frames, function(value) {
  numeric <- value[vapply(value, is.numeric, logical(1L))]
  !length(numeric) || all(is.finite(as.matrix(numeric)))
}, logical(1L)))
record_check("finite_outputs", finite_ok,
             "all numeric evidence values are finite")

validation <- do.call(rbind, checks)
if (!all(validation$passed)) {
  stop("MV6-C independent validation failed: ",
       paste(validation$category[!validation$passed], collapse = ", "),
       call. = FALSE)
}
dir.create(dirname(validation_path), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(validation, validation_path, row.names = FALSE, na = "",
                 quote = TRUE)
message("MV6-C independent validation passed ", nrow(validation),
        " categories.")
