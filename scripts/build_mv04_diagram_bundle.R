#!/usr/bin/env Rscript

options(warn = 2)

if (!requireNamespace("digest", quietly = TRUE)) {
  stop("digest is required for MV-04 bundle identities.", call. = FALSE)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop(
    "Usage: build_mv04_diagram_bundle.R <stage-a-metrics> ",
    "<stage-b-metrics> <stage-c-metrics> <manifest-csv> ",
    "<intervals-csv> <dedup-audit-csv>", call. = FALSE
  )
}

metric_tables <- lapply(args[1:3], function(path) {
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
})
columns <- unique(unlist(lapply(metric_tables, names)))
metric_tables <- lapply(metric_tables, function(value) {
  for (column in setdiff(columns, names(value))) value[[column]] <- NA
  value[columns]
})
metrics <- do.call(rbind, metric_tables)
metrics <- metrics[
  metrics$seed == 20260805L & metrics$disposition == "completed", ,
  drop = FALSE
]
identity_columns <- c(
  "cohort", "representation", "sample_id", "seed", "view_id"
)
identity <- do.call(paste, c(metrics[identity_columns], sep = "::"))
groups <- split(seq_len(nrow(metrics)), identity)

dedup_rows <- list()
manifest_rows <- list()
interval_rows <- list()
for (key in sort(names(groups), method = "radix")) {
  indices <- groups[[key]]
  candidates <- metrics[indices, , drop = FALSE]
  if (length(unique(candidates$diagram_sha256)) != 1L ||
      length(unique(candidates$view_cache_key)) != 1L) {
    stop("Duplicate MV-03 identity has conflicting hashes: ", key,
         call. = FALSE)
  }
  priority <- match(candidates$stage, c("B", "C", "A"))
  selected <- candidates[order(priority)[[1L]], , drop = FALSE]
  result <- readRDS(selected$result_file)
  if (!isTRUE(result$provenance$scientific_eligible) ||
      !identical(result$provenance$diagram_sha256,
                 selected$diagram_sha256)) {
    stop("MV-03 result failed scientific/hash preflight: ", key,
         call. = FALSE)
  }
  diagram_id <- paste0(
    "mv04_diagram_v1:",
    digest::digest(
      list(
        cohort = selected$cohort,
        representation = selected$representation,
        sample_id = selected$sample_id,
        seed = as.integer(selected$seed),
        view_id = selected$view_id,
        diagram_sha256 = selected$diagram_sha256
      ),
      algo = "sha256", serialize = TRUE
    )
  )
  stratum_id <- paste(
    selected$cohort, selected$representation, selected$view_id, sep = "__"
  )
  finite_counts <- integer(2L)
  for (dimension in 0:1) {
    keep <- result$diagram[, "dimension"] == dimension &
      is.finite(result$diagram[, "birth"]) &
      is.finite(result$diagram[, "death"]) &
      result$diagram[, "birth"] < result$diagram[, "death"]
    values <- result$diagram[keep, c("birth", "death"), drop = FALSE]
    finite_counts[[dimension + 1L]] <- nrow(values)
    if (nrow(values)) {
      interval_rows[[length(interval_rows) + 1L]] <- data.frame(
        diagram_id = diagram_id,
        homology_dimension = dimension,
        interval_order = seq_len(nrow(values)),
        birth = values[, "birth"],
        death = values[, "death"],
        stringsAsFactors = FALSE
      )
    }
  }
  manifest_rows[[length(manifest_rows) + 1L]] <- data.frame(
    bundle_contract_id = "mv04_eligible_diagram_bundle_v1",
    diagram_id = diagram_id,
    stratum_id = stratum_id,
    cohort = selected$cohort,
    representation = selected$representation,
    sample_id = selected$sample_id,
    seed = as.integer(selected$seed),
    view_id = selected$view_id,
    point_count = as.integer(selected$point_count),
    h0_finite_intervals = finite_counts[[1L]],
    h1_finite_intervals = finite_counts[[2L]],
    diagram_sha256 = selected$diagram_sha256,
    result_file_sha256 = selected$result_file_sha256,
    source_stage = selected$stage,
    result_file = selected$result_file,
    stringsAsFactors = FALSE
  )
  dedup_rows[[length(dedup_rows) + 1L]] <- data.frame(
    identity = key,
    source_rows = nrow(candidates),
    source_stages = paste(sort(unique(candidates$stage)), collapse = ";"),
    unique_diagram_hashes = length(unique(candidates$diagram_sha256)),
    unique_view_cache_keys = length(unique(candidates$view_cache_key)),
    unique_result_hashes = length(unique(candidates$result_file_sha256)),
    selected_stage = selected$stage,
    selected_result_file = selected$result_file,
    stringsAsFactors = FALSE
  )
}

manifest <- do.call(rbind, manifest_rows)
expected <- c(
  large__sct_whole__cell_topology_v1 = 10L,
  large__sct_whole__gene_topology_v1 = 10L,
  large__seurat_integration__cell_topology_v1 = 10L,
  large__seurat_integration__gene_topology_v1 = 10L,
  bone__sct_whole__cell_topology_v1 = 4L,
  bone__sct_whole__gene_topology_v1 = 4L,
  bone__integrated__cell_topology_v1 = 4L,
  bone__integrated__gene_topology_v1 = 4L
)
observed <- table(manifest$stratum_id)
if (!identical(as.integer(observed[names(expected)]), as.integer(expected))) {
  stop("Deduplicated MV-04 stratum sizes do not match the frozen pilot.",
       call. = FALSE)
}
if (anyDuplicated(manifest$diagram_id) || nrow(manifest) != 56L) {
  stop("MV-04 manifest must contain exactly 56 unique diagrams.",
       call. = FALSE)
}

dir.create(dirname(args[[4L]]), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(args[[5L]]), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(args[[6L]]), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(manifest, args[[4L]], row.names = FALSE)
utils::write.csv(do.call(rbind, interval_rows), args[[5L]], row.names = FALSE)
utils::write.csv(do.call(rbind, dedup_rows), args[[6L]], row.names = FALSE)
message("Built immutable MV-04 input bundle with ", nrow(manifest),
        " diagrams across ", length(unique(manifest$stratum_id)),
        " strata.")
