#!/usr/bin/env Rscript

Sys.setenv(
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1"
)
options(warn = 1)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 10L) {
  stop(
    "usage: run_mv05f_mapping_group.R PILOT_MANIFEST GROUP_ID D1_CACHE_DIR ",
    "RAW_RESOURCE_CSV RAW_CACHE_DIR SCT_RESOURCE_CSV SCT_CACHE_DIR ",
    "RESULT_DIR GROUP_AUDIT_CSV RUN_MODE", call. = FALSE
  )
}
manifest_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
group_id <- args[[2L]]
d1_cache_dir <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
raw_resource_path <- normalizePath(args[[4L]], winslash = "/", mustWork = TRUE)
raw_cache_dir <- normalizePath(args[[5L]], winslash = "/", mustWork = TRUE)
sct_resource_path <- normalizePath(args[[6L]], winslash = "/", mustWork = TRUE)
sct_cache_dir <- normalizePath(args[[7L]], winslash = "/", mustWork = TRUE)
result_dir <- args[[8L]]
audit_path <- args[[9L]]
run_mode <- match.arg(args[[10L]], c("build_or_resume", "validate_resume"))
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(audit_path), recursive = TRUE, showWarnings = FALSE)

source("R/provenance_utils.R")
source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/mv03_pilot.R")
source("R/mv05_resource_safe_execution.R")
source("R/mv05_benchmark_execution.R")
source("R/mv05_inductive_mapping.R")
source("R/mv05f_integration_gate.R")
if (file.exists("R/mv05g_coordinate_production.R")) {
  source("R/mv05g_coordinate_production.R")
}

if (!requireNamespace("future", quietly = TRUE)) {
  stop("future is required to freeze sequential MV5-F execution.", call. = FALSE)
}
future::plan(future::sequential)
previous_future_max <- getOption("future.globals.maxSize")
options(future.globals.maxSize = 4 * 1024^3)
on.exit(options(future.globals.maxSize = previous_future_max), add = TRUE)

file_sha <- .mv05f_file_sha256
implementation_files <- c(
  "R/provenance_utils.R", "R/toy_baseline.R", "R/dual_view_topology.R",
  "R/mv03_pilot.R", "R/mv05_resource_safe_execution.R",
  "R/mv05_benchmark_execution.R", "R/mv05_inductive_mapping.R",
  "R/mv05f_integration_gate.R", "R/mv05g_coordinate_production.R",
  "scripts/run_mv05f_mapping_group.R"
)
implementation_sha <- .mv05f_sha256(stats::setNames(
  vapply(implementation_files, file_sha, character(1L)), implementation_files
))

manifest <- utils::read.csv(
  manifest_path, stringsAsFactors = FALSE, check.names = FALSE
)
contract <- unique(manifest$contract_id)
if (identical(contract, "mv05f_mapping_pilot_manifest_v1")) {
  mv05f_validate_pilot_manifest_v1(manifest)
} else if (identical(contract, "mv05g_integrated_coordinate_manifest_v1")) {
  mv05g_validate_full_manifest_v1(manifest)
} else {
  stop("Unsupported mapping manifest contract.", call. = FALSE)
}
job <- manifest[manifest$group_id == group_id, , drop = FALSE]
if (nrow(job) != 1L) stop("Requested MV5-F group is absent.", call. = FALSE)

d1_path <- file.path(d1_cache_dir, job$private_cache_file)
if (!file.exists(d1_path) || file_sha(d1_path) != job$private_cache_sha256) {
  stop("Accepted MV5-D1 fold cache is missing or stale.", call. = FALSE)
}
d1_record <- readRDS(d1_path)
mv05d1_validate_cell_fold_record_v1(d1_record)
if (d1_record$cache_key != job$fold_cache_key) {
  stop("MV5-D1 fold-cache identity differs from the pilot.", call. = FALSE)
}
sample_ids <- sort(c(
  d1_record$identity$training_ids, d1_record$identity$query_ids
), method = "radix")
training_ids <- d1_record$identity$training_ids
query_ids <- d1_record$identity$query_ids
panel <- d1_record$payload$panel
if (length(sample_ids) != 90L || nrow(panel) != 500L ||
    length(query_ids) != job$held_out_samples ||
    length(training_ids) != job$training_samples) {
  stop("MV5-F source axes differ from the frozen 90/500 contract.", call. = FALSE)
}

raw_resources <- utils::read.csv(
  raw_resource_path, stringsAsFactors = FALSE, check.names = FALSE
)
sct_resources <- utils::read.csv(
  sct_resource_path, stringsAsFactors = FALSE, check.names = FALSE
)
.mv05f_assert_label_closed(raw_resources, "MV5-D0 raw resources")
.mv05f_assert_label_closed(sct_resources, "MV5-D0 SCT resources")
raw_resources <- raw_resources[match(sample_ids, raw_resources$sample_id),,
                               drop = FALSE]
sct_resources <- sct_resources[
  sct_resources$seed == job$seed & sct_resources$sample_id %in% sample_ids,,
  drop = FALSE
]
sct_resources <- sct_resources[match(sample_ids, sct_resources$sample_id),,
                               drop = FALSE]
if (anyNA(raw_resources$sample_id) || anyNA(sct_resources$sample_id) ||
    !identical(raw_resources$sample_id, sample_ids) ||
    !identical(sct_resources$sample_id, sample_ids) ||
    any(raw_resources$disposition != "built_atomic") ||
    any(sct_resources$disposition != "built_atomic")) {
  stop("Accepted D0 source records do not cover the exact 90 samples.",
       call. = FALSE)
}

input_started <- proc.time()[["elapsed"]]
selected_cells <- vector("list", length(sample_ids))
names(selected_cells) <- sample_ids
sct_source_keys <- stats::setNames(
  sct_resources$normalization_cache_key, sample_ids
)
for (index in seq_along(sample_ids)) {
  path <- file.path(sct_cache_dir, sct_resources$private_cache_file[[index]])
  if (!file.exists(path) || file_sha(path) !=
      sct_resources$private_cache_sha256[[index]]) {
    stop("Accepted D0 SCT cache is missing or stale: ", sample_ids[[index]],
         call. = FALSE)
  }
  record <- readRDS(path)
  mv05d0_validate_normalization_cache_record_v2(record)
  if (record$cache_key != sct_source_keys[[sample_ids[[index]]]]) {
    stop("D0 SCT cache key is stale: ", sample_ids[[index]], call. = FALSE)
  }
  selected_cells[[index]] <- as.character(record$payload$selected_cell_ids)
  rm(record)
}
raw_hashes <- stats::setNames(raw_resources$private_raw_cache_sha256, sample_ids)
counts <- vector("list", length(sample_ids))
names(counts) <- sample_ids
present_features <- vector("list", length(sample_ids))
names(present_features) <- sample_ids
for (index in seq_along(sample_ids)) {
  sample_id <- sample_ids[[index]]
  path <- file.path(raw_cache_dir, raw_resources$private_raw_cache_file[[index]])
  if (!file.exists(path) || file_sha(path) != raw_hashes[[sample_id]]) {
    stop("Accepted D0 raw shard is missing or stale: ", sample_id, call. = FALSE)
  }
  raw <- readRDS(path)
  mv05d0_validate_raw_sample_cache_v2(raw)
  if (!identical(raw$sample_id, sample_id) ||
      !all(selected_cells[[sample_id]] %in% colnames(raw$counts))) {
    stop("Raw shard identity or selected cells are stale: ", sample_id,
         call. = FALSE)
  }
  present <- panel$feature_id %in% rownames(raw$counts)
  present_features[[sample_id]] <- panel$gene[present]
  value <- raw$counts[
    panel$feature_id[present], selected_cells[[sample_id]], drop = FALSE
  ]
  rownames(value) <- panel$gene[present]
  colnames(value) <- paste(sample_id, selected_cells[[sample_id]], sep = "__")
  counts[[sample_id]] <- Matrix::Matrix(value, sparse = TRUE)
  rm(raw, value)
}
input_seconds <- proc.time()[["elapsed"]] - input_started
if (any(vapply(counts[training_ids], nrow, integer(1L)) != 500L) ||
    any(vapply(counts, ncol, integer(1L)) != 384L)) {
  stop("Training panel or frozen selected-cell dimensions are incomplete.",
       call. = FALSE)
}

identity <- mv05f_group_identity_v1(
  job, d1_record, raw_hashes, sct_source_keys, implementation_sha
)
result_file <- paste0(gsub("[^A-Za-z0-9_.-]", "_", group_id), ".rds")
result_path <- file.path(result_dir, result_file)
if (file.exists(result_path)) {
  record <- readRDS(result_path)
  mv05f_validate_group_record_v1(record)
  if (record$cache_key != identity$cache_key) {
    stop("Existing MV5-F result has a stale identity.", call. = FALSE)
  }
  if (run_mode == "build_or_resume") disposition <- "reuse_validated" else
    disposition <- "validated_resume"
  stage <- list(
    input_seconds = input_seconds, reference_sct_pca_seconds = NA_real_,
    query_sct_seconds = NA_real_, mapping_seconds = NA_real_,
    assembly_seconds = NA_real_
  )
} else {
  if (run_mode == "validate_resume") {
    stop("MV5-F resume validation found a missing result.", call. = FALSE)
  }
  reference_started <- proc.time()[["elapsed"]]
  reference_counts <- do.call(cbind, counts[training_ids])
  reference <- Seurat::CreateSeuratObject(
    counts = reference_counts, project = paste0("mv05f_", job$held_out_study),
    min.cells = 0L, min.features = 0L
  )
  reference$sample_id <- sub("__.*$", "", colnames(reference))
  reference <- Seurat::SCTransform(
    reference, variable.features.n = 3000L, return.only.var.genes = FALSE,
    verbose = FALSE, seed.use = as.integer(job$seed)
  )
  reference <- Seurat::RunPCA(
    reference, assay = "SCT", features = panel$gene, npcs = 30L,
    approx = FALSE, seed.use = as.integer(job$seed), verbose = FALSE
  )
  reference_seconds <- proc.time()[["elapsed"]] - reference_started
  reference_before <- .mv05f_sha256(.mv05_reference_mapping_identity(
    reference, panel$gene, 1:30
  ))
  reference_embeddings <- Seurat::Embeddings(reference, "pca")[, 1:30,
                                                                  drop = FALSE]
  colnames(reference_embeddings) <- paste0("PC", 1:30)

  query_sct_seconds <- 0
  mapping_seconds <- 0
  mappings <- vector("list", length(query_ids))
  names(mappings) <- query_ids
  active_features <- vector("list", length(query_ids))
  names(active_features) <- query_ids
  for (sample_id in query_ids) {
    query_started <- proc.time()[["elapsed"]]
    query <- Seurat::CreateSeuratObject(
      counts = counts[[sample_id]], project = sample_id,
      min.cells = 0L, min.features = 0L
    )
    query$sample_id <- sample_id
    query <- Seurat::SCTransform(
      query, variable.features.n = 3000L, return.only.var.genes = FALSE,
      verbose = FALSE, seed.use = as.integer(job$seed)
    )
    query_sct_seconds <- query_sct_seconds +
      (proc.time()[["elapsed"]] - query_started)
    query_features <- rownames(Seurat::GetAssayData(
      query, assay = "SCT", layer = "data"
    ))
    active <- panel$gene[panel$gene %in% query_features]
    if (length(active) < 31L) {
      stop("Too few fixed-panel active query features: ", sample_id,
           call. = FALSE)
    }
    active_features[[sample_id]] <- active
    mapping_started <- proc.time()[["elapsed"]]
    attempt <- mv05_try_inductive_mapping_v1(
      reference = reference, query = query, features = active,
      dimensions = 1:30, fold_id = identity$fold_id,
      fit_scope_id = identity$fit_scope_id,
      reference_sample_ids = training_ids,
      held_out_sample_id = sample_id, seed = as.integer(job$seed),
      k_anchor = 3L, k_score = 10L, k_weight = 20L, verbose = FALSE
    )
    mapping_seconds <- mapping_seconds +
      (proc.time()[["elapsed"]] - mapping_started)
    if (attempt$status != "completed") {
      stop("Held-out mapping failed for ", sample_id, ": ", attempt$reason,
           call. = FALSE)
    }
    mappings[[sample_id]] <- attempt$result
    rm(query, attempt)
    invisible(gc())
  }
  reference_after <- .mv05f_sha256(.mv05_reference_mapping_identity(
    reference, panel$gene, 1:30
  ))
  if (!identical(reference_before, reference_after)) {
    stop("The full fixed-panel reference changed during query mapping.",
         call. = FALSE)
  }
  assembly_started <- proc.time()[["elapsed"]]
  coordinates <- lapply(sample_ids, function(sample_id) {
    value <- if (sample_id %in% training_ids) {
      reference_embeddings[
        startsWith(rownames(reference_embeddings), paste0(sample_id, "__")),
        , drop = FALSE
      ]
    } else {
      mappings[[sample_id]]$query_embeddings
    }
    expected <- paste(sample_id, selected_cells[[sample_id]], sep = "__")
    value <- value[expected, , drop = FALSE]
    if (!identical(rownames(value), expected) || nrow(value) != 384L) {
      stop("Ordered coordinate assembly failed: ", sample_id, call. = FALSE)
    }
    value
  })
  names(coordinates) <- sample_ids
  assembly_seconds <- proc.time()[["elapsed"]] - assembly_started
  record <- mv05f_new_group_record_v1(
    identity, coordinates, mappings, active_features,
    reference_before, reference_after
  )
  stage <- list(
    input_seconds = input_seconds,
    reference_sct_pca_seconds = reference_seconds,
    query_sct_seconds = query_sct_seconds,
    mapping_seconds = mapping_seconds,
    assembly_seconds = assembly_seconds
  )
  mv05f_validate_group_record_v1(record)
  partial <- tempfile(pattern = paste0(result_file, "."), tmpdir = result_dir)
  saveRDS(record, partial, compress = FALSE, version = 3)
  if (file.exists(result_path) || !file.rename(partial, result_path)) {
    unlink(partial)
    stop("Failed to atomically publish the MV5-F result.", call. = FALSE)
  }
  disposition <- "built_atomic"
}

record <- readRDS(result_path)
mv05f_validate_group_record_v1(record)
active_counts <- vapply(record$payload$active_features, length, integer(1L))
anchor_counts <- vapply(
  record$payload$query_mappings, `[[`, integer(1L), "anchor_count"
)
audit <- data.frame(
  contract_id = "mv05f_mapping_group_audit_v1", group_id = group_id,
  group_order = job$group_order,
  execution_role = if ("pilot_role" %in% names(job)) job$pilot_role else
    job$production_role,
  fold_id = identity$fold_id, fit_scope_id = identity$fit_scope_id,
  held_out_study = identity$held_out_study, seed = identity$seed,
  training_samples = length(training_ids), held_out_samples = length(query_ids),
  completed_query_mappings = length(record$payload$query_mappings),
  completed_coordinate_views = length(record$payload$coordinates),
  panel_genes = nrow(panel), minimum_active_query_features = min(active_counts),
  maximum_dropped_query_features = max(500L - active_counts),
  minimum_anchor_count = min(anchor_counts), maximum_anchor_count = max(anchor_counts),
  reference_identity_sha256_before =
    record$payload$reference_identity_sha256_before,
  reference_identity_sha256_after =
    record$payload$reference_identity_sha256_after,
  reference_immutable = identical(
    record$payload$reference_identity_sha256_before,
    record$payload$reference_identity_sha256_after
  ),
  coordinate_set_sha256 = record$payload$coordinate_set_sha256,
  cache_key = record$cache_key, payload_sha256 = record$payload_sha256,
  input_seconds = stage$input_seconds,
  reference_sct_pca_seconds = stage$reference_sct_pca_seconds,
  query_sct_seconds = stage$query_sct_seconds,
  mapping_seconds = stage$mapping_seconds,
  assembly_seconds = stage$assembly_seconds,
  private_result_file = result_file,
  private_result_bytes = unname(file.info(result_path)$size),
  private_result_sha256 = file_sha(result_path), disposition = disposition,
  implementation_sha256 = implementation_sha,
  label_transfer_jobs_executed = 0L, ph_jobs_executed = 0L,
  landscape_jobs_executed = 0L, distance_jobs_executed = 0L,
  clustering_jobs_executed = 0L, gene_view_jobs_executed = 0L,
  fusion_jobs_executed = 0L, new_data_jobs_executed = 0L,
  biological_outcomes_computed = FALSE, outcome_label_state = "closed",
  stringsAsFactors = FALSE
)
write_provenance_csv(audit, audit_path)
message("Completed MV5-F mapping group: ", group_id, " (", disposition, ")")
