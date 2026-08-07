#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) {
  stop(
    "usage: run_mv05c_existing_data_job.R CACHE SAMPLE_MANIFEST ",
    "CELL_MANIFEST FOLD_ID SEED PRIVATE_DIR AUDIT_DIR", call. = FALSE
  )
}

cache_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
sample_manifest_path <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
cell_manifest_path <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
fold_id <- args[[4L]]
seed <- as.integer(args[[5L]])
private_dir <- args[[6L]]
audit_dir <- args[[7L]]
dir.create(private_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

source("R/provenance_utils.R")
source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/topological_distance_contract.R")
source("R/mv03_pilot.R")
source("R/mv05_benchmark_execution.R")
source("R/mv05_inductive_mapping.R")

if (!seed %in% 20260805:20260809 ||
    !grepl("^mv05c_loso_v1:SRA[0-9]+$", fold_id)) {
  stop("Fold or seed is outside the frozen MV5-C manifest.", call. = FALSE)
}
held_out_study <- sub("^mv05c_loso_v1:", "", fold_id)
fit_scope_id <- paste0(fold_id, ":training_reference")
job_id <- paste(fold_id, seed, sep = "__")
safe_job_id <- gsub("[^A-Za-z0-9_.-]", "_", job_id)
bundle_path <- file.path(private_dir, paste0(safe_job_id, ".rds"))
status_path <- file.path(
  audit_dir, paste0("mv05c-job-status-", safe_job_id, ".csv")
)
if (file.exists(bundle_path) || file.exists(status_path)) {
  stop("Refusing to overwrite an existing MV5-C job artifact.", call. = FALSE)
}

sample_manifest <- utils::read.csv(
  sample_manifest_path, stringsAsFactors = FALSE, check.names = FALSE
)
cell_manifest <- utils::read.csv(
  cell_manifest_path, stringsAsFactors = FALSE, check.names = FALSE
)
required_sample <- c(
  "sample_id", "study", "filtered_cells", "outcome_label_state",
  "biological_outcomes_computed"
)
required_cell <- c(
  "sample_id", "seed", "cell_id", "selected_cell_sha256",
  "outcome_label_state", "biological_outcomes_computed"
)
if (!all(required_sample %in% names(sample_manifest)) ||
    !all(required_cell %in% names(cell_manifest)) ||
    any(c("tissue", "approach") %in% names(sample_manifest)) ||
    any(c("tissue", "approach") %in% names(cell_manifest)) ||
    any(sample_manifest$outcome_label_state != "closed") ||
    any(cell_manifest$outcome_label_state != "closed") ||
    any(as.logical(sample_manifest$biological_outcomes_computed)) ||
    any(as.logical(cell_manifest$biological_outcomes_computed))) {
  stop("MV5-C manifests violate the frozen label boundary.", call. = FALSE)
}
if (!held_out_study %in% sample_manifest$study) {
  stop("Held-out study is absent from the closed sample manifest.", call. = FALSE)
}

cache_started <- proc.time()[["elapsed"]]
cache <- readRDS(cache_path)
cache_seconds <- proc.time()[["elapsed"]] - cache_started
if (!identical(cache$contract_id, "mv05c_existing_data_raw_cache_v1") ||
    !identical(cache$sample_ids, sample_manifest$sample_id)) {
  stop("Private source cache and closed sample manifest disagree.", call. = FALSE)
}
raw_samples <- cache$samples
source_preflight_path <- file.path(
  dirname(audit_dir), "mv05c-source-preflight-2026-08-06.csv"
)
source_preflight <- utils::read.csv(
  source_preflight_path, stringsAsFactors = FALSE, check.names = FALSE
)
if (nrow(source_preflight) != 1L ||
    source_preflight$private_cache_size_bytes[[1L]] !=
      unname(file.info(cache_path)$size) ||
    !grepl("^[0-9a-f]{64}$",
           source_preflight$private_cache_sha256[[1L]])) {
  stop("Private cache does not match the frozen source preflight.")
}
selected_rows <- cell_manifest[cell_manifest$seed == seed, , drop = FALSE]
selected_cells <- split(selected_rows$cell_id, selected_rows$sample_id)
if (!identical(sort(names(selected_cells)), sort(sample_manifest$sample_id)) ||
    any(lengths(selected_cells) != 384L)) {
  stop("Matched-cell manifest is incomplete for this seed.", call. = FALSE)
}

prefix_cells <- function(sample_id, ids) paste(sample_id, ids, sep = "__")
make_counts <- function(sample_id) {
  value <- raw_samples[[sample_id]]
  ids <- selected_cells[[sample_id]]
  if (!all(ids %in% colnames(value))) {
    stop("Selected cells are absent from raw counts for ", sample_id)
  }
  value <- value[, ids, drop = FALSE]
  colnames(value) <- prefix_cells(sample_id, ids)
  Matrix::Matrix(value, sparse = TRUE)
}
counts <- lapply(sample_manifest$sample_id, make_counts)
names(counts) <- sample_manifest$sample_id
rm(cache, raw_samples)
invisible(gc())

make_sct <- function(sample_id) {
  object <- Seurat::CreateSeuratObject(
    counts = counts[[sample_id]], project = sample_id, min.cells = 0L,
    min.features = 0L
  )
  object$sample_id <- sample_id
  Seurat::SCTransform(
    object, variable.features.n = 3000L, return.only.var.genes = FALSE,
    verbose = FALSE, seed.use = seed
  )
}

select_training_panel <- function(matrices, training_ids, panel_size = 500L) {
  common <- Reduce(intersect, lapply(matrices, rownames))
  common <- sort(common, method = "radix")
  canonical <- canonical_mv03_gene_ids(common)
  category <- mv03_feature_category(common)
  training_variances <- lapply(training_ids, function(sample_id) {
    training_matrix <- as.matrix(
      matrices[[sample_id]][common, , drop = FALSE]
    )
    .mv03_row_variance(training_matrix)
  })
  ranks <- lapply(training_variances, function(variance) {
    rank(-variance, ties.method = "min", na.last = "keep")
  })
  rank_matrix <- do.call(cbind, ranks)
  median_rank <- apply(rank_matrix, 1L, stats::median, na.rm = TRUE)
  variance_matrix <- do.call(cbind, training_variances)
  finite_training <- apply(variance_matrix, 1L, function(x) {
    all(is.finite(x) & x > .Machine$double.eps)
  })
  candidates <- data.frame(
    feature_id = common, gene = canonical, category = category,
    median_training_variance_rank = median_rank,
    finite_training = finite_training, stringsAsFactors = FALSE
  )
  candidates <- candidates[
    candidates$category == "retained_candidate" & candidates$finite_training,
    , drop = FALSE
  ]
  candidates <- candidates[order(
    candidates$median_training_variance_rank, candidates$gene,
    candidates$feature_id, method = "radix"
  ), , drop = FALSE]
  candidates <- candidates[!duplicated(candidates$gene), , drop = FALSE]
  if (nrow(candidates) < panel_size) {
    stop("Fewer than 500 unique training-ranked genes are available.")
  }
  candidates[seq_len(panel_size), , drop = FALSE]
}

prepare_sources <- function(matrices, panel, representation, training_ids) {
  selected <- lapply(matrices, function(value) {
    result <- as.matrix(value[panel$feature_id, , drop = FALSE])
    rownames(result) <- panel$gene
    result
  })
  training_pool <- do.call(cbind, selected[training_ids])
  center <- rowMeans(training_pool)
  scale <- apply(training_pool, 1L, stats::sd)
  if (any(!is.finite(center)) || any(!is.finite(scale)) ||
      any(scale <= sqrt(.Machine$double.eps))) {
    stop("Training-only standardization produced invalid parameters.")
  }
  standardized <- lapply(selected, function(value) {
    sweep(sweep(value, 1L, center, "-"), 1L, scale, "/")
  })
  if (any(vapply(standardized, function(value) any(!is.finite(value)),
                 logical(1L)))) {
    stop("One or more held-out sources contain nonfinite panel values.")
  }
  standardization_identity <- list(
    contract_id = "mv05c_training_standardization_v1",
    fold_id = fold_id, fit_scope_id = fit_scope_id,
    representation = representation, seed = seed,
    training_ids = training_ids, panel = panel,
    center = center, scale = scale
  )
  standardization_id <- paste0(
    "mv05c_training_standardization_v1:",
    digest::digest(standardization_identity, algo = "sha256", serialize = TRUE)
  )
  cell_sources <- lapply(names(standardized), function(sample_id) {
    new_cell_projection_source(
      standardized[[sample_id]], sample_id = sample_id,
      cohort = "mv05c_one_tissue_feasibility_v1",
      representation = representation, fit_scope_id = fit_scope_id,
      subsample_seed = seed, standardization_id = standardization_id
    )
  })
  names(cell_sources) <- names(standardized)
  gene_eligible <- vapply(standardized, function(value) {
    all(apply(value, 1L, stats::sd) > sqrt(.Machine$double.eps))
  }, logical(1L))
  dual_sources <- lapply(names(standardized), function(sample_id) {
    if (!gene_eligible[[sample_id]]) return(NULL)
    new_dual_view_source(
      standardized[[sample_id]], sample_id = sample_id,
      cohort = "mv05c_one_tissue_feasibility_v1",
      representation = representation, fit_scope_id = fit_scope_id,
      subsample_seed = seed, standardization_id = standardization_id
    )
  })
  names(dual_sources) <- names(standardized)
  list(
    sources = cell_sources, dual_sources = dual_sources,
    gene_eligible = gene_eligible, selected = selected, center = center,
    scale = scale, standardization_id = standardization_id
  )
}

run_ph_views <- function(views) {
  lapply(views, function(view) {
    started <- proc.time()[["elapsed"]]
    result <- run_topology_view_ph(view, max_dim = 1L, threshold = -1, field = 2L)
    result$execution <- list(
      elapsed_seconds = proc.time()[["elapsed"]] - started,
      completed_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
    )
    result
  })
}

status_rows <- list()
add_status <- function(representation, view_id, status, reason = "",
                       elapsed_seconds = NA_real_, sample_count = 6L,
                       diagram_count = NA_integer_, h0_intervals = NA_integer_,
                       h1_intervals = NA_integer_) {
  status_rows[[length(status_rows) + 1L]] <<- data.frame(
    contract_id = "mv05c_existing_data_job_v1",
    job_id = job_id, fold_id = fold_id, fit_scope_id = fit_scope_id,
    held_out_study = held_out_study, seed = seed,
    representation = representation, view_id = view_id,
    status = status, reason = reason, sample_count = sample_count,
    diagram_count = diagram_count, h0_finite_intervals = h0_intervals,
    h1_finite_intervals = h1_intervals, elapsed_seconds = elapsed_seconds,
    source_cache_read_seconds = cache_seconds,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
}

summarize_diagrams <- function(diagrams, dimension) {
  sum(vapply(diagrams, function(result) {
    sum(result$diagram[, "dimension"] == dimension &
          is.finite(result$diagram[, "birth"]) &
          is.finite(result$diagram[, "death"]) &
          result$diagram[, "birth"] < result$diagram[, "death"])
  }, integer(1L)))
}

training_ids <- sample_manifest$sample_id[sample_manifest$study != held_out_study]
query_ids <- sample_manifest$sample_id[sample_manifest$study == held_out_study]
bundle <- list(
  contract_id = "mv05c_existing_data_job_v1", job_id = job_id,
  fold_id = fold_id, fit_scope_id = fit_scope_id,
  held_out_study = held_out_study, seed = seed,
  training_sample_ids = training_ids, held_out_sample_ids = query_ids,
  source_cache_sha256 = source_preflight$private_cache_sha256[[1L]],
  sample_manifest_sha256 = digest::digest(
    file = sample_manifest_path, algo = "sha256", serialize = FALSE
  ),
  cell_manifest_sha256 = digest::digest(
    file = cell_manifest_path, algo = "sha256", serialize = FALSE
  ),
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  sct_fold = NULL, inductive_integrated = NULL, failures = list()
)

sct_started <- proc.time()[["elapsed"]]
sct_step <- "initialize"
sct_phase <- tryCatch({
  sct_step <- "sctransform_samples"
  sct_objects <- lapply(sample_manifest$sample_id, make_sct)
  names(sct_objects) <- sample_manifest$sample_id
  sct_step <- "extract_sct_data"
  sct_matrices <- lapply(sct_objects, function(object) {
    Seurat::GetAssayData(object, assay = "SCT", layer = "data")
  })
  sct_step <- "select_training_panel"
  panel <- select_training_panel(sct_matrices, training_ids)
  sct_step <- "prepare_training_scaled_sources"
  prepared <- prepare_sources(sct_matrices, panel, "sct_fold", training_ids)
  sct_step <- "fit_training_pca"
  pca_model <- fit_cell_topology_pca(
    prepared$dual_sources[training_ids], n_components = 30L, pca_seed = seed
  )
  sct_step <- "construct_cell_views"
  cell_views <- lapply(prepared$sources, function(source) {
    coordinates <- t(source$matrix) %*% pca_model$rotation
    construct_frozen_cell_topology_view(
      source, coordinates, "mv05c_training_fitted_pca_v1",
      pca_model$cache_key
    )
  })
  gene_views <- NULL
  gene_diagrams <- NULL
  gene_reason <- ""
  if (all(prepared$gene_eligible)) {
    sct_step <- "construct_gene_views"
    gene_views <- lapply(prepared$dual_sources, construct_gene_topology_view)
  } else {
    failed_ids <- names(prepared$gene_eligible)[!prepared$gene_eligible]
    gene_reason <- paste0(
      "training_only_panel_has_constant_genes_in_held_out_samples: ",
      paste(failed_ids, collapse = ";")
    )
  }
  sct_step <- "cell_ph"
  cell_diagrams <- run_ph_views(cell_views)
  if (!is.null(gene_views)) {
    sct_step <- "gene_ph"
    gene_diagrams <- run_ph_views(gene_views)
  }
  sct_step <- "matched_baselines"
  baselines <- list(
    cell_energy = mv05_cell_energy_baseline_v1(cell_views),
    pseudobulk = mv05_pseudobulk_baseline_v1(prepared$sources)
  )
  if (!is.null(gene_views)) {
    baselines$gene_correlation <- mv05_gene_correlation_baseline_v1(
      prepared$dual_sources
    )
  }
  list(
    status = "completed", panel = panel, prepared = prepared,
    pca_model = pca_model, cell_views = cell_views, gene_views = gene_views,
    cell_diagrams = cell_diagrams, gene_diagrams = gene_diagrams,
    gene_reason = gene_reason, baselines = baselines,
    sct_objects = sct_objects
  )
}, error = identity)
sct_elapsed <- proc.time()[["elapsed"]] - sct_started
if (inherits(sct_phase, "error")) {
  reason <- paste0(sct_step, ": ", conditionMessage(sct_phase))
  bundle$failures$sct_fold <- reason
  add_status("sct_fold", "cell_topology_v1", "failed", reason, sct_elapsed)
  add_status("sct_fold", "gene_topology_v1", "failed", reason, sct_elapsed)
  sct_objects <- NULL
} else {
  sct_objects <- sct_phase$sct_objects
  sct_phase$sct_objects <- NULL
  bundle$sct_fold <- sct_phase
  add_status(
    "sct_fold", "cell_topology_v1", "completed", elapsed_seconds = sct_elapsed,
    diagram_count = 6L,
    h0_intervals = summarize_diagrams(sct_phase$cell_diagrams, 0L),
    h1_intervals = summarize_diagrams(sct_phase$cell_diagrams, 1L)
  )
  if (is.null(sct_phase$gene_diagrams)) {
    add_status(
      "sct_fold", "gene_topology_v1",
      "structured_unavailable_held_out_constant_genes",
      sct_phase$gene_reason, sct_elapsed, diagram_count = 0L
    )
  } else {
    add_status(
      "sct_fold", "gene_topology_v1", "completed",
      elapsed_seconds = sct_elapsed, diagram_count = 6L,
      h0_intervals = summarize_diagrams(sct_phase$gene_diagrams, 0L),
      h1_intervals = summarize_diagrams(sct_phase$gene_diagrams, 1L)
    )
  }
}

integrated_started <- proc.time()[["elapsed"]]
integrated_step <- "initialize"
integrated_phase <- tryCatch({
  integrated_step <- "ensure_sample_sct"
  if (is.null(sct_objects)) {
    sct_objects <- lapply(sample_manifest$sample_id, make_sct)
    names(sct_objects) <- sample_manifest$sample_id
  }
  integrated_step <- "align_reference_counts"
  reference_features <- Reduce(
    intersect, lapply(counts[training_ids], rownames)
  )
  reference_counts <- do.call(cbind, lapply(counts[training_ids], function(value) {
    value[reference_features, , drop = FALSE]
  }))
  integrated_step <- "reference_sctransform"
  reference <- Seurat::CreateSeuratObject(
    counts = reference_counts, project = paste0("reference_", safe_job_id),
    min.cells = 0L, min.features = 0L
  )
  reference$sample_id <- sub("__.*$", "", colnames(reference))
  reference <- Seurat::SCTransform(
    reference, variable.features.n = 3000L, return.only.var.genes = FALSE,
    verbose = FALSE, seed.use = seed
  )
  integrated_step <- "select_mapping_features"
  features <- Seurat::VariableFeatures(reference)
  features <- features[mv03_feature_category(features) == "retained_candidate"]
  feature_genes <- canonical_mv03_gene_ids(features)
  features <- features[!duplicated(feature_genes)]
  features <- features[vapply(features, function(feature) {
    all(vapply(sct_objects[query_ids], function(object) {
      feature %in% rownames(Seurat::GetAssayData(
        object, assay = "SCT", layer = "data"
      ))
    }, logical(1L)))
  }, logical(1L))]
  if (length(features) < 500L) {
    stop("Fewer than 500 reference-selected mapping features are query-compatible.")
  }
  features <- features[seq_len(500L)]
  integrated_step <- "reference_pca"
  reference <- Seurat::RunPCA(
    reference, assay = "SCT", features = features, npcs = 30L,
    approx = FALSE, seed.use = seed, verbose = FALSE
  )
  integrated_step <- "held_out_mapping"
  previous_future_max <- getOption("future.globals.maxSize")
  options(future.globals.maxSize = 2 * 1024^3)
  if (requireNamespace("future", quietly = TRUE)) {
    future::plan("sequential")
  }
  mapping <- lapply(query_ids, function(sample_id) {
    mv05_try_inductive_mapping_v1(
      reference = reference, query = sct_objects[[sample_id]],
      features = features, dimensions = 1:30, fold_id = fold_id,
      fit_scope_id = fit_scope_id, reference_sample_ids = training_ids,
      held_out_sample_id = sample_id, seed = seed, k_anchor = 3L,
      k_score = 10L, k_weight = 20L, verbose = FALSE
    )
  })
  options(future.globals.maxSize = previous_future_max)
  names(mapping) <- query_ids
  failed <- names(mapping)[vapply(mapping, `[[`, character(1L), "status") !=
                                   "completed"]
  if (length(failed)) {
    reasons <- vapply(mapping[failed], `[[`, character(1L), "reason")
    stop("Held-out mapping failed: ", paste(failed, reasons, collapse = "; "))
  }
  integrated_step <- "assemble_shared_coordinates"
  reference_embeddings <- Seurat::Embeddings(reference, "pca")[, 1:30,
                                                                   drop = FALSE]
  colnames(reference_embeddings) <- paste0("PC", 1:30)
  coordinates <- lapply(sample_manifest$sample_id, function(sample_id) {
    if (sample_id %in% training_ids) {
      value <- reference_embeddings[
        startsWith(rownames(reference_embeddings), paste0(sample_id, "__")),
        , drop = FALSE
      ]
    } else {
      value <- mapping[[sample_id]]$result$query_embeddings
    }
    expected <- prefix_cells(sample_id, selected_cells[[sample_id]])
    value <- value[expected, , drop = FALSE]
    if (!identical(rownames(value), expected) || nrow(value) != 384L) {
      stop("Integrated coordinates failed ordered-cell validation: ", sample_id)
    }
    value
  })
  names(coordinates) <- sample_manifest$sample_id
  if (is.null(bundle$sct_fold)) {
    stop("Integrated cell coordinates require validated matched-cell identity sources.")
  }
  integrated_step <- "prepare_integrated_identity_sources"
  identity_sources <- lapply(bundle$sct_fold$prepared$sources, function(source) {
    new_cell_projection_source(
      source$matrix, sample_id = source$sample_id, cohort = source$cohort,
      representation = "inductive_integrated", fit_scope_id = fit_scope_id,
      subsample_seed = seed,
      standardization_id = paste0(
        "mv05c_coordinate_identity_carrier_v1:", source$cache_key
      )
    )
  })
  names(identity_sources) <- sample_manifest$sample_id
  prepared <- list(
    contract_id = "mv05c_coordinate_identity_carrier_v1",
    sources = identity_sources,
    mathematical_coordinate_input = FALSE,
    role = "ordered_cell_identity_only"
  )
  reference_fit_key <- paste0(
    "mv05c_reference_pca_v1:", digest::digest(
      list(
        fold_id = fold_id, seed = seed, training_ids = training_ids,
        features = features, embeddings = reference_embeddings,
        loadings = Seurat::Loadings(reference, "pca")[features, 1:30, drop = FALSE]
      ), algo = "sha256", serialize = TRUE
    )
  )
  integrated_step <- "construct_integrated_cell_views"
  cell_views <- lapply(sample_manifest$sample_id, function(sample_id) {
    fit_key <- if (sample_id %in% training_ids) reference_fit_key else
      mapping[[sample_id]]$result$cache_key
    construct_frozen_cell_topology_view(
      identity_sources[[sample_id]], coordinates[[sample_id]],
      "seurat_5_inductive_reference_pca_v1", fit_key
    )
  })
  names(cell_views) <- sample_manifest$sample_id
  integrated_step <- "integrated_cell_ph"
  cell_diagrams <- run_ph_views(cell_views)
  integrated_step <- "integrated_cell_baseline"
  list(
    status = "completed", features = features, prepared = prepared,
    reference_fit_key = reference_fit_key, mapping = mapping,
    coordinates = coordinates, cell_views = cell_views,
    cell_diagrams = cell_diagrams,
    baselines = list(cell_energy = mv05_cell_energy_baseline_v1(cell_views))
  )
}, error = identity)
integrated_elapsed <- proc.time()[["elapsed"]] - integrated_started
if (inherits(integrated_phase, "error")) {
  reason <- paste0(integrated_step, ": ", conditionMessage(integrated_phase))
  bundle$failures$inductive_integrated <- reason
  add_status(
    "inductive_integrated", "cell_topology_v1", "held_out_mapping_unavailable",
    reason, integrated_elapsed
  )
} else {
  bundle$inductive_integrated <- integrated_phase
  add_status(
    "inductive_integrated", "cell_topology_v1", "completed",
    elapsed_seconds = integrated_elapsed, diagram_count = 6L,
    h0_intervals = summarize_diagrams(integrated_phase$cell_diagrams, 0L),
    h1_intervals = summarize_diagrams(integrated_phase$cell_diagrams, 1L)
  )
}
add_status(
  "inductive_integrated", "gene_topology_v1",
  "structured_unavailable_no_inductive_gene_correction",
  "IntegrateEmbeddings supplies shared held-out cell coordinates but no corrected held-out gene-expression matrix.",
  integrated_elapsed, diagram_count = 0L
)

bundle$completed_utc <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
bundle$status <- do.call(rbind, status_rows)
partial <- paste0(bundle_path, ".partial")
saveRDS(bundle, partial, compress = FALSE, version = 3)
if (!file.rename(partial, bundle_path)) {
  stop("Failed to atomically publish the MV5-C job bundle.")
}
bundle_sha256 <- digest::digest(
  file = bundle_path, algo = "sha256", serialize = FALSE
)
status <- bundle$status
status$private_bundle_file <- basename(bundle_path)
status$private_bundle_size_bytes <- unname(file.info(bundle_path)$size)
status$private_bundle_sha256 <- bundle_sha256
utils::write.csv(status, status_path, row.names = FALSE, na = "")
message("MV5-C job completed with explicit representation dispositions: ", job_id)
