# Internal MV5-B label-closed inductive Seurat mapping helpers.

.mv05_mapping_identity <- function(result) {
  list(
    contract_id = result$contract_id,
    engine_id = result$engine_id,
    fold_id = result$fold_id,
    fit_scope_id = result$fit_scope_id,
    reference_sample_ids = result$reference_sample_ids,
    held_out_sample_id = result$held_out_sample_id,
    features = result$features,
    dimensions = result$dimensions,
    seed = result$seed,
    outcome_label_state = result$outcome_label_state,
    biological_outcomes_computed = result$biological_outcomes_computed,
    reference_identity_sha256 = result$reference_identity_sha256,
    anchor_count = result$anchor_count,
    query_embedding_sha256 = result$query_embedding_sha256
  )
}

.mv05_reference_mapping_identity <- function(reference, features, dimensions) {
  list(
    cells = colnames(reference),
    features = features,
    pca_embeddings = Seurat::Embeddings(reference, "pca")[, dimensions, drop = FALSE],
    pca_loadings = Seurat::Loadings(reference, "pca")[features, dimensions,
                                                            drop = FALSE],
    sct_model_count = length(reference[["SCT"]]@SCTModel.list)
  )
}

.mv05_new_mapping_result <- function(fold_id, fit_scope_id,
                                     reference_sample_ids,
                                     held_out_sample_id, features, dimensions,
                                     seed, reference_identity_sha256,
                                     anchor_count, query_embeddings) {
  query_embeddings <- .validate_named_numeric_matrix(
    query_embeddings, "query_embeddings"
  )
  result <- structure(
    list(
      contract_id = "mv05_inductive_mapping_v1",
      engine_id = "seurat_5_transfer_pcaproject_integrate_embeddings_v1",
      fold_id = .one_nonempty_string(fold_id, "fold_id"),
      fit_scope_id = .one_nonempty_string(fit_scope_id, "fit_scope_id"),
      reference_sample_ids = sort(unique(as.character(reference_sample_ids)),
                                  method = "radix"),
      held_out_sample_id = .one_nonempty_string(
        held_out_sample_id, "held_out_sample_id"
      ),
      features = as.character(features),
      dimensions = as.integer(dimensions),
      seed = .one_integer(seed, "seed", 0L),
      reference_identity_sha256 = .one_nonempty_string(
        reference_identity_sha256, "reference_identity_sha256"
      ),
      anchor_count = .one_integer(anchor_count, "anchor_count", 1L),
      query_embeddings = query_embeddings,
      query_embedding_sha256 = .mv05_execution_digest(query_embeddings),
      outcome_label_state = "closed",
      biological_outcomes_computed = FALSE,
      cache_key = NULL
    ),
    class = "scph_mv05_inductive_mapping_v1"
  )
  if (!length(result$reference_sample_ids) ||
      anyNA(result$reference_sample_ids) ||
      any(!nzchar(result$reference_sample_ids)) ||
      result$held_out_sample_id %in% result$reference_sample_ids ||
      anyNA(result$features) || any(!nzchar(result$features)) ||
      anyDuplicated(result$features) || !length(result$dimensions) ||
      anyNA(result$dimensions) || any(result$dimensions < 1L) ||
      anyDuplicated(result$dimensions)) {
    stop("Inductive mapping identifiers or dimensions are invalid.", call. = FALSE)
  }
  result$cache_key <- paste0(
    "mv05_inductive_mapping_v1:",
    .mv05_execution_digest(.mv05_mapping_identity(result))
  )
  mv05_validate_inductive_mapping_v1(result)
}

mv05_validate_inductive_mapping_v1 <- function(result) {
  if (!inherits(result, "scph_mv05_inductive_mapping_v1") ||
      !is.list(result)) {
    stop("result must be an scph_mv05_inductive_mapping_v1 object.",
         call. = FALSE)
  }
  embeddings <- .validate_named_numeric_matrix(
    result$query_embeddings, "result$query_embeddings"
  )
  if (!identical(colnames(embeddings), paste0("PC", result$dimensions)) ||
      !identical(.mv05_execution_digest(embeddings),
                 result$query_embedding_sha256) ||
      !identical(result$outcome_label_state, "closed") ||
      !identical(result$biological_outcomes_computed, FALSE)) {
    stop("Inductive mapping payload or label boundary is stale.", call. = FALSE)
  }
  expected <- paste0(
    "mv05_inductive_mapping_v1:",
    .mv05_execution_digest(.mv05_mapping_identity(result))
  )
  if (!identical(result$cache_key, expected)) {
    stop("Inductive mapping provenance or cache identity is stale.",
         call. = FALSE)
  }
  invisible(result)
}

mv05_run_inductive_mapping_v1 <- function(reference, query, features,
                                           dimensions = 1:10,
                                           fold_id, fit_scope_id,
                                           reference_sample_ids,
                                           held_out_sample_id,
                                           seed = 20260805L,
                                           k_anchor = 3L,
                                           k_score = 10L,
                                           k_weight = 20L,
                                           verbose = FALSE) {
  if (!inherits(reference, "Seurat") || !inherits(query, "Seurat")) {
    stop("reference and query must be Seurat objects.", call. = FALSE)
  }
  if (!"SCT" %in% names(reference@assays) || !"SCT" %in% names(query@assays) ||
      !"pca" %in% names(reference@reductions)) {
    stop("Mapping requires reference/query SCT assays and a reference PCA.",
         call. = FALSE)
  }
  if (any(colnames(reference) %in% colnames(query))) {
    stop("Reference and held-out query cells must be disjoint.", call. = FALSE)
  }
  features <- as.character(features)
  dimensions <- as.integer(dimensions)
  if (!length(features) || anyNA(features) || any(!nzchar(features)) ||
      anyDuplicated(features) || !length(dimensions) || anyNA(dimensions) ||
      any(dimensions < 1L) || anyDuplicated(dimensions)) {
    stop("features and dimensions must be non-empty and unique.", call. = FALSE)
  }
  reference_features <- rownames(Seurat::GetAssayData(
    reference, assay = "SCT", layer = "data"
  ))
  query_features <- rownames(Seurat::GetAssayData(
    query, assay = "SCT", layer = "data"
  ))
  pca_loadings <- rownames(Seurat::Loadings(reference, "pca"))
  if (!all(features %in% reference_features) ||
      !all(features %in% query_features) || !all(features %in% pca_loadings) ||
      max(dimensions) > ncol(Seurat::Embeddings(reference, "pca"))) {
    stop("Mapping features or dimensions are absent from the frozen reference.",
         call. = FALSE)
  }
  seed <- .one_integer(seed, "seed", 0L)
  k_anchor <- .one_integer(k_anchor, "k_anchor", 1L)
  k_score <- .one_integer(k_score, "k_score", 1L)
  k_weight <- .one_integer(k_weight, "k_weight", 1L)
  before <- .mv05_execution_digest(.mv05_reference_mapping_identity(
    reference, features, dimensions
  ))
  anchors <- tryCatch(
    .with_preserved_seed(seed, Seurat::FindTransferAnchors(
      reference = reference,
      query = query,
      normalization.method = "SCT",
      recompute.residuals = TRUE,
      reference.assay = "SCT",
      query.assay = "SCT",
      reduction = "pcaproject",
      reference.reduction = "pca",
      project.query = FALSE,
      features = features,
      scale = TRUE,
      npcs = max(dimensions),
      l2.norm = TRUE,
      dims = dimensions,
      k.anchor = k_anchor,
      k.filter = NA,
      k.score = k_score,
      nn.method = "rann",
      eps = 0,
      approx.pca = FALSE,
      verbose = verbose
    )),
    error = function(error) stop(
      "find_transfer_anchors_failed: ", conditionMessage(error), call. = FALSE
    )
  )
  anchor_count <- nrow(methods::slot(anchors, "anchors"))
  if (!is.finite(anchor_count) || anchor_count < 1L) {
    stop("No transfer anchors were found.", call. = FALSE)
  }
  effective_k_weight <- min(k_weight, anchor_count)
  mapped <- tryCatch(
    .with_preserved_seed(seed, Seurat::IntegrateEmbeddings(
      anchorset = anchors,
      reference = reference,
      query = query,
      new.reduction.name = "ref.pca",
      k.weight = effective_k_weight,
      reuse.weights.matrix = FALSE,
      verbose = verbose
    )),
    error = function(error) stop(
      "integrate_embeddings_failed: ", conditionMessage(error), call. = FALSE
    )
  )
  after <- .mv05_execution_digest(.mv05_reference_mapping_identity(
    reference, features, dimensions
  ))
  if (!identical(before, after)) {
    stop("Reference identity changed during held-out mapping.", call. = FALSE)
  }
  embeddings <- Seurat::Embeddings(mapped, "ref.pca")
  embeddings <- embeddings[, seq_along(dimensions), drop = FALSE]
  colnames(embeddings) <- paste0("PC", dimensions)
  .mv05_new_mapping_result(
    fold_id = fold_id, fit_scope_id = fit_scope_id,
    reference_sample_ids = reference_sample_ids,
    held_out_sample_id = held_out_sample_id, features = features,
    dimensions = dimensions, seed = seed,
    reference_identity_sha256 = before, anchor_count = anchor_count,
    query_embeddings = embeddings
  )
}

mv05_try_inductive_mapping_v1 <- function(...) {
  started <- proc.time()[["elapsed"]]
  result <- tryCatch(mv05_run_inductive_mapping_v1(...), error = identity)
  elapsed <- proc.time()[["elapsed"]] - started
  if (inherits(result, "error")) {
    return(list(
      status = "held_out_mapping_unavailable",
      reason = conditionMessage(result),
      elapsed_seconds = elapsed,
      result = NULL
    ))
  }
  list(
    status = "completed", reason = "", elapsed_seconds = elapsed,
    result = result
  )
}
