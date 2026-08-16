# MV8-G exact common-475 reference sensitivity contracts.

.mv08g_seeds <- 20260805:20260809
.mv08g_views <- c("cell_topology_v1", "gene_topology_v1")
.mv08g_dimensions <- c("H0", "H1")
.mv08g_components <- c("cell_H0", "cell_H1", "gene_H0", "gene_H1")
.mv08g_forbidden_fields <- c(
  "tissue", "approach", "study", "endpoint", "outcome", "class", "label"
)

mv08g_validate_common475_panel_v1 <- function(panel) {
  required <- c("panel_order_475", "panel_order_500", "feature_id", "gene",
                "parent_panel_sha256", "common475_axis_sha256",
                "selection_rule", "hca_expression_used_for_selection",
                "biological_outcomes_used_for_selection")
  truth <- function(value) tolower(as.character(value)) %in% c("true", "t", "1")
  if (!is.data.frame(panel) || nrow(panel) != 475L ||
      !all(required %in% names(panel)) ||
      !identical(as.integer(panel$panel_order_475), seq_len(475L)) ||
      any(diff(as.integer(panel$panel_order_500)) <= 0L) ||
      anyDuplicated(panel$feature_id) || anyDuplicated(panel$panel_order_500) ||
      length(unique(panel$parent_panel_sha256)) != 1L ||
      length(unique(panel$common475_axis_sha256)) != 1L ||
      unique(panel$common475_axis_sha256) !=
        "b7b802ca862a63d7a4bbcaeab5af1192577663992a5ebde831371b6efafbc0ba" ||
      any(truth(panel$hca_expression_used_for_selection)) ||
      any(truth(panel$biological_outcomes_used_for_selection)) ||
      any(tolower(names(panel)) %in% .mv08g_forbidden_fields)) {
    stop("MV8-G common-475 panel differs from the frozen exact-ID axis.",
         call. = FALSE)
  }
  invisible(TRUE)
}

mv08g_resource_caps_v1 <- function() {
  data.frame(
    contract_id = "mv08g_resource_caps_v1",
    stage = c("source_views", "cell_ph", "gene_ph_primary",
              "gene_ph_fallback", "landscape_group", "comparison", "aggregate"),
    elapsed_cap_seconds = c(1800, 600, 1800, 1800, 3600, 3600, 172800),
    rss_cap_bytes = c(8, 4, 8, 12, 12, 8, 12) * 1024^3,
    storage_cap_bytes = c(NA, NA, NA, NA, NA, NA, 4 * 1024^3),
    workers = 1L, retries = 0L, atomic_write = TRUE,
    stringsAsFactors = FALSE
  )
}

mv08g_sample_seed_axis_v1 <- function(cache_manifest) {
  mv07h_validate_cache_manifest_v1(cache_manifest)
  axis <- cache_manifest[
    order(cache_manifest$seed, cache_manifest$sample_id, method = "radix"),
    c("sample_id", "seed", "source_tier"), drop = FALSE]
  rownames(axis) <- NULL
  axis$contract_id <- "mv08g_sample_seed_axis_v1"
  axis$selected_cells <- 384L
  axis$panel_genes <- 475L
  axis$outcome_label_state <- "closed"
  axis$biological_outcomes_computed <- FALSE
  axis[c("contract_id", "sample_id", "seed", "source_tier",
         "selected_cells", "panel_genes", "outcome_label_state",
         "biological_outcomes_computed")]
}

mv08g_source_queue_v1 <- function(axis) {
  if (!is.data.frame(axis) || nrow(axis) != 620L ||
      any(table(axis$seed) != 124L) || any(axis$panel_genes != 475L)) {
    stop("MV8-G source queue requires 124 samples by five seeds.", call. = FALSE)
  }
  cap <- mv08g_resource_caps_v1()
  cap <- cap[cap$stage == "source_views", , drop = FALSE]
  seeds <- sort(unique(as.integer(axis$seed)))
  data.frame(
    contract_id = "mv08g_source_queue_v1", job_order = seq_along(seeds),
    job_id = paste0("source__", seeds), stage = "source_views", seed = seeds,
    sample_count = 124L, typed_view_count = 248L,
    output_file = paste0("source/mv08g__", seeds, "__source.rds"),
    elapsed_cap_seconds = cap$elapsed_cap_seconds,
    rss_cap_bytes = cap$rss_cap_bytes, workers = 1L, retries = 0L,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
}

mv08g_ph_queue_v1 <- function(axis) {
  if (!is.data.frame(axis) || nrow(axis) != 620L ||
      any(table(axis$seed) != 124L)) {
    stop("MV8-G PH queue requires 124 samples by five seeds.", call. = FALSE)
  }
  queue <- merge(axis[c("sample_id", "seed")],
                 data.frame(view_id = .mv08g_views, stringsAsFactors = FALSE),
                 all = TRUE, sort = FALSE)
  queue$stage <- ifelse(queue$view_id == "cell_topology_v1",
                        "cell_ph", "gene_ph_primary")
  queue <- queue[order(queue$seed, queue$sample_id, queue$view_id,
                       method = "radix"), , drop = FALSE]
  rownames(queue) <- NULL
  caps <- mv08g_resource_caps_v1()
  matched <- match(queue$stage, caps$stage)
  queue$contract_id <- "mv08g_ph_queue_v1"
  queue$job_order <- seq_len(nrow(queue))
  queue$job_id <- paste("ph", queue$seed, queue$sample_id, queue$view_id,
                        sep = "__")
  queue$output_file <- paste0("ph/mv08g__", queue$seed, "__",
    queue$sample_id, "__", queue$view_id, "__ph.rds")
  queue$elapsed_cap_seconds <- caps$elapsed_cap_seconds[matched]
  queue$rss_cap_bytes <- caps$rss_cap_bytes[matched]
  queue$workers <- 1L
  queue$retries <- 0L
  queue$outcome_label_state <- "closed"
  queue$biological_outcomes_computed <- FALSE
  queue[c("contract_id", "job_order", "job_id", "stage", "seed",
          "sample_id", "view_id", "output_file", "elapsed_cap_seconds",
          "rss_cap_bytes", "workers", "retries", "outcome_label_state",
          "biological_outcomes_computed")]
}

mv08g_landscape_queue_v1 <- function() {
  groups <- expand.grid(
    seed = .mv08g_seeds, view_id = .mv08g_views,
    homology_dimension = .mv08g_dimensions, stringsAsFactors = FALSE)
  groups <- groups[order(groups$seed, groups$view_id,
                         groups$homology_dimension, method = "radix"), ]
  rownames(groups) <- NULL
  groups$contract_id <- "mv08g_landscape_queue_v1"
  groups$group_order <- seq_len(nrow(groups))
  groups$group_id <- paste("mv08g", groups$seed, groups$view_id,
                           groups$homology_dimension, sep = "__")
  groups$samples <- 124L
  groups$unordered_pairs <- as.integer(choose(124L, 2L))
  groups$component_rows <- groups$unordered_pairs
  groups$output_directory <- paste0("landscape/", groups$group_id)
  cap <- mv08g_resource_caps_v1()
  cap <- cap[cap$stage == "landscape_group", , drop = FALSE]
  groups$elapsed_cap_seconds <- cap$elapsed_cap_seconds
  groups$rss_cap_bytes <- cap$rss_cap_bytes
  groups$workers <- 1L
  groups$retries <- 0L
  groups$outcome_label_state <- "closed"
  groups$biological_outcomes_computed <- FALSE
  groups[c("contract_id", "group_order", "group_id", "seed", "view_id",
           "homology_dimension", "samples", "unordered_pairs",
           "component_rows", "output_directory", "elapsed_cap_seconds",
           "rss_cap_bytes", "workers", "retries", "outcome_label_state",
           "biological_outcomes_computed")]
}

mv08g_matched_shift_queue_v1 <- function() {
  groups <- mv08g_landscape_queue_v1()
  groups$contract_id <- "mv08g_matched_shift_queue_v1"
  groups$group_id <- sub("^mv08g__", "mv08g_shift__", groups$group_id)
  groups$component_rows <- 124L
  groups$unordered_pairs <- NULL
  groups$output_directory <- paste0("matched-shift/", groups$group_id)
  groups
}

.mv08g_source_payload_v1 <- function(record) {
  list(identity = record$identity, panel = record$panel,
       center = record$center, scale = record$scale,
       pca_model = record$pca_model, views = record$views)
}

mv08g_new_source_record_v1 <- function(identity, panel, prepared, pca_model,
                                       views) {
  record <- list(
    contract_id = "mv08g_common475_source_record_v1", identity = identity,
    panel = panel, center = prepared$center, scale = prepared$scale,
    pca_model = pca_model, views = views, payload_sha256 = NULL,
    cache_key = NULL, downstream_execution = list(
      ph_jobs = 0L, landscape_jobs = 0L, comparison_jobs = 0L,
      clustering_jobs = 0L, label_jobs = 0L, biological_outcome_jobs = 0L))
  record$payload_sha256 <- digest::digest(
    .mv08g_source_payload_v1(record), algo = "sha256", serialize = TRUE)
  record$cache_key <- paste0("mv08g_common475_source_record_v1:",
                             record$payload_sha256)
  class(record) <- c("scph_mv08g_source_record_v1", "list")
  mv08g_validate_source_record_v1(record)
  record
}

mv08g_validate_source_record_v1 <- function(record) {
  required <- c("contract_id", "identity", "panel", "center", "scale",
                "pca_model", "views", "payload_sha256", "cache_key",
                "downstream_execution")
  if (!is.list(record) || !all(required %in% names(record)) ||
      !inherits(record, "scph_mv08g_source_record_v1") ||
      record$contract_id != "mv08g_common475_source_record_v1" ||
      record$identity$fit_samples != 124L ||
      record$identity$fit_cells != 47616L ||
      record$identity$panel_size != 475L ||
      !(record$identity$seed %in% .mv08g_seeds) ||
      record$identity$common475_axis_sha256 !=
        "b7b802ca862a63d7a4bbcaeab5af1192577663992a5ebde831371b6efafbc0ba" ||
      length(record$identity$input_cache_keys) != 124L ||
      length(record$identity$selected_cell_sha256) != 124L ||
      length(record$identity$sample_ids) != 124L ||
      !identical(names(record$views), record$identity$sample_ids) ||
      nrow(record$panel) != 475L || length(record$center) != 475L ||
      length(record$scale) != 475L || any(!is.finite(record$center)) ||
      any(!is.finite(record$scale)) ||
      any(record$scale <= sqrt(.Machine$double.eps)) ||
      record$identity$outcome_label_state != "closed" ||
      isTRUE(record$identity$biological_outcomes_computed)) {
    stop("MV8-G common-475 source record is incomplete or stale.",
         call. = FALSE)
  }
  .validate_cell_pca_model(record$pca_model)
  if (record$pca_model$n_components != 30L ||
      length(record$pca_model$fit_sample_ids) != 124L ||
      record$pca_model$cache_key != record$identity$pca_model_cache_key ||
      record$pca_model$standardization_id != record$identity$standardization_id) {
    stop("MV8-G PCA identity is stale.", call. = FALSE)
  }
  for (sample_id in names(record$views)) {
    pair <- record$views[[sample_id]]
    if (!identical(sort(names(pair)), sort(.mv08g_views))) {
      stop("MV8-G source lacks one typed-view pair.", call. = FALSE)
    }
    lapply(pair, validate_topology_view)
    if (any(vapply(pair, `[[`, character(1L), "sample_id") != sample_id)) {
      stop("MV8-G typed-view identity is stale.", call. = FALSE)
    }
  }
  expected <- digest::digest(.mv08g_source_payload_v1(record), algo = "sha256",
                             serialize = TRUE)
  if (record$payload_sha256 != expected ||
      record$cache_key != paste0("mv08g_common475_source_record_v1:", expected) ||
      any(unlist(record$downstream_execution, use.names = FALSE) != 0)) {
    stop("MV8-G source payload or stop boundary is stale.", call. = FALSE)
  }
  invisible(TRUE)
}

.mv08g_ph_payload_v1 <- function(record) {
  list(identity = record$identity, topology_result = record$topology_result,
       h0_mst_oracle = record$h0_mst_oracle)
}

mv08g_new_ph_record_v1 <- function(source_record, sample_id, view_id, result) {
  mv08g_validate_source_record_v1(source_record)
  view <- source_record$views[[sample_id]][[view_id]]
  oracle <- mv07g_validate_ph_against_view_v1(result, view)
  identity <- list(
    contract_id = "mv08g_ph_identity_v1",
    source_cache_key = source_record$cache_key, sample_id = sample_id,
    seed = source_record$identity$seed,
    fit_scope_id = source_record$identity$fit_scope_id,
    common475_axis_sha256 = source_record$identity$common475_axis_sha256,
    view_id = view_id, view_cache_key = view$cache_key,
    point_axis_role = view$point_axis_role,
    coordinate_axis_role = view$coordinate_axis_role,
    point_metric = view$point_metric, point_count = length(view$point_ids),
    max_dim = 1L, threshold = -1, field = 2L,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE)
  record <- list(
    contract_id = "mv08g_ph_record_v1", identity = identity,
    topology_result = result, h0_mst_oracle = oracle,
    payload_sha256 = NULL, cache_key = NULL,
    downstream_execution = list(landscape_jobs = 0L, comparison_jobs = 0L,
      clustering_jobs = 0L, label_jobs = 0L, biological_outcome_jobs = 0L))
  record$payload_sha256 <- digest::digest(
    .mv08g_ph_payload_v1(record), algo = "sha256", serialize = TRUE)
  record$cache_key <- paste0("mv08g_ph_record_v1:", record$payload_sha256)
  class(record) <- c("scph_mv08g_ph_record_v1", "list")
  mv08g_validate_ph_record_v1(record, view)
  record
}

mv08g_validate_ph_record_v1 <- function(record, view = NULL) {
  if (!is.list(record) || !inherits(record, "scph_mv08g_ph_record_v1") ||
      record$contract_id != "mv08g_ph_record_v1" ||
      !(record$identity$seed %in% .mv08g_seeds) ||
      !(record$identity$view_id %in% .mv08g_views) ||
      record$identity$outcome_label_state != "closed" ||
      isTRUE(record$identity$biological_outcomes_computed) ||
      !isTRUE(record$h0_mst_oracle$passed)) {
    stop("MV8-G PH record is incomplete.", call. = FALSE)
  }
  expected <- digest::digest(.mv08g_ph_payload_v1(record), algo = "sha256",
                             serialize = TRUE)
  if (record$payload_sha256 != expected ||
      record$cache_key != paste0("mv08g_ph_record_v1:", expected) ||
      any(unlist(record$downstream_execution, use.names = FALSE) != 0)) {
    stop("MV8-G PH payload or stop boundary is stale.", call. = FALSE)
  }
  if (!is.null(view)) {
    if (record$identity$view_cache_key != view$cache_key) {
      stop("MV8-G PH record belongs to another typed view.", call. = FALSE)
    }
    mv07g_validate_ph_against_view_v1(record$topology_result, view)
  }
  invisible(TRUE)
}

mv08g_finite_intervals_v1 <- function(record, homology_dimension) {
  mv08g_validate_ph_record_v1(record)
  dimension <- match.arg(homology_dimension, .mv08g_dimensions)
  number <- as.integer(sub("H", "", dimension, fixed = TRUE))
  diagram <- record$topology_result$diagram
  value <- diagram[diagram[, "dimension"] == number &
    is.finite(diagram[, "death"]) & diagram[, "death"] > diagram[, "birth"],
    c("birth", "death"), drop = FALSE]
  storage.mode(value) <- "double"
  value[order(value[, "birth"], -value[, "death"], method = "radix"),,
        drop = FALSE]
}

mv08g_nonnegative_scale_stress_v1 <- function(reference, candidate) {
  reference <- as.numeric(reference)
  candidate <- as.numeric(candidate)
  if (length(reference) < 2L || length(reference) != length(candidate) ||
      any(!is.finite(reference)) || any(!is.finite(candidate)) ||
      any(reference < 0) || any(candidate < 0)) {
    stop("MV8-G stress requires paired finite nonnegative distances.",
         call. = FALSE)
  }
  denominator <- sum(candidate^2)
  scale <- if (denominator > 0) max(0, sum(reference * candidate) / denominator) else 0
  reference_norm <- sqrt(sum(reference^2))
  stress <- if (reference_norm > 0) {
    sqrt(sum((reference - scale * candidate)^2)) / reference_norm
  } else if (all(candidate == 0)) 0 else Inf
  list(scale = scale, normalized_stress = stress)
}

mv08g_top_k_neighbor_overlap_v1 <- function(reference, candidate, k = 10L) {
  validate <- function(value) {
    value <- as.matrix(value)
    if (nrow(value) != ncol(value) || is.null(rownames(value)) ||
        !identical(rownames(value), colnames(value)) ||
        any(!is.finite(value)) || any(value < 0) ||
        max(abs(value - t(value))) > 1e-10 ||
        any(abs(diag(value)) > 1e-10)) {
      stop("MV8-G neighbor overlap requires a named distance matrix.",
           call. = FALSE)
    }
    value
  }
  reference <- validate(reference)
  candidate <- validate(candidate)
  if (!identical(rownames(reference), rownames(candidate))) {
    stop("MV8-G neighbor matrices have different sample axes.", call. = FALSE)
  }
  k <- as.integer(k)
  if (k < 1L || k >= nrow(reference)) stop("MV8-G neighbor k is invalid.")
  ids <- rownames(reference)
  overlap <- vapply(seq_along(ids), function(index) {
    ordered <- function(value) {
      others <- seq_along(ids)[-index]
      others[order(value[index, others], ids[others], method = "radix")][seq_len(k)]
    }
    length(intersect(ordered(reference), ordered(candidate))) / k
  }, numeric(1L))
  data.frame(sample_id = ids, k = k, overlap = overlap,
             stringsAsFactors = FALSE)
}

mv08g_harmonization_class_v1 <- function(component_summary) {
  required <- c("component_id", "median_spearman", "median_top10_overlap",
                "median_fixed_k_pam_ari")
  if (!is.data.frame(component_summary) ||
      !all(required %in% names(component_summary)) ||
      !setequal(component_summary$component_id, .mv08g_components) ||
      nrow(component_summary) != 4L ||
      any(!is.finite(as.matrix(component_summary[required[-1L]])))) {
    stop("MV8-G classification requires four complete component summaries.",
         call. = FALSE)
  }
  misses <- (component_summary$median_spearman < 0.95) +
    (component_summary$median_top10_overlap < 0.80) +
    (component_summary$median_fixed_k_pam_ari < 0.80)
  if (all(misses == 0L)) "high_harmonization_stability" else
    if (any(misses >= 2L)) "material_panel_sensitivity" else
      "mixed_harmonization_stability"
}
