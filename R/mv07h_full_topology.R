# MV7-H full-124 descriptive dual-view topology and landscape contracts.

.mv07h_seeds <- 20260805:20260809
.mv07h_views <- c("cell_topology_v1", "gene_topology_v1")
.mv07h_dimensions <- c("H0", "H1")

mv07h_ordered_axis_identical_v1 <- function(left, right) {
  identical(unname(left), unname(right))
}

mv07h_canonical_sample_axis_v1 <- function(sample_ids) {
  unname(sort(sample_ids, method = "radix"))
}
.mv07h_forbidden_fields <- c(
  "tissue", "approach", "endpoint", "outcome", "class", "label",
  "cluster", "ari", "nmi"
)

.mv07h_sha256 <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

mv07h_validate_cache_manifest_v1 <- function(manifest) {
  required <- c(
    "sample_id", "seed", "source_tier", "private_cache_file",
    "private_cache_sha256", "normalization_cache_key", "payload_sha256",
    "outcome_label_state", "biological_outcomes_computed"
  )
  if (!is.data.frame(manifest) || nrow(manifest) != 620L ||
      !all(required %in% names(manifest)) ||
      any(tolower(names(manifest)) %in% .mv07h_forbidden_fields) ||
      length(unique(manifest$sample_id)) != 124L ||
      !identical(sort(unique(as.integer(manifest$seed))), .mv07h_seeds) ||
      any(table(manifest$seed) != 124L) ||
      anyDuplicated(paste(manifest$sample_id, manifest$seed, sep = "\r")) ||
      any(manifest$outcome_label_state != "closed") ||
      any(as.logical(manifest$biological_outcomes_computed))) {
    stop("MV7-H cache manifest differs from the frozen 124 by five axis.",
         call. = FALSE)
  }
  invisible(TRUE)
}

mv07h_sample_seed_axis_v1 <- function(manifest) {
  mv07h_validate_cache_manifest_v1(manifest)
  axis <- manifest[order(manifest$seed, manifest$sample_id, method = "radix"),
                   c("sample_id", "seed", "source_tier"), drop = FALSE]
  rownames(axis) <- NULL
  axis$contract_id <- "mv07h_sample_seed_axis_v1"
  axis$selected_cells <- 384L
  axis$panel_genes <- 500L
  axis$outcome_label_state <- "closed"
  axis$biological_outcomes_computed <- FALSE
  axis[c("contract_id", "sample_id", "seed", "source_tier",
         "selected_cells", "panel_genes", "outcome_label_state",
         "biological_outcomes_computed")]
}

mv07h_resource_caps_v1 <- function() {
  data.frame(
    contract_id = "mv07h_resource_caps_v1",
    stage = c("source_views", "cell_ph", "gene_ph", "landscape_group",
              "aggregate"),
    elapsed_cap_seconds = c(1800, 600, 1800, 3600, 172800),
    rss_cap_bytes = c(8, 4, 8, 12, 12) * 1024^3,
    storage_cap_bytes = c(NA, NA, NA, NA, 4 * 1024^3),
    workers = 1L, retries = 0L, atomic_write = TRUE,
    stringsAsFactors = FALSE
  )
}

mv07h_ph_fallback_policy_v1 <- function() {
  data.frame(
    contract_id = "mv07h_exact_ph_resource_fallback_v1",
    primary_engine = "ripserr",
    eligible_view_id = "gene_topology_v1",
    eligible_primary_disposition = "rss_cap_exceeded",
    fallback_engine = "TDA_ripsDiag_GUDHI",
    mathematical_estimand =
      "complete_vietoris_rips_H0_H1_field_2_threshold_minus_1",
    capped_essential_h0_normalization = TRUE,
    fallback_elapsed_cap_seconds = 1800,
    fallback_rss_cap_bytes = 12 * 1024^3,
    fallback_workers = 1L,
    fallback_attempts = 1L,
    fallback_repeat_required = TRUE,
    full_interval_equivalence_required_where_both_complete = TRUE,
    outcome_label_state = "closed",
    biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
}

mv07h_validate_ph_fallback_policy_v1 <- function(policy) {
  expected <- mv07h_ph_fallback_policy_v1()
  if (!is.data.frame(policy) || nrow(policy) != 1L ||
      !identical(names(policy), names(expected)) ||
      !isTRUE(all.equal(policy, expected, check.attributes = FALSE)) ||
      any(tolower(names(policy)) %in% .mv07h_forbidden_fields)) {
    stop("MV7-H PH fallback policy differs from its exact resource contract.",
         call. = FALSE)
  }
  invisible(TRUE)
}

mv07h_ph_fallback_eligible_v1 <- function(stage, view_id, disposition,
                                          policy) {
  mv07h_validate_ph_fallback_policy_v1(policy)
  identical(as.character(stage), "gene_ph") &&
    identical(as.character(view_id), policy$eligible_view_id) &&
    identical(as.character(disposition),
              policy$eligible_primary_disposition)
}

mv07h_source_queue_v1 <- function(axis) {
  if (!is.data.frame(axis) || nrow(axis) != 620L ||
      any(table(axis$seed) != 124L)) {
    stop("MV7-H source queue requires the complete sample-seed axis.",
         call. = FALSE)
  }
  caps <- mv07h_resource_caps_v1()
  cap <- caps[caps$stage == "source_views", , drop = FALSE]
  seeds <- sort(unique(as.integer(axis$seed)))
  data.frame(
    contract_id = "mv07h_source_queue_v1", job_order = seq_along(seeds),
    job_id = paste0("source__", seeds), stage = "source_views", seed = seeds,
    sample_count = 124L, typed_view_count = 248L,
    output_file = paste0("source/mv07h__", seeds, "__source.rds"),
    elapsed_cap_seconds = cap$elapsed_cap_seconds,
    rss_cap_bytes = cap$rss_cap_bytes, workers = 1L, retries = 0L,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
}

mv07h_ph_queue_v1 <- function(axis) {
  if (!is.data.frame(axis) || nrow(axis) != 620L ||
      any(table(axis$seed) != 124L)) {
    stop("MV7-H PH queue requires the complete sample-seed axis.",
         call. = FALSE)
  }
  queue <- merge(
    axis[c("sample_id", "seed")],
    data.frame(view_id = .mv07h_views, stringsAsFactors = FALSE),
    all = TRUE, sort = FALSE
  )
  queue$stage <- ifelse(queue$view_id == "cell_topology_v1",
                        "cell_ph", "gene_ph")
  queue <- queue[order(queue$seed, queue$sample_id, queue$view_id,
                       method = "radix"), , drop = FALSE]
  rownames(queue) <- NULL
  caps <- mv07h_resource_caps_v1()
  matched <- match(queue$stage, caps$stage)
  queue$contract_id <- "mv07h_ph_queue_v1"
  queue$job_order <- seq_len(nrow(queue))
  queue$job_id <- paste("ph", queue$seed, queue$sample_id, queue$view_id,
                        sep = "__")
  queue$output_file <- paste0(
    "ph/mv07h__", queue$seed, "__", queue$sample_id, "__", queue$view_id,
    "__ph.rds"
  )
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

mv07h_landscape_queue_v1 <- function(sentinel_metrics) {
  required <- c("seed", "view_id", "finite_h0_intervals",
                "finite_h1_intervals")
  if (!is.data.frame(sentinel_metrics) || nrow(sentinel_metrics) != 60L ||
      !all(required %in% names(sentinel_metrics)) ||
      any(table(sentinel_metrics$seed) != 12L) ||
      !setequal(unique(sentinel_metrics$view_id), .mv07h_views)) {
    stop("MV7-H landscape ordering requires the accepted MV7-G metrics.",
         call. = FALSE)
  }
  groups <- expand.grid(
    seed = .mv07h_seeds, view_id = .mv07h_views,
    homology_dimension = .mv07h_dimensions,
    stringsAsFactors = FALSE
  )
  groups$sentinel_interval_sum <- vapply(seq_len(nrow(groups)), function(i) {
    row <- groups[i, , drop = FALSE]
    part <- sentinel_metrics[
      sentinel_metrics$seed == row$seed &
        sentinel_metrics$view_id == row$view_id, , drop = FALSE
    ]
    field <- if (row$homology_dimension == "H0") {
      "finite_h0_intervals"
    } else "finite_h1_intervals"
    sum(as.integer(part[[field]]))
  }, integer(1L))
  groups <- groups[order(-groups$sentinel_interval_sum, groups$seed,
                         groups$view_id, groups$homology_dimension,
                         method = "radix"), , drop = FALSE]
  rownames(groups) <- NULL
  identity <- lapply(seq_len(nrow(groups)), function(i) list(
    contract_id = "mv07h_landscape_group_identity_v1",
    seed = as.integer(groups$seed[[i]]),
    view_id = groups$view_id[[i]],
    homology_dimension = groups$homology_dimension[[i]],
    samples = 124L, unordered_pairs = choose(124L, 2L),
    pair_scope = "within_seed_all_124_unordered"
  ))
  groups$contract_id <- "mv07h_landscape_queue_v1"
  groups$group_order <- seq_len(nrow(groups))
  groups$group_id <- paste0("mv07h_landscape_group_v1:", vapply(
    identity, digest::digest, character(1L), algo = "sha256", serialize = TRUE
  ))
  groups$stage <- ifelse(groups$group_order == 1L,
                         "stage_1_stress", "stage_2")
  groups$samples <- 124L
  groups$unordered_pairs <- as.integer(choose(124L, 2L))
  groups$component_rows <- groups$unordered_pairs
  groups$output_directory <- paste0(
    "landscape/", gsub(":", "_", groups$group_id, fixed = TRUE)
  )
  cap <- mv07h_resource_caps_v1()
  cap <- cap[cap$stage == "landscape_group", , drop = FALSE]
  groups$elapsed_cap_seconds <- cap$elapsed_cap_seconds
  groups$rss_cap_bytes <- cap$rss_cap_bytes
  groups$workers <- 1L
  groups$retries <- 0L
  groups$outcome_label_state <- "closed"
  groups$biological_outcomes_computed <- FALSE
  groups$clustering_jobs <- 0L
  groups$label_jobs <- 0L
  groups$outcome_jobs <- 0L
  groups[c("contract_id", "group_order", "group_id", "stage", "seed",
           "view_id", "homology_dimension", "samples", "unordered_pairs",
           "component_rows", "sentinel_interval_sum", "output_directory",
           "elapsed_cap_seconds", "rss_cap_bytes", "workers", "retries",
           "outcome_label_state", "biological_outcomes_computed",
           "clustering_jobs", "label_jobs", "outcome_jobs")]
}

mv07h_pair_id_v1 <- function(seed, first_sample_id, second_sample_id, view_id,
                             homology_dimension) {
  ids <- sort(c(as.character(first_sample_id), as.character(second_sample_id)),
              method = "radix")
  identity <- list(
    contract_id = "mv07h_landscape_pair_identity_v1", seed = as.integer(seed),
    first_sample_id = ids[[1L]], second_sample_id = ids[[2L]],
    view_id = match.arg(view_id, .mv07h_views),
    homology_dimension = match.arg(homology_dimension, .mv07h_dimensions),
    pair_scope = "within_seed_all_124_unordered"
  )
  paste0("mv07h_pair_v1:", digest::digest(
    identity, algo = "sha256", serialize = TRUE
  ))
}

.mv07h_source_payload <- function(record) {
  list(identity = record$identity, panel = record$panel,
       center = record$center, scale = record$scale,
       pca_model = record$pca_model, views = record$views)
}

mv07h_new_source_record_v1 <- function(parent, views) {
  mv07g_validate_source_record_v1(parent)
  sample_ids <- mv07h_canonical_sample_axis_v1(
    parent$pca_model$fit_sample_ids)
  if (!is.list(views) || !identical(names(views), sample_ids)) {
    stop("MV7-H requires all 124 ordered typed-view pairs.", call. = FALSE)
  }
  identity <- list(
    contract_id = "mv07h_source_identity_v1", seed = parent$identity$seed,
    fit_scope_id = parent$identity$fit_scope_id, fit_samples = 124L,
    fit_cells = 47616L, panel_size = 500L,
    panel_sha256 = parent$identity$panel_sha256,
    parent_mv07g_source_cache_key = parent$cache_key,
    input_cache_keys = parent$identity$input_cache_keys,
    selected_cell_sha256 = parent$identity$selected_cell_sha256,
    sample_ids = sample_ids,
    standardization_id = parent$identity$standardization_id,
    pca_model_cache_key = parent$pca_model$cache_key,
    view_cache_keys = lapply(views, function(pair) {
      vapply(pair, `[[`, character(1L), "cache_key")
    }), outcome_label_state = "closed", biological_outcomes_computed = FALSE
  )
  record <- list(
    contract_id = "mv07h_full124_typed_views_v1", identity = identity,
    panel = parent$panel, center = parent$center, scale = parent$scale,
    pca_model = parent$pca_model, views = views, payload_sha256 = NULL,
    cache_key = NULL, downstream_execution = list(
      ph_jobs = 0L, landscape_jobs = 0L, clustering_jobs = 0L,
      label_jobs = 0L, biological_outcome_jobs = 0L)
  )
  record$payload_sha256 <- digest::digest(
    .mv07h_source_payload(record), algo = "sha256", serialize = TRUE)
  record$cache_key <- paste0("mv07h_full124_typed_views_v1:",
                             record$payload_sha256)
  class(record) <- c("scph_mv07h_source_record_v1", "list")
  mv07h_validate_source_record_v1(record)
  record
}

mv07h_validate_source_record_v1 <- function(record) {
  required <- c("contract_id", "identity", "panel", "center", "scale",
                "pca_model", "views", "payload_sha256", "cache_key",
                "downstream_execution")
  if (!is.list(record) || !all(required %in% names(record)) ||
      !inherits(record, "scph_mv07h_source_record_v1") ||
      record$contract_id != "mv07h_full124_typed_views_v1" ||
      record$identity$fit_samples != 124L ||
      record$identity$fit_cells != 47616L ||
      record$identity$panel_size != 500L ||
      !(record$identity$seed %in% .mv07h_seeds) ||
      length(record$identity$input_cache_keys) != 124L ||
      length(record$identity$selected_cell_sha256) != 124L ||
      length(record$identity$sample_ids) != 124L ||
      !identical(names(record$views), record$identity$sample_ids) ||
      nrow(record$panel) != 500L || length(record$center) != 500L ||
      length(record$scale) != 500L || any(!is.finite(record$center)) ||
      any(!is.finite(record$scale)) ||
      any(record$scale <= sqrt(.Machine$double.eps)) ||
      record$identity$outcome_label_state != "closed" ||
      isTRUE(record$identity$biological_outcomes_computed)) {
    stop("MV7-H source record is incomplete or stale.", call. = FALSE)
  }
  .validate_cell_pca_model(record$pca_model)
  if (record$pca_model$n_components != 30L ||
      length(record$pca_model$fit_sample_ids) != 124L ||
      record$pca_model$cache_key != record$identity$pca_model_cache_key ||
      record$pca_model$standardization_id !=
        record$identity$standardization_id) {
    stop("MV7-H source record has a stale PCA identity.", call. = FALSE)
  }
  for (sample_id in names(record$views)) {
    pair <- record$views[[sample_id]]
    if (!identical(sort(names(pair)), sort(.mv07h_views))) {
      stop("MV7-H source record lacks one typed-view pair.", call. = FALSE)
    }
    lapply(pair, validate_topology_view)
    if (any(vapply(pair, `[[`, character(1L), "sample_id") != sample_id)) {
      stop("MV7-H typed-view identity is stale.", call. = FALSE)
    }
  }
  expected <- digest::digest(.mv07h_source_payload(record), algo = "sha256",
                             serialize = TRUE)
  if (record$payload_sha256 != expected ||
      record$cache_key != paste0("mv07h_full124_typed_views_v1:", expected) ||
      any(unlist(record$downstream_execution, use.names = FALSE) != 0)) {
    stop("MV7-H source payload or stop boundary is stale.", call. = FALSE)
  }
  invisible(TRUE)
}

.mv07h_ph_payload <- function(record) {
  list(identity = record$identity, topology_result = record$topology_result,
       h0_mst_oracle = record$h0_mst_oracle)
}

mv07h_new_ph_record_v1 <- function(source_record, sample_id, view_id, result) {
  mv07h_validate_source_record_v1(source_record)
  view <- source_record$views[[sample_id]][[view_id]]
  oracle <- mv07g_validate_ph_against_view_v1(result, view)
  identity <- list(
    contract_id = "mv07h_ph_identity_v1",
    source_cache_key = source_record$cache_key, sample_id = sample_id,
    seed = source_record$identity$seed,
    fit_scope_id = source_record$identity$fit_scope_id,
    panel_sha256 = source_record$identity$panel_sha256, view_id = view_id,
    view_cache_key = view$cache_key, point_axis_role = view$point_axis_role,
    coordinate_axis_role = view$coordinate_axis_role,
    point_metric = view$point_metric, point_count = length(view$point_ids),
    max_dim = 1L, threshold = -1, field = 2L,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE
  )
  record <- list(
    contract_id = "mv07h_ph_record_v1", identity = identity,
    topology_result = result, h0_mst_oracle = oracle,
    payload_sha256 = NULL, cache_key = NULL,
    downstream_execution = list(landscape_jobs = 0L, clustering_jobs = 0L,
      label_jobs = 0L, biological_outcome_jobs = 0L)
  )
  record$payload_sha256 <- digest::digest(
    .mv07h_ph_payload(record), algo = "sha256", serialize = TRUE)
  record$cache_key <- paste0("mv07h_ph_record_v1:", record$payload_sha256)
  class(record) <- c("scph_mv07h_ph_record_v1", "list")
  mv07h_validate_ph_record_v1(record, view)
  record
}

mv07h_run_topology_view_ph_gudhi_v1 <- function(
    view, max_dim = 1L, threshold = -1, field = 2L) {
  validate_topology_view(view)
  if (view$view_id != "gene_topology_v1" || !identical(max_dim, 1L) ||
      !identical(as.numeric(threshold), -1) || !identical(field, 2L)) {
    stop("MV7-H GUDHI fallback is restricted to the exact gene H0/H1 contract.",
         call. = FALSE)
  }
  if (!requireNamespace("TDA", quietly = TRUE)) {
    stop("TDA is required for the MV7-H GUDHI fallback.", call. = FALSE)
  }
  distances <- as.matrix(view$payload)
  maximum_scale <- max(distances)
  raw <- TDA::ripsDiag(
    X = distances, maxdimension = max_dim, maxscale = maximum_scale,
    dist = "arbitrary", library = "GUDHI", location = FALSE,
    printProgress = FALSE
  )$diagram
  diagram <- .as_diagram_matrix(raw)
  h0 <- which(diagram[, "dimension"] == 0)
  if (length(h0) != length(view$point_ids)) {
    stop("MV7-H GUDHI did not expose one capped essential H0 interval.",
         call. = FALSE)
  }
  capped <- h0[[which.max(diagram[h0, "death"])]]
  capped_death <- diagram[capped, "death"]
  diagram <- diagram[-capped, , drop = FALSE]
  diagram <- rbind(diagram, c(dimension = 0, birth = 0, death = Inf))
  invalid_intervals <- !is.finite(diagram[, "dimension"]) |
    !is.finite(diagram[, "birth"]) | is.na(diagram[, "death"]) |
    diagram[, "death"] < diagram[, "birth"]
  zero_persistence <- is.finite(diagram[, "death"]) &
    diagram[, "death"] == diagram[, "birth"]
  provenance <- list(
    result_contract_id = "corrected_topology_result_v1",
    view_id = view$view_id, view_cache_key = view$cache_key,
    contract_version = view$contract_version,
    contract_profile = view$contract_profile,
    scientific_eligible = view$scientific_eligible,
    sample_id = view$sample_id, point_axis_role = view$point_axis_role,
    coordinate_axis_role = view$coordinate_axis_role,
    point_metric = view$point_metric, point_count = length(view$point_ids),
    max_dim = max_dim, threshold = as.numeric(threshold), field = field,
    ph_engine = "TDA_ripsDiag_GUDHI",
    ph_engine_version = as.character(utils::packageVersion("TDA")),
    engine_maxscale = maximum_scale,
    capped_essential_h0_removed = TRUE,
    capped_essential_h0_death = capped_death,
    essential_h0_added = TRUE, essential_h0_count = 1L,
    finite_interval_count = sum(is.finite(diagram[, "death"])),
    infinite_interval_count = sum(is.infinite(diagram[, "death"])),
    zero_persistence_count = sum(zero_persistence),
    invalid_interval_count = sum(invalid_intervals),
    diagram_sha256 = .scientific_digest(diagram)
  )
  structure(
    list(
      diagram = diagram, provenance = provenance,
      cache_key = paste0("corrected_topology_result_v1:",
                         .scientific_digest(provenance))
    ),
    class = c("scph_topology_result_v1", "scph_topology_result")
  )
}

mv07h_validate_ph_record_v1 <- function(record, view = NULL) {
  if (!is.list(record) || !inherits(record, "scph_mv07h_ph_record_v1") ||
      record$contract_id != "mv07h_ph_record_v1" ||
      !(record$identity$seed %in% .mv07h_seeds) ||
      !(record$identity$view_id %in% .mv07h_views) ||
      record$identity$outcome_label_state != "closed" ||
      isTRUE(record$identity$biological_outcomes_computed) ||
      !isTRUE(record$h0_mst_oracle$passed)) {
    stop("MV7-H PH record is incomplete.", call. = FALSE)
  }
  expected <- digest::digest(.mv07h_ph_payload(record), algo = "sha256",
                             serialize = TRUE)
  if (record$payload_sha256 != expected ||
      record$cache_key != paste0("mv07h_ph_record_v1:", expected) ||
      any(unlist(record$downstream_execution, use.names = FALSE) != 0)) {
    stop("MV7-H PH payload or stop boundary is stale.", call. = FALSE)
  }
  if (!is.null(view)) {
    if (record$identity$view_cache_key != view$cache_key) {
      stop("MV7-H PH record belongs to another typed view.", call. = FALSE)
    }
    mv07g_validate_ph_against_view_v1(record$topology_result, view)
  }
  invisible(TRUE)
}

mv07h_finite_intervals_v1 <- function(record, homology_dimension) {
  mv07h_validate_ph_record_v1(record)
  dimension <- match.arg(homology_dimension, .mv07h_dimensions)
  number <- as.integer(sub("H", "", dimension, fixed = TRUE))
  diagram <- record$topology_result$diagram
  value <- diagram[
    diagram[, "dimension"] == number & is.finite(diagram[, "death"]) &
      diagram[, "death"] > diagram[, "birth"],
    c("birth", "death"), drop = FALSE
  ]
  storage.mode(value) <- "double"
  value[order(value[, "birth"], -value[, "death"], method = "radix"),,
        drop = FALSE]
}
