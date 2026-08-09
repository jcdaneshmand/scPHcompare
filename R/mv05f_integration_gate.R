# Internal MV5-F label-closed integration-induction resource-gate helpers.

mv05f_forbidden_columns_v1 <- function() {
  c("tissue", "approach", "outcome", "label", "cluster", "ari", "mrr")
}

.mv05f_sha256 <- function(value) {
  digest::digest(value, algo = "sha256", serialize = TRUE)
}

.mv05f_file_sha256 <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

.mv05f_assert_label_closed <- function(value, name = "input") {
  forbidden <- mv05f_forbidden_columns_v1()
  if (!is.data.frame(value) ||
      any(tolower(names(value)) %in% forbidden) ||
      !all(c("outcome_label_state", "biological_outcomes_computed") %in%
           names(value)) ||
      any(value$outcome_label_state != "closed") ||
      any(as.logical(value$biological_outcomes_computed))) {
    stop(name, " violates the MV5-F label-closed boundary.", call. = FALSE)
  }
  invisible(value)
}

mv05f_select_pilot_groups_v1 <- function(resources, pilot_seed = 20260805L) {
  required <- c(
    "fold_id", "fit_scope_id", "held_out_study", "seed",
    "training_samples", "held_out_samples",
    "held_out_missing_feature_instances",
    "maximum_missing_features_per_view", "private_cache_file",
    "private_cache_sha256", "fold_cache_key", "outcome_label_state",
    "biological_outcomes_computed"
  )
  if (!is.data.frame(resources) || !all(required %in% names(resources))) {
    stop("MV5-F requires the accepted MV5-D1 resource schema.", call. = FALSE)
  }
  .mv05f_assert_label_closed(resources, "MV5-D1 resources")
  pilot_seed <- as.integer(pilot_seed)
  if (length(pilot_seed) != 1L || is.na(pilot_seed) ||
      !pilot_seed %in% sort(unique(as.integer(resources$seed)))) {
    stop("pilot_seed is outside the accepted MV5-D1 seed axis.", call. = FALSE)
  }
  aggregate_key <- stats::aggregate(
    cbind(
      held_out_missing_feature_instances,
      maximum_missing_features_per_view
    ) ~ held_out_study + held_out_samples + training_samples,
    data = resources,
    FUN = sum
  )
  names(aggregate_key)[names(aggregate_key) ==
    "held_out_missing_feature_instances"] <-
    "aggregate_missing_feature_instances"
  names(aggregate_key)[names(aggregate_key) ==
    "maximum_missing_features_per_view"] <-
    "aggregate_maximum_missing_features"
  base <- resources[resources$seed == pilot_seed, , drop = FALSE]
  base <- merge(
    base, aggregate_key,
    by = c("held_out_study", "held_out_samples", "training_samples"),
    sort = FALSE
  )
  pick_one <- function(rows, order_columns, decreasing = FALSE) {
    order_args <- lapply(order_columns, function(column) rows[[column]])
    if (decreasing) {
      order_args <- lapply(order_args, function(value) {
        if (is.numeric(value)) -value else xtfrm(value) * -1
      })
    }
    rows[do.call(order, c(order_args, list(method = "radix")))[1L],,
         drop = FALSE]
  }
  minimum_query <- pick_one(base, c("held_out_samples", "held_out_study"))
  maximum_query <- base[order(
    -base$held_out_samples, base$held_out_study, method = "radix"
  )[1L], , drop = FALSE]
  remaining <- base[!base$held_out_study %in%
    c(minimum_query$held_out_study, maximum_query$held_out_study),,
    drop = FALSE]
  maximum_missing <- remaining[order(
    -remaining$aggregate_missing_feature_instances,
    -remaining$aggregate_maximum_missing_features,
    remaining$held_out_study, method = "radix"
  )[1L], , drop = FALSE]
  remaining <- remaining[remaining$held_out_study !=
    maximum_missing$held_out_study, , drop = FALSE]
  missing <- remaining[remaining$aggregate_missing_feature_instances > 0L,,
                       drop = FALSE]
  median_query <- stats::median(base$held_out_samples)
  median_missing <- missing[order(
    abs(missing$held_out_samples - median_query),
    missing$held_out_samples, missing$held_out_study, method = "radix"
  )[1L], , drop = FALSE]
  selected <- rbind(minimum_query, maximum_query, maximum_missing, median_missing)
  selected$pilot_role <- c(
    "minimum_query_maximum_reference", "maximum_query_minimum_reference",
    "maximum_aggregate_missing_burden", "median_query_nonzero_missing"
  )
  selected$group_order <- seq_len(nrow(selected))
  selected$group_id <- paste0(
    "mv05f_group__", selected$held_out_study, "__", selected$seed
  )
  selected$contract_id <- "mv05f_mapping_pilot_manifest_v1"
  selected$panel_genes <- 500L
  selected$cells_per_sample <- 384L
  selected$pca_components <- 30L
  selected$group_timeout_seconds <- 1800
  selected$rss_cap_bytes <- 8 * 1024^3
  selected$stage_cap_seconds <- 7200
  selected$storage_cap_bytes <- 10 * 1000^3
  selected$maximum_heavy_workers <- 1L
  selected$ph_jobs_executed <- 0L
  selected$landscape_jobs_executed <- 0L
  selected$distance_jobs_executed <- 0L
  selected$clustering_jobs_executed <- 0L
  selected$gene_view_jobs_executed <- 0L
  selected$fusion_jobs_executed <- 0L
  selected$new_data_jobs_executed <- 0L
  selected$biological_outcomes_computed <- FALSE
  selected$outcome_label_state <- "closed"
  selected <- selected[order(selected$group_order), , drop = FALSE]
  rownames(selected) <- NULL
  mv05f_validate_pilot_manifest_v1(selected)
  selected
}

mv05f_validate_pilot_manifest_v1 <- function(manifest) {
  required <- c(
    "contract_id", "group_id", "group_order", "pilot_role", "fold_id", "fit_scope_id",
    "held_out_study", "seed", "training_samples", "held_out_samples",
    "panel_genes", "cells_per_sample", "pca_components",
    "private_cache_file", "private_cache_sha256", "fold_cache_key",
    "group_timeout_seconds", "rss_cap_bytes", "stage_cap_seconds",
    "storage_cap_bytes", "maximum_heavy_workers", "outcome_label_state",
    "biological_outcomes_computed", "ph_jobs_executed",
    "landscape_jobs_executed", "distance_jobs_executed",
    "clustering_jobs_executed", "gene_view_jobs_executed",
    "fusion_jobs_executed", "new_data_jobs_executed"
  )
  zero <- c(
    "ph_jobs_executed", "landscape_jobs_executed", "distance_jobs_executed",
    "clustering_jobs_executed", "gene_view_jobs_executed",
    "fusion_jobs_executed", "new_data_jobs_executed"
  )
  roles <- c(
    "minimum_query_maximum_reference", "maximum_query_minimum_reference",
    "maximum_aggregate_missing_burden", "median_query_nonzero_missing"
  )
  if (!is.data.frame(manifest) || !all(required %in% names(manifest)) ||
      nrow(manifest) != 4L ||
      !identical(as.integer(manifest$group_order), 1:4) ||
      !setequal(manifest$pilot_role, roles) ||
      anyDuplicated(manifest$group_id) ||
      anyDuplicated(manifest$held_out_study) ||
      length(unique(manifest$seed)) != 1L ||
      any(manifest$training_samples + manifest$held_out_samples != 90L) ||
      any(manifest$panel_genes != 500L) ||
      any(manifest$cells_per_sample != 384L) ||
      any(manifest$pca_components != 30L) ||
      any(manifest$maximum_heavy_workers != 1L) ||
      any(as.matrix(manifest[zero]) != 0) ||
      any(!grepl("^[0-9a-f]{64}$", manifest$private_cache_sha256)) ||
      any(!grepl("^mv05d1_sct_cell_fold_v1:[0-9a-f]{64}$",
                 manifest$fold_cache_key))) {
    stop("MV5-F pilot manifest violates its structural contract.", call. = FALSE)
  }
  .mv05f_assert_label_closed(manifest, "MV5-F pilot manifest")
  invisible(manifest)
}

mv05f_runtime_v1 <- function() {
  package_version <- function(package) {
    if (!requireNamespace(package, quietly = TRUE)) {
      stop("Required MV5-F package is unavailable: ", package, call. = FALSE)
    }
    as.character(utils::packageVersion(package))
  }
  list(
    contract_id = "mv05f_runtime_v1", r_version = R.version.string,
    seurat_version = package_version("Seurat"),
    seurat_object_version = package_version("SeuratObject"),
    sctransform_version = package_version("sctransform"),
    matrix_version = package_version("Matrix"),
    digest_version = package_version("digest"),
    rng_kind = unname(RNGkind()),
    omp_num_threads = Sys.getenv("OMP_NUM_THREADS", unset = ""),
    openblas_num_threads = Sys.getenv("OPENBLAS_NUM_THREADS", unset = ""),
    mkl_num_threads = Sys.getenv("MKL_NUM_THREADS", unset = "")
  )
}

mv05f_group_identity_v1 <- function(
    manifest_row, d1_record, raw_source_hashes, sct_source_keys,
    implementation_sha256, runtime = mv05f_runtime_v1()) {
  required_row <- c(
    "fold_id", "fit_scope_id", "held_out_study", "seed",
    "cells_per_sample", "fold_cache_key", "outcome_label_state",
    "biological_outcomes_computed"
  )
  if (!is.data.frame(manifest_row) || nrow(manifest_row) != 1L ||
      !all(required_row %in% names(manifest_row)) ||
      manifest_row$outcome_label_state != "closed" ||
      as.logical(manifest_row$biological_outcomes_computed)) {
    stop("manifest_row must contain exactly one MV5-F group.", call. = FALSE)
  }
  mv05d1_validate_cell_fold_record_v1(d1_record)
  ids <- sort(c(d1_record$identity$training_ids, d1_record$identity$query_ids),
              method = "radix")
  if (!identical(d1_record$cache_key, manifest_row$fold_cache_key) ||
      !identical(d1_record$identity$seed, as.integer(manifest_row$seed)) ||
      !identical(d1_record$identity$held_out_study,
                 manifest_row$held_out_study) ||
      !identical(sort(names(raw_source_hashes), method = "radix"), ids) ||
      !identical(sort(names(sct_source_keys), method = "radix"), ids)) {
    stop("MV5-F group inputs differ from the accepted fold identity.",
         call. = FALSE)
  }
  identity <- list(
    contract_id = "mv05f_mapping_group_identity_v1",
    mapping_contract_id = "mv05_inductive_mapping_v1",
    reference_reconstruction_contract_id =
      "mv05f_fixed_d1_panel_joint_reference_sct_v1",
    query_reconstruction_contract_id =
      "mv05f_fixed_cells_individual_query_sct_v1",
    missing_feature_policy =
      "fixed_panel_active_intersection_no_replacement_v1",
    fold_id = d1_record$identity$fold_id,
    fit_scope_id = d1_record$identity$fit_scope_id,
    held_out_study = d1_record$identity$held_out_study,
    seed = d1_record$identity$seed,
    training_ids = d1_record$identity$training_ids,
    query_ids = d1_record$identity$query_ids,
    panel = d1_record$payload$panel,
    dimensions = seq_len(d1_record$identity$n_components),
    cells_per_sample = as.integer(manifest_row$cells_per_sample),
    d1_fold_cache_key = d1_record$cache_key,
    d1_payload_sha256 = d1_record$payload_sha256,
    raw_source_hashes = raw_source_hashes[ids],
    sct_source_keys = sct_source_keys[ids],
    implementation_sha256 = implementation_sha256,
    runtime = runtime, outcome_label_state = "closed",
    biological_outcomes_computed = FALSE
  )
  identity$cache_key <- paste0(
    "mv05f_mapping_group_v1:", .mv05f_sha256(identity)
  )
  identity
}

mv05f_validate_group_record_v1 <- function(record) {
  forbidden <- c(
    "diagrams", "landscapes", "distances", "rankings", "clustering",
    "gene_views", "fusion", "outcomes", "labels", "tissue"
  )
  if (!inherits(record, "scph_mv05f_mapping_group_v1") || !is.list(record) ||
      !identical(record$contract_id, "mv05f_mapping_group_record_v1") ||
      any(forbidden %in% names(record$payload)) ||
      !identical(record$outcome_label_state, "closed") ||
      !identical(record$biological_outcomes_computed, FALSE) ||
      !identical(record$cache_key, record$identity$cache_key) ||
      !identical(record$payload_sha256, .mv05f_sha256(record$payload))) {
    stop("MV5-F mapping record violates identity, payload, or scope.",
         call. = FALSE)
  }
  identity <- record$identity
  payload <- record$payload
  identity_without_key <- identity
  identity_without_key$cache_key <- NULL
  ids <- sort(c(identity$training_ids, identity$query_ids), method = "radix")
  coordinates <- payload$coordinates
  mappings <- payload$query_mappings
  active <- payload$active_features
  if (!identical(identity$cache_key, paste0(
        "mv05f_mapping_group_v1:", .mv05f_sha256(identity_without_key)
      )) ||
      !is.list(coordinates) ||
      !identical(sort(names(coordinates), method = "radix"), ids) ||
      any(vapply(coordinates, function(value) {
        !is.matrix(value) || !is.numeric(value) || nrow(value) != 384L ||
          ncol(value) != 30L || anyNA(value) || any(!is.finite(value)) ||
          !identical(colnames(value), paste0("PC", 1:30))
      }, logical(1L))) ||
      !is.list(mappings) ||
      !identical(sort(names(mappings), method = "radix"),
                 sort(identity$query_ids, method = "radix")) ||
      !is.list(active) ||
      !identical(sort(names(active), method = "radix"),
                 sort(identity$query_ids, method = "radix")) ||
      any(vapply(names(mappings), function(sample_id) {
        value <- mappings[[sample_id]]
        invalid <- tryCatch({mv05_validate_inductive_mapping_v1(value); FALSE},
                 error = function(error) TRUE)
        invalid || value$held_out_sample_id != sample_id ||
          !identical(value$reference_sample_ids,
                     sort(identity$training_ids, method = "radix")) ||
          value$fold_id != identity$fold_id ||
          value$fit_scope_id != identity$fit_scope_id ||
          value$seed != identity$seed ||
          !identical(value$dimensions, identity$dimensions) ||
          !identical(value$features, active[[sample_id]]) ||
          !identical(value$query_embeddings, coordinates[[sample_id]]) ||
          !all(active[[sample_id]] %in% identity$panel$gene) ||
          anyDuplicated(active[[sample_id]])
      }, logical(1L))) ||
      any(vapply(names(coordinates), function(sample_id) {
        !all(startsWith(
          rownames(coordinates[[sample_id]]), paste0(sample_id, "__")
        ))
      }, logical(1L))) ||
      !identical(payload$reference_identity_sha256_before,
                 payload$reference_identity_sha256_after) ||
      !identical(payload$coordinate_set_sha256,
                 .mv05f_sha256(coordinates)) ||
      !identical(payload$downstream_execution,
                 list(ph_jobs = 0L, landscape_jobs = 0L,
                      distance_jobs = 0L, clustering_jobs = 0L,
                      gene_view_jobs = 0L, fusion_jobs = 0L,
                      biological_outcome_jobs = 0L))) {
    stop("MV5-F mapping payload is incomplete, nonfinite, or crossed scope.",
         call. = FALSE)
  }
  invisible(record)
}

mv05f_new_group_record_v1 <- function(identity, coordinates, mappings,
                                       active_features,
                                       reference_before,
                                       reference_after) {
  ids <- sort(names(coordinates), method = "radix")
  coordinates <- coordinates[ids]
  mappings <- mappings[sort(names(mappings), method = "radix")]
  payload <- list(
    contract_id = "mv05f_mapping_group_payload_v1",
    coordinates = coordinates, query_mappings = mappings,
    active_features = active_features[sort(names(active_features),
                                           method = "radix")],
    reference_identity_sha256_before = reference_before,
    reference_identity_sha256_after = reference_after,
    coordinate_set_sha256 = .mv05f_sha256(coordinates),
    downstream_execution = list(
      ph_jobs = 0L, landscape_jobs = 0L, distance_jobs = 0L,
      clustering_jobs = 0L, gene_view_jobs = 0L, fusion_jobs = 0L,
      biological_outcome_jobs = 0L
    )
  )
  record <- structure(
    list(
      contract_id = "mv05f_mapping_group_record_v1", identity = identity,
      payload = payload, payload_sha256 = .mv05f_sha256(payload),
      cache_key = identity$cache_key, outcome_label_state = "closed",
      biological_outcomes_computed = FALSE
    ),
    class = "scph_mv05f_mapping_group_v1"
  )
  mv05f_validate_group_record_v1(record)
  record
}

mv05f_validate_resource_metrics_v1 <- function(
    metrics, expected_groups = 4L, elapsed_cap_seconds = 1800,
    rss_cap_bytes = 8 * 1024^3, storage_cap_bytes = 10 * 1000^3) {
  required <- c(
    "group_id", "held_out_study", "seed", "disposition", "exit_status",
    "completed_query_mappings", "completed_coordinate_views",
    "elapsed_seconds", "peak_process_tree_rss_bytes", "private_result_bytes",
    "reference_immutable", "label_transfer_jobs_executed",
    "ph_jobs_executed", "landscape_jobs_executed", "distance_jobs_executed",
    "clustering_jobs_executed", "gene_view_jobs_executed",
    "fusion_jobs_executed", "new_data_jobs_executed",
    "biological_outcomes_computed", "outcome_label_state"
  )
  zero <- c(
    "label_transfer_jobs_executed", "ph_jobs_executed",
    "landscape_jobs_executed", "distance_jobs_executed",
    "clustering_jobs_executed", "gene_view_jobs_executed",
    "fusion_jobs_executed", "new_data_jobs_executed"
  )
  if (!is.data.frame(metrics) || !all(required %in% names(metrics)) ||
      nrow(metrics) != as.integer(expected_groups) ||
      anyDuplicated(metrics$group_id) ||
      any(metrics$disposition != "completed") || any(metrics$exit_status != 0L) ||
      any(metrics$completed_coordinate_views != 90L) ||
      any(metrics$completed_query_mappings < 1L) ||
      any(metrics$elapsed_seconds > elapsed_cap_seconds) ||
      any(metrics$peak_process_tree_rss_bytes > rss_cap_bytes) ||
      sum(metrics$private_result_bytes) > storage_cap_bytes ||
      any(!as.logical(metrics$reference_immutable)) ||
      any(as.matrix(metrics[zero]) != 0) ||
      any(metrics$outcome_label_state != "closed") ||
      any(as.logical(metrics$biological_outcomes_computed))) {
    stop("MV5-F resource metrics violate completion, scope, or caps.",
         call. = FALSE)
  }
  invisible(metrics)
}

mv05f_project_full_workload_v1 <- function(
    pilot_metrics, d3_metrics, d4_metrics, d5_metrics,
    full_groups = 75L, reserve_fraction = 0.25) {
  mv05f_validate_resource_metrics_v1(
    pilot_metrics, expected_groups = nrow(pilot_metrics)
  )
  for (value in list(d3_metrics, d4_metrics, d5_metrics)) {
    if (!is.data.frame(value) || any(value$outcome_label_state != "closed") ||
        any(as.logical(value$biological_outcomes_computed))) {
      stop("Projection sources must remain label closed.", call. = FALSE)
    }
  }
  full_groups <- as.integer(full_groups)
  reserve_fraction <- as.numeric(reserve_fraction)
  mapping_seconds <- stats::median(pilot_metrics$elapsed_seconds) * full_groups
  downstream_seconds <- c(
    ph = sum(d3_metrics$elapsed_seconds),
    landscape = sum(d4_metrics$elapsed_seconds),
    retrieval = sum(d5_metrics$elapsed_seconds)
  )
  integrated_multiplier <- 1 + reserve_fraction
  rows <- data.frame(
    component = c("mapping", names(downstream_seconds)),
    basis = c(
      "median_four_group_real_mapping_pilot_times_75",
      "accepted_sct_full_stage_with_25_percent_integrated_reserve",
      "accepted_sct_full_stage_with_25_percent_integrated_reserve",
      "accepted_sct_full_stage_with_25_percent_integrated_reserve"
    ),
    projected_worker_seconds = c(
      mapping_seconds, downstream_seconds * integrated_multiplier
    ),
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
  rows$projected_worker_hours <- rows$projected_worker_seconds / 3600
  rows
}

mv05f_project_full_workload_v2 <- function(
    pilot_metrics, fold_resources, d3_metrics, d4_metrics, d5_metrics,
    mv05c_status, reserve_fraction = 0.25,
    worker_hour_cap = 21.6, storage_cap_bytes = 10 * 1000^3,
    rss_cap_bytes = 8 * 1024^3, group_cap_seconds = 1800) {
  mv05f_validate_resource_metrics_v1(
    pilot_metrics, expected_groups = nrow(pilot_metrics)
  )
  .mv05f_assert_label_closed(fold_resources, "MV5-D1 fold resources")
  for (value in list(d3_metrics, d4_metrics, d5_metrics, mv05c_status)) {
    if (!is.data.frame(value) || any(value$outcome_label_state != "closed") ||
        any(as.logical(value$biological_outcomes_computed))) {
      stop("Projection sources must remain label closed.", call. = FALSE)
    }
  }
  if (nrow(fold_resources) != 75L ||
      length(unique(fold_resources$held_out_study)) != 15L ||
      length(unique(fold_resources$seed)) != 5L) {
    stop("MV5-F projection requires the complete 75-group D1 axis.",
         call. = FALSE)
  }
  total_training <- sum(fold_resources$training_samples)
  total_query <- sum(fold_resources$held_out_samples)
  input_rate <- stats::median(pilot_metrics$input_seconds)
  reference_rate <- stats::median(
    pilot_metrics$reference_sct_pca_seconds / pilot_metrics$training_samples
  )
  query_sct_rate <- stats::median(
    pilot_metrics$query_sct_seconds / pilot_metrics$held_out_samples
  )
  mapping_rate <- stats::median(
    pilot_metrics$mapping_seconds / pilot_metrics$held_out_samples
  )
  assembly_rate <- stats::median(pilot_metrics$assembly_seconds)
  mapping_components <- c(
    input_validation = input_rate * nrow(fold_resources),
    reference_sct_pca = reference_rate * total_training,
    query_sct = query_sct_rate * total_query,
    held_out_mapping = mapping_rate * total_query,
    coordinate_assembly = assembly_rate * nrow(fold_resources)
  )
  cell <- mv05c_status[
    mv05c_status$status == "completed" &
      mv05c_status$view_id == "cell_topology_v1" &
      mv05c_status$representation %in% c("sct_fold", "inductive_integrated"),,
    drop = FALSE
  ]
  split_status <- split(cell, cell$job_id)
  interval_ratios <- vapply(split_status, function(rows) {
    sct <- rows[rows$representation == "sct_fold", , drop = FALSE]
    integrated <- rows[rows$representation == "inductive_integrated", ,
                       drop = FALSE]
    if (nrow(sct) != 1L || nrow(integrated) != 1L) return(NA_real_)
    (integrated$h0_finite_intervals + integrated$h1_finite_intervals) /
      (sct$h0_finite_intervals + sct$h1_finite_intervals)
  }, numeric(1L))
  geometry_multiplier <- max(1, interval_ratios, na.rm = TRUE)
  reserve_multiplier <- 1 + as.numeric(reserve_fraction)
  downstream <- c(
    cell_ph = sum(d3_metrics$elapsed_seconds) * geometry_multiplier,
    landscape_distances = sum(d4_metrics$elapsed_seconds) * geometry_multiplier,
    baseline_retrieval = sum(d5_metrics$elapsed_seconds)
  )
  rows <- data.frame(
    contract_id = "mv05f_full_integrated_resource_projection_v2",
    component = c(names(mapping_components), names(downstream)),
    evidence_basis = c(
      "pilot_median_per_group", "pilot_median_per_training_sample",
      "pilot_median_per_query_sample", "pilot_median_per_query_sample",
      "pilot_median_per_group", "accepted_D3_scaled_by_MV5C_interval_ratio",
      "accepted_D4_scaled_by_MV5C_interval_ratio", "accepted_D5_same_axes"
    ),
    base_projected_worker_seconds = c(mapping_components, downstream),
    reserve_multiplier = reserve_multiplier,
    projected_worker_seconds = c(mapping_components, downstream) *
      reserve_multiplier,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
  rows$projected_worker_hours <- rows$projected_worker_seconds / 3600
  mapping_storage <- stats::median(pilot_metrics$private_result_bytes) * 75L
  downstream_storage <- sum(c(
    d3_metrics$private_result_bytes, d4_metrics$private_result_bytes,
    d5_metrics$private_result_bytes
  ))
  total_hours <- sum(rows$projected_worker_hours)
  projected_storage <- (mapping_storage + downstream_storage) *
    reserve_multiplier
  projected_rss <- max(pilot_metrics$peak_process_tree_rss_bytes) *
    reserve_multiplier
  projected_group <- max(pilot_metrics$elapsed_seconds) * reserve_multiplier
  decision <- if (
    total_hours <= worker_hour_cap && projected_storage <= storage_cap_bytes &&
      projected_rss <= rss_cap_bytes && projected_group <= group_cap_seconds
  ) "go_separately_authorized_full_label_closed_integrated_execution" else
    "no_go_or_narrow_scope_required"
  summary <- data.frame(
    contract_id = "mv05f_full_integrated_resource_decision_v2",
    pilot_groups = nrow(pilot_metrics), projected_groups = 75L,
    projected_training_sample_instances = total_training,
    projected_query_sample_instances = total_query,
    mv05c_max_integrated_to_sct_interval_ratio = max(
      interval_ratios, na.rm = TRUE
    ),
    conservative_geometry_multiplier = geometry_multiplier,
    reserve_multiplier = reserve_multiplier,
    projected_worker_hours = total_hours, worker_hour_cap = worker_hour_cap,
    projected_storage_bytes = projected_storage,
    storage_cap_bytes = storage_cap_bytes,
    projected_peak_rss_bytes = projected_rss, rss_cap_bytes = rss_cap_bytes,
    projected_max_group_seconds = projected_group,
    group_cap_seconds = group_cap_seconds,
    decision = decision, full_execution_authorized = FALSE,
    biological_outcomes_computed = FALSE, outcome_label_state = "closed",
    stringsAsFactors = FALSE
  )
  list(components = rows, summary = summary)
}
