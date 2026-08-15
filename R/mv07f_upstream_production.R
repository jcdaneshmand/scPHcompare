# MV7-F omitted-34 raw/SCT upstream production helpers.

.mv07f_truth <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  tolower(trimws(as.character(x))) == "true"
}

#' Build the exact raw and SCT queue from the frozen MV7-E axis.
mv07f_upstream_queue_v1 <- function(reconciliation, sample_seed_axis) {
  required_reconciliation <- c("sample_id", "post_qc_cells", "corpus_class")
  required_axis <- c("sample_id", "seed", "selected_cells",
                     "outcome_label_state", "biological_outcomes_computed")
  if (!is.data.frame(reconciliation) ||
      !all(required_reconciliation %in% names(reconciliation)) ||
      !is.data.frame(sample_seed_axis) ||
      !all(required_axis %in% names(sample_seed_axis))) {
    stop("MV7-F queue inputs are incomplete.", call. = FALSE)
  }
  omitted <- reconciliation[
    reconciliation$corpus_class == "single_study_tissue_descriptive_only",
    c("sample_id", "post_qc_cells"), drop = FALSE]
  omitted <- omitted[order(omitted$sample_id, method = "radix"), , drop = FALSE]
  axis <- sample_seed_axis[sample_seed_axis$sample_id %in% omitted$sample_id,
                           , drop = FALSE]
  axis <- axis[order(axis$seed, axis$sample_id, method = "radix"), , drop = FALSE]
  if (nrow(omitted) != 34L || nrow(axis) != 170L ||
      anyDuplicated(omitted$sample_id) ||
      anyDuplicated(paste(axis$sample_id, axis$seed)) ||
      !identical(sort(unique(axis$seed)), 20260805:20260809) ||
      any(axis$selected_cells != 384L) ||
      any(axis$outcome_label_state != "closed") ||
      any(.mv07f_truth(axis$biological_outcomes_computed))) {
    stop("MV7-F requires the fixed 34-sample, five-seed label-closed axis.",
         call. = FALSE)
  }
  raw <- data.frame(
    contract_id = "mv07f_upstream_queue_v1", stage = "raw",
    stage_order = seq_len(nrow(omitted)), sample_id = omitted$sample_id,
    seed = NA_integer_, expected_post_qc_cells = omitted$post_qc_cells,
    selected_cells = NA_integer_, queue_role = "one_atomic_raw_shard",
    stringsAsFactors = FALSE)
  cells <- omitted$post_qc_cells[match(axis$sample_id, omitted$sample_id)]
  sct <- data.frame(
    contract_id = "mv07f_upstream_queue_v1", stage = "sct",
    stage_order = seq_len(nrow(axis)), sample_id = axis$sample_id,
    seed = as.integer(axis$seed), expected_post_qc_cells = cells,
    selected_cells = 384L, queue_role = "one_atomic_sct_cache",
    stringsAsFactors = FALSE)
  result <- rbind(raw, sct)
  result$queue_order <- seq_len(nrow(result))
  result$outcome_label_state <- "closed"
  result$biological_outcomes_computed <- FALSE
  result$panel_fit <- FALSE; result$pca <- FALSE; result$ph <- FALSE
  result$landscape <- FALSE
  result[c("contract_id", "queue_order", "stage", "stage_order", "sample_id",
    "seed", "expected_post_qc_cells", "selected_cells", "queue_role",
    "outcome_label_state", "biological_outcomes_computed", "panel_fit", "pca",
    "ph", "landscape")]
}

mv07f_resource_caps_v1 <- function() {
  data.frame(
    contract_id = "mv07f_resource_caps_v1",
    scope = c("raw_child", "sct_child", "aggregate_worker", "aggregate_storage"),
    elapsed_cap_seconds = c(1800, 1800, 14400, NA),
    rss_cap_bytes = c(8, 8, 8, NA) * 1024^3,
    storage_cap_bytes = c(NA, NA, NA, 4) * 1024^3,
    stop_policy = c("kill_child_and_stop", "kill_child_and_stop",
      "stop_before_publication", "stop_before_publication"),
    safety_basis = c("MV7-D depth-extreme source profile",
      "MV7-D depth-extreme SCT profile", "greater_than_3x_linear_projection",
      "greater_than_2x_linear_projection"),
    stringsAsFactors = FALSE)
}

mv07f_repeat_target_v1 <- function(queue) {
  sct <- queue[queue$stage == "sct", , drop = FALSE]
  if (nrow(sct) != 170L) stop("MV7-F SCT queue is incomplete.", call. = FALSE)
  sample_id <- sort(unique(sct$sample_id), method = "radix")[[1L]]
  row <- sct[sct$sample_id == sample_id & sct$seed == max(sct$seed), , drop = FALSE]
  if (nrow(row) != 1L) stop("MV7-F repeat target is not unique.", call. = FALSE)
  data.frame(
    contract_id = "mv07f_representative_repeat_target_v1",
    sample_id = row$sample_id, seed = row$seed,
    selection_basis = "lexicographically_first_added_sample_last_frozen_seed",
    required_comparisons =
      "raw_counts_identity;selected_cell_identity;sct_payload_identity;cache_bytes",
    resource_selection_role = "determinism_not_worst_case",
    stringsAsFactors = FALSE)
}

mv07f_prefreeze_decision_v1 <- function(queue, caps) {
  if (sum(queue$stage == "raw") != 34L || sum(queue$stage == "sct") != 170L ||
      nrow(caps) != 4L) stop("MV7-F prefreeze decision inputs changed.", call. = FALSE)
  data.frame(
    contract_id = "mv07f_prefreeze_decision_v1",
    decision = "authorize_serial_atomic_upstream_production_only",
    raw_jobs = 34L, sct_jobs = 170L, workers = 1L, retries = 0L,
    validated_reuse_allowed = TRUE, partial_publication_allowed = FALSE,
    panel_fit_authorized = FALSE, pca_authorized = FALSE,
    ph_authorized = FALSE, landscape_authorized = FALSE,
    labels_authorized = FALSE, outcomes_authorized = FALSE,
    next_gate = "post_MV7F_exact_124_panel_lock",
    stringsAsFactors = FALSE)
}
