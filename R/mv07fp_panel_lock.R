# MV7-FP full-124 global-core panel-lock helpers.

.mv07fp_truth <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  tolower(trimws(as.character(x))) == "true"
}

mv07fp_combined_cache_manifest_v1 <- function(primary, added) {
  primary_required <- c(
    "sample_id", "seed", "selected_cell_sha256", "normalization_cache_key",
    "payload_contract_id", "payload_sha256", "private_cache_file",
    "private_cache_size_bytes", "private_cache_sha256", "outcome_label_state",
    "biological_outcomes_computed")
  added_required <- c(
    "sample_id", "seed", "selected_cell_sha256", "normalization_cache_key",
    "payload_sha256", "private_cache_file", "private_cache_bytes",
    "private_cache_sha256", "finite_sct_payload", "outcome_label_state",
    "biological_outcomes_computed")
  if (!all(primary_required %in% names(primary)) ||
      !all(added_required %in% names(added))) {
    stop("MV7-FP cache ledgers are incomplete.", call. = FALSE)
  }
  standardize <- function(value, source_tier, bytes_name) data.frame(
    contract_id = "mv07fp_combined_cache_manifest_v1",
    source_tier = source_tier, sample_id = as.character(value$sample_id),
    seed = as.integer(value$seed), selected_cells = 384L,
    selected_cell_sha256 = as.character(value$selected_cell_sha256),
    normalization_cache_key = as.character(value$normalization_cache_key),
    payload_contract_id = "mv05d0_sct_data_matrix_v1",
    payload_sha256 = as.character(value$payload_sha256),
    private_cache_file = as.character(value$private_cache_file),
    private_cache_bytes = as.numeric(value[[bytes_name]]),
    private_cache_sha256 = as.character(value$private_cache_sha256),
    outcome_label_state = as.character(value$outcome_label_state),
    biological_outcomes_computed =
      .mv07fp_truth(value$biological_outcomes_computed),
    stringsAsFactors = FALSE)
  p <- standardize(primary, "primary90", "private_cache_size_bytes")
  a <- standardize(added, "added34", "private_cache_bytes")
  result <- rbind(p, a)
  result <- result[order(result$seed, result$sample_id, method = "radix"),,
                   drop = FALSE]
  rownames(result) <- NULL
  if (nrow(p) != 450L || nrow(a) != 170L || nrow(result) != 620L ||
      length(unique(p$sample_id)) != 90L || length(unique(a$sample_id)) != 34L ||
      length(unique(result$sample_id)) != 124L ||
      anyDuplicated(paste(result$sample_id, result$seed)) ||
      !identical(sort(unique(result$seed)), 20260805:20260809) ||
      any(table(result$seed) != 124L) ||
      any(result$outcome_label_state != "closed") ||
      any(result$biological_outcomes_computed)) {
    stop("MV7-FP requires the exact label-closed 124 by five cache axis.",
         call. = FALSE)
  }
  result
}

mv07fp_resource_caps_v1 <- function() data.frame(
  contract_id = "mv07fp_resource_caps_v1",
  elapsed_cap_seconds = 2700,
  rss_cap_bytes = 8 * 1024^3,
  panel_size = 500L, cache_records = 620L,
  stop_policy = "stop_without_panel_authorization",
  basis = "MV6-C measured 840s scaled by source bytes with greater than 2x margin",
  stringsAsFactors = FALSE)
