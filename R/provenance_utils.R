# Internal provenance helpers for sample-flow and persistent-homology runs.

.validate_provenance_ids <- function(ids, label) {
  ids <- as.character(ids)
  if (anyNA(ids) || any(!nzchar(ids))) {
    stop(label, " contains missing or empty sample IDs.", call. = FALSE)
  }
  duplicated_ids <- unique(ids[duplicated(ids)])
  if (length(duplicated_ids) > 0) {
    stop(
      label, " contains duplicated sample IDs: ",
      paste(utils::head(duplicated_ids, 10), collapse = ", "),
      call. = FALSE
    )
  }
  ids
}

new_sample_flow <- function(sample_ids, cohort, source_counts = NA_real_,
                            loaded_features = NA_integer_, loaded_cells = NA_integer_) {
  sample_ids <- .validate_provenance_ids(sample_ids, "sample_ids")
  n <- length(sample_ids)
  recycle_or_stop <- function(x, label) {
    if (length(x) == 1L) {
      return(rep(x, n))
    }
    if (length(x) != n) {
      stop(label, " must have length 1 or match sample_ids.", call. = FALSE)
    }
    x
  }

  data.frame(
    cohort = rep(as.character(cohort), n),
    input_order = seq_len(n),
    sample_id = sample_ids,
    source_count = recycle_or_stop(source_counts, "source_counts"),
    loaded_features = as.integer(recycle_or_stop(loaded_features, "loaded_features")),
    loaded_cells = as.integer(recycle_or_stop(loaded_cells, "loaded_cells")),
    pre_qc_cells = rep(NA_integer_, n),
    post_qc_cells = rep(NA_integer_, n),
    min_cells_required = rep(NA_integer_, n),
    ph_eligible = rep(NA, n),
    disposition = rep("loaded", n),
    reason_code = rep("", n),
    stringsAsFactors = FALSE
  )
}

record_pre_qc_sample_flow <- function(flow, sample_ids, pre_qc_cells) {
  sample_ids <- .validate_provenance_ids(sample_ids, "pre-QC sample_ids")
  if (length(pre_qc_cells) != length(sample_ids)) {
    stop("pre_qc_cells must match pre-QC sample_ids.", call. = FALSE)
  }
  idx <- match(sample_ids, flow$sample_id)
  if (anyNA(idx)) {
    stop("Pre-QC IDs are absent from the input sample flow: ",
         paste(sample_ids[is.na(idx)], collapse = ", "), call. = FALSE)
  }
  flow$pre_qc_cells[idx] <- as.integer(pre_qc_cells)
  flow
}

record_post_qc_sample_flow <- function(flow, sample_ids, post_qc_cells, min_cells) {
  sample_ids <- .validate_provenance_ids(sample_ids, "post-QC sample_ids")
  if (length(post_qc_cells) != length(sample_ids)) {
    stop("post_qc_cells must match post-QC sample_ids.", call. = FALSE)
  }
  idx <- match(sample_ids, flow$sample_id)
  if (anyNA(idx)) {
    stop("Post-QC IDs are absent from the input sample flow: ",
         paste(sample_ids[is.na(idx)], collapse = ", "), call. = FALSE)
  }
  eligible <- as.integer(post_qc_cells) >= as.integer(min_cells)
  flow$post_qc_cells[idx] <- as.integer(post_qc_cells)
  flow$min_cells_required[idx] <- as.integer(min_cells)
  flow$ph_eligible[idx] <- eligible
  flow$disposition[idx] <- ifelse(eligible, "eligible_for_ph", "excluded_before_ph")
  flow$reason_code[idx] <- ifelse(eligible, "", "excluded_post_qc_min_cells")
  flow
}

write_provenance_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = paste0(basename(path), "."), tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  utils::write.csv(x, temporary, row.names = FALSE, na = "")
  if (file.exists(path)) {
    unlink(path)
  }
  if (!file.rename(temporary, path)) {
    stop("Could not atomically write provenance file: ", path, call. = FALSE)
  }
  invisible(path)
}

append_ph_attempts <- function(attempts, path) {
  if (is.null(attempts) || nrow(attempts) == 0L) {
    return(invisible(path))
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  existing <- if (file.exists(path)) {
    utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    attempts[0, , drop = FALSE]
  }
  required <- c("cohort", "representation", "sample_id", "job_id", "attempt", "status")
  if (!all(required %in% names(existing)) || !all(required %in% names(attempts))) {
    stop("Existing PH attempt log has an incompatible schema: ", path, call. = FALSE)
  }
  columns <- union(names(existing), names(attempts))
  for (column in setdiff(columns, names(existing))) existing[[column]] <- NA
  for (column in setdiff(columns, names(attempts))) attempts[[column]] <- NA
  existing <- existing[, columns, drop = FALSE]
  attempts <- attempts[, columns, drop = FALSE]
  write_provenance_csv(rbind(existing, attempts), path)
}

new_ph_attempt <- function(cohort, representation, sample_id, job_id, attempt,
                           threshold, input_dimensions, started_at, finished_at,
                           exit_status = NA_integer_, timed_out = FALSE,
                           pd_written = FALSE, status, error_message = "",
                           poll_interval_seconds = NA_real_,
                           memory_samples = 0L,
                           monitor_peak_rss_bytes = NA_real_,
                           child_peak_rss_bytes = NA_real_,
                           descendant_peak_rss_bytes = NA_real_,
                           process_tree_peak_rss_bytes = NA_real_,
                           process_tree_peak_count = 0L) {
  elapsed <- as.numeric(difftime(finished_at, started_at, units = "secs"))
  data.frame(
    cohort = as.character(cohort),
    representation = as.character(representation),
    sample_id = as.character(sample_id),
    job_id = as.integer(job_id),
    attempt = as.integer(attempt),
    threshold = as.numeric(threshold),
    input_rows = as.integer(input_dimensions[[1]]),
    input_columns = as.integer(input_dimensions[[2]]),
    started_at = format(started_at, "%Y-%m-%dT%H:%M:%OS3%z"),
    finished_at = format(finished_at, "%Y-%m-%dT%H:%M:%OS3%z"),
    elapsed_seconds = elapsed,
    exit_status = as.integer(exit_status),
    timed_out = isTRUE(timed_out),
    pd_written = isTRUE(pd_written),
    poll_interval_seconds = as.numeric(poll_interval_seconds),
    memory_samples = as.integer(memory_samples),
    monitor_peak_rss_bytes = as.numeric(monitor_peak_rss_bytes),
    child_peak_rss_bytes = as.numeric(child_peak_rss_bytes),
    descendant_peak_rss_bytes = as.numeric(descendant_peak_rss_bytes),
    process_tree_peak_rss_bytes = as.numeric(process_tree_peak_rss_bytes),
    process_tree_peak_count = as.integer(process_tree_peak_count),
    status = as.character(status),
    error_message = substr(paste(error_message, collapse = "\n"), 1L, 4000L),
    stringsAsFactors = FALSE
  )
}

next_ph_attempt_number <- function(log_file, job_id) {
  if (is.null(log_file) || !file.exists(log_file)) {
    return(1L)
  }
  prior_log <- tryCatch(
    utils::read.csv(log_file, stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  if (is.null(prior_log) || !"job_id" %in% names(prior_log)) {
    return(1L)
  }
  as.integer(sum(prior_log$job_id == job_id, na.rm = TRUE) + 1L)
}

summarize_ph_attempts <- function(attempts) {
  if (is.null(attempts) || nrow(attempts) == 0L) {
    return(attempts)
  }
  required <- c("cohort", "representation", "sample_id", "attempt")
  missing <- setdiff(required, names(attempts))
  if (length(missing) > 0L) {
    stop("Attempt data are missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  key <- interaction(attempts[required[1:3]], drop = TRUE, lex.order = TRUE)
  ordered <- attempts[order(key, attempts$attempt), , drop = FALSE]
  final_key <- interaction(ordered[required[1:3]], drop = TRUE, lex.order = TRUE)
  ordered[!duplicated(final_key, fromLast = TRUE), , drop = FALSE]
}

assert_sample_reconciliation <- function(input_ids, excluded_ids, eligible_ids,
                                         completed_ids = character(), failed_ids = character(),
                                         strict_completion = TRUE) {
  input_ids <- .validate_provenance_ids(input_ids, "input_ids")
  excluded_ids <- .validate_provenance_ids(excluded_ids, "excluded_ids")
  eligible_ids <- .validate_provenance_ids(eligible_ids, "eligible_ids")
  completed_ids <- .validate_provenance_ids(completed_ids, "completed_ids")
  failed_ids <- .validate_provenance_ids(failed_ids, "failed_ids")

  overlap <- intersect(excluded_ids, eligible_ids)
  if (length(overlap) > 0L) {
    stop("Excluded and eligible sample sets overlap: ", paste(overlap, collapse = ", "), call. = FALSE)
  }
  partition <- c(excluded_ids, eligible_ids)
  if (!setequal(input_ids, partition) || length(partition) != length(input_ids)) {
    stop("Input IDs do not reconcile exactly with excluded and eligible IDs.", call. = FALSE)
  }
  final_overlap <- intersect(completed_ids, failed_ids)
  if (length(final_overlap) > 0L) {
    stop("Completed and failed sample sets overlap: ", paste(final_overlap, collapse = ", "), call. = FALSE)
  }
  if (length(c(completed_ids, failed_ids)) > 0L &&
      !setequal(eligible_ids, c(completed_ids, failed_ids))) {
    stop("Eligible IDs do not reconcile exactly with completed and failed IDs.", call. = FALSE)
  }
  if (isTRUE(strict_completion) && length(failed_ids) > 0L) {
    stop("PH did not complete for eligible samples: ", paste(failed_ids, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

bounded_knn_k <- function(n_observations, requested_k = 100L) {
  n_observations <- as.integer(n_observations)
  requested_k <- as.integer(requested_k)
  if (length(n_observations) != 1L || is.na(n_observations) || n_observations < 2L) {
    stop("At least two observations are required for k-nearest-neighbor estimation.", call. = FALSE)
  }
  if (length(requested_k) != 1L || is.na(requested_k) || requested_k < 1L) {
    stop("requested_k must be a positive integer.", call. = FALSE)
  }
  min(requested_k, n_observations - 1L)
}
