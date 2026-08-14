# Internal helpers for structured stage timing.

new_stage_timing <- function(stage, sample_id = "", started_at, finished_at,
                             status = "completed", error_message = "",
                             rss_before_bytes = NA_real_,
                             rss_after_bytes = NA_real_) {
  data.frame(
    stage = as.character(stage),
    sample_id = as.character(sample_id),
    started_at = format(started_at, "%Y-%m-%dT%H:%M:%OS3%z"),
    finished_at = format(finished_at, "%Y-%m-%dT%H:%M:%OS3%z"),
    elapsed_seconds = as.numeric(difftime(finished_at, started_at, units = "secs")),
    rss_before_bytes = as.numeric(rss_before_bytes),
    rss_after_bytes = as.numeric(rss_after_bytes),
    rss_delta_bytes = as.numeric(rss_after_bytes - rss_before_bytes),
    status = as.character(status),
    error_message = substr(paste(error_message, collapse = "\n"), 1L, 4000L),
    stringsAsFactors = FALSE
  )
}

current_process_rss <- function() {
  tryCatch(
    as.numeric(ps::ps_memory_info(ps::ps_handle())[["rss"]]),
    error = function(e) NA_real_
  )
}

sample_process_tree_rss <- function(pid) {
  empty <- list(
    root_rss_bytes = NA_real_,
    descendant_rss_bytes = NA_real_,
    tree_rss_bytes = NA_real_,
    process_count = 0L
  )
  if (length(pid) != 1L || is.na(pid) || !is.numeric(pid) || pid <= 0) {
    return(empty)
  }
  root <- tryCatch(ps::ps_handle(as.integer(pid)), error = function(e) NULL)
  if (is.null(root)) {
    return(empty)
  }
  descendants <- tryCatch(
    ps::ps_children(root, recursive = TRUE),
    error = function(e) list()
  )
  handles <- c(list(root), descendants)
  rss <- vapply(handles, function(handle) {
    tryCatch(
      as.numeric(ps::ps_memory_info(handle)[["rss"]]),
      error = function(e) NA_real_
    )
  }, numeric(1))
  valid <- is.finite(rss)
  if (!any(valid)) {
    return(empty)
  }
  root_rss <- rss[[1]]
  descendant_rss <- if (length(rss) > 1L && any(valid[-1L])) {
    sum(rss[-1L][valid[-1L]])
  } else {
    0
  }
  list(
    root_rss_bytes = if (is.finite(root_rss)) root_rss else NA_real_,
    descendant_rss_bytes = descendant_rss,
    tree_rss_bytes = sum(rss[valid]),
    process_count = as.integer(sum(valid))
  )
}

time_stage <- function(stage, code, sample_id = "") {
  started_at <- Sys.time()
  rss_before <- current_process_rss()
  value <- tryCatch(
    force(code),
    error = function(e) e
  )
  finished_at <- Sys.time()
  rss_after <- current_process_rss()
  failed <- inherits(value, "error")
  timing <- new_stage_timing(
    stage = stage,
    sample_id = sample_id,
    started_at = started_at,
    finished_at = finished_at,
    status = if (failed) "failed" else "completed",
    error_message = if (failed) conditionMessage(value) else "",
    rss_before_bytes = rss_before,
    rss_after_bytes = rss_after
  )
  if (failed) {
    attr(value, "stage_timing") <- timing
    stop(value)
  }
  list(value = value, timing = timing)
}
