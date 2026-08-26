mv17g_parallel_recovery_contract_v1 <- function(workers = 8L) {
  workers <- as.integer(workers)
  if (length(workers) != 1L || is.na(workers) || workers != 8L) {
    stop("MV17-G recovery requires exactly eight workers", call. = FALSE)
  }
  data.frame(
    contract_id = "mv17g_parallel_recovery_v1",
    workers = workers,
    threads_per_child = 1L,
    retries = 0L,
    child_timeout_seconds = 1800L,
    child_RSS_cap_bytes = 8589934592,
    concurrent_child_RSS_cap_bytes = workers * 8589934592,
    aggregate_timeout_seconds = 604800L,
    private_cap_bytes = 12884901888,
    public_cap_bytes = 67108864,
    stringsAsFactors = FALSE
  )
}

mv17g_parallel_thread_environment_v1 <- function() {
  c(
    OMP_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    BLIS_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1",
    RCPP_PARALLEL_NUM_THREADS = "1"
  )
}

mv17g_job_stem_v1 <- function(q) {
  sprintf(
    "%04d__%s__%03d__%s",
    as.integer(q$job_order),
    as.character(q$view),
    as.integer(q$unit_order),
    as.character(q$null_family)
  )
}

mv17g_job_artifacts_v1 <- function(q, private_root) {
  stem <- mv17g_job_stem_v1(q)
  c(
    result = file.path(private_root, "jobs", paste0(stem, ".rds")),
    time = file.path(private_root, "logs", paste0(stem, ".time.txt")),
    stdout = file.path(private_root, "logs", paste0(stem, ".stdout.txt")),
    stderr = file.path(private_root, "logs", paste0(stem, ".stderr.txt"))
  )
}

mv17g_checkpoint_scan_v1 <- function(queue, private_root) {
  required <- c("job_order", "view", "unit_order", "null_family")
  if (!all(required %in% names(queue)) || anyDuplicated(queue$job_order)) {
    stop("invalid MV17-G queue for checkpoint scan", call. = FALSE)
  }
  rows <- lapply(seq_len(nrow(queue)), function(i) {
    paths <- mv17g_job_artifacts_v1(queue[i, , drop = FALSE], private_root)
    present <- file.exists(paths)
    state <- if (all(present)) "complete" else if (any(present)) "partial" else "absent"
    data.frame(
      job_order = as.integer(queue$job_order[i]),
      state = state,
      present_artifacts = sum(present),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out[order(out$job_order, method = "radix"), , drop = FALSE]
}

mv17g_complete_prefix_v1 <- function(scan, require_incomplete = TRUE) {
  if (!all(c("job_order", "state") %in% names(scan)) || any(scan$state == "partial")) {
    stop("MV17-G checkpoint contains partial evidence", call. = FALSE)
  }
  complete <- sort(scan$job_order[scan$state == "complete"])
  prefix <- if (length(complete)) max(complete) else 0L
  if (length(complete) && !identical(as.integer(complete), seq_len(prefix))) {
    stop("MV17-G checkpoint is not a consecutive prefix", call. = FALSE)
  }
  if (isTRUE(require_incomplete) && prefix >= nrow(scan)) {
    stop("MV17-G recovery prefix must leave pending work", call. = FALSE)
  }
  as.integer(prefix)
}

mv17g_parallel_batches_v1 <- function(job_orders, workers = 8L) {
  workers <- as.integer(workers)
  job_orders <- as.integer(job_orders)
  if (!length(job_orders)) return(list())
  if (is.na(workers) || workers < 1L || workers > 8L || anyNA(job_orders) || anyDuplicated(job_orders)) {
    stop("invalid MV17-G parallel batch request", call. = FALSE)
  }
  split(job_orders, ceiling(seq_along(job_orders) / workers))
}

mv17g_run_parallel_wave_v1 <- function(queue_rows, private_root, matrix_root, worker_path, contract) {
  if (.Platform$OS.type != "unix") stop("MV17-G parallel recovery requires Unix/WSL", call. = FALSE)
  if (nrow(queue_rows) < 1L || nrow(queue_rows) > contract$workers) {
    stop("MV17-G wave cardinality exceeds contract", call. = FALSE)
  }
  env <- mv17g_parallel_thread_environment_v1()
  old <- Sys.getenv(names(env), unset = NA_character_)
  on.exit({
    restore <- !is.na(old)
    if (any(restore)) do.call(Sys.setenv, as.list(setNames(old[restore], names(old)[restore])))
    if (any(!restore)) Sys.unsetenv(names(old)[!restore])
  }, add = TRUE)
  do.call(Sys.setenv, as.list(env))
  run_one <- function(i) {
    q <- queue_rows[i, , drop = FALSE]
    paths <- mv17g_job_artifacts_v1(q, private_root)
    if (any(file.exists(paths))) {
      return(data.frame(job_order = q$job_order, exit_status = 125L, artifacts = sum(file.exists(paths))))
    }
    matrix_path <- file.path(matrix_root, "matrices", sprintf("%s__%03d.rds", q$view, q$unit_order))
    status <- system2(
      "/usr/bin/time",
      c(
        "-v", "-o", shQuote(paths[["time"]]),
        "timeout", as.character(contract$child_timeout_seconds),
        "Rscript", "--vanilla", shQuote(worker_path), shQuote(matrix_path),
        q$null_family, as.character(q$seed_first), as.character(q$replicate_count),
        shQuote(paths[["result"]])
      ),
      stdout = paths[["stdout"]],
      stderr = paths[["stderr"]]
    )
    data.frame(job_order = q$job_order, exit_status = as.integer(status), artifacts = sum(file.exists(paths)))
  }
  result <- parallel::mclapply(
    seq_len(nrow(queue_rows)),
    run_one,
    mc.cores = nrow(queue_rows),
    mc.preschedule = FALSE,
    mc.set.seed = FALSE
  )
  out <- do.call(rbind, result)
  out[order(out$job_order, method = "radix"), , drop = FALSE]
}
