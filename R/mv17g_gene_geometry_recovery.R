# Runtime helpers for the prospectively gated MV17-G geometry recovery.

mv17g_recovery_thread_environment_v2 <- function() {
  c(
    OMP_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    BLIS_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1",
    RCPP_PARALLEL_NUM_THREADS = "1",
    NUMEXPR_NUM_THREADS = "1"
  )
}

mv17g_complete_queue_prefix_v2 <- function(scan, queue) {
  if (!all(c("job_order", "state") %in% names(scan)) ||
      !"job_order" %in% names(queue) || anyDuplicated(queue$job_order) ||
      any(scan$state == "partial")) {
    stop("invalid MV17-G recovery checkpoint", call. = FALSE)
  }
  ordered <- scan[match(queue$job_order, scan$job_order), , drop = FALSE]
  if (anyNA(ordered$job_order) || !identical(ordered$job_order, queue$job_order)) {
    stop("MV17-G recovery checkpoint/queue mismatch", call. = FALSE)
  }
  complete <- which(ordered$state == "complete")
  prefix <- if (length(complete)) max(complete) else 0L
  if (length(complete) && !identical(complete, seq_len(prefix))) {
    stop("MV17-G recovery checkpoint is not a queue prefix", call. = FALSE)
  }
  if (prefix < nrow(queue) && any(ordered$state[seq.int(prefix + 1L, nrow(queue))] != "absent")) {
    stop("MV17-G recovery checkpoint has evidence after its prefix", call. = FALSE)
  }
  as.integer(prefix)
}

mv17g_run_parallel_wave_v2 <- function(queue_rows, private_root, matrix_root,
                                       worker_path, contract) {
  if (.Platform$OS.type != "unix") {
    stop("MV17-G geometry recovery requires Unix/WSL", call. = FALSE)
  }
  if (nrow(queue_rows) < 1L || nrow(queue_rows) > contract$workers ||
      !all(queue_rows$view %in% c("cell", "gene"))) {
    stop("MV17-G geometry-recovery wave violates its contract", call. = FALSE)
  }
  environment <- mv17g_recovery_thread_environment_v2()
  previous <- Sys.getenv(names(environment), unset = NA_character_)
  on.exit({
    restore <- !is.na(previous)
    if (any(restore)) {
      do.call(Sys.setenv, as.list(setNames(
        previous[restore], names(previous)[restore]
      )))
    }
    if (any(!restore)) Sys.unsetenv(names(previous)[!restore])
  }, add = TRUE)
  do.call(Sys.setenv, as.list(environment))

  run_one <- function(i) {
    q <- queue_rows[i, , drop = FALSE]
    paths <- mv17g_job_artifacts_v1(q, private_root)
    if (any(file.exists(paths))) {
      return(data.frame(
        job_order = q$job_order, exit_status = 125L,
        artifacts = sum(file.exists(paths))
      ))
    }
    matrix_path <- file.path(
      matrix_root, "matrices", sprintf("%s__%03d.rds", q$view, q$unit_order)
    )
    status <- system2(
      "/usr/bin/time",
      c(
        "-v", "-o", shQuote(paths[["time"]]),
        "timeout", as.character(contract$child_timeout_seconds),
        "Rscript", "--vanilla", shQuote(worker_path), shQuote(matrix_path),
        q$view, q$null_family, as.character(q$seed_first),
        as.character(q$replicate_count), shQuote(paths[["result"]])
      ),
      stdout = paths[["stdout"]], stderr = paths[["stderr"]]
    )
    data.frame(
      job_order = q$job_order,
      exit_status = as.integer(status),
      artifacts = sum(file.exists(paths))
    )
  }

  result <- parallel::mclapply(
    seq_len(nrow(queue_rows)), run_one,
    mc.cores = nrow(queue_rows), mc.preschedule = FALSE, mc.set.seed = FALSE
  )
  output <- do.call(rbind, result)
  output[order(output$job_order, method = "radix"), , drop = FALSE]
}
