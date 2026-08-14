.with_preserved_seed <- function(seed, code) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) {
    previous_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    if (had_seed) {
      assign(".Random.seed", previous_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(as.integer(seed))
  force(code)
}

.generate_toy_point_clouds <- function(seed) {
  .with_preserved_seed(seed, {
    angle <- stats::runif(1L, min = 0, max = 2 * pi)
    rotation <- matrix(
      c(cos(angle), -sin(angle), sin(angle), cos(angle)),
      nrow = 2L,
      byrow = TRUE
    )
    square <- matrix(
      c(0, 0, 1, 0, 0, 1, 1, 1),
      ncol = 2L,
      byrow = TRUE
    ) %*% rotation
    line <- matrix(
      c(0, 0, 1, 0, 2, 0, 3, 0),
      ncol = 2L,
      byrow = TRUE
    ) %*% rotation
    list(square = unname(square), line = unname(line))
  })
}

#' Run the deterministic analytical toy baseline
#'
#' Runs persistent homology on a rotated unit square and a rotated four-point
#' line. Rotation is selected from an explicit seed and preserves the known
#' topology: the square has three finite H0 features and one H1 feature, while
#' the line has three finite H0 features and no H1 feature. The caller's random
#' number generator state is restored on exit.
#'
#' The function writes the point clouds, persistence diagrams, PH-attempt log,
#' stable scientific manifest, and stage timings into `output_dir`. Wall-clock
#' timestamps and elapsed times are deliberately kept out of the scientific
#' manifest so repeated runs can compare the stable fields directly. The
#' manifest records the maximum error against the known birth/death values and
#' enforces a numerical tolerance of `1e-10`.
#'
#' @param output_dir Empty or new directory in which to write baseline files.
#' @param seed One finite integer-compatible random seed.
#' @param max_time_per_sample Maximum seconds allowed for each PH child process.
#' @param poll_interval Seconds between child-process status polls.
#'
#' @return A list containing `inputs`, `persistence_diagrams`, `manifest`,
#'   `timings`, `attempts`, and normalized output `files`.
#' @export
run_toy_baseline <- function(output_dir,
                             seed = 20260805L,
                             max_time_per_sample = 180,
                             poll_interval = 0.05) {
  if (missing(output_dir) || length(output_dir) != 1L || is.na(output_dir) ||
      !nzchar(output_dir)) {
    stop("output_dir must be one non-empty path.", call. = FALSE)
  }
  if (length(seed) != 1L || is.na(seed) || !is.numeric(seed) ||
      !is.finite(seed) || seed != as.integer(seed)) {
    stop("seed must be one finite integer-compatible value.", call. = FALSE)
  }
  if (length(max_time_per_sample) != 1L || is.na(max_time_per_sample) ||
      !is.numeric(max_time_per_sample) || max_time_per_sample <= 0) {
    stop("max_time_per_sample must be one positive number.", call. = FALSE)
  }
  if (length(poll_interval) != 1L || is.na(poll_interval) ||
      !is.numeric(poll_interval) || poll_interval <= 0) {
    stop("poll_interval must be one positive number.", call. = FALSE)
  }
  if (dir.exists(output_dir) && length(list.files(output_dir, all.files = TRUE,
                                                   no.. = TRUE)) > 0L) {
    stop("output_dir must be empty to prevent stale baseline artifacts.", call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)

  generated <- time_stage("generate_inputs", .generate_toy_point_clouds(seed))
  inputs <- generated$value
  timings <- generated$timing

  ph_results <- list()
  attempts <- list()
  quiet_log <- function(...) invisible(NULL)
  for (i in seq_along(inputs)) {
    sample_id <- names(inputs)[[i]]
    timed <- time_stage(
      "persistent_homology",
      process_and_monitor(
        expr_matrix = inputs[[i]],
        i = i,
        DIM = 1L,
        log_message = quiet_log,
        memory_threshold = 0.25,
        max_time_per_iteration = max_time_per_sample,
        poll_interval = poll_interval,
        results_file = file.path(output_dir, paste0(sample_id, "_intermediate.rds")),
        log_file = file.path(output_dir, paste0(sample_id, "_progress.csv")),
        temp_dataset_dir = file.path(output_dir, "subprocess"),
        sample_id = sample_id,
        cohort = "analytical_toy",
        representation = "point_cloud"
      ),
      sample_id = sample_id
    )
    ph_results[[sample_id]] <- timed$value$PD
    attempts[[sample_id]] <- timed$value$attempt
    timings <- rbind(timings, timed$timing)
  }
  attempts <- do.call(rbind, attempts)

  feature_counts <- lapply(ph_results, function(pd) {
    c(h0 = sum(pd[, 1] == 0), h1 = sum(pd[, 1] == 1))
  })
  observed <- c(
    square_h0 = unname(feature_counts$square[["h0"]]),
    square_h1 = unname(feature_counts$square[["h1"]]),
    line_h0 = unname(feature_counts$line[["h0"]]),
    line_h1 = unname(feature_counts$line[["h1"]])
  )
  expected <- c(square_h0 = 3L, square_h1 = 1L, line_h0 = 3L, line_h1 = 0L)
  if (!identical(as.integer(observed), as.integer(expected))) {
    stop(
      "Toy baseline topology changed: ",
      paste(names(observed), observed, sep = "=", collapse = ", "),
      call. = FALSE
    )
  }
  if (any(attempts$status != "completed") || any(attempts$exit_status != 0L)) {
    stop("One or more toy PH child processes did not complete.", call. = FALSE)
  }

  numeric_tolerance <- 1e-10
  reference_values <- c(
    square_h0_deaths = sort(ph_results$square[ph_results$square[, 1] == 0, 3]),
    square_h1_birth = ph_results$square[ph_results$square[, 1] == 1, 2],
    square_h1_death = ph_results$square[ph_results$square[, 1] == 1, 3],
    line_h0_deaths = sort(ph_results$line[ph_results$line[, 1] == 0, 3])
  )
  expected_values <- c(
    square_h0_deaths = rep(1, 3),
    square_h1_birth = 1,
    square_h1_death = sqrt(2),
    line_h0_deaths = rep(1, 3)
  )
  max_abs_reference_error <- max(abs(reference_values - expected_values))
  if (max_abs_reference_error > numeric_tolerance) {
    stop(
      "Toy baseline birth/death values exceeded the numerical tolerance: ",
      format(max_abs_reference_error, scientific = TRUE),
      call. = FALSE
    )
  }

  manifest <- data.frame(
    schema_version = "1",
    seed = as.integer(seed),
    max_dimension = 1L,
    threshold = -1,
    input_sha256 = digest::digest(inputs, algo = "sha256", serialize = TRUE),
    persistence_sha256 = digest::digest(
      ph_results, algo = "sha256", serialize = TRUE
    ),
    square_h0 = observed[["square_h0"]],
    square_h1 = observed[["square_h1"]],
    line_h0 = observed[["line_h0"]],
    line_h1 = observed[["line_h1"]],
    numeric_tolerance = numeric_tolerance,
    max_abs_reference_error = max_abs_reference_error,
    scPHcompare_version = as.character(utils::packageVersion("scPHcompare")),
    ripserr_version = as.character(utils::packageVersion("ripserr")),
    r_version = paste(R.version$major, R.version$minor, sep = "."),
    stringsAsFactors = FALSE
  )

  inputs_file <- file.path(output_dir, "toy_inputs.rds")
  pd_file <- file.path(output_dir, "persistence_diagrams.rds")
  attempts_file <- file.path(output_dir, "ph_attempt_log.csv")
  manifest_file <- file.path(output_dir, "baseline_manifest.csv")
  timings_file <- file.path(output_dir, "stage_timings.csv")

  written <- time_stage("write_artifacts", {
    saveRDS(inputs, inputs_file, version = 3)
    saveRDS(ph_results, pd_file, version = 3)
    utils::write.csv(attempts, attempts_file, row.names = FALSE, na = "")
    utils::write.csv(manifest, manifest_file, row.names = FALSE, na = "")
    invisible(NULL)
  })
  timings <- rbind(timings, written$timing)
  utils::write.csv(timings, timings_file, row.names = FALSE, na = "")

  files <- c(
    inputs = inputs_file,
    persistence_diagrams = pd_file,
    attempts = attempts_file,
    manifest = manifest_file,
    timings = timings_file
  )
  file_names <- names(files)
  files <- normalizePath(files, winslash = "/", mustWork = TRUE)
  names(files) <- file_names

  list(
    inputs = inputs,
    persistence_diagrams = ph_results,
    manifest = manifest,
    timings = timings,
    attempts = attempts,
    files = files
  )
}
