.generate_realistic_fixture_inputs <- function(input_dir, seed) {
  dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
  input_dir <- normalizePath(input_dir, winslash = "/", mustWork = TRUE)
  sample_ids <- c("fixture_batch_a", "fixture_batch_b")
  genes <- c(
    paste0("RPL", seq_len(30L)),
    paste0("RPS", seq_len(30L)),
    paste0("MT-FIXTURE", seq_len(5L)),
    c("HBA1", "HBA2", "HBB", "HBD", "HBG1"),
    paste0("GENE", seq_len(450L))
  )

  matrices <- .with_preserved_seed(seed, {
    lapply(seq_along(sample_ids), function(i) {
      counts <- matrix(
        stats::rpois(length(genes) * 20L, lambda = 2) + 1L,
        nrow = length(genes),
        dimnames = list(genes, paste0(sample_ids[[i]], "_cell", seq_len(20L)))
      )
      if (i == 2L) {
        counts[71:110, ] <- counts[71:110, ] + 2L
      }
      Matrix::Matrix(counts, sparse = TRUE)
    })
  })
  names(matrices) <- sample_ids

  matrix_paths <- character(length(matrices))
  for (i in seq_along(matrices)) {
    synthetic_counts <- matrices[[i]]
    matrix_paths[[i]] <- file.path(input_dir, paste0(sample_ids[[i]], ".sparse.RData"))
    save(synthetic_counts, file = matrix_paths[[i]], version = 3)
  }
  metadata <- data.frame(
    `File Path` = matrix_paths,
    `Sample Name` = sample_ids,
    SRA = c("FIXTURE_SRA_A", "FIXTURE_SRA_B"),
    Tissue = c("fixture_control", "fixture_treatment"),
    Approach = "synthetic_scRNA-seq",
    `Number of Cells` = 20L,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  metadata_path <- file.path(input_dir, "metadata.csv")
  utils::write.csv(metadata, metadata_path, row.names = FALSE)
  list(metadata = metadata, metadata_path = metadata_path, matrices = matrices)
}

#' Run the realistic production-route fixture
#'
#' Generates two deterministic, redistributable synthetic count matrices and
#' runs them through the real sparse-RData loader, Seurat construction, fixed
#' QC gates, SCTransform normalization, Harmony integration, and persistent
#' homology path. Each sample has 520 detected genes and 20 cells, including 60
#' ribosomal features and five mitochondrial features, so the fixture exercises
#' rather than bypasses the production QC criteria. Scientific defaults remain
#' unchanged; the fixture explicitly requests only the Harmony PH
#' representation to keep a CI smoke run bounded.
#'
#' @param output_dir Empty or new directory for fixture inputs and run outputs.
#' @param seed One finite integer-compatible random seed.
#' @param num_cores Number of cores supplied to preprocessing and integration.
#' @param max_dimension Maximum PH homology dimension. H0 is the bounded default;
#'   H1 can be requested explicitly for performance and topology audits.
#' @param ph_poll_interval Seconds between child-process status and sampled
#'   process-tree memory observations.
#' @param ph_max_time_per_sample Maximum runtime in seconds for one PH child.
#'
#' @return A list with the pipeline `results`, stable `manifest`, generated
#'   `metadata`, and normalized artifact `files`.
#' @export
run_realistic_fixture <- function(output_dir,
                                  seed = 20260805L,
                                  num_cores = 1L,
                                  max_dimension = 0L,
                                  ph_poll_interval = 0.25,
                                  ph_max_time_per_sample = 20 * 60) {
  if (missing(output_dir) || length(output_dir) != 1L || is.na(output_dir) ||
      !nzchar(output_dir)) {
    stop("output_dir must be one non-empty path.", call. = FALSE)
  }
  if (length(seed) != 1L || is.na(seed) || !is.numeric(seed) ||
      !is.finite(seed) || seed != as.integer(seed)) {
    stop("seed must be one finite integer-compatible value.", call. = FALSE)
  }
  if (length(num_cores) != 1L || is.na(num_cores) || num_cores < 1 ||
      num_cores != as.integer(num_cores)) {
    stop("num_cores must be one positive integer-compatible value.", call. = FALSE)
  }
  if (length(max_dimension) != 1L || is.na(max_dimension) ||
      !is.numeric(max_dimension) || max_dimension < 0 ||
      max_dimension != as.integer(max_dimension)) {
    stop("max_dimension must be one non-negative integer-compatible value.",
         call. = FALSE)
  }
  if (length(ph_poll_interval) != 1L || is.na(ph_poll_interval) ||
      !is.finite(ph_poll_interval) || ph_poll_interval <= 0) {
    stop("ph_poll_interval must be one positive finite value.", call. = FALSE)
  }
  if (length(ph_max_time_per_sample) != 1L || is.na(ph_max_time_per_sample) ||
      !is.finite(ph_max_time_per_sample) || ph_max_time_per_sample <= 0) {
    stop("ph_max_time_per_sample must be one positive finite value.", call. = FALSE)
  }
  if (dir.exists(output_dir) && length(list.files(
    output_dir, all.files = TRUE, no.. = TRUE
  )) > 0L) {
    stop("output_dir must be empty to prevent stale fixture artifacts.", call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
  input_dir <- file.path(output_dir, "inputs")
  run_dir <- file.path(output_dir, "run")
  results_dir <- file.path(run_dir, "results")
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

  generated <- .generate_realistic_fixture_inputs(input_dir, seed)
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(run_dir)
  results <- .with_preserved_seed(seed, run_unified_pipeline(
    metadata_path = generated$metadata_path,
    results_dir = results_dir,
    num_cores = as.integer(num_cores),
    integration_methods = "harmony",
    MIN_CELLS = 20L,
    DIM = as.integer(max_dimension),
    THRESHOLD = -1,
    dataset_tag = "realistic_fixture",
    provenance_dir = file.path(run_dir, "provenance"),
    strict_reconciliation = TRUE,
    ph_representations = "harmony_integration",
    ph_poll_interval = ph_poll_interval,
    ph_max_time_per_sample = ph_max_time_per_sample
  ))

  flow <- results$provenance$sample_flow
  iteration_names <- vapply(results$data_iterations, `[[`, character(1), "name")
  harmony_iteration <- results$data_iterations[[which(
    iteration_names == HARMONY_INTEGRATION_LABEL
  )[[1]]]]
  expression_dimensions <- vapply(
    harmony_iteration$expr_list,
    function(x) paste(dim(x), collapse = "x"),
    character(1)
  )
  finite_fraction <- vapply(
    harmony_iteration$expr_list,
    function(x) mean(is.finite(as.matrix(x))),
    numeric(1)
  )
  attempts <- utils::read.csv(
    results$provenance$attempt_log_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  pd_file <- file.path(
    run_dir,
    paste0(
      "PD_list_dim", as.integer(max_dimension),
      "_th-1_harmony_integration_realistic_fixture.Rds"
    )
  )
  pd_list <- readRDS(pd_file)
  count_dimension <- function(dimension) {
    paste(vapply(pd_list, function(pd) {
      sum(pd[, 1] == dimension)
    }, integer(1)), collapse = ";")
  }
  manifest <- data.frame(
    schema_version = "2",
    seed = as.integer(seed),
    ph_max_dimension = as.integer(max_dimension),
    samples_loaded = nrow(flow),
    samples_eligible = sum(flow$ph_eligible),
    loaded_features = paste(flow$loaded_features, collapse = ";"),
    post_qc_cells = paste(flow$post_qc_cells, collapse = ";"),
    integration_iteration = paste(iteration_names, collapse = ";"),
    expression_dimensions = paste(expression_dimensions, collapse = ";"),
    finite_fraction_min = min(finite_fraction),
    ph_completed = sum(attempts$status == "completed"),
    h0_feature_counts = count_dimension(0L),
    h1_feature_counts = if (max_dimension >= 1L) count_dimension(1L) else "not_requested",
    ph_elapsed_seconds = paste(round(attempts$elapsed_seconds, 6L), collapse = ";"),
    ph_process_tree_peak_rss_bytes = paste(
      attempts$process_tree_peak_rss_bytes, collapse = ";"
    ),
    ph_process_tree_peak_count = paste(
      attempts$process_tree_peak_count, collapse = ";"
    ),
    ph_sha256_rounded_8 = digest::digest(
      lapply(pd_list, function(pd) round(as.matrix(pd), 8L)),
      algo = "sha256",
      serialize = TRUE
    ),
    input_sha256 = digest::digest(generated$matrices, algo = "sha256", serialize = TRUE),
    harmony_sha256_rounded_8 = digest::digest(
      lapply(harmony_iteration$expr_list, function(x) round(as.matrix(x), 8L)),
      algo = "sha256",
      serialize = TRUE
    ),
    scPHcompare_version = as.character(utils::packageVersion("scPHcompare")),
    seurat_version = as.character(utils::packageVersion("Seurat")),
    harmony_version = as.character(utils::packageVersion("harmony")),
    ripserr_version = as.character(utils::packageVersion("ripserr")),
    r_version = paste(R.version$major, R.version$minor, sep = "."),
    stringsAsFactors = FALSE
  )
  reference_path <- system.file(
    "extdata", "realistic_fixture_reference.csv", package = "scPHcompare"
  )
  reference <- NULL
  if (seed == 20260805L && nzchar(reference_path)) {
    reference <- utils::read.csv(
      reference_path, stringsAsFactors = FALSE, check.names = FALSE
    )
    reference <- reference[
      reference$ph_max_dimension == as.integer(max_dimension), , drop = FALSE
    ]
  }
  if (!is.null(reference) && nrow(reference) == 1L) {
    exact_fields <- c(
      "schema_version", "seed", "ph_max_dimension", "samples_loaded", "samples_eligible",
      "loaded_features", "post_qc_cells", "integration_iteration",
      "expression_dimensions", "ph_completed", "input_sha256",
      "harmony_sha256_rounded_8", "h0_feature_counts", "h1_feature_counts",
      "ph_sha256_rounded_8"
    )
    mismatched <- exact_fields[vapply(exact_fields, function(field) {
      !identical(as.character(manifest[[field]]), as.character(reference[[field]]))
    }, logical(1))]
    tolerance <- as.numeric(reference$numeric_tolerance[[1]])
    finite_ok <- manifest$finite_fraction_min + tolerance >=
      as.numeric(reference$finite_fraction_minimum[[1]])
    if (length(mismatched) > 0L || !finite_ok) {
      stop(
        "Realistic fixture reference contract changed: ",
        paste(c(mismatched, if (!finite_ok) "finite_fraction_min"), collapse = ", "),
        call. = FALSE
      )
    }
  }
  manifest_file <- file.path(output_dir, "realistic_fixture_manifest.csv")
  utils::write.csv(manifest, manifest_file, row.names = FALSE, na = "")
  files <- c(
    metadata = generated$metadata_path,
    manifest = manifest_file,
    sample_flow = results$provenance$sample_flow_file,
    attempts = results$provenance$attempt_log_file,
    persistence_diagrams = pd_file,
    pipeline_metrics = results$provenance$pipeline_metrics_file
  )
  file_names <- names(files)
  files <- normalizePath(files, winslash = "/", mustWork = TRUE)
  names(files) <- file_names
  list(
    results = results,
    manifest = manifest,
    metadata = generated$metadata,
    files = files
  )
}
