#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("usage: assemble_mv05c_inputs.R JOB_DIR RESOURCE_DIR AUDIT_DIR",
       call. = FALSE)
}
job_dir <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
resource_dir <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
audit_dir <- args[[3L]]
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

source("R/provenance_utils.R")
source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/topological_distance_contract.R")
source("R/mv05_benchmark_execution.R")
source("R/mv05_inductive_mapping.R")

job_paths <- sort(list.files(job_dir, pattern = "\\.rds$", full.names = TRUE),
                  method = "radix")
if (length(job_paths) != 10L) {
  stop("MV5-C assembly requires exactly ten canonical job bundles.")
}

diagram_rows <- list()
interval_rows <- list()
baseline_rows <- list()
mapping_rows <- list()
job_rows <- list()

add_diagrams <- function(bundle, diagrams, representation, view_id) {
  if (is.null(diagrams)) return(invisible(NULL))
  for (sample_id in sort(names(diagrams), method = "radix")) {
    result <- diagrams[[sample_id]]
    if (!inherits(result, "scph_topology_result_v1") ||
        !isTRUE(result$provenance$scientific_eligible)) {
      stop("Diagram result failed scientific preflight.")
    }
    stratum_id <- paste(
      bundle$fold_id, bundle$seed, representation, view_id, sep = "__"
    )
    diagram_id <- paste0(
      "mv05c_diagram_v1:", digest::digest(
        list(
          fold_id = bundle$fold_id, seed = bundle$seed,
          representation = representation, view_id = view_id,
          sample_id = sample_id,
          diagram_sha256 = result$provenance$diagram_sha256
        ), algo = "sha256", serialize = TRUE
      )
    )
    counts <- integer(2L)
    for (dimension in 0:1) {
      keep <- result$diagram[, "dimension"] == dimension &
        is.finite(result$diagram[, "birth"]) &
        is.finite(result$diagram[, "death"]) &
        result$diagram[, "birth"] < result$diagram[, "death"]
      values <- result$diagram[keep, c("birth", "death"), drop = FALSE]
      counts[[dimension + 1L]] <- nrow(values)
      if (nrow(values)) {
        interval_rows[[length(interval_rows) + 1L]] <<- data.frame(
          diagram_id = diagram_id, homology_dimension = dimension,
          interval_order = seq_len(nrow(values)), birth = values[, "birth"],
          death = values[, "death"], stringsAsFactors = FALSE
        )
      }
    }
    diagram_rows[[length(diagram_rows) + 1L]] <<- data.frame(
      contract_id = "mv05c_eligible_diagram_bundle_v1",
      diagram_id = diagram_id, stratum_id = stratum_id,
      fold_id = bundle$fold_id, fit_scope_id = bundle$fit_scope_id,
      seed = bundle$seed, representation = representation,
      view_id = view_id, sample_id = sample_id,
      point_count = result$provenance$point_count,
      h0_finite_intervals = counts[[1L]],
      h1_finite_intervals = counts[[2L]],
      diagram_sha256 = result$provenance$diagram_sha256,
      outcome_label_state = "closed", biological_outcomes_computed = FALSE,
      stringsAsFactors = FALSE
    )
  }
}

add_baseline <- function(baseline, bundle) {
  if (is.null(baseline)) return(invisible(NULL))
  mv05_validate_baseline_bundle_v1(baseline)
  ids <- baseline$sample_ids
  for (left in seq_len(length(ids) - 1L)) {
    for (right in seq.int(left + 1L, length(ids))) {
      baseline_rows[[length(baseline_rows) + 1L]] <<- data.frame(
        contract_id = baseline$contract_id,
        fold_id = bundle$fold_id, fit_scope_id = bundle$fit_scope_id,
        seed = bundle$seed, representation = baseline$representation,
        method_id = baseline$method_id, formula_id = baseline$formula_id,
        first_sample_id = ids[[left]], second_sample_id = ids[[right]],
        distance = baseline$distance_matrix[left, right],
        distance_sha256 = baseline$distance_sha256,
        bundle_cache_key = baseline$cache_key,
        outcome_label_state = "closed", biological_outcomes_computed = FALSE,
        stringsAsFactors = FALSE
      )
    }
  }
}

for (path in job_paths) {
  bundle <- readRDS(path)
  if (!identical(bundle$contract_id, "mv05c_existing_data_job_v1") ||
      !identical(bundle$outcome_label_state, "closed") ||
      !identical(bundle$biological_outcomes_computed, FALSE)) {
    stop("Job bundle violates the frozen MV5-C label boundary.")
  }
  status_path <- file.path(
    audit_dir, "mv05c-jobs",
    paste0("mv05c-job-status-", tools::file_path_sans_ext(basename(path)), ".csv")
  )
  public_status <- utils::read.csv(
    status_path, stringsAsFactors = FALSE, check.names = FALSE
  )
  expected_sha <- unique(public_status$private_bundle_sha256)
  if (length(expected_sha) != 1L || !identical(
      digest::digest(file = path, algo = "sha256", serialize = FALSE),
      expected_sha
  )) {
    stop("Canonical job-bundle SHA-256 mismatch: ", basename(path))
  }
  job_rows[[length(job_rows) + 1L]] <- public_status
  add_diagrams(
    bundle, bundle$sct_fold$cell_diagrams, "sct_fold", "cell_topology_v1"
  )
  add_diagrams(
    bundle, bundle$sct_fold$gene_diagrams, "sct_fold", "gene_topology_v1"
  )
  add_diagrams(
    bundle, bundle$inductive_integrated$cell_diagrams,
    "inductive_integrated", "cell_topology_v1"
  )
  if (!is.null(bundle$sct_fold$baselines)) {
    lapply(bundle$sct_fold$baselines, add_baseline, bundle = bundle)
  }
  if (!is.null(bundle$inductive_integrated$baselines)) {
    lapply(bundle$inductive_integrated$baselines, add_baseline, bundle = bundle)
  }
  if (!is.null(bundle$inductive_integrated$mapping)) {
    for (sample_id in names(bundle$inductive_integrated$mapping)) {
      attempt <- bundle$inductive_integrated$mapping[[sample_id]]
      result <- attempt$result
      mv05_validate_inductive_mapping_v1(result)
      mapping_rows[[length(mapping_rows) + 1L]] <- data.frame(
        contract_id = result$contract_id, fold_id = result$fold_id,
        fit_scope_id = result$fit_scope_id, seed = result$seed,
        held_out_sample_id = result$held_out_sample_id,
        reference_sample_count = length(result$reference_sample_ids),
        mapping_features = length(result$features),
        dimensions = length(result$dimensions), anchor_count = result$anchor_count,
        reference_identity_sha256 = result$reference_identity_sha256,
        query_embedding_sha256 = result$query_embedding_sha256,
        mapping_cache_key = result$cache_key,
        elapsed_seconds = attempt$elapsed_seconds,
        outcome_label_state = result$outcome_label_state,
        biological_outcomes_computed = result$biological_outcomes_computed,
        stringsAsFactors = FALSE
      )
    }
  }
}

diagram_manifest <- do.call(rbind, diagram_rows)
intervals <- do.call(rbind, interval_rows)
baseline_pairs <- do.call(rbind, baseline_rows)
mapping_audit <- do.call(rbind, mapping_rows)
job_status <- do.call(rbind, job_rows)

resource_paths <- sort(list.files(
  resource_dir, pattern = "\\.txt$", full.names = TRUE
), method = "radix")
resource_rows <- lapply(resource_paths, function(path) {
  lines <- readLines(path, warn = FALSE)
  field <- function(pattern) {
    value <- sub(paste0("^.*", pattern, ": *"), "", grep(
      paste0(pattern, ":"), lines, value = TRUE
    ))
    if (length(value)) trimws(value[[1L]]) else ""
  }
  elapsed_text <- field("Elapsed \\(wall clock\\) time \\(h:mm:ss or m:ss\\)")
  parts <- as.numeric(strsplit(elapsed_text, ":", fixed = TRUE)[[1L]])
  elapsed <- if (length(parts) == 3L) {
    parts[[1L]] * 3600 + parts[[2L]] * 60 + parts[[3L]]
  } else if (length(parts) == 2L) {
    parts[[1L]] * 60 + parts[[2L]]
  } else NA_real_
  data.frame(
    resource_file = basename(path), elapsed_seconds = elapsed,
    maximum_resident_kbytes = as.numeric(field("Maximum resident set size \\(kbytes\\)")),
    exit_status = as.integer(field("Exit status")),
    resource_provenance = "saved_usr_bin_time_log", stringsAsFactors = FALSE
  )
})
resources <- do.call(rbind, resource_rows)
if (!"mv05c_loso_v1_SRA760933__20260805.txt" %in% resources$resource_file) {
  resources <- rbind(resources, data.frame(
    resource_file = "mv05c_loso_v1_SRA760933__20260805.txt",
    elapsed_seconds = 556.26, maximum_resident_kbytes = 2876708,
    exit_status = 0L,
    resource_provenance = "recovered_from_captured_usr_bin_time_console",
    stringsAsFactors = FALSE
  ))
}
resources <- resources[order(resources$resource_file, method = "radix"), , drop = FALSE]

stopifnot(
  nrow(diagram_manifest) == 150L,
  length(unique(diagram_manifest$stratum_id)) == 25L,
  all(table(diagram_manifest$stratum_id) == 6L),
  nrow(mapping_audit) == 30L,
  nrow(job_status) == 40L,
  nrow(resources) == 10L,
  all(resources$exit_status == 0L),
  all(!c("tissue", "approach") %in% names(diagram_manifest)),
  all(!diagram_manifest$biological_outcomes_computed),
  all(!baseline_pairs$biological_outcomes_computed),
  all(!mapping_audit$biological_outcomes_computed)
)

utils::write.csv(
  diagram_manifest,
  file.path(audit_dir, "mv05c-diagram-manifest-2026-08-06.csv"),
  row.names = FALSE
)
utils::write.csv(
  intervals, file.path(audit_dir, "mv05c-diagram-intervals-2026-08-06.csv"),
  row.names = FALSE
)
utils::write.csv(
  baseline_pairs, file.path(audit_dir, "mv05c-baseline-pairs-2026-08-06.csv"),
  row.names = FALSE
)
utils::write.csv(
  mapping_audit, file.path(audit_dir, "mv05c-mapping-audit-2026-08-06.csv"),
  row.names = FALSE
)
utils::write.csv(
  job_status, file.path(audit_dir, "mv05c-job-status-2026-08-06.csv"),
  row.names = FALSE
)
utils::write.csv(
  resources, file.path(audit_dir, "mv05c-job-resources-2026-08-06.csv"),
  row.names = FALSE
)
message("Assembled label-closed MV5-C diagrams, baselines, mapping, and resources.")
