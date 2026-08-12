#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || !args[[2L]] %in% c("sentinel", "pressure")) {
  stop("Usage: run_mv05aq_numerical_remediation.R OUTPUT_DIR sentinel|pressure")
}
output_dir <- normalizePath(args[[1L]], mustWork = FALSE)
mode <- args[[2L]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
setwd("/mnt/e/Repositories/Jonah/PH Pipeline Repo/scPHcompare")
pkgload::load_all(".", quiet = TRUE)

write_csv <- function(value, id) utils::write.csv(
  value, file.path(output_dir, paste0("mv05aq-", id, "-2026-08-12.csv")),
  row.names = FALSE, na = "", quote = TRUE
)

paths <- if (mode == "sentinel") c(
  "results/mv03/diagrams/stage-b/B__bone__integrated__GSM3396163_C1.RData__20260805__cell_topology_v1.rds",
  "results/mv03/diagrams/stage-b/B__bone__integrated__GSM3396179_S1.RData__20260805__cell_topology_v1.rds"
) else c(
  "results/mv03/diagrams/stage-b/B__large__sct_whole__SRA645804_SRS2823407__20260805__gene_topology_v1.rds",
  "results/mv03/diagrams/stage-b/B__large__sct_whole__SRA728025_SRS3454425__20260805__gene_topology_v1.rds"
)
objects <- lapply(paths, readRDS)
diagrams <- lapply(objects, `[[`, "diagram")
ids <- vapply(objects, function(value) value$provenance$sample_id, character(1))
verified <- vapply(seq_along(objects), function(index) {
  identical(
    objects[[index]]$provenance$diagram_sha256,
    digest::digest(diagrams[[index]], algo = "sha256")
  ) && isTRUE(objects[[index]]$provenance$scientific_eligible)
}, logical(1))
if (!all(verified)) stop("MV5-AQ input provenance failed.")

run_one <- function(method) {
  gc()
  started <- proc.time()[["elapsed"]]
  value <- persistence_landscape_distance(
    diagrams[[1]], diagrams[[2]], method = method,
    first_id = ids[[1]], second_id = ids[[2]]
  )
  elapsed <- unname(proc.time()[["elapsed"]] - started)
  list(value = value, elapsed = elapsed)
}

if (mode == "sentinel") {
  exact <- run_one("exact")
  adaptive <- run_one("adaptive")
  auto <- run_one("auto")
  results <- list(exact = exact, adaptive = adaptive, auto = auto)
} else {
  auto <- run_one("auto")
  results <- list(auto = auto)
}

method_rows <- do.call(rbind, lapply(names(results), function(method) {
  run <- results[[method]]
  value <- run$value
  data.frame(
    contract_id = paste0("mv05aq_", mode, "_execution_v1"),
    mode = mode, method_requested = method,
    first_sample_id = ids[[1]], second_sample_id = ids[[2]],
    first_h0_intervals = value$dimensions$H0$first_finite_intervals,
    second_h0_intervals = value$dimensions$H0$second_finite_intervals,
    first_h1_intervals = value$dimensions$H1$first_finite_intervals,
    second_h1_intervals = value$dimensions$H1$second_finite_intervals,
    h0_method = value$dimensions$H0$method,
    h1_method = value$dimensions$H1$method,
    h0_distance = unname(value$distances[["H0"]]),
    h1_distance = unname(value$distances[["H1"]]),
    combined_distance = unname(value$distances[["combined"]]),
    h0_error_estimate = value$dimensions$H0$achieved_absolute_error_estimate,
    h1_error_estimate = value$dimensions$H1$achieved_absolute_error_estimate,
    h0_within_tolerance = value$dimensions$H0$within_requested_tolerance,
    h1_within_tolerance = value$dimensions$H1$within_requested_tolerance,
    elapsed_seconds = run$elapsed, result_bytes = as.numeric(object.size(value)),
    cache_key = value$cache_key, engine_version = value$provenance$engine_version,
    exact_guard = value$provenance$exact_max_intervals,
    stringsAsFactors = FALSE
  )
}))
write_csv(method_rows, paste0(mode, "-execution"))

if (mode == "sentinel") {
  exact_value <- exact$value
  adaptive_value <- adaptive$value
  dimensions <- c("H0", "H1", "combined")
  agreement <- data.frame(
    contract_id = "mv05aq_strict_exact_adaptive_agreement_v1",
    dimension = dimensions,
    exact_distance = unname(exact_value$distances[dimensions]),
    adaptive_distance = unname(adaptive_value$distances[dimensions]),
    absolute_difference = abs(unname(
      exact_value$distances[dimensions] - adaptive_value$distances[dimensions]
    )),
    comparison_limit = 1e-8,
    passed = abs(unname(
      exact_value$distances[dimensions] - adaptive_value$distances[dimensions]
    )) <= 1e-8,
    stringsAsFactors = FALSE
  )
  if (!all(agreement$passed) ||
      !all(vapply(adaptive_value$dimensions, `[[`, logical(1),
                  "within_requested_tolerance"))) {
    stop("MV5-AQ sentinel did not certify strict exact/adaptive agreement.")
  }
  write_csv(agreement, "strict-exact-adaptive-agreement")

  diagnostics <- do.call(rbind, lapply(names(adaptive_value$dimensions), function(d) {
    value <- adaptive_value$dimensions[[d]]
    data.frame(
      contract_id = "mv05aq_adaptive_error_certificate_v1",
      dimension = d, method = value$method,
      tolerance_allocation = value$tolerance_allocation,
      error_estimate_policy = value$error_estimate_policy,
      requested_absolute_tolerance = value$requested_absolute_tolerance,
      requested_relative_tolerance = value$requested_relative_tolerance,
      fine_summed_quadrature_error = value$fine_summed_quadrature_error,
      refinement_delta = value$refinement_delta,
      achieved_absolute_error_estimate = value$achieved_absolute_error_estimate,
      final_threshold = max(
        value$requested_absolute_tolerance,
        value$requested_relative_tolerance * abs(value$squared_distance)
      ),
      within_requested_tolerance = value$within_requested_tolerance,
      integration_nodes = value$integration_nodes,
      integration_subdivisions = value$integration_subdivisions,
      stringsAsFactors = FALSE
    )
  }))
  write_csv(diagnostics, "adaptive-error-certificate")

  serial_path <- file.path(output_dir, "mv05aq-sentinel-private.rds")
  saveRDS(adaptive_value, serial_path, version = 3)
  reloaded <- readRDS(serial_path)
  serialization <- data.frame(
    contract_id = "mv05aq_serialization_v1",
    object_identical = identical(adaptive_value, reloaded),
    distances_identical = identical(adaptive_value$distances, reloaded$distances),
    cache_key_identical = identical(adaptive_value$cache_key, reloaded$cache_key),
    cache_key = reloaded$cache_key,
    serialized_bytes = file.info(serial_path)$size,
    stringsAsFactors = FALSE
  )
  write_csv(serialization, "serialization")

  policy <- data.frame(
    contract_id = "mv05aq_engine_routing_policy_v1",
    public_default_method = "auto",
    exact_max_intervals = 500L,
    evidence_corpus = "accepted_mv04_56_diagram_manifest",
    observed_h0_min = 383L, observed_h0_max = 499L,
    observed_h1_min = 79L, observed_h1_max = 2802L,
    guard_semantics = "numerical_engine_routing_only",
    interval_removal = FALSE, landscape_level_cap = FALSE,
    grid_fallback = FALSE, silent_tolerance_relaxation = FALSE,
    stringsAsFactors = FALSE
  )
  write_csv(policy, "engine-routing-policy")
}

cat("MV5-AQ", mode, "execution passed; methods:",
    paste(method_rows$h0_method, method_rows$h1_method, sep = "/",
          collapse = ", "), "\n")
