#!/usr/bin/env Rscript
options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 10L) {
  stop(
    paste(
      "usage: build_mv17gr_resource_profile_prefreeze.R",
      "<recovery-public> <recovery-private> <failed-private>",
      "<outer-time> <outer-stdout> <outer-stderr> <matrix-root>",
      "<private-output> <public-output> <implementation-head>"
    ), call. = FALSE
  )
}
recovery_public <- normalizePath(args[[1L]], mustWork = TRUE)
recovery_private <- normalizePath(args[[2L]], mustWork = TRUE)
failed_private <- normalizePath(args[[3L]], mustWork = TRUE)
outer_time <- normalizePath(args[[4L]], mustWork = TRUE)
outer_stdout <- normalizePath(args[[5L]], mustWork = TRUE)
outer_stderr <- normalizePath(args[[6L]], mustWork = TRUE)
matrix_root <- normalizePath(args[[7L]], mustWork = TRUE)
private <- args[[8L]]
public <- args[[9L]]
implementation_head <- tolower(args[[10L]])
if (dir.exists(private) || dir.exists(public) ||
    !grepl("^[0-9a-f]{40}$", implementation_head)) {
  stop("invalid MV17-GR prefreeze target/head", call. = FALSE)
}

source("R/mv08z_landscape_production.R")
source("R/mv17_null_models.R")
source("R/mv17_calibration.R")
source("R/mv17_localization.R")
source("R/mv17_full_calibration.R")
source("R/mv17_full_calibration_geometry_v2.R")
source("R/mv17g_parallel_recovery.R")
read_csv <- .mv08z_read_csv
write_csv <- .mv08z_atomic_csv
sha256 <- .mv08z_sha256_file

plan_text <- paste(readLines("PROJECT_PLAN.md", warn = FALSE), collapse = "\n")
if (!grepl("mv17gr_resource_profile_authorized_v1", plan_text, fixed = TRUE)) {
  stop("MV17-GR resource profile lacks owner authorization", call. = FALSE)
}
processes <- system2("ps", c("-eo", "args="), stdout = TRUE)
if (any(grepl("run_mv17g_gene_geometry_recovery[.]R", processes)) ||
    any(grepl("run_mv17g_calibration_group_worker_v2[.]R", processes))) {
  stop("MV17-G corrected production process still active", call. = FALSE)
}

manifest <- read_csv(file.path(
  recovery_public, "mv17g-geometry-recovery-manifest.csv"
))
manifest_paths <- file.path(recovery_public, manifest$artifact)
if (!all(file.exists(manifest_paths)) ||
    !identical(unname(as.numeric(file.info(manifest_paths)$size)),
               unname(as.numeric(manifest$bytes))) ||
    !identical(unname(vapply(manifest_paths, sha256, character(1L))),
               unname(tolower(manifest$sha256)))) {
  stop("MV17-GR recovery manifest drift", call. = FALSE)
}
private_binding_source <- read_csv(file.path(
  recovery_public, "mv17g-geometry-recovery-private-binding.csv"
))
private_paths <- file.path(recovery_private, private_binding_source$artifact)
if (!all(file.exists(private_paths)) ||
    !identical(unname(as.numeric(file.info(private_paths)$size)),
               unname(as.numeric(private_binding_source$bytes))) ||
    !identical(unname(vapply(private_paths, sha256, character(1L))),
               unname(tolower(private_binding_source$sha256)))) {
  stop("MV17-GR recovery private binding drift", call. = FALSE)
}

queue <- read_csv(file.path(
  recovery_private, "mv17g-corrected-gene-primary-queue.csv"
))
matrix_catalog <- read_csv(file.path(
  recovery_private, "mv17g-geometry-recovery-matrix-catalog.csv"
))
matrix_paths <- file.path(
  matrix_root, "matrices",
  sprintf("%s__%03d.rds", matrix_catalog$view, matrix_catalog$unit_order)
)
if (nrow(queue) != 660L ||
    !identical(as.integer(queue$job_order), 529:1188) ||
    nrow(matrix_catalog) != 264L || !all(file.exists(matrix_paths)) ||
    !identical(unname(as.numeric(file.info(matrix_paths)$size)),
               unname(as.numeric(matrix_catalog$bytes))) ||
    !identical(unname(vapply(matrix_paths, sha256, character(1L))),
               unname(tolower(matrix_catalog$sha256)))) {
  stop("MV17-GR queue/matrix drift", call. = FALSE)
}

expected_orders <- 529:536
expected_success <- c(TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE)
inventory_parts <- vector("list", length(expected_orders))
failure_rows <- vector("list", length(expected_orders))
for (j in seq_along(expected_orders)) {
  order <- expected_orders[[j]]
  q <- queue[queue$job_order == order, , drop = FALSE]
  paths <- mv17g_job_artifacts_v1(q, failed_private)
  present <- file.exists(paths)
  required <- if (expected_success[[j]]) rep(TRUE, 4L) else c(
    result = FALSE, time = TRUE, stdout = TRUE, stderr = TRUE
  )
  if (!identical(unname(present), unname(required)) ||
      file.info(paths[["stdout"]])$size != 0L ||
      file.info(paths[["stderr"]])$size != 0L) {
    stop("MV17-GR first-wave artifact drift at order ", order, call. = FALSE)
  }
  resource <- mv17c_parse_gnu_time_v1(paths[["time"]])
  if (expected_success[[j]]) {
    result <- readRDS(paths[["result"]])
    if (!identical(result$contract_id, "mv17g_group_result_v2") ||
        !identical(result$geometry, "euclidean_correlation_chord_v1") ||
        !identical(result$view, "gene") || !isTRUE(result$finite) ||
        result$matrix_sha256 != sha256(file.path(
          matrix_root, "matrices", sprintf("gene__%03d.rds", q$unit_order)
        )) || result$labels_opened || result$outcomes_opened ||
        resource$exit_status != 0L) {
      stop("MV17-GR successful payload drift at order ", order, call. = FALSE)
    }
  } else if (resource$exit_status != 124L) {
    stop("MV17-GR expected timeout receipt drift at order ", order, call. = FALSE)
  }
  existing <- paths[present]
  inventory_parts[[j]] <- do.call(rbind, lapply(names(existing), function(role) {
    data.frame(
      contract_id = "mv17gr_failure_inventory_v1", job_order = order,
      role = role, artifact = normalizePath(existing[[role]]),
      bytes = as.numeric(file.info(existing[[role]])$size),
      sha256 = sha256(existing[[role]]),
      disposition = if (expected_success[[j]])
        "preserved_corrected_profile_success" else "preserved_timeout_failure",
      stringsAsFactors = FALSE
    )
  }))
  failure_rows[[j]] <- data.frame(
    contract_id = "mv17gr_failure_summary_private_v1",
    job_order = order, unit_order = q$unit_order,
    null_family = q$null_family,
    success = expected_success[[j]], exit_status = resource$exit_status,
    wall_seconds = resource$wall_seconds,
    maximum_RSS_bytes = resource$maximum_RSS_bytes,
    stringsAsFactors = FALSE
  )
}
failure_inventory <- do.call(rbind, inventory_parts)
failure_private <- do.call(rbind, failure_rows)
if (nrow(failure_inventory) != 29L || sum(failure_private$success) != 5L ||
    sum(!failure_private$success) != 3L ||
    any(failure_private$maximum_RSS_bytes[!failure_private$success] <= 8 * 1024 ^ 3)) {
  stop("MV17-GR failure cardinality/resource drift", call. = FALSE)
}
outer <- mv17c_parse_gnu_time_v1(outer_time)
if (outer$exit_status != 1L || file.info(outer_stderr)$size < 1L ||
    file.info(outer_stdout)$size < 1L ||
    !any(grepl("wave failed", readLines(outer_stderr, warn = FALSE),
               fixed = TRUE))) {
  stop("MV17-GR outer failure receipt drift", call. = FALSE)
}

add_case <- function(source_order, case_role, engines) {
  q <- queue[queue$job_order == source_order, , drop = FALSE]
  do.call(rbind, lapply(engines, function(engine) data.frame(
    case_role = case_role, source_job_order = source_order,
    unit_order = q$unit_order, null_family = q$null_family,
    seed = if (q$null_family == "observed") 0L else q$seed_first,
    engine = engine, stringsAsFactors = FALSE
  )))
}
profile_queue <- rbind(
  add_case(529L, "observed_control", c(
    "ripserr_infinite_v1", "ripserr_cone_exact_v1", "gudhi_cone_exact_v1"
  )),
  add_case(535L, "completed_null_control", c(
    "ripserr_infinite_v1", "ripserr_cone_exact_v1", "gudhi_cone_exact_v1"
  )),
  add_case(530L, "failed_coordinate", "ripserr_cone_exact_v1"),
  add_case(532L, "failed_radial", "ripserr_cone_exact_v1"),
  add_case(533L, "failed_shuffle", c(
    "ripserr_cone_exact_v1", "gudhi_cone_exact_v1"
  ))
)
profile_queue$profile_order <- seq_len(nrow(profile_queue))
profile_queue$contract_id <- "mv17gr_profile_queue_v1"
profile_queue <- profile_queue[c(
  "contract_id", "profile_order", "case_role", "source_job_order",
  "unit_order", "null_family", "seed", "engine"
)]
if (nrow(profile_queue) != 10L || anyDuplicated(profile_queue$profile_order) ||
    sum(grepl("control$", profile_queue$case_role)) != 6L ||
    sum(grepl("^failed_", profile_queue$case_role)) != 4L) {
  stop("MV17-GR profile queue drift", call. = FALSE)
}

implementation_files <- c(
  "R/mv17gr_exact_h1_resource_profile.R",
  "scripts/run_mv17gr_exact_h1_profile_worker.R",
  "scripts/run_mv17gr_exact_h1_resource_profile.R",
  "scripts/build_mv17gr_resource_profile_prefreeze.R",
  "scripts/build_mv17gr_resource_profile_closure.R",
  "tests/testthat/test-mv17gr-exact-h1-resource-profile.R"
)
implementation <- data.frame(
  contract_id = "mv17gr_implementation_binding_v1",
  file = implementation_files,
  bytes = as.numeric(file.info(implementation_files)$size),
  sha256 = vapply(implementation_files, sha256, character(1L)),
  stringsAsFactors = FALSE
)
if (anyNA(implementation$bytes)) stop("MV17-GR implementation missing", call. = FALSE)

dir.create(private, recursive = TRUE)
dir.create(public, recursive = TRUE)
private_items <- list(
  "mv17gr-profile-queue.csv" = profile_queue,
  "mv17gr-failure-inventory.csv" = failure_inventory,
  "mv17gr-failure-summary-private.csv" = failure_private,
  "mv17gr-matrix-catalog.csv" = matrix_catalog
)
for (name in names(private_items)) {
  write_csv(private_items[[name]], file.path(private, name))
}
private_names <- names(private_items)
private_binding <- data.frame(
  contract_id = "mv17gr_private_binding_v1",
  artifact = private_names,
  rows = vapply(private_items, nrow, integer(1L)),
  bytes = as.numeric(file.info(file.path(private, private_names))$size),
  sha256 = vapply(file.path(private, private_names), sha256, character(1L)),
  tracking_state = "private_not_tracked", stringsAsFactors = FALSE
)
failure_public <- do.call(rbind, lapply(
  split(failure_private, failure_private$success),
  function(z) data.frame(
    contract_id = "mv17gr_failure_summary_public_v1",
    status = if (z$success[[1L]]) "success" else "timeout",
    children = nrow(z), minimum_wall_seconds = min(z$wall_seconds),
    maximum_wall_seconds = max(z$wall_seconds),
    maximum_RSS_bytes = max(z$maximum_RSS_bytes),
    stringsAsFactors = FALSE
  )
))
source_paths <- c(
  file.path(recovery_public, "mv17g-geometry-recovery-manifest.csv"),
  outer_time, outer_stdout, outer_stderr, "PROJECT_PLAN.md"
)
source_binding <- data.frame(
  contract_id = "mv17gr_source_binding_v1",
  role = c("recovery_prefreeze_manifest", "failed_outer_time",
           "failed_outer_stdout", "failed_outer_stderr",
           "owner_authorization_plan"),
  bytes = as.numeric(file.info(source_paths)$size),
  sha256 = vapply(source_paths, sha256, character(1L)),
  stringsAsFactors = FALSE
)
contract <- data.frame(
  contract_id = "mv17gr_exact_h1_resource_profile_v1",
  implementation_head = implementation_head,
  execution_authorized_after_commit = TRUE,
  failed_production_retry_authorized = FALSE,
  profile_attempts = 10L, workers = 1L, threads_per_child = 1L,
  retries = 0L, attempt_timeout_seconds = 3600L,
  attempt_address_space_cap_bytes = 80 * 1024 ^ 3,
  aggregate_timeout_seconds = 10 * 3600L,
  private_cap_bytes = 4 * 1024 ^ 3,
  public_cap_bytes = 16 * 1024 ^ 2,
  geometry = "euclidean_correlation_chord_v1",
  filtration = "complete_vietoris_rips_field_2_H1",
  exact_cone_stop = TRUE,
  exact_cone_reason = "one_vertex_adjacent_to_all_points_makes_complex_a_cone",
  labels_opened = FALSE, outcomes_opened = FALSE,
  downstream_surfaces = "closed", stringsAsFactors = FALSE
)
validation <- data.frame(
  contract_id = "mv17gr_prefreeze_validation_v1",
  check_id = c(
    "recovery_manifest_bound", "recovery_private_bound", "queue_660",
    "queue_orders_529_1188", "matrices_264", "all_matrices_bound",
    "controller_absent", "worker_absent", "first_wave_eight_children",
    "first_wave_five_success", "first_wave_three_timeout",
    "failure_artifacts_29", "success_payloads_typed",
    "failure_streams_empty", "outer_exit_one", "outer_failure_bound",
    "profile_queue_ten", "control_attempts_six", "failed_attempts_four",
    "full_controls_only", "cone_ripser_failed_three",
    "gudhi_worst_failed_one", "one_worker", "one_thread",
    "zero_retries", "one_hour_attempt_cap", "80_GiB_address_cap",
    "production_retry_closed", "exact_cone_argument_frozen",
    "implementation_bound", "owner_authorization_recorded",
    "labels_closed", "outcomes_closed", "downstream_closed",
    "aggregate_only_public"
  ),
  passed = c(
    nrow(manifest) >= 1L, nrow(private_binding_source) >= 1L,
    nrow(queue) == 660L, identical(as.integer(queue$job_order), 529:1188),
    nrow(matrix_catalog) == 264L, all(file.exists(matrix_paths)),
    !any(grepl("run_mv17g_gene_geometry_recovery[.]R", processes)),
    !any(grepl("run_mv17g_calibration_group_worker_v2[.]R", processes)),
    nrow(failure_private) == 8L, sum(failure_private$success) == 5L,
    sum(!failure_private$success) == 3L, nrow(failure_inventory) == 29L,
    TRUE, all(file.info(failure_inventory$artifact[
      failure_inventory$role %in% c("stdout", "stderr")
    ])$size == 0L), outer$exit_status == 1L,
    file.info(outer_stderr)$size >= 1L, nrow(profile_queue) == 10L,
    sum(grepl("control$", profile_queue$case_role)) == 6L,
    sum(grepl("^failed_", profile_queue$case_role)) == 4L,
    all(profile_queue$engine[profile_queue$engine == "ripserr_infinite_v1"] %in%
      "ripserr_infinite_v1") && sum(profile_queue$engine == "ripserr_infinite_v1") == 2L,
    sum(profile_queue$engine == "ripserr_cone_exact_v1" &
          grepl("^failed_", profile_queue$case_role)) == 3L,
    sum(profile_queue$engine == "gudhi_cone_exact_v1" &
          profile_queue$case_role == "failed_shuffle") == 1L,
    contract$workers == 1L, contract$threads_per_child == 1L,
    contract$retries == 0L, contract$attempt_timeout_seconds == 3600L,
    contract$attempt_address_space_cap_bytes == 80 * 1024 ^ 3,
    !contract$failed_production_retry_authorized, contract$exact_cone_stop,
    nrow(implementation) == length(implementation_files),
    grepl("mv17gr_resource_profile_authorized_v1", plan_text, fixed = TRUE),
    !contract$labels_opened, !contract$outcomes_opened,
    contract$downstream_surfaces == "closed",
    !any(c("job_order", "unit_order", "seed", "artifact") %in%
      names(failure_public))
  ), stringsAsFactors = FALSE
)
if (!all(validation$passed)) stop("MV17-GR prefreeze failed", call. = FALSE)
decision <- data.frame(
  contract_id = "mv17gr_prefreeze_decision_v1",
  decision_token = "mv17gr_resource_profile_authorized_v1",
  preserve_failed_evidence = TRUE,
  production_retry_authorized = FALSE,
  exact_profile_authorized_after_commit = TRUE,
  result_dependent_method_selection = FALSE,
  localization_authorized = FALSE,
  downstream_surfaces = "closed", stringsAsFactors = FALSE
)
items <- list(
  "mv17gr-contract.csv" = contract,
  "mv17gr-failure-summary.csv" = failure_public,
  "mv17gr-private-binding.csv" = private_binding,
  "mv17gr-source-binding.csv" = source_binding,
  "mv17gr-implementation-binding.csv" = implementation,
  "mv17gr-validation.csv" = validation,
  "mv17gr-decision.csv" = decision
)
for (name in names(items)) write_csv(items[[name]], file.path(public, name))
writeLines(c(
  "# MV17-GR exact H1 resource-profile prefreeze", "",
  "The corrected eight-worker primary stopped after five successful children and",
  "three 30-minute timeouts. All evidence is preserved; production retry is closed.",
  "", "The profile is single-worker and exact. A finite cutoff is permitted only at",
  "the minimum vertex eccentricity: that vertex is then adjacent to every point,",
  "so the Rips complex is a cone and every H1 class has died. This changes neither",
  "the correlation-chord geometry nor any H1 interval or landscape definition.",
  "", "Ten prospectively fixed attempts compare full and cone Ripserr plus exact",
  "TDA/GUDHI controls, then profile the three failed cases. Each attempt has a",
  "one-hour timeout and enforced 80-GiB address-space cap. No retry, label, outcome,",
  "clustering, fusion, localization, biology, or manuscript claim is authorized."
), file.path(public, "MV17GR_EXACT_H1_RESOURCE_PROFILE_PREFREEZE_2026-08-27.md"))
files <- sort(list.files(public))
write_csv(data.frame(
  contract_id = "mv17gr_prefreeze_manifest_v1", artifact = files,
  bytes = as.numeric(file.info(file.path(public, files))$size),
  sha256 = vapply(file.path(public, files), sha256, character(1L)),
  stringsAsFactors = FALSE
), file.path(public, "mv17gr-prefreeze-manifest.csv"))
message("Built MV17-GR resource profile prefreeze; checks=", nrow(validation),
        "; attempts=10; failed_production_retry=closed")
