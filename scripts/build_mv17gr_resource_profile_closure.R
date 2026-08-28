#!/usr/bin/env Rscript
options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8L) {
  stop(
    paste(
      "usage: build_mv17gr_resource_profile_closure.R",
      "<public-prefreeze> <private-prefreeze> <private-execution>",
      "<public-execution> <outer-time> <outer-stdout> <outer-stderr>",
      "<output>"
    ), call. = FALSE
  )
}
public_prefreeze <- normalizePath(args[[1L]], mustWork = TRUE)
private_prefreeze <- normalizePath(args[[2L]], mustWork = TRUE)
private_execution <- normalizePath(args[[3L]], mustWork = TRUE)
public_execution <- normalizePath(args[[4L]], mustWork = TRUE)
outer_time <- normalizePath(args[[5L]], mustWork = TRUE)
outer_stdout <- normalizePath(args[[6L]], mustWork = TRUE)
outer_stderr <- normalizePath(args[[7L]], mustWork = TRUE)
output <- args[[8L]]
if (dir.exists(output)) stop("MV17-GR closure exists", call. = FALSE)

source("R/mv08z_landscape_production.R")
source("R/mv17_calibration.R")
source("R/mv17gr_exact_h1_resource_profile.R")
read_csv <- .mv08z_read_csv
write_csv <- .mv08z_atomic_csv
sha256 <- .mv08z_sha256_file
verify_manifest <- function(root, name) {
  manifest <- read_csv(file.path(root, name))
  paths <- file.path(root, manifest$artifact)
  if (!all(file.exists(paths)) ||
      !identical(unname(as.numeric(file.info(paths)$size)),
                 unname(as.numeric(manifest$bytes))) ||
      !identical(unname(vapply(paths, sha256, character(1L))),
                 unname(tolower(manifest$sha256)))) {
    stop("MV17-GR manifest drift", call. = FALSE)
  }
  manifest
}
prefreeze_manifest <- verify_manifest(
  public_prefreeze, "mv17gr-prefreeze-manifest.csv"
)
execution_manifest <- verify_manifest(
  public_execution, "mv17gr-artifact-manifest.csv"
)
private_binding <- read_csv(file.path(
  public_prefreeze, "mv17gr-private-binding.csv"
))
private_paths <- file.path(private_prefreeze, private_binding$artifact)
if (!all(file.exists(private_paths)) ||
    !identical(unname(as.numeric(file.info(private_paths)$size)),
               unname(as.numeric(private_binding$bytes))) ||
    !identical(unname(vapply(private_paths, sha256, character(1L))),
               unname(tolower(private_binding$sha256)))) {
  stop("MV17-GR closure private prefreeze drift", call. = FALSE)
}
queue <- read_csv(file.path(private_prefreeze, "mv17gr-profile-queue.csv"))
ledger <- read_csv(file.path(private_execution, "mv17gr-private-ledger.csv"))
contract <- read_csv(file.path(public_prefreeze, "mv17gr-contract.csv"))
if (nrow(queue) != 10L || nrow(ledger) != 10L ||
    !identical(as.integer(queue$profile_order), as.integer(ledger$profile_order)) ||
    any(ledger$wall_seconds > contract$attempt_timeout_seconds + 10) ||
    any(ledger$maximum_RSS_bytes > contract$attempt_address_space_cap_bytes)) {
  stop("MV17-GR closure queue/ledger drift", call. = FALSE)
}
outer <- mv17c_parse_gnu_time_v1(outer_time)
if (outer$exit_status != 0L || file.info(outer_stderr)$size != 0L ||
    outer$wall_seconds > contract$aggregate_timeout_seconds ||
    !file.exists(outer_stdout)) {
  stop("MV17-GR closure outer receipt drift", call. = FALSE)
}
result_path <- function(i) {
  q <- queue[i, , drop = FALSE]
  file.path(
    private_execution, "results",
    sprintf("%02d__%s__%s.rds", q$profile_order, q$case_role, q$engine)
  )
}
results <- lapply(seq_len(nrow(queue)), function(i) {
  path <- result_path(i)
  if (!ledger$success[[i]]) return(NULL)
  if (!file.exists(path)) stop("MV17-GR successful result missing", call. = FALSE)
  result <- readRDS(path)
  if (!identical(result$contract_id, "mv17gr_exact_h1_result_v1") ||
      result$engine != queue$engine[[i]] ||
      result$null_family != queue$null_family[[i]] ||
      result$seed != queue$seed[[i]] || !result$finite ||
      result$geometry != "euclidean_correlation_chord_v1" ||
      result$labels_opened || result$outcomes_opened) {
    stop("MV17-GR successful result drift", call. = FALSE)
  }
  result
})
control_roles <- c("observed_control", "completed_null_control")
control_rows <- lapply(control_roles, function(role) {
  idx <- which(queue$case_role == role)
  success <- all(ledger$success[idx])
  difference_full_cone <- difference_cone_gudhi <- Inf
  if (success) {
    by_engine <- setNames(results[idx], queue$engine[idx])
    difference_full_cone <- mv17gr_maximum_h1_difference_v1(
      by_engine[["ripserr_infinite_v1"]]$h1,
      by_engine[["ripserr_cone_exact_v1"]]$h1
    )
    difference_cone_gudhi <- mv17gr_maximum_h1_difference_v1(
      by_engine[["ripserr_cone_exact_v1"]]$h1,
      by_engine[["gudhi_cone_exact_v1"]]$h1
    )
  }
  data.frame(
    contract_id = "mv17gr_exactness_control_v1", case_role = role,
    all_engines_succeeded = success,
    maximum_difference_full_vs_cone = difference_full_cone,
    maximum_difference_cone_vs_gudhi = difference_cone_gudhi,
    exactness_passed = success && difference_full_cone <= 1e-8 &&
      difference_cone_gudhi <= 1e-6,
    stringsAsFactors = FALSE
  )
})
controls <- do.call(rbind, control_rows)
failed <- queue[grepl("^failed_", queue$case_role), , drop = FALSE]
failed$success <- ledger$success[match(failed$profile_order, ledger$profile_order)]
failed$wall_seconds <- ledger$wall_seconds[
  match(failed$profile_order, ledger$profile_order)
]
failed$maximum_RSS_bytes <- ledger$maximum_RSS_bytes[
  match(failed$profile_order, ledger$profile_order)
]
ripser_failed <- failed$engine == "ripserr_cone_exact_v1"
all_ripser_failed_cases_succeed <- all(failed$success[ripser_failed])
worst_gudhi_succeeds <- failed$success[
  failed$case_role == "failed_shuffle" & failed$engine == "gudhi_cone_exact_v1"
]
if (length(worst_gudhi_succeeds) != 1L) {
  stop("MV17-GR worst-case GUDHI row drift", call. = FALSE)
}
recommendation <- if (all(controls$exactness_passed) &&
    all_ripser_failed_cases_succeed) {
  "implement_exact_cone_recovery_with_profiled_concurrency_and_small_chunks"
} else if (all(controls$exactness_passed) && worst_gudhi_succeeds) {
  "profile_exact_gudhi_on_remaining_failed_cases_before_recovery"
} else {
  "exact_H1_infeasible_under_profile_caps_owner_estimand_decision_required"
}
decision <- data.frame(
  contract_id = "mv17gr_resource_profile_decision_v1",
  controls_exact = all(controls$exactness_passed),
  all_failed_ripser_cone_succeeded = all_ripser_failed_cases_succeed,
  worst_failed_gudhi_cone_succeeded = worst_gudhi_succeeds,
  recommendation = recommendation,
  production_retry_authorized = FALSE,
  localization_authorized = FALSE,
  downstream_surfaces = "closed", stringsAsFactors = FALSE
)
failed_public <- failed[c(
  "case_role", "engine", "success", "wall_seconds", "maximum_RSS_bytes"
)]
failed_public$contract_id <- "mv17gr_failed_case_profile_v1"
validation <- data.frame(
  contract_id = "mv17gr_closure_validation_v1",
  check_id = c(
    "prefreeze_manifest_bound", "execution_manifest_bound",
    "private_prefreeze_bound", "ten_attempts", "queue_ledger_exact",
    "outer_exit_zero", "outer_stderr_empty", "all_attempt_caps",
    "all_success_results_valid", "two_controls", "full_cone_compared",
    "cone_gudhi_compared", "three_failed_ripser_cases",
    "one_failed_gudhi_case", "one_worker", "zero_retries",
    "no_partials", "labels_closed", "outcomes_closed",
    "production_retry_closed", "localization_closed",
    "downstream_closed", "aggregate_only_public"
  ),
  passed = c(
    nrow(prefreeze_manifest) >= 1L, nrow(execution_manifest) >= 1L,
    nrow(private_binding) == 4L, nrow(ledger) == 10L,
    identical(as.integer(queue$profile_order), as.integer(ledger$profile_order)),
    outer$exit_status == 0L, file.info(outer_stderr)$size == 0L,
    all(ledger$wall_seconds <= contract$attempt_timeout_seconds + 10) &&
      all(ledger$maximum_RSS_bytes <= contract$attempt_address_space_cap_bytes),
    all(vapply(seq_len(nrow(queue)), function(i) {
      !ledger$success[[i]] || !is.null(results[[i]])
    }, logical(1L))), nrow(controls) == 2L,
    all(is.finite(controls$maximum_difference_full_vs_cone) |
          !controls$all_engines_succeeded),
    all(is.finite(controls$maximum_difference_cone_vs_gudhi) |
          !controls$all_engines_succeeded), sum(ripser_failed) == 3L,
    sum(failed$engine == "gudhi_cone_exact_v1") == 1L,
    contract$workers == 1L, contract$retries == 0L,
    length(list.files(private_execution, pattern = "[.]partial$",
                      recursive = TRUE)) == 0L,
    !any(ledger$labels_opened), !any(ledger$outcomes_opened),
    !decision$production_retry_authorized, !decision$localization_authorized,
    decision$downstream_surfaces == "closed",
    !any(c("source_job_order", "unit_order", "seed", "artifact") %in%
      names(failed_public))
  ), stringsAsFactors = FALSE
)
if (!all(validation$passed)) stop("MV17-GR closure failed", call. = FALSE)
resource <- read_csv(file.path(public_execution, "mv17gr-resource-summary.csv"))
source_binding <- data.frame(
  contract_id = "mv17gr_closure_source_binding_v1",
  role = c("prefreeze_manifest", "execution_manifest", "outer_time",
           "outer_stdout", "outer_stderr"),
  bytes = as.numeric(file.info(c(
    file.path(public_prefreeze, "mv17gr-prefreeze-manifest.csv"),
    file.path(public_execution, "mv17gr-artifact-manifest.csv"),
    outer_time, outer_stdout, outer_stderr
  ))$size),
  sha256 = vapply(c(
    file.path(public_prefreeze, "mv17gr-prefreeze-manifest.csv"),
    file.path(public_execution, "mv17gr-artifact-manifest.csv"),
    outer_time, outer_stdout, outer_stderr
  ), sha256, character(1L)), stringsAsFactors = FALSE
)
dir.create(output, recursive = TRUE)
items <- list(
  "mv17gr-exactness-controls.csv" = controls,
  "mv17gr-failed-case-profile.csv" = failed_public,
  "mv17gr-resource-summary.csv" = resource,
  "mv17gr-source-binding.csv" = source_binding,
  "mv17gr-validation.csv" = validation,
  "mv17gr-decision.csv" = decision
)
for (name in names(items)) write_csv(items[[name]], file.path(output, name))
writeLines(c(
  "# MV17-GR exact H1 resource-profile closure", "",
  paste0("Exactness controls pass: ", decision$controls_exact, "."),
  paste0("All three failed Ripserr cone cases succeed: ",
         decision$all_failed_ripser_cone_succeeded, "."),
  paste0("Worst failed GUDHI cone case succeeds: ",
         decision$worst_failed_gudhi_cone_succeeded, "."),
  paste0("Recommendation: `", decision$recommendation, "`."), "",
  "This closure authorizes no production retry, localization, labels, outcomes,",
  "clustering, fusion, biology, manuscript claims, cleanup, or deletion."
), file.path(output, "MV17GR_EXACT_H1_RESOURCE_PROFILE_CLOSURE_2026-08-27.md"))
files <- sort(list.files(output))
write_csv(data.frame(
  contract_id = "mv17gr_closure_manifest_v1", artifact = files,
  bytes = as.numeric(file.info(file.path(output, files))$size),
  sha256 = vapply(file.path(output, files), sha256, character(1L)),
  stringsAsFactors = FALSE
), file.path(output, "mv17gr-artifact-manifest.csv"))
message("Built MV17-GR resource profile closure; recommendation=", recommendation)
