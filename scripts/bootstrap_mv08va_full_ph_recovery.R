#!/usr/bin/env Rscript

# Create a fresh, one-record completed prefix from MV8-VA's independently
# admitted job-1 output. The PH record is byte-copied, never recomputed.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 13L) stop(paste(
  "usage: bootstrap_mv08va_full_ph_recovery.R <mv08va-prefreeze>",
  "<mv08u-prefreeze> <mv08p-private> <mv08pr-private> <mv08ps-private>",
  "<mv08o-internal-private> <common-panel> <original-private>",
  "<original-public> <new-private> <new-public> <expected-head>",
  "<launch-stderr>"
), call. = FALSE)
recovery_prefreeze <- normalizePath(args[[1L]], mustWork = TRUE)
prefreeze <- normalizePath(args[[2L]], mustWork = TRUE)
source_roots <- list(
  mv08p_original_v1 = normalizePath(args[[3L]], mustWork = TRUE),
  mv08pr_overlay_v1 = normalizePath(args[[4L]], mustWork = TRUE),
  mv08ps_overlay_v1 = normalizePath(args[[5L]], mustWork = TRUE),
  mv08o_internal_primary_v8 = normalizePath(args[[6L]], mustWork = TRUE)
)
panel_path <- normalizePath(args[[7L]], mustWork = TRUE)
original_private <- normalizePath(args[[8L]], mustWork = TRUE)
original_public <- normalizePath(args[[9L]], mustWork = TRUE)
new_private <- normalizePath(args[[10L]], mustWork = FALSE)
new_public <- normalizePath(args[[11L]], mustWork = FALSE)
expected_head <- tolower(trimws(args[[12L]]))
launch_stderr <- normalizePath(args[[13L]], mustWork = TRUE)
environment_head <- tolower(trimws(Sys.getenv("MV08VA_GIT_HEAD", unset = "")))
amendment_prefreeze_value <- trimws(Sys.getenv(
  "MV08V_BOOTSTRAP_AMENDMENT_PREFREEZE", unset = ""
))
if (!grepl("^[0-9a-f]{40}$", expected_head) ||
    !identical(expected_head, environment_head)) {
  stop("MV8-VA exact recovery HEAD binding failed", call. = FALSE)
}
if (!nzchar(amendment_prefreeze_value)) {
  stop("MV8-V bootstrap amendment prefreeze absent", call. = FALSE)
}
amendment_prefreeze <- normalizePath(amendment_prefreeze_value, mustWork = TRUE)
if (dir.exists(new_private) || dir.exists(new_public)) {
  stop("MV8-VA recovery roots must be fresh", call. = FALSE)
}
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required", call. = FALSE)
source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/mv07g_sentinel.R")
source("R/mv08s_ph_sentinel.R")
read_csv <- function(path) utils::read.csv(
  path, check.names = FALSE, stringsAsFactors = FALSE
)
sha_file <- function(path) digest::digest(
  file = path, algo = "sha256", serialize = FALSE
)
atomic_csv <- function(value, path) {
  partial <- paste0(path, ".partial")
  utils::write.csv(value, partial, row.names = FALSE, quote = TRUE, na = "")
  if (!file.rename(partial, path)) stop("failed to publish ", basename(path), call. = FALSE)
}
resolve_source <- function(row) file.path(
  source_roots[[row$source_cache_root_role]], row$source_cache_relative_file
)

decision <- read_csv(file.path(recovery_prefreeze, "mv08va-decision.csv"))
failure <- read_csv(file.path(recovery_prefreeze, "mv08va-failure-receipt.csv"))
evidence <- read_csv(file.path(recovery_prefreeze, "mv08va-stopped-evidence.csv"))
implementation <- read_csv(file.path(
  recovery_prefreeze, "mv08va-implementation-bindings.csv"
))
amendment_decision <- read_csv(file.path(
  amendment_prefreeze, "mv08vc-decision.csv"
))
amendment_binding <- read_csv(file.path(
  amendment_prefreeze, "mv08vc-implementation-binding.csv"
))
amendment_manifest <- read_csv(file.path(
  amendment_prefreeze, "mv08vc-artifact-manifest.csv"
))
bootstrap_file <- "scripts/bootstrap_mv08va_full_ph_recovery.R"
bootstrap_binding <- implementation$file == bootstrap_file
if (nrow(decision) != 1L || decision$decision !=
      "authorize_no_retry_job1_bootstrap_and_resume_at_job2" ||
    nrow(failure) != 1L || !failure$job1_independently_validated ||
    sum(bootstrap_binding) != 1L ||
    !all(file.exists(implementation$file[!bootstrap_binding])) ||
    !all(vapply(implementation$file[!bootstrap_binding], sha_file, character(1L)) ==
           implementation$sha256[!bootstrap_binding]) ||
    nrow(amendment_decision) != 1L ||
    amendment_decision$decision !=
      "authorize_hash_bound_zero_retry_MV8VA_bootstrap_after_commit" ||
    !amendment_decision$corrected_bootstrap_authorized ||
    nrow(amendment_binding) != 1L ||
    amendment_binding$file != bootstrap_file ||
    amendment_binding$mv08va_sha256 != implementation$sha256[bootstrap_binding] ||
    sha_file(bootstrap_file) != amendment_binding$sha256 ||
    !all(vapply(file.path(amendment_prefreeze, amendment_manifest$artifact),
                sha_file, character(1L)) == amendment_manifest$sha256)) {
  stop("MV8-VA committed recovery prefreeze drift", call. = FALSE)
}
original_paths <- c(
  file.path(original_public, "mv08v-progress.csv"),
  file.path(original_public, "mv08v-resource-ledger.csv"),
  file.path(original_private, "logs", "ph_primary__1__stdout.txt"),
  file.path(original_private, "logs", "ph_primary__1__stderr.txt"),
  file.path(original_private, failure$output_file), launch_stderr
)
if (length(original_paths) != nrow(evidence) || !all(file.exists(original_paths)) ||
    !all(as.numeric(file.info(original_paths)$size) ==
           as.numeric(evidence$bytes)) ||
    !identical(unname(vapply(original_paths, sha_file, character(1L))),
               evidence$sha256)) {
  stop("MV8-VA stopped evidence drift", call. = FALSE)
}
queue <- read_csv(file.path(prefreeze, "mv08u-full-ph-queue.csv"))
runtime_input <- read_csv(file.path(prefreeze, "mv08u-runtime-input-bindings.csv"))
ledger <- read_csv(original_paths[[2L]])
row <- queue[queue$production_order == 1L, , drop = FALSE]
source_path <- resolve_source(row)
if (nrow(row) != 1L || nrow(ledger) != 1L ||
    ledger$attempt_id != "ph_primary__1" || ledger$disposition != "completed" ||
    sha_file(panel_path) != runtime_input$file_sha256 ||
    sha_file(source_path) != row$source_cache_sha256) {
  stop("MV8-VA accepted job-1 binding drift", call. = FALSE)
}
cache <- readRDS(source_path)
view <- mv08s_residual_gene_view_v1(cache, row, read_csv(panel_path))
record <- readRDS(original_paths[[5L]])
mv08s_validate_ph_record_v1(record, row, view)

dir.create(file.path(new_private, "ph"), recursive = TRUE)
dir.create(file.path(new_private, "logs"), recursive = TRUE)
dir.create(new_public, recursive = TRUE)
new_output <- file.path(new_private, row$output_file)
dir.create(dirname(new_output), recursive = TRUE, showWarnings = FALSE)
if (!file.copy(original_paths[[5L]], new_output, overwrite = FALSE,
               copy.mode = FALSE, copy.date = FALSE) ||
    sha_file(new_output) != ledger$output_sha256 ||
    as.numeric(file.info(new_output)$size) != ledger$output_bytes) {
  stop("MV8-VA byte-preserving job-1 copy failed", call. = FALSE)
}
atomic_csv(ledger, file.path(new_public, "mv08v-resource-ledger.csv"))
metric <- data.frame(
  contract_id = "mv08v_selected_ph_metric_v1",
  execution_head = failure$failed_execution_head,
  production_order = 1L, mv08r_job_order = row$mv08r_job_order,
  job_id = row$job_id, unit_id = row$unit_id, seed = row$seed,
  representation_id = row$representation_id, panel_id = row$panel_id,
  view_kind = row$view_kind, selected_engine = "ripserr",
  disposition = "completed", elapsed_seconds = ledger$elapsed_seconds,
  peak_process_tree_rss_bytes = ledger$peak_process_tree_rss_bytes,
  output_file = row$output_file, output_bytes = ledger$output_bytes,
  output_sha256 = ledger$output_sha256,
  point_count = record$h0_mst_oracle$point_count,
  finite_h0_intervals = record$h0_mst_oracle$finite_h0_intervals,
  finite_h1_intervals = record$h0_mst_oracle$finite_h1_intervals,
  h0_mst_passed = record$h0_mst_oracle$passed,
  diagram_sha256 = record$topology_result$provenance$diagram_sha256,
  ph_cache_key = record$cache_key, landscape_records = 0L,
  comparison_records = 0L, clustering_records = 0L, fusion_records = 0L,
  label_records = 0L, biological_outcome_records = 0L,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
atomic_csv(metric, file.path(new_public, "mv08v-selected-ph-metrics.csv"))
progress <- data.frame(
  contract_id = "mv08v_progress_v1", execution_head = failure$failed_execution_head,
  state = "mv08va_bootstrap_complete_resume_pending", completed_records = 1L,
  total_records = 1257L, last_production_order = 1L,
  last_job_id = row$job_id, aggregate_attempt_seconds = ledger$elapsed_seconds,
  private_bytes = as.numeric(file.info(new_output)$size), workers = 1L,
  retries = 0L, landscape_records = 0L, comparison_records = 0L,
  clustering_records = 0L, fusion_records = 0L, label_records = 0L,
  biological_outcome_records = 0L, outcome_label_state = "closed",
  biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
)
atomic_csv(progress, file.path(new_public, "mv08v-progress.csv"))
receipt <- data.frame(
  contract_id = "mv08va_bootstrap_receipt_v1",
  recovery_execution_head = expected_head,
  original_execution_head = failure$failed_execution_head,
  accepted_records = 1L, recomputed_records = 0L, retry_records = 0L,
  resume_at_production_order = 2L, output_bytes = ledger$output_bytes,
  output_sha256 = ledger$output_sha256,
  byte_identical_to_original = sha_file(new_output) == sha_file(original_paths[[5L]]),
  h0_mst_passed = record$h0_mst_oracle$passed,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
atomic_csv(receipt, file.path(new_public, "mv08va-bootstrap-receipt.csv"))
cat("MV8-VA bootstrap accepted 1/1 valid record; retries=0; resume_at=2\n")
