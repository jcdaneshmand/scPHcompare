#!/usr/bin/env Rscript

# Bind MV8-W's independent 1,280-record closure to the audited MV8-VA through
# MV8-VD recovery history. Outputs are aggregate-only: no unit identifiers,
# accessions, donor metadata, source paths, or private artifact names.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 9L) stop(paste(
  "usage: build_mv08wx_recovery_provenance_closure.R <mv08w-root>",
  "<mv08va-root> <mv08vb-root> <mv08vc-root> <mv08vd-root>",
  "<mv08v-public> <launch-stdout> <launch-stderr> <output-dir>"
), call. = FALSE)
w_root <- normalizePath(args[[1L]], mustWork = TRUE)
recovery_roots <- c(
  mv08va = normalizePath(args[[2L]], mustWork = TRUE),
  mv08vb = normalizePath(args[[3L]], mustWork = TRUE),
  mv08vc = normalizePath(args[[4L]], mustWork = TRUE),
  mv08vd = normalizePath(args[[5L]], mustWork = TRUE)
)
public_root <- normalizePath(args[[6L]], mustWork = TRUE)
launch_stdout <- normalizePath(args[[7L]], mustWork = TRUE)
launch_stderr <- normalizePath(args[[8L]], mustWork = TRUE)
output_dir <- normalizePath(args[[9L]], mustWork = FALSE)
if (dir.exists(output_dir)) stop("refusing to overwrite MV8-WX output", call. = FALSE)
runner_head <- tolower(trimws(Sys.getenv("MV08WX_RUNNER_HEAD", unset = "")))
if (!grepl("^[0-9a-f]{40}$", runner_head)) {
  stop("MV8-WX exact runner HEAD absent", call. = FALSE)
}
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required", call. = FALSE)
dir.create(output_dir, recursive = TRUE)
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
atomic_text <- function(value, path) {
  partial <- paste0(path, ".partial")
  writeLines(value, partial, useBytes = TRUE)
  if (!file.rename(partial, path)) stop("failed to publish ", basename(path), call. = FALSE)
}
verify_manifest <- function(root, name) {
  path <- file.path(root, name)
  rows <- read_csv(path)
  files <- file.path(root, rows$artifact)
  if (!all(file.exists(files)) ||
      !all(vapply(files, sha_file, character(1L)) == rows$sha256)) {
    stop("MV8-WX manifest drift: ", name, call. = FALSE)
  }
  data.frame(
    audit_stage = sub("-artifact-manifest\\.csv$", "", name),
    artifact_count = nrow(rows),
    aggregate_public_bytes = sum(as.numeric(rows$bytes)),
    manifest_sha256 = sha_file(path), stringsAsFactors = FALSE
  )
}

manifest_names <- c(
  mv08va = "mv08va-artifact-manifest.csv",
  mv08vb = "mv08vb-artifact-manifest.csv",
  mv08vc = "mv08vc-artifact-manifest.csv",
  mv08vd = "mv08vd-artifact-manifest.csv"
)
audit_summary <- do.call(rbind, lapply(names(recovery_roots), function(stage) {
  verify_manifest(recovery_roots[[stage]], manifest_names[[stage]])
}))
w_summary <- verify_manifest(w_root, "mv08w-artifact-manifest.csv")
w_summary$audit_stage <- "mv08w"
audit_summary <- rbind(audit_summary, w_summary)

va_decision <- read_csv(file.path(recovery_roots[["mv08va"]], "mv08va-decision.csv"))
vb_decision <- read_csv(file.path(recovery_roots[["mv08vb"]], "mv08vb-decision.csv"))
vc_decision <- read_csv(file.path(recovery_roots[["mv08vc"]], "mv08vc-decision.csv"))
vd_decision <- read_csv(file.path(recovery_roots[["mv08vd"]], "mv08vd-decision.csv"))
vd_binding <- read_csv(file.path(
  recovery_roots[["mv08vd"]], "mv08vd-implementation-bindings.csv"
))
w_decision <- read_csv(file.path(w_root, "mv08w-decision.csv"))
progress <- read_csv(file.path(public_root, "mv08v-progress.csv"))
ledger <- read_csv(file.path(public_root, "mv08v-resource-ledger.csv"))
metrics <- read_csv(file.path(public_root, "mv08v-selected-ph-metrics.csv"))
receipt <- read_csv(file.path(public_root, "mv08va-bootstrap-receipt.csv"))

if (nrow(progress) != 1L || nrow(receipt) != 1L || nrow(w_decision) != 1L ||
    nrow(va_decision) != 1L || nrow(vb_decision) != 1L ||
    nrow(vc_decision) != 1L || nrow(vd_decision) != 1L) {
  stop("MV8-WX singleton contract drift", call. = FALSE)
}
original_head <- tolower(receipt$original_execution_head)
bootstrap_head <- tolower(receipt$recovery_execution_head)
metric_heads <- tolower(metrics$execution_head)
ledger_heads <- tolower(ledger$execution_head)
downstream_columns <- c(
  "landscape_records", "comparison_records", "clustering_records",
  "fusion_records", "label_records", "biological_outcome_records"
)
decision_firewall_columns <- c(
  "landscape_groups_authorized", "comparison_strata_authorized",
  "clustering_jobs_authorized", "fusion_jobs_authorized",
  "label_jobs_authorized", "outcome_jobs_authorized"
)
recovery_decisions <- list(va_decision, vb_decision, vc_decision, vd_decision)
recovery_decision_names <- vapply(
  recovery_decisions, function(value) value$decision, character(1L)
)
expected_recovery_decisions <- c(
  "authorize_no_retry_job1_bootstrap_and_resume_at_job2",
  "authorize_corrected_zero_retry_MV8VA_bootstrap_after_commit",
  "authorize_hash_bound_zero_retry_MV8VA_bootstrap_after_commit",
  "authorize_amendment_bound_resume_at_job2"
)
recovery_firewalls <- vapply(recovery_decisions, function(value) {
  all(as.numeric(unlist(value[decision_firewall_columns], use.names = FALSE)) == 0) &&
    identical(value$outcome_label_state, "closed") &&
    !value$biological_outcomes_computed
}, logical(1L))
current_binding_hashes <- vapply(vd_binding$file, sha_file, character(1L))
stdout_lines <- trimws(readLines(launch_stdout, warn = FALSE))
stdout_lines <- stdout_lines[nzchar(stdout_lines)]
completion_line <- "MV8-V full PH production completed 1257/1257; MV8-W closure required"

head_summary <- data.frame(
  contract_id = "mv08wx_execution_head_summary_v1",
  role = c("original_job1_science", "bootstrap_byte_copy", "recovery_science"),
  execution_head = c(original_head, bootstrap_head, runner_head),
  ph_metric_rows = c(sum(metric_heads == original_head), 0L,
                     sum(metric_heads == runner_head)),
  resource_ledger_rows = c(sum(ledger_heads == original_head), 0L,
                           sum(ledger_heads == runner_head)),
  ph_recomputations = c(0L, receipt$recomputed_records, 0L),
  retries = c(0L, receipt$retry_records, progress$retries),
  stringsAsFactors = FALSE
)

validation <- data.frame(
  check_id = c(
    "recovery_manifests", "recovery_decisions", "mv08w_manifest",
    "mv08w_1280_closure",
    "terminal_progress", "production_cardinality", "strict_metric_order",
    "completed_attempt_cardinality", "original_head_job1_only",
    "recovery_head_remaining_metrics", "ledger_head_distribution",
    "bootstrap_head_distinct", "byte_identical_bootstrap", "h0_mst_bootstrap",
    "zero_recomputation", "zero_retry", "one_worker", "aggregate_cap",
    "private_cap", "runner_completion_line", "empty_launch_stderr",
    "final_implementations",
    "recovery_firewalls", "progress_firewalls", "mv08w_firewalls",
    "aggregate_only_schema"
  ),
  passed = c(
    nrow(audit_summary[audit_summary$audit_stage != "mv08w", , drop = FALSE]) == 4L,
    identical(recovery_decision_names, expected_recovery_decisions) &&
      vd_decision$accepted_completed_records == 1L &&
      vd_decision$remaining_ph_records_authorized == 1256L &&
      vd_decision$resume_at_production_order == 2L &&
      vd_decision$workers == 1L && vd_decision$automatic_retries == 0L,
    nrow(audit_summary[audit_summary$audit_stage == "mv08w", , drop = FALSE]) == 1L,
    w_decision$full_ph_records == 1280L &&
      w_decision$validations_passed == w_decision$validations_total,
    progress$state == "ph_production_complete_closure_pending" &&
      tolower(progress$execution_head) == runner_head,
    nrow(metrics) == 1257L && progress$completed_records == 1257L,
    identical(as.integer(metrics$production_order), 1:1257),
    sum(ledger$disposition == "completed") == 1257L,
    metric_heads[[1L]] == original_head && sum(metric_heads == original_head) == 1L,
    all(metric_heads[-1L] == runner_head),
    ledger_heads[[1L]] == original_head && all(ledger_heads[-1L] == runner_head),
    length(unique(c(original_head, bootstrap_head, runner_head))) == 3L,
    receipt$accepted_records == 1L && receipt$byte_identical_to_original,
    receipt$h0_mst_passed,
    receipt$recomputed_records == 0L,
    receipt$retry_records == 0L && progress$retries == 0L &&
      all(ledger$retries == 0L),
    progress$workers == 1L && all(ledger$workers == 1L),
    as.numeric(progress$aggregate_attempt_seconds) <= 72 * 3600 + 5,
    as.numeric(progress$private_bytes) <= 1024^3,
    length(stdout_lines) >= 1L && identical(tail(stdout_lines, 1L), completion_line),
    as.numeric(file.info(launch_stderr)$size) == 0,
    identical(unname(current_binding_hashes), vd_binding$sha256),
    all(recovery_firewalls),
    all(as.numeric(unlist(progress[downstream_columns], use.names = FALSE)) == 0) &&
      identical(progress$outcome_label_state, "closed") &&
      !progress$biological_outcomes_computed,
    all(as.numeric(unlist(
      w_decision[decision_firewall_columns], use.names = FALSE
    )) == 0) &&
      identical(w_decision$outcome_label_state, "closed") &&
      !w_decision$biological_outcomes_computed,
    !any(c(
      "unit_id", "job_id", "accession", "source_path", "private_path",
      "private_file"
    ) %in% c(names(audit_summary), names(head_summary)))
  ),
  evidence = c(
    "MV8-VA through MV8-VD manifests independently rehash",
    "all recovery decisions retain the exact admitted bootstrap/resume sequence",
    "MV8-W closure manifest independently rehashes",
    "MV8-W reports exact independently validated 1,280-record closure",
    "terminal production state is bound to the exact recovery runner head",
    "1,257 selected PH metrics and completed production records",
    "selected metrics retain production order 1 through 1,257",
    "exactly 1,257 primary-or-fallback attempts completed",
    "the original science head appears only for admitted job 1",
    "orders 2 through 1,257 bind to the recovery science head",
    "resource ledger contains admitted job 1 then recovery attempts only",
    "original science, bootstrap copy, and recovery science heads are distinct",
    "job 1 is a byte-identical copy of the independently admitted record",
    "job-1 H0 MST oracle passed before bootstrap admission",
    "bootstrap performed zero PH recomputations",
    "bootstrap, progress, and ledger record zero retries",
    "one worker retained throughout",
    "aggregate attempt time remains within 72 hours",
    "private production evidence remains within one GiB",
    "runner stdout ends with the exact 1,257/1,257 completion line",
    "recovery runner launch stderr is empty",
    "current bootstrap and recovery runner match MV8-VD bindings",
    "all four recovery decisions retain downstream and outcome firewalls",
    "terminal progress retains downstream and outcome firewalls",
    "MV8-W retains downstream and outcome firewalls",
    "published summaries contain no unit/private metadata columns"
  ), stringsAsFactors = FALSE
)
if (!all(validation$passed)) stop("MV8-WX recovery closure failed", call. = FALSE)

decision <- data.frame(
  contract_id = "mv08wx_decision_v1",
  decision = "recovery_provenance_bound_to_full_1280_PH_closure",
  recovery_audits = 4L, mv08w_records = w_decision$full_ph_records,
  production_records = nrow(metrics), accepted_bootstrap_records = 1L,
  recomputed_records = receipt$recomputed_records,
  retry_records = progress$retries, execution_heads = nrow(head_summary),
  validations_passed = sum(validation$passed),
  validations_total = nrow(validation),
  implementation_sha256 = sha_file("scripts/build_mv08wx_recovery_provenance_closure.R"),
  landscape_groups_authorized = 0L, comparison_strata_authorized = 0L,
  clustering_jobs_authorized = 0L, fusion_jobs_authorized = 0L,
  label_jobs_authorized = 0L, outcome_jobs_authorized = 0L,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
atomic_csv(audit_summary, file.path(output_dir, "mv08wx-audit-chain-summary.csv"))
atomic_csv(head_summary, file.path(output_dir, "mv08wx-execution-head-summary.csv"))
atomic_csv(validation, file.path(output_dir, "mv08wx-validation.csv"))
atomic_csv(decision, file.path(output_dir, "mv08wx-decision.csv"))
atomic_text(c(
  "# MV8-WX recovery-provenance closure", "",
  paste0("**Result:** ", sum(validation$passed), "/", nrow(validation),
         " checks pass; recovery provenance is bound to 1,280/1,280 PH closure."),
  "",
  paste0(
    "MV8-WX binds the independently validated MV8-W inventory to the complete ",
    "MV8-VA through MV8-VD recovery chain. Job 1 remains an independently ",
    "admitted byte copy with zero PH recomputation; production orders 2 through ",
    "1,257 bind to the exact committed recovery runner head."
  ), "",
  paste0(
    "Only aggregate audit and execution-head counts are published. No unit ",
    "identifier, accession, donor metadata, source path, or private filename is ",
    "included. Landscapes and all downstream analyses remain closed."
  )
), file.path(output_dir, "MV08WX_RECOVERY_PROVENANCE_CLOSURE.md"))
artifacts <- list.files(output_dir, full.names = TRUE)
artifacts <- artifacts[basename(artifacts) != "mv08wx-artifact-manifest.csv"]
manifest <- data.frame(
  artifact = basename(artifacts), bytes = as.numeric(file.info(artifacts)$size),
  sha256 = vapply(artifacts, sha_file, character(1L)), stringsAsFactors = FALSE
)
atomic_csv(manifest, file.path(output_dir, "mv08wx-artifact-manifest.csv"))
cat("MV8-WX checks=", sum(validation$passed), "/", nrow(validation),
    "; PH=", w_decision$full_ph_records, "/1280; retries=0\n", sep = "")
