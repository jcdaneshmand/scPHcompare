#!/usr/bin/env Rscript

# Freeze the MV8-ZF order-503 scientific failure before any Rust repair.
# This builder is deliberately unable to resume production or mutate the
# stopped production roots.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 10L) stop(paste(
  "usage: build_mv08zo_landscape_kernel_failure_prefreeze.R",
  "<mv08z-root> <mv08zf-root> <private-bindings> <mv08s-private>",
  "<mv08v-private> <failed-private> <failed-public> <rust-library>",
  "<output> <parent-head>"
), call. = FALSE)
if (!requireNamespace("digest", quietly = TRUE)) {
  stop("digest required", call. = FALSE)
}

z_root <- normalizePath(args[[1L]], mustWork = TRUE)
zf_root <- normalizePath(args[[2L]], mustWork = TRUE)
bindings_path <- normalizePath(args[[3L]], mustWork = TRUE)
s_root <- normalizePath(args[[4L]], mustWork = TRUE)
v_root <- normalizePath(args[[5L]], mustWork = TRUE)
failed_private <- normalizePath(args[[6L]], mustWork = TRUE)
failed_public <- normalizePath(args[[7L]], mustWork = TRUE)
rust_library <- normalizePath(args[[8L]], mustWork = TRUE)
output <- normalizePath(args[[9L]], mustWork = FALSE)
parent_head <- tolower(args[[10L]])
if (dir.exists(output) || !grepl("^[0-9a-f]{40}$", parent_head)) {
  stop("MV8-ZO requires fresh output and exact parent HEAD", call. = FALSE)
}

source("R/landscape_contract.R")
source("R/landscape_reference.R")
source("R/landscape_rust_prototype.R")
source("R/mv08z_landscape_production.R")

.mv08z_verify_manifest(z_root, "mv08z-artifact-manifest.csv")
.mv08z_verify_manifest(zf_root, "mv08zf-artifact-manifest.csv")
read_csv <- function(path) .mv08z_read_csv(path)
sha_file <- function(path) .mv08z_sha256_file(path)
truth <- function(value) .mv08z_truth(value)

progress_path <- file.path(failed_public, "mv08zf-progress.csv")
completion_path <- file.path(failed_public, "mv08zf-chunk-completions.csv")
ledger_path <- file.path(failed_public, "mv08zf-resource-ledger.csv")
progress <- read_csv(progress_path)
completion <- read_csv(completion_path)
ledger <- read_csv(ledger_path)
queue <- read_csv(file.path(zf_root, "mv08zf-production-queue.csv"))
groups <- read_csv(file.path(z_root, "mv08z-group-queue.csv"))
chunks <- read_csv(file.path(z_root, "mv08z-chunk-queue.csv"))
bindings <- read_csv(bindings_path)

failed_order <- 503L
failed_ledger <- ledger[as.integer(ledger$production_order) == failed_order, , drop = FALSE]
failed_queue <- queue[as.integer(queue$production_order) == failed_order, , drop = FALSE]
failed_chunk <- chunks[
  as.integer(chunks$group_order) == as.integer(failed_queue$group_order) &
    as.integer(chunks$chunk_order) == as.integer(failed_queue$chunk_order),
  , drop = FALSE
]
failed_group <- groups[
  as.integer(groups$group_order) == as.integer(failed_queue$group_order),
  , drop = FALSE
]
stderr_path <- file.path(failed_private, "logs", "chunk_0503.stderr")
stdout_path <- file.path(failed_private, "logs", "chunk_0503.stdout")
stderr_text <- paste(readLines(stderr_path, warn = FALSE), collapse = "\n")
failed_pair_ordinal <- as.integer(sub(
  ".*pair ordinal ([0-9]+).*", "\\1", stderr_text
))

group_bindings <- bindings[
  as.integer(bindings$group_order) == as.integer(failed_group$group_order),
  , drop = FALSE
]
pairs <- .mv08z_add_pair_identities(
  .mv08z_group_pairs(group_bindings), failed_group$group_id
)
failed_pair <- pairs[
  as.integer(pairs$pair_ordinal) == failed_pair_ordinal, , drop = FALSE
]

root_for <- function(role) {
  if (role == "mv08s_private_v3") return(s_root)
  if (role == "mv08v_recovery_private_v2") return(v_root)
  stop("MV8-ZO unknown PH source role", call. = FALSE)
}
intervals_for <- function(axis) {
  row <- group_bindings[
    as.integer(group_bindings$axis_order) == as.integer(axis), , drop = FALSE
  ]
  path <- file.path(root_for(row$source_role), row$output_file)
  if (nrow(row) != 1L || !file.exists(path) ||
      as.numeric(file.info(path)$size) != as.numeric(row$output_bytes) ||
      sha_file(path) != row$output_sha256) {
    stop("MV8-ZO private PH binding drift", call. = FALSE)
  }
  record <- readRDS(path)
  .canonical_rust_intervals(.mv08z_finite_intervals(
    record, failed_group$homology_dimension
  ))
}
first <- intervals_for(failed_pair$first_axis_order)
second <- intervals_for(failed_pair$second_axis_order)

if (!.load_scph_rust_prototype(rust_library)) {
  stop("MV8-ZO admitted Rust library unavailable", call. = FALSE)
}
rust_call <- function(first, second) {
  base::.C(
    .scph_rust_prototype_state$symbol,
    first_births = as.double(first[, 1L]),
    first_deaths = as.double(first[, 2L]),
    first_len = as.integer(nrow(first)),
    second_births = as.double(second[, 1L]),
    second_deaths = as.double(second[, 2L]),
    second_len = as.integer(nrow(second)),
    dimension = 1L,
    squared_distance = double(1L), active_levels = double(1L),
    event_segments = double(1L), first_finite_intervals = double(1L),
    second_finite_intervals = double(1L), engine_version = integer(1L),
    status = integer(1L), NAOK = TRUE
  )
}
empty <- matrix(numeric(), nrow = 0L, ncol = 2L,
                dimnames = list(NULL, c("birth", "death")))
failed_result <- rust_call(first, second)
first_norm <- rust_call(first, empty)
second_norm <- rust_call(second, empty)
canonical_norm <- function(intervals) {
  sum((intervals[, "death"] - intervals[, "birth"]) ^ 3 / 12)
}

synthetic <- matrix(c(
  1, 4,
  3, 4,
  0, 2,
  1, 2
), ncol = 2L, byrow = TRUE,
dimnames = list(NULL, c("birth", "death")))
synthetic_result <- rust_call(synthetic, empty)
synthetic_exact <- canonical_norm(synthetic)

# Axis order is private runtime structure. Publish only exposure counts.
prior_pairs <- pairs[
  as.integer(pairs$pair_ordinal) < failed_pair_ordinal &
    (as.integer(pairs$first_axis_order) == as.integer(failed_pair$second_axis_order) |
       as.integer(pairs$second_axis_order) == as.integer(failed_pair$second_axis_order)),
  , drop = FALSE
]
prior_partner_axes <- ifelse(
  as.integer(prior_pairs$first_axis_order) == as.integer(failed_pair$second_axis_order),
  prior_pairs$second_axis_order, prior_pairs$first_axis_order
)
prior_partner_depths <- vapply(prior_partner_axes, function(axis) {
  .mv08z_active_depth(intervals_for(axis))
}, integer(1L))

production_dir <- file.path(
  failed_private, "production",
  .mv08z_safe_group(failed_queue$group_order),
  .mv08z_safe_chunk(failed_queue$chunk_order)
)
private_partials <- list.files(
  failed_private, pattern = "[.]partial$|__partial__|[.]tmp$|[.]part$",
  recursive = TRUE, full.names = TRUE
)

failure <- data.frame(
  contract_id = "mv08zo_failure_evidence_v1",
  execution_head = progress$execution_head,
  stopped_production_order = failed_order,
  completed_chunks = nrow(completion),
  completed_pairs = sum(as.integer(completion$pair_count)),
  failed_group_order = as.integer(failed_queue$group_order),
  failed_chunk_order = as.integer(failed_queue$chunk_order),
  failed_pair_ordinal = failed_pair_ordinal,
  homology_dimension = failed_group$homology_dimension,
  exit_status = as.integer(failed_ledger$exit_status),
  elapsed_seconds = as.numeric(failed_ledger$elapsed_seconds),
  peak_process_tree_rss_bytes = as.numeric(failed_ledger$peak_process_tree_rss_bytes),
  private_partial_files = length(private_partials),
  failed_chunk_published = dir.exists(production_dir),
  retries = sum(as.integer(ledger$retries)),
  outcome_label_state = "closed",
  biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)

diagnostic <- data.frame(
  contract_id = "mv08zo_kernel_diagnostic_v1",
  diagnostic = c(
    "failed_pair", "failed_first_norm", "failed_second_norm",
    "public_synthetic_norm"
  ),
  rust_status = c(
    failed_result$status, first_norm$status, second_norm$status,
    synthetic_result$status
  ),
  rust_active_levels = c(
    failed_result$active_levels, first_norm$active_levels,
    second_norm$active_levels, synthetic_result$active_levels
  ),
  canonical_active_levels = c(
    max(.mv08z_active_depth(first), .mv08z_active_depth(second)),
    .mv08z_active_depth(first), .mv08z_active_depth(second),
    .mv08z_active_depth(synthetic)
  ),
  rust_squared_norm = c(
    NA_real_, first_norm$squared_distance, second_norm$squared_distance,
    synthetic_result$squared_distance
  ),
  canonical_squared_norm = c(
    NA_real_, canonical_norm(first), canonical_norm(second), synthetic_exact
  ),
  absolute_norm_error = c(
    NA_real_, abs(first_norm$squared_distance - canonical_norm(first)),
    abs(second_norm$squared_distance - canonical_norm(second)),
    abs(synthetic_result$squared_distance - synthetic_exact)
  ),
  relative_norm_error = c(
    NA_real_, abs(first_norm$squared_distance - canonical_norm(first)) /
      canonical_norm(first),
    abs(second_norm$squared_distance - canonical_norm(second)) /
      canonical_norm(second),
    abs(synthetic_result$squared_distance - synthetic_exact) / synthetic_exact
  ),
  interpretation = c(
    "success_status_but_active_level_mismatch",
    "exact_control",
    "active_level_and_exact_norm_failure",
    "exact_norm_failure_despite_level_count_match"
  ),
  stringsAsFactors = FALSE
)

exposure <- data.frame(
  contract_id = "mv08zo_prefix_exposure_v1",
  stopped_completed_chunks = nrow(completion),
  stopped_completed_pairs = sum(as.integer(completion$pair_count)),
  known_prior_pairs_using_failed_diagram = nrow(prior_pairs),
  known_prior_partners_masking_level_mismatch = sum(prior_partner_depths >=
    as.integer(second_norm$active_levels)),
  active_level_guard_sufficient = FALSE,
  stopped_prefix_scientifically_accepted = FALSE,
  stopped_prefix_preservation = "immutable_rejected_evidence",
  production_resume_authorized = FALSE,
  stringsAsFactors = FALSE
)

fixture <- data.frame(
  contract_id = "mv08zo_public_synthetic_fixture_v1",
  interval_order = seq_len(nrow(synthetic)),
  birth = synthetic[, "birth"], death = synthetic[, "death"],
  expected_active_levels = .mv08z_active_depth(synthetic),
  stringsAsFactors = FALSE
)

bound_paths <- c(
  progress_path, completion_path, ledger_path, stdout_path, stderr_path,
  bindings_path, rust_library, "rust/scph_landscape_kernel/src/lib.rs",
  "R/landscape_rust_prototype.R", "R/mv08z_landscape_production.R",
  "scripts/run_mv08z_landscape_chunk.R",
  "scripts/run_mv08zf_full_landscape_production.R",
  "scripts/build_mv08zo_landscape_kernel_failure_prefreeze.R"
)
bindings_out <- data.frame(
  contract_id = "mv08zo_evidence_binding_v1",
  role = c(
    "stopped_progress", "stopped_completion", "stopped_ledger",
    "failed_stdout", "failed_stderr", "private_unit_axis",
    "admitted_private_rust", "faulted_rust_source", "rust_R_shim",
    "landscape_helpers", "chunk_worker", "serial_runner", "prefreeze_builder"
  ),
  bytes = as.numeric(file.info(bound_paths)$size),
  sha256 = vapply(bound_paths, sha_file, character(1L)),
  public_content = c(rep(TRUE, 3L), rep(FALSE, 4L), rep(TRUE, 6L)),
  stringsAsFactors = FALSE
)

checks <- data.frame(
  check_id = c(
    "z_manifest", "zf_manifest", "one_progress", "clean_502_prefix",
    "failed_503_ledger", "failure_log_bound", "failure_not_resource_cap",
    "zero_retry", "zero_private_partials", "failed_chunk_unpublished",
    "failed_pair_in_chunk", "h1_failure", "rust_status_success",
    "failed_level_mismatch", "control_exact_norm", "failed_exact_norm_error",
    "synthetic_level_guard_passes", "synthetic_exact_norm_fails",
    "known_prefix_exposure", "prefix_not_accepted", "resume_closed",
    "candidate_only_authorized", "downstream_closed", "parent_bound"
  ),
  passed = c(
    TRUE, TRUE, nrow(progress) == 1L,
    nrow(completion) == 502L &&
      identical(as.integer(completion$production_order), 1:502),
    nrow(failed_ledger) == 1L && failed_ledger$disposition == "failed" &&
      failed_ledger$exit_status == 1L,
    grepl("Rust calculation failed closed at pair ordinal 1531", stderr_text,
          fixed = TRUE),
    failed_ledger$elapsed_seconds < failed_ledger$elapsed_cap_seconds &&
      failed_ledger$peak_process_tree_rss_bytes < failed_ledger$rss_cap_bytes,
    sum(as.integer(ledger$retries)) == 0L,
    length(private_partials) == 0L, !dir.exists(production_dir),
    nrow(failed_pair) == 1L && failed_pair_ordinal >= failed_chunk$pair_start &&
      failed_pair_ordinal <= failed_chunk$pair_end,
    failed_group$homology_dimension == "H1", failed_result$status == 0L,
    as.integer(failed_result$active_levels) !=
      max(.mv08z_active_depth(first), .mv08z_active_depth(second)),
    abs(first_norm$squared_distance - canonical_norm(first)) <= 1e-15,
    abs(second_norm$squared_distance - canonical_norm(second)) > 1e-12,
    as.integer(synthetic_result$active_levels) == .mv08z_active_depth(synthetic),
    abs(synthetic_result$squared_distance - synthetic_exact) > 1e-12,
    nrow(prior_pairs) == 13L && all(prior_partner_depths >=
      as.integer(second_norm$active_levels)),
    !exposure$stopped_prefix_scientifically_accepted,
    !exposure$production_resume_authorized, TRUE,
    all(c("comparison_jobs", "clustering_jobs", "fusion_jobs", "label_jobs",
          "outcome_jobs", "adoption_jobs", "manuscript_claim_jobs") %in%
        names(progress)) &&
      all(unlist(progress[c("comparison_jobs", "clustering_jobs",
        "fusion_jobs", "label_jobs", "outcome_jobs", "adoption_jobs",
        "manuscript_claim_jobs")]) == 0) &&
      progress$outcome_label_state == "closed" &&
      !truth(progress$biological_outcomes_computed),
    nzchar(parent_head)
  ),
  stringsAsFactors = FALSE
)
if (!all(checks$passed)) {
  stop("MV8-ZO prefreeze failed: ",
       paste(checks$check_id[!checks$passed], collapse = ", "), call. = FALSE)
}

decision <- data.frame(
  contract_id = "mv08zo_recovery_decision_v1",
  decision = "preserve_and_reject_stopped_prefix_then_repair_privately",
  stopped_prefix_chunks_preserved = 502L,
  stopped_prefix_scientifically_accepted = FALSE,
  old_root_resume_authorized = FALSE,
  old_root_delete_authorized = FALSE,
  candidate_rust_source_change_authorized = TRUE,
  private_candidate_builds_authorized = 1L,
  public_synthetic_fixture_runs_authorized = 1L,
  private_failed_pair_runs_authorized = 1L,
  canonical_norm_checks_required = TRUE,
  prior_output_reuse_authorized = FALSE,
  fresh_production_authorized = FALSE,
  comparison_jobs = 0L, clustering_jobs = 0L, fusion_jobs = 0L,
  label_jobs = 0L, outcome_jobs = 0L, manuscript_claim_jobs = 0L,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  next_gate = "MV8_ZP_kernel_repair_admission",
  stringsAsFactors = FALSE
)

dir.create(output, recursive = TRUE)
.mv08z_atomic_csv(failure, file.path(output, "mv08zo-failure-evidence.csv"))
.mv08z_atomic_csv(diagnostic, file.path(output, "mv08zo-kernel-diagnostic.csv"))
.mv08z_atomic_csv(exposure, file.path(output, "mv08zo-prefix-exposure.csv"))
.mv08z_atomic_csv(fixture, file.path(output, "mv08zo-public-synthetic-fixture.csv"))
.mv08z_atomic_csv(bindings_out, file.path(output, "mv08zo-evidence-bindings.csv"))
.mv08z_atomic_csv(checks, file.path(output, "mv08zo-validation.csv"))
.mv08z_atomic_csv(decision, file.path(output, "mv08zo-decision.csv"))
writeLines(c(
  "# MV8-ZO landscape-kernel failure recovery prefreeze", "",
  "**Result:** 24/24 checks pass; the stopped 502-chunk prefix is preserved but scientifically rejected.",
  "",
  "- Order 503 failed on a valid H1 pair without reaching time or memory caps.",
  "- Rust returned status 0 but produced 656 levels for canonical depth 655.",
  "- The affected diagram violates the exact landscape-norm identity by about 2.09e-4 relatively.",
  "- A four-interval public fixture violates the same exact identity even though its level count matches, proving that the old guard is insufficient.",
  "- Thirteen completed pairs used the affected diagram behind higher-depth partners; the entire stopped prefix is therefore evidence, not accepted science.",
  "- MV8-ZO authorizes one private candidate repair and bounded diagnostics only. It does not authorize resuming, overwriting, deleting, comparing, clustering, labeling, or drawing biological claims."
), file.path(output, "MV08ZO_LANDSCAPE_KERNEL_FAILURE_PREFREEZE.md"))
artifacts <- list.files(output, full.names = TRUE)
manifest <- data.frame(
  artifact = basename(artifacts),
  bytes = as.numeric(file.info(artifacts)$size),
  sha256 = vapply(artifacts, sha_file, character(1L)),
  stringsAsFactors = FALSE
)
.mv08z_atomic_csv(manifest, file.path(output, "mv08zo-artifact-manifest.csv"))
cat("MV8-ZO failure recovery prefreeze passed 24/24; old_resume=FALSE; next=MV8-ZP\n")
