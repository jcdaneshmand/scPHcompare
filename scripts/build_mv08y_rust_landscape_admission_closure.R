#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8L) stop(paste(
  "usage: build_mv08y_rust_landscape_admission_closure.R <prefreeze-root>",
  "<amendment-root> <toolchain-root> <build-evidence> <oracle-run-a> <oracle-run-b>",
  "<private-selection.csv> <output-dir>"
), call. = FALSE)

prefreeze_root <- normalizePath(args[[1L]], mustWork = TRUE)
amendment_root <- normalizePath(args[[2L]], mustWork = TRUE)
toolchain_root <- normalizePath(args[[3L]], mustWork = TRUE)
build_root <- normalizePath(args[[4L]], mustWork = TRUE)
run_a_root <- normalizePath(args[[5L]], mustWork = TRUE)
run_b_root <- normalizePath(args[[6L]], mustWork = TRUE)
selection_path <- normalizePath(args[[7L]], mustWork = TRUE)
output_dir <- normalizePath(args[[8L]], mustWork = FALSE)
if (dir.exists(output_dir)) stop("refusing to overwrite MV8-Y output", call. = FALSE)
oracle_execution_head <- tolower(trimws(
  Sys.getenv("MV08Y_EXECUTION_HEAD", unset = "")
))
if (!grepl("^[0-9a-f]{40}$", oracle_execution_head)) {
  stop("MV8-Y exact execution HEAD absent", call. = FALSE)
}
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required", call. = FALSE)

read_csv <- function(path) utils::read.csv(
  path, check.names = FALSE, stringsAsFactors = FALSE
)
sha_file <- function(path) digest::digest(
  file = path, algo = "sha256", serialize = FALSE
)
truth <- function(value) {
  if (is.logical(value)) return(!is.na(value) & value)
  tolower(trimws(as.character(value))) %in% c("true", "t", "1", "yes")
}
atomic_csv <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  partial <- paste0(path, ".partial")
  utils::write.csv(value, partial, row.names = FALSE, quote = TRUE, na = "")
  if (!file.rename(partial, path)) stop("failed to publish ", basename(path), call. = FALSE)
}
atomic_text <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  partial <- paste0(path, ".partial")
  writeLines(value, partial, useBytes = TRUE)
  if (!file.rename(partial, path)) stop("failed to publish ", basename(path), call. = FALSE)
}
verify_manifest <- function(root, name) {
  path <- file.path(root, name)
  manifest <- read_csv(path)
  files <- file.path(root, manifest$artifact)
  if (!all(file.exists(files)) ||
      !all(vapply(files, sha_file, character(1L)) == manifest$sha256) ||
      !all(as.numeric(file.info(files)$size) == as.numeric(manifest$bytes))) {
    stop("MV8-Y manifest drift: ", name, call. = FALSE)
  }
  manifest
}
parse_elapsed <- function(path) {
  lines <- readLines(path, warn = FALSE)
  line <- grep("Elapsed (wall clock) time", lines, value = TRUE, fixed = TRUE)
  if (length(line) != 1L) return(NA_real_)
  value <- tail(strsplit(line, "): ", fixed = TRUE)[[1L]], 1L)
  parts <- as.numeric(strsplit(value, ":", fixed = TRUE)[[1L]])
  if (anyNA(parts)) return(NA_real_)
  if (length(parts) == 2L) parts[[1L]] * 60 + parts[[2L]] else
    if (length(parts) == 3L) parts[[1L]] * 3600 + parts[[2L]] * 60 + parts[[3L]] else
      NA_real_
}
parse_rss <- function(path) {
  lines <- readLines(path, warn = FALSE)
  line <- grep("Maximum resident set size (kbytes):", lines,
               value = TRUE, fixed = TRUE)
  if (length(line) != 1L) return(NA_real_)
  as.numeric(sub("^.*:[[:space:]]*", "", line)) * 1024
}
tree_bytes <- function(root) {
  files <- list.files(root, recursive = TRUE, full.names = TRUE,
                      all.files = TRUE, no.. = TRUE)
  files <- files[!file.info(files)$isdir]
  sum(as.numeric(file.info(files)$size))
}

prefreeze_manifest <- verify_manifest(
  prefreeze_root, "mv08x-artifact-manifest.csv"
)
amendment_manifest <- verify_manifest(
  amendment_root, "mv08xa-artifact-manifest.csv"
)
build_manifest <- verify_manifest(build_root, "artifact-manifest.csv")
run_a_manifest <- verify_manifest(run_a_root, "artifact-manifest.csv")
run_b_manifest <- verify_manifest(run_b_root, "artifact-manifest.csv")

contract <- read_csv(file.path(prefreeze_root, "mv08x-contract.csv"))
selection_public <- read_csv(file.path(prefreeze_root, "mv08x-oracle-selection.csv"))
resource_policy <- read_csv(file.path(prefreeze_root, "mv08x-resource-policy.csv"))
acceptance <- read_csv(file.path(prefreeze_root, "mv08x-acceptance.csv"))
implementation <- read_csv(file.path(prefreeze_root, "mv08x-implementation-bindings.csv"))
input_manifest <- read_csv(file.path(prefreeze_root, "mv08x-input-manifest.csv"))
prefreeze_validation <- read_csv(file.path(prefreeze_root, "mv08x-validation.csv"))
prefreeze_decision <- read_csv(file.path(prefreeze_root, "mv08x-decision.csv"))
amendment_bindings <- read_csv(file.path(
  amendment_root, "mv08xa-amendment-bindings.csv"
))
failure <- read_csv(file.path(amendment_root, "mv08xa-failure.csv"))
diagnostic <- read_csv(file.path(amendment_root, "mv08xa-diagnostic.csv"))
amendment_decision <- read_csv(file.path(amendment_root, "mv08xa-decision.csv"))
build <- read_csv(file.path(build_root, "build-validation.csv"))
source_bindings <- read_csv(file.path(build_root, "source-bindings.csv"))
results_a <- read_csv(file.path(run_a_root, "oracle-results.csv"))
results_b <- read_csv(file.path(run_b_root, "oracle-results.csv"))
fixtures_a <- read_csv(file.path(run_a_root, "fixture-results.csv"))
fixtures_b <- read_csv(file.path(run_b_root, "fixture-results.csv"))
resource_a <- read_csv(file.path(run_a_root, "resource.csv"))
resource_b <- read_csv(file.path(run_b_root, "resource.csv"))
private_selection <- read_csv(selection_path)

if (nrow(contract) != 1L || nrow(prefreeze_decision) != 1L ||
    nrow(build) != 1L || nrow(resource_a) != 1L || nrow(resource_b) != 1L ||
    nrow(selection_public) != 28L || nrow(private_selection) != 28L ||
    nrow(results_a) != 28L || nrow(results_b) != 28L ||
    nrow(fixtures_a) != 9L || nrow(fixtures_b) != 9L) {
  stop("MV8-Y singleton/cardinality drift", call. = FALSE)
}

has_old <- !is.na(amendment_bindings$old_sha256) &
  nzchar(amendment_bindings$old_sha256)
amended_roles <- amendment_bindings$role[has_old]
unamended <- !implementation$role %in% amended_roles
current_hashes <- vapply(implementation$file[unamended], sha_file, character(1L))
implementation_match <- all(current_hashes == implementation$sha256[unamended])
amendment_current_hashes <- vapply(
  amendment_bindings$file, sha_file, character(1L)
)
amendment_old <- implementation[
  match(amendment_bindings$role[has_old], implementation$role), , drop = FALSE
]
amendment_binding_ok <-
  !anyNA(amendment_old$role) &&
  all(amendment_old$bytes == amendment_bindings$old_bytes[has_old]) &&
  all(amendment_old$sha256 == amendment_bindings$old_sha256[has_old]) &&
  all(as.numeric(file.info(amendment_bindings$file)$size) ==
        as.numeric(amendment_bindings$new_bytes)) &&
  all(amendment_current_hashes == amendment_bindings$new_sha256)
failure_preserved <- nrow(failure) == 1L && nrow(diagnostic) == 1L &&
  nrow(amendment_decision) == 1L &&
  failure$rebuild_execution_head == build$execution_head &&
  !truth(failure$output_root_created) &&
  failure$terminal_error == "MV8-X canonical R oracle gate failed" &&
  diagnostic$pairs == 28L && diagnostic$engine_passes == 28L &&
  diagnostic$reverse_bit_passes == 28L &&
  diagnostic$reverse_count_passes == 28L &&
  diagnostic$reverse_diagnostic_passes == 28L &&
  diagnostic$first_self_zero_passes == 28L &&
  diagnostic$second_self_zero_passes == 28L &&
  diagnostic$all_active_level_passes == 14L &&
  diagnostic$production_landscape_jobs == 0L &&
  truth(amendment_decision$replacement_authorized) &&
  amendment_decision$replacement_attempts == 1L
selection_binding <- input_manifest[
  input_manifest$role == "private_oracle_selection", , drop = FALSE
]
installer_path <- file.path(toolchain_root, "rustup-init")
receipt_path <- file.path(toolchain_root, "toolchain-acquisition-receipt.txt")
installer_ok <- file.exists(installer_path) &&
  sha_file(installer_path) == contract$rustup_installer_sha256
receipt <- if (file.exists(receipt_path)) readLines(receipt_path, warn = FALSE) else character()
receipt_ok <- all(c(
  paste0("installer_sha256=", contract$rustup_installer_sha256),
  "rustc_release=1.97.1", "rustc_host=x86_64-unknown-linux-gnu",
  "profile=minimal_plus_rustfmt_clippy", "home_mutation=false",
  "system_package_mutation=false"
) %in% receipt)

source_roles <- c("cargo_manifest", "cargo_lock", "rust_kernel", "c_abi_header")
expected_source <- implementation[
  match(source_roles, implementation$role), c("role", "bytes", "sha256")
]
source_binding_ok <- identical(source_bindings$role, source_roles) &&
  all(as.numeric(source_bindings$bytes) == as.numeric(expected_source$bytes)) &&
  all(source_bindings$sha256 == expected_source$sha256)

scientific_a_sha <- sha_file(file.path(run_a_root, "oracle-results.csv"))
scientific_b_sha <- sha_file(file.path(run_b_root, "oracle-results.csv"))
fixture_a_sha <- sha_file(file.path(run_a_root, "fixture-results.csv"))
fixture_b_sha <- sha_file(file.path(run_b_root, "fixture-results.csv"))
selection_public_match <- identical(results_a$oracle_id, selection_public$oracle_id) &&
  identical(results_a$pair_identity_sha256, selection_public$pair_identity_sha256)
selection_private_match <- identical(
  selection_public$oracle_id, private_selection$oracle_id
) && nrow(selection_binding) == 1L &&
  sha_file(selection_path) == selection_binding$sha256 &&
  as.numeric(file.info(selection_path)$size) == as.numeric(selection_binding$bytes)

build_times <- c(
  parse_elapsed(file.path(build_root, "build-a-time.txt")),
  parse_elapsed(file.path(build_root, "build-b-time.txt"))
)
build_rss <- c(
  parse_rss(file.path(build_root, "build-a-time.txt")),
  parse_rss(file.path(build_root, "build-b-time.txt"))
)
build_cap <- resource_policy[resource_policy$stage %in% c(
  "clean_build_a", "clean_build_b"
), , drop = FALSE]
oracle_cap <- resource_policy[resource_policy$stage %in% c(
  "oracle_run_a", "oracle_run_b"
), , drop = FALSE]
build_resource_pass <- length(build_times) == 2L && all(is.finite(build_times)) &&
  all(is.finite(build_rss)) && all(build_times <= build_cap$elapsed_cap_seconds) &&
  all(build_rss <= build_cap$rss_cap_bytes)
oracle_resource_pass <- all(c(resource_a$total_seconds, resource_b$total_seconds) <=
                              oracle_cap$elapsed_cap_seconds) &&
  all(c(resource_a$peak_process_rss_bytes, resource_b$peak_process_rss_bytes) <=
        oracle_cap$rss_cap_bytes) &&
  all(c(resource_a$workers, resource_b$workers) == 1L) &&
  all(c(resource_a$retries, resource_b$retries) == 0L)
private_bytes <- tree_bytes(toolchain_root) + sum(as.numeric(build_manifest$bytes)) +
  sum(as.numeric(run_a_manifest$bytes)) + sum(as.numeric(run_b_manifest$bytes)) +
  as.numeric(file.info(selection_path)$size)
private_cap <- 4 * 1024^3

dependency_lines <- readLines(
  file.path(build_root, "native-dependencies.txt"), warn = FALSE
)
dependency_ok <- any(grepl("libc\\.so", dependency_lines)) &&
  !any(grepl("libR|libpython|libcurl|libssl", dependency_lines, ignore.case = TRUE))

firewall_columns <- c(
  "production_landscape_jobs", "comparison_jobs", "clustering_jobs",
  "fusion_jobs", "label_jobs", "outcome_jobs"
)
resource_firewalls <- all(as.numeric(unlist(
  resource_a[firewall_columns], use.names = FALSE
)) == 0) && all(as.numeric(unlist(
  resource_b[firewall_columns], use.names = FALSE
)) == 0)
result_firewalls <- all(results_a$outcome_label_state == "closed") &&
  !any(truth(results_a$biological_outcomes_computed))

validation <- data.frame(
  check_id = c(
    "prefreeze_manifest", "amendment_manifest", "failed_attempt_preserved",
    "amendment_bindings", "prefreeze_checks", "implementation_hashes",
    "private_selection_binding", "installer_hash", "toolchain_receipt",
    "source_binding", "execution_head", "rust_identity", "one_build_job",
    "zero_external_crates", "format_unit_clippy", "two_clean_builds",
    "prior_binary_comparison_recorded", "native_c_abi", "dependencies",
    "build_resource_caps", "private_storage_cap", "oracle_cardinality",
    "all_group_coverage", "dimension_coverage", "scope_panel_representation",
    "reference_certificates", "reverse_invariants", "self_zero_invariants",
    "all_active_levels", "two_run_scientific_identity", "fixture_identity",
    "analytical_and_fallback_fixtures", "oracle_resource_caps",
    "aggregate_public_schema", "production_landscapes_closed",
    "downstream_firewalls", "labels_outcomes_closed"
  ),
  passed = c(
    nrow(prefreeze_manifest) == 9L, nrow(amendment_manifest) == 6L,
    failure_preserved, amendment_binding_ok,
    nrow(prefreeze_validation) == 20L && all(truth(prefreeze_validation$passed)),
    implementation_match, selection_private_match, installer_ok, receipt_ok,
    source_binding_ok, build$execution_head == failure$rebuild_execution_head &&
      resource_a$execution_head == oracle_execution_head &&
      resource_b$execution_head == oracle_execution_head,
    build$rustc_release == "1.97.1" &&
      build$rustc_host == "x86_64-unknown-linux-gnu",
    build$build_jobs == 1L, build$external_crates == 0L,
    truth(build$format_passed) && truth(build$unit_tests_passed) &&
      truth(build$clippy_passed),
    truth(build$clean_builds_byte_identical) &&
      build$build_a_sha256 == build$build_b_sha256 &&
      build$build_b_sha256 == build$candidate_sha256,
    nchar(build$prior_accepted_sha256) == 64L &&
      !is.na(build$matches_prior_accepted_binary),
    truth(build$native_c_abi_passed), dependency_ok, build_resource_pass,
    private_bytes <= private_cap, nrow(results_a) == 28L && all(truth(results_a$passed)),
    selection_public_match && length(unique(results_a$oracle_id)) == 28L,
    identical(sort(unique(results_a$homology_dimension)), c("H0", "H1")) &&
      all(table(results_a$homology_dimension) == 14L),
    length(unique(results_a$dataset_scope)) == 2L &&
      length(unique(results_a$panel_id)) == 2L &&
      length(unique(results_a$representation_id)) == 3L,
    all(results_a$absolute_error <= results_a$acceptance_threshold),
    all(truth(results_a$reverse_bit_identical)) &&
      all(truth(results_a$reverse_counts_swap)) &&
      all(truth(results_a$reverse_diagnostics_match)),
    all(truth(results_a$first_self_exact_zero)) &&
      all(truth(results_a$second_self_exact_zero)),
    all(truth(results_a$all_active_levels)) &&
      all(results_a$active_levels == results_a$expected_active_levels) &&
      all(results_a$expected_active_levels == pmax(
        results_a$expected_first_active_levels,
        results_a$expected_second_active_levels
      )),
    scientific_a_sha == scientific_b_sha,
    fixture_a_sha == fixture_b_sha,
    all(truth(fixtures_a$passed)) && nrow(fixtures_a) == 9L,
    oracle_resource_pass,
    !any(c("job_id", "unit_id", "output_file", "source_role", "accession",
           "private_path") %in% names(results_a)),
    all(c(resource_a$production_landscape_jobs,
          resource_b$production_landscape_jobs) == 0L),
    resource_firewalls,
    result_firewalls && all(c(resource_a$outcome_label_state,
                              resource_b$outcome_label_state) == "closed") &&
      !any(truth(c(resource_a$biological_outcomes_computed,
                   resource_b$biological_outcomes_computed)))
  ),
  evidence = c(
    "MV8-X manifest rehashes", "MV8-XA manifest rehashes",
    "failed opaque run and aggregate diagnosis are preserved",
    "old/new observability amendment bytes rehash",
    "20/20 prospective checks remain accepted",
    "all prefrozen implementation bytes rehash", "private locator selection rehashes",
    "official rustup installer rehashes", "isolated toolchain receipt is exact",
    "build source equals prefreeze", "build binds exact committed prefreeze head",
    "Rust 1.97.1 native Linux host", "one Cargo build job", "dependency-free crate",
    "format, four unit tests, and strict Clippy pass", "two clean libraries are identical",
    "comparison with prior accepted binary is explicit", "native C ABI harness passes",
    "native dependency inventory is bounded", "both clean builds stay under caps",
    "all ignored evidence stays below four GiB", "28/28 oracle pairs pass",
    "one immutable pair per future group", "14 H0 and 14 H1",
    "both scopes/panels and all representations", "all Rust errors fit R certificates",
    "every reverse is invariant", "both self distances are exact zero",
    "every nonzero active landscape level is represented", "A/B oracle tables are byte-identical",
    "A/B fixture tables are byte-identical", "9/9 analytical/fallback fixtures pass",
    "both oracle runs stay under caps", "published result schema is aggregate-only",
    "production landscape jobs remain zero", "all downstream job counts remain zero",
    "labels and biological outcomes remain closed"
  ), stringsAsFactors = FALSE
)
if (!all(validation$passed)) {
  stop("MV8-Y Rust admission closure failed at: ",
       paste(validation$check_id[!validation$passed], collapse = ", "), call. = FALSE)
}

build_summary <- data.frame(
  contract_id = "mv08y_build_summary_v1",
  execution_head = build$execution_head, rustc_release = build$rustc_release,
  oracle_execution_head = oracle_execution_head,
  rustc_host = build$rustc_host, build_jobs = build$build_jobs,
  external_crates = build$external_crates,
  candidate_sha256 = build$candidate_sha256,
  candidate_bytes = build$candidate_bytes,
  prior_accepted_sha256 = build$prior_accepted_sha256,
  matches_prior_accepted_binary = truth(build$matches_prior_accepted_binary),
  clean_builds_byte_identical = truth(build$clean_builds_byte_identical),
  native_c_abi_passed = truth(build$native_c_abi_passed),
  dependency_inventory_sha256 = sha_file(file.path(build_root, "native-dependencies.txt")),
  published_release = FALSE, public_default_changed = FALSE,
  stringsAsFactors = FALSE
)
oracle_summary <- results_a
oracle_summary$contract_id <- "mv08y_oracle_summary_v1"
resource_summary <- data.frame(
  contract_id = "mv08y_resource_summary_v1",
  stage = c("clean_build_a", "clean_build_b", "oracle_run_a", "oracle_run_b"),
  elapsed_seconds = c(build_times, resource_a$total_seconds, resource_b$total_seconds),
  peak_process_rss_bytes = c(build_rss, resource_a$peak_process_rss_bytes,
                             resource_b$peak_process_rss_bytes),
  elapsed_cap_seconds = c(build_cap$elapsed_cap_seconds,
                          oracle_cap$elapsed_cap_seconds),
  rss_cap_bytes = c(build_cap$rss_cap_bytes, oracle_cap$rss_cap_bytes),
  workers = 1L, retries = 0L, cap_passed = TRUE,
  stringsAsFactors = FALSE
)
fixture_summary <- fixtures_a
fixture_summary$contract_id <- "mv08y_fixture_summary_v1"
decision <- data.frame(
  contract_id = "mv08y_decision_v1",
  decision = "admit_hash_bound_Rust_kernel_for_separate_landscape_execution_prefreeze",
  full_ph_records = 1280L, landscape_groups_covered = 28L,
  oracle_pairs_passed = sum(truth(results_a$passed)), oracle_pairs_total = nrow(results_a),
  fixtures_passed = sum(truth(fixtures_a$passed)), fixtures_total = nrow(fixtures_a),
  candidate_sha256 = build$candidate_sha256,
  candidate_bytes = build$candidate_bytes,
  validations_passed = sum(validation$passed), validations_total = nrow(validation),
  private_wsl_candidate_admitted = TRUE, canonical_r_oracle_retained = TRUE,
  grouped_persim_fallback_retained = TRUE, public_default_changed = FALSE,
  binary_published = FALSE, landscape_execution_authorized = FALSE,
  production_landscape_jobs = 0L, comparison_jobs_authorized = 0L,
  clustering_jobs_authorized = 0L, fusion_jobs_authorized = 0L,
  label_jobs_authorized = 0L, outcome_jobs_authorized = 0L,
  next_gate = "separate_MV8_landscape_execution_prefreeze",
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)

dir.create(output_dir, recursive = TRUE)
atomic_csv(build_summary, file.path(output_dir, "mv08y-build-summary.csv"))
atomic_csv(oracle_summary, file.path(output_dir, "mv08y-oracle-summary.csv"))
atomic_csv(fixture_summary, file.path(output_dir, "mv08y-fixture-summary.csv"))
atomic_csv(resource_summary, file.path(output_dir, "mv08y-resource-summary.csv"))
atomic_csv(validation, file.path(output_dir, "mv08y-validation.csv"))
atomic_csv(decision, file.path(output_dir, "mv08y-decision.csv"))
atomic_text(c(
  "# MV8-Y Rust landscape-kernel admission closure", "",
  paste0("**Result:** ", sum(validation$passed), "/", nrow(validation),
         " checks pass; 28/28 MV8 group oracles and 9/9 fixtures pass."), "",
  paste0(
    "The dependency-free Rust kernel was rebuilt twice with the isolated ",
    "official Rust 1.97.1 Linux toolchain. Clean candidate libraries are ",
    "byte-identical and bound by SHA-256; the native C ABI, strict Rust checks, ",
    "failure fallback, determinism, and resource gates pass."
  ), "",
  paste0(
    "One prospectively selected stress pair from every future MV8 landscape ",
    "group agrees with the canonical R exact or error-certified oracle. H0 and ",
    "H1 remain separate, all active levels are retained, and no grid or level ",
    "cap is introduced."
  ), "",
  paste0(
    "Admission is limited to an explicit hash-verified private WSL candidate. ",
    "No binary is published or made default. Production landscapes and every ",
    "comparison, clustering, fusion, label, outcome, adoption, and manuscript ",
    "claim remain closed pending a separate execution prefreeze."
  )
), file.path(output_dir, "MV08Y_RUST_LANDSCAPE_ADMISSION_CLOSURE.md"))
artifacts <- list.files(output_dir, full.names = TRUE)
manifest <- data.frame(
  artifact = basename(artifacts), bytes = as.numeric(file.info(artifacts)$size),
  sha256 = vapply(artifacts, sha_file, character(1L)), stringsAsFactors = FALSE
)
atomic_csv(manifest, file.path(output_dir, "mv08y-artifact-manifest.csv"))
cat("MV8-Y checks=", sum(validation$passed), "/", nrow(validation),
    "; oracles=", sum(truth(results_a$passed)), "/28; fixtures=",
    sum(truth(fixtures_a$passed)), "/9; production_landscapes=0\n", sep = "")
