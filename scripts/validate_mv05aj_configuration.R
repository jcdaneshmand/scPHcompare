#!/usr/bin/env Rscript

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1")
options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 11L) {
  stop(paste(
    "usage: validate_mv05aj_configuration.R QUEUE INVENTORY RESULT_ROOT",
    "RESOURCE_CSV REPEAT_ROOT RESUME_BEFORE RESUME_AFTER RESUME_RESOURCE",
    "OUTPUT_DIR EXECUTION_HEAD PYTHON_EXECUTABLE"
  ), call. = FALSE)
}
for (package in c("digest")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(package, " is required for MV5-AJ validation.", call. = FALSE)
  }
}
source("R/provenance_utils.R")
source("R/mv05v_robustness_prefreeze.R")
source("R/mv05aj_nested_execution.R")

read_csv <- function(path) utils::read.csv(
  path, stringsAsFactors = FALSE, check.names = FALSE
)
file_sha <- function(path) digest::digest(
  file = path, algo = "sha256", serialize = FALSE
)
safe_name <- function(value) gsub("[^A-Za-z0-9_.-]", "_", value)
strict <- function(values) tolower(as.character(values)) == "true"
close_numeric <- function(first, second, tolerance = 1e-10) {
  length(first) == length(second) &&
    all(is.finite(first)) && all(is.finite(second)) &&
    all(abs(first - second) <= tolerance * pmax(1, abs(first), abs(second)))
}
closed_frame <- function(value) {
  state <- grep("outcome_label_state$", names(value), value = TRUE)
  flags <- grep("(biological_)?outcomes_computed$|rankings_computed$",
                names(value), value = TRUE)
  !any(c("tissue", "approach") %in% names(value)) &&
    all(vapply(state, function(column) {
      all(as.character(value[[column]]) == "closed")
    }, logical(1L))) &&
    all(vapply(flags, function(column) {
      !any(strict(value[[column]]))
    }, logical(1L)))
}
write_once <- function(value, path) {
  if (file.exists(path)) stop("Refusing to overwrite: ", path, call. = FALSE)
  write_provenance_csv(value, path)
}

queue_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
inventory_path <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
result_root <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
resource_path <- normalizePath(args[[4L]], winslash = "/", mustWork = TRUE)
repeat_root <- normalizePath(args[[5L]], winslash = "/", mustWork = TRUE)
before_path <- normalizePath(args[[6L]], winslash = "/", mustWork = TRUE)
after_path <- normalizePath(args[[7L]], winslash = "/", mustWork = TRUE)
resume_resource_path <- normalizePath(args[[8L]], winslash = "/",
                                      mustWork = TRUE)
output_dir <- args[[9L]]
execution_head <- args[[10L]]
python <- normalizePath(args[[11L]], winslash = "/", mustWork = TRUE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

queue <- read_csv(queue_path)
mv05aj_validate_nested_queue_v1(queue)
authorized <- queue[strict(queue$execution_authorized), , drop = FALSE]
authorized <- authorized[order(authorized$configuration_execution_order),
                         , drop = FALSE]
inventory <- read_csv(inventory_path)
resources <- read_csv(resource_path)
before <- read_csv(before_path)
after <- read_csv(after_path)
resume_resources <- read_csv(resume_resource_path)
current_head <- trimws(system2("git", c("rev-parse", "HEAD"), stdout = TRUE))
queue_sha <- file_sha(queue_path)
python_sha <- file_sha(python)

if (!identical(current_head, execution_head) || nrow(authorized) != 150L ||
    nrow(inventory) != 150L || !closed_frame(inventory) ||
    any(strict(inventory$admission_executed)) ||
    nrow(resources) != 150L ||
    !identical(resources$robustness_group_id, authorized$robustness_group_id) ||
    any(resources$disposition != "completed") ||
    any(resources$exit_status != 0L) || !closed_frame(resources) ||
    any(resources$maximum_workers != 1L) ||
    any(resources$elapsed_seconds > 600) ||
    sum(resources$elapsed_seconds) > 6 * 60 * 60 ||
    any(resources$peak_process_tree_rss_bytes > 4 * 1024^3) ||
    max(resources$cumulative_private_bytes) > 4 * 1024^3 ||
    nrow(before) != 1650L || nrow(after) != 1650L ||
    !identical(before, after) || !closed_frame(before) ||
    nrow(resume_resources) != 150L ||
    !identical(resume_resources$robustness_group_id,
               authorized$robustness_group_id) ||
    any(resume_resources$disposition != "reused_validated") ||
    any(resume_resources$exit_status != 0L) ||
    !closed_frame(resume_resources)) {
  stop("MV5-AJ aggregate execution/resume contract failed.", call. = FALSE)
}

method_ids <- c(
  "cell_landscape_h0_v1", "cell_landscape_h1_v1",
  "cell_landscape_h0_h1_raw_euclidean_v1",
  "cell_distribution_energy_shared_pca_v1"
)
completion <- vector("list", nrow(authorized))
total_summaries <- total_landscapes <- total_energy <- total_methods <- 0
total_pairs <- total_views <- total_files <- 0

for (index in seq_len(nrow(authorized))) {
  unit <- authorized[index, , drop = FALSE]
  path <- file.path(result_root, safe_name(unit$robustness_group_id))
  files <- sort(list.files(path, full.names = FALSE, all.files = TRUE,
                           no.. = TRUE), method = "radix")
  if (!dir.exists(path) || length(files) != 11L || any(file.info(
      file.path(path, files)
    )$isdir)) stop("MV5-AJ group file cardinality failed.", call. = FALSE)
  status <- read_csv(file.path(path, "status.csv"))
  source_identity <- read_csv(file.path(path, "source_identity.csv"))
  manifest <- read_csv(file.path(path, "artifact_manifest.csv"))
  if (nrow(status) != 1L || nrow(source_identity) != 1L ||
      nrow(manifest) != 9L || anyDuplicated(manifest$artifact_file) ||
      status$status != "completed" ||
      status$robustness_group_id != unit$robustness_group_id ||
      status$execution_order != unit$execution_order ||
      status$execution_queue_sha256 != queue_sha ||
      status$source_freeze_sha256 != unit$source_freeze_sha256 ||
      status$implementation_sha256 != unit$implementation_sha256 ||
      status$python_executable_sha256 != python_sha ||
      status$prospective_head != execution_head ||
      source_identity$robustness_group_id != unit$robustness_group_id ||
      source_identity$fold_id != unit$fold_id ||
      source_identity$seed != unit$seed ||
      source_identity$representation != unit$representation ||
      source_identity$configuration_id != unit$configuration_id ||
      source_identity$execution_queue_sha256 != queue_sha ||
      source_identity$source_freeze_sha256 != unit$source_freeze_sha256 ||
      source_identity$implementation_sha256 != unit$implementation_sha256 ||
      source_identity$python_executable_sha256 != python_sha ||
      source_identity$prospective_head != execution_head ||
      !closed_frame(status) || !closed_frame(source_identity) ||
      !closed_frame(manifest) || any(!strict(manifest$private_artifact))) {
    stop("MV5-AJ group identity/manifest failed at row ", index,
         call. = FALSE)
  }
  manifest_ok <- vapply(seq_len(nrow(manifest)), function(row) {
    artifact <- file.path(path, manifest$artifact_file[[row]])
    file.exists(artifact) && file_sha(artifact) == manifest$sha256[[row]] &&
      as.numeric(file.info(artifact)$size) == manifest$bytes[[row]]
  }, logical(1L))
  if (!all(manifest_ok) || file_sha(file.path(path, "artifact_manifest.csv")) !=
        status$artifact_manifest_sha256) {
    stop("MV5-AJ artifact hash validation failed at row ", index,
         call. = FALSE)
  }

  fold_study <- sub("^large_loso_v1:", "", unit$fold_id)
  expected_d1 <- inventory[
    inventory$source_type == "sct" & inventory$fold_study == fold_study &
      inventory$seed == unit$seed, , drop = FALSE
  ]
  expected_g <- inventory[
    inventory$source_type == "integrated" &
      inventory$fold_study == fold_study & inventory$seed == unit$seed,
    , drop = FALSE
  ]
  if (nrow(expected_d1) != 1L || nrow(expected_g) != 1L ||
      source_identity$d1_source_sha256 != expected_d1$sha256 ||
      source_identity$integrated_source_sha256 != expected_g$sha256) {
    stop("MV5-AJ private source identity failed at row ", index,
         call. = FALSE)
  }

  views <- read_csv(file.path(path, "view_metrics.csv"))
  pairs <- read_csv(file.path(path, "pair_scope.csv"))
  summaries <- read_csv(file.path(path, "landscape_summary.csv"))
  landscapes <- read_csv(file.path(path, "landscape_pairs.csv"))
  energy <- read_csv(file.path(path, "energy_pairs.csv"))
  methods <- read_csv(file.path(path, "method_rows.csv"))
  internal <- read_csv(file.path(path, "internal_resources.csv"))
  if (nrow(views) != 90L || any(views$point_count != unit$cells) ||
      any(views$coordinate_count != unit$coordinates) ||
      any(views$point_metric != unit$point_metric) ||
      any(!is.finite(views$minimum_row_norm) |
            views$minimum_row_norm < 0) ||
      any(!is.finite(views$maximum_row_norm) |
            views$maximum_row_norm < views$minimum_row_norm) ||
      any(views$point_selection !=
            "sha256_sample_seed_cell_nested_v1") ||
      any(!grepl("^[0-9a-f]{64}$", views$selected_point_ids_sha256)) ||
      any(!grepl("^[0-9a-f]{64}$", views$nested192_point_ids_sha256)) ||
      any(!strict(views$nested192_is_exact_prefix)) ||
      any(views$h0_intervals != views$point_count) ||
      any(views$essential_h0_intervals != 1L) ||
      any(views$finite_intervals !=
            views$h0_intervals + views$h1_intervals - 1L) ||
      any(!strict(views$h0_mst_oracle_passed)) || !closed_frame(views) ||
      nrow(pairs) != unit$biological_pairs ||
      anyDuplicated(pairs$pair_request_id) ||
      any(pairs$pair_scope != "held_out_query_to_training_reference") ||
      !closed_frame(pairs) || nrow(summaries) != 180L ||
      any(table(summaries$sample_id) != 2L) ||
      !setequal(summaries$homology_dimension, c("H0", "H1")) ||
      any(summaries$level_policy != "all_consecutive_active_levels") ||
      any(summaries$integration_policy !=
            "exact_linear_critical_pair_segments") ||
      !closed_frame(summaries) ||
      nrow(landscapes) != unit$landscape_request_rows ||
      anyDuplicated(paste(landscapes$pair_request_id,
                          landscapes$homology_dimension)) ||
      any(table(landscapes$pair_request_id) != 2L) ||
      !setequal(landscapes$homology_dimension, c("H0", "H1")) ||
      any(!strict(landscapes$exact)) ||
      any(!strict(landscapes$all_active_levels)) ||
      any(strict(landscapes$level_cap_applied)) ||
      any(!is.finite(landscapes$distance) | landscapes$distance < 0) ||
      !close_numeric(landscapes$squared_distance,
                     landscapes$distance^2) ||
      length(unique(landscapes$subchunk_id)) != unit$landscape_subchunks ||
      any(table(landscapes$subchunk_id) > 250L) ||
      !closed_frame(landscapes) ||
      nrow(energy) != unit$energy_request_rows ||
      anyDuplicated(energy$pair_request_id) ||
      any(!is.finite(energy$distance) | energy$distance < 0) ||
      !closed_frame(energy) || nrow(methods) != unit$assembled_method_rows ||
      !setequal(methods$method_id, method_ids) ||
      any(table(methods$pair_request_id) != 4L) ||
      any(table(methods$method_id) != unit$biological_pairs) ||
      any(!is.finite(methods$distance) | methods$distance < 0) ||
      !closed_frame(methods) || nrow(internal) != 1L ||
      internal$transformed_views != 90L ||
      internal$pair_coverage_rows != unit$biological_pairs ||
      internal$landscape_summary_rows != 180L ||
      internal$landscape_pair_rows != unit$landscape_request_rows ||
      internal$energy_pair_rows != unit$energy_request_rows ||
      internal$method_rows != unit$assembled_method_rows ||
      !closed_frame(internal)) {
    stop("MV5-AJ scientific/cardinality contract failed at row ", index,
         call. = FALSE)
  }

  h0 <- landscapes[landscapes$homology_dimension == "H0", , drop = FALSE]
  h1 <- landscapes[landscapes$homology_dimension == "H1", , drop = FALSE]
  h1 <- h1[match(h0$pair_request_id, h1$pair_request_id), , drop = FALSE]
  energy_match <- energy[match(h0$pair_request_id, energy$pair_request_id),
                         , drop = FALSE]
  expected <- list(
    cell_landscape_h0_v1 = h0$distance,
    cell_landscape_h1_v1 = h1$distance,
    cell_landscape_h0_h1_raw_euclidean_v1 = sqrt(
      h0$distance^2 + h1$distance^2
    ),
    cell_distribution_energy_shared_pca_v1 = energy_match$distance
  )
  for (method_id in names(expected)) {
    observed <- methods[methods$method_id == method_id, , drop = FALSE]
    observed <- observed[match(h0$pair_request_id, observed$pair_request_id),
                         , drop = FALSE]
    if (anyNA(observed$pair_request_id) ||
        !close_numeric(observed$distance, expected[[method_id]])) {
      stop("MV5-AJ four-method reconstruction failed at row ", index,
           call. = FALSE)
    }
  }

  completion[[index]] <- data.frame(
    contract_id = "mv05aj_nested256_unit_completion_v1",
    robustness_group_id = unit$robustness_group_id,
    configuration_execution_order = unit$configuration_execution_order,
    fold_id = unit$fold_id, seed = unit$seed,
    representation = unit$representation,
    configuration_id = unit$configuration_id,
    files = length(files), manifest_artifacts = nrow(manifest),
    views = nrow(views), biological_pairs = nrow(pairs),
    landscape_summary_rows = nrow(summaries),
    landscape_rows = nrow(landscapes), energy_rows = nrow(energy),
    method_rows = nrow(methods),
    private_bytes = sum(file.info(file.path(path, files))$size),
    execution_elapsed_seconds = resources$elapsed_seconds[[index]],
    peak_process_tree_rss_bytes =
      resources$peak_process_tree_rss_bytes[[index]],
    manifest_hashes_valid = TRUE, method_reconstruction_valid = TRUE,
    labels_opened = FALSE, outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
  total_views <- total_views + nrow(views)
  total_pairs <- total_pairs + nrow(pairs)
  total_summaries <- total_summaries + nrow(summaries)
  total_landscapes <- total_landscapes + nrow(landscapes)
  total_energy <- total_energy + nrow(energy)
  total_methods <- total_methods + nrow(methods)
  total_files <- total_files + length(files)
}
completion <- do.call(rbind, completion)
if (total_views != 13500L || total_pairs != 70700L ||
    total_summaries != 27000L || total_landscapes != 141400L ||
    total_energy != 70700L || total_methods != 282800L ||
    total_files != 1650L || !closed_frame(completion)) {
  stop("MV5-AJ aggregate scientific cardinalities failed.", call. = FALSE)
}

repeat_units <- authorized[strict(authorized$deterministic_repeat_required),
                           , drop = FALSE]
if (nrow(repeat_units) != 2L) {
  stop("MV5-AJ requires two prospectively marked repeat groups.", call. = FALSE)
}
repeat_rows <- list()
cursor <- 0L
for (index in seq_len(nrow(repeat_units))) {
  unit <- repeat_units[index, , drop = FALSE]
  stem <- safe_name(unit$robustness_group_id)
  first_manifest <- read_csv(file.path(result_root, stem,
                                       "artifact_manifest.csv"))
  second_manifest <- read_csv(file.path(repeat_root, stem,
                                        "artifact_manifest.csv"))
  deterministic <- first_manifest[
    strict(first_manifest$deterministic_repeat_required), , drop = FALSE
  ]
  matched <- merge(
    deterministic[c("artifact_file", "sha256", "bytes")],
    second_manifest[c("artifact_file", "sha256", "bytes")],
    by = "artifact_file", suffixes = c("_first", "_repeat"), sort = TRUE
  )
  if (nrow(matched) != 8L ||
      any(matched$sha256_first != matched$sha256_repeat) ||
      any(matched$bytes_first != matched$bytes_repeat)) {
    stop("MV5-AJ deterministic repeat failed for ",
         unit$robustness_group_id, call. = FALSE)
  }
  for (row in seq_len(nrow(matched))) {
    cursor <- cursor + 1L
    repeat_rows[[cursor]] <- data.frame(
      contract_id = "mv05aj_deterministic_repeat_v1",
      robustness_group_id = unit$robustness_group_id,
      representation = unit$representation,
      artifact_file = matched$artifact_file[[row]],
      sha256 = matched$sha256_first[[row]],
      bytes = matched$bytes_first[[row]], exact_match = TRUE,
      labels_opened = FALSE, outcomes_computed = FALSE,
      stringsAsFactors = FALSE
    )
  }
}
repeat_evidence <- do.call(rbind, repeat_rows)

validation <- data.frame(
  contract_id = "mv05aj_independent_validation_v1",
  category = c(
    "queue_authorization", "execution_identity", "private_sources",
    "group_publications", "manifest_hashes", "strict_nested_point_identity",
    "view_ph_mst",
    "biological_pair_scope", "exact_all_active_landscapes",
    "energy_rows", "four_method_reconstruction", "aggregate_cardinalities",
    "resource_caps", "deterministic_repeat", "immutable_resume",
    "label_outcome_closure"
  ),
  passed = TRUE,
  observed = c(
    "150_nested256_only", execution_head, "150_hash_bound",
    "150_complete_1650_files", "1350_of_1350_valid",
    "13500_192_prefix_256_subset_384",
    "13500_views_all_mst_passed", "70700_heldout_training_pairs",
    "141400_rows_27000_summaries", "70700_rows",
    "282800_rows_reconstructed", "all_exact",
    paste0(round(sum(resources$elapsed_seconds), 3), "s_",
           max(resources$peak_process_tree_rss_bytes), "B_",
           max(resources$cumulative_private_bytes), "B"),
    "16_of_16_exact", "1650_of_1650_unchanged_150_reused",
    "labels_0_outcomes_0"
  ),
  labels_opened = FALSE, outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
resume_evidence <- data.frame(
  contract_id = "mv05aj_resume_validation_v1",
  snapshot_rows_before = nrow(before), snapshot_rows_after = nrow(after),
  unchanged_rows = sum(before$relative_path == after$relative_path &
                         before$bytes == after$bytes &
                         before$sha256 == after$sha256),
  reused_validated_groups = sum(
    resume_resources$disposition == "reused_validated"
  ),
  exact_snapshot_match = identical(before, after),
  labels_opened = FALSE, outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
decision <- data.frame(
  contract_id = "mv05aj_nested256_configuration_decision_v1",
  configuration_id = "nested_cells_256_pc30_euclidean_v1",
  calculation_complete = TRUE, independent_validation_passed = TRUE,
  configuration_calculation_accepted = TRUE,
  robustness_outcome_execution_authorized = FALSE,
  other_configuration_execution_authorized = FALSE,
  next_action = "prefreeze_nested256_prediction_locked_outcome_evaluation",
  labels_opened = FALSE, outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
resources$contract_id <- "mv05aj_nested256_production_resource_metric_v1"
write_once(resources, file.path(
  output_dir, "mv05aj-nested256-production-resources-2026-08-11.csv"
))
write_once(completion, file.path(
  output_dir, "mv05aj-nested256-unit-completion-2026-08-11.csv"
))
write_once(validation, file.path(
  output_dir, "mv05aj-independent-validation-2026-08-11.csv"
))
write_once(repeat_evidence, file.path(
  output_dir, "mv05aj-deterministic-repeat-2026-08-11.csv"
))
write_once(resume_evidence, file.path(
  output_dir, "mv05aj-resume-validation-2026-08-11.csv"
))
write_once(decision, file.path(
  output_dir, "mv05aj-configuration-decision-2026-08-11.csv"
))
message("MV5-AJ independent artifact validation passed: groups=150 pairs=70700 outcomes=0")
