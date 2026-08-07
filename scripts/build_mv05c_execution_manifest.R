#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("usage: build_mv05c_execution_manifest.R HISTORICAL_DIR AUDIT_DIR PRIVATE_DIR",
       call. = FALSE)
}

historical_dir <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
audit_dir <- args[[2L]]
private_dir <- args[[3L]]
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(private_dir, recursive = TRUE, showWarnings = FALSE)

source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/mv05_benchmark_execution.R")

metadata_path <- file.path(historical_dir, "joined_metadata_cellcounts.csv")
raw_list_path <- file.path(historical_dir, "expr_list_raw.Rds")
metadata_raw <- utils::read.csv(
  metadata_path, stringsAsFactors = FALSE, check.names = TRUE
)
metadata <- data.frame(
  sample_id = metadata_raw$orig.ident,
  study = metadata_raw$SRA,
  tissue = tolower(metadata_raw$Tissue.x),
  approach = metadata_raw$Approach.x,
  filtered_cells = as.integer(metadata_raw$Number_of_Cells_After_Filtering),
  source_path = metadata_raw$File.Path,
  stringsAsFactors = FALSE
)

tissue_counts <- aggregate(
  sample_id ~ tissue, metadata, function(x) length(unique(x))
)
names(tissue_counts)[[2L]] <- "sample_count"
study_counts <- aggregate(
  study ~ tissue, metadata, function(x) length(unique(x))
)
names(study_counts)[[2L]] <- "study_count"
selection <- merge(tissue_counts, study_counts, by = "tissue", sort = TRUE)
selection$all_samples_cell_eligible <- vapply(selection$tissue, function(value) {
  all(metadata$filtered_cells[metadata$tissue == value] >= 384L)
}, logical(1L))
selection$eligible <- selection$study_count >= 2L &
  selection$sample_count >= 2L & selection$all_samples_cell_eligible
eligible <- selection[selection$eligible, , drop = FALSE]
eligible <- eligible[order(eligible$sample_count, eligible$tissue,
                           method = "radix"), , drop = FALSE]
if (!nrow(eligible)) stop("No tissue satisfies the frozen bounded-feasibility rule.")
selected_tissue <- eligible$tissue[[1L]]
selection$selected <- selection$tissue == selected_tissue
selection$selection_rule <- paste(
  "at_least_two_studies;all_samples_at_least_384_cells;",
  "minimum_sample_count;lexical_tissue_tie_break", sep = ""
)
selection$outcome_labels_used_after_selection <- FALSE

pilot <- metadata[metadata$tissue == selected_tissue, , drop = FALSE]
pilot <- pilot[order(pilot$study, pilot$sample_id, method = "radix"), , drop = FALSE]
if (nrow(pilot) != 6L || length(unique(pilot$study)) != 2L ||
    any(pilot$filtered_cells < 384L)) {
  stop("Frozen bounded-feasibility selection did not reproduce the six-sample case.")
}

private_cache_path <- file.path(private_dir, "mv05c-existing-data-raw-cache.rds")
started <- proc.time()[["elapsed"]]
if (file.exists(private_cache_path)) {
  private_cache <- readRDS(private_cache_path)
  read_mode <- "resume_private_cache"
  if (!identical(private_cache$contract_id,
                 "mv05c_existing_data_raw_cache_v1") ||
      !identical(private_cache$sample_ids, pilot$sample_id) ||
      !identical(private_cache$source_size_bytes,
                 unname(file.info(raw_list_path)$size))) {
    stop("Existing private MV5-C cache has stale identity.")
  }
  raw_samples <- private_cache$samples
} else {
  raw_list <- readRDS(raw_list_path)
  read_mode <- "full_historical_source"
  if (!is.list(raw_list) || is.null(names(raw_list)) ||
      !all(pilot$sample_id %in% names(raw_list))) {
    stop("Historical raw-count list is missing one or more selected samples.")
  }
  raw_samples <- raw_list[pilot$sample_id]
  rm(raw_list)
  invisible(gc())
  private_cache <- list(
    contract_id = "mv05c_existing_data_raw_cache_v1",
    pilot_id = "mv05c_one_tissue_feasibility_v1",
    sample_ids = pilot$sample_id,
    source_path = normalizePath(raw_list_path, winslash = "/", mustWork = TRUE),
    source_size_bytes = unname(file.info(raw_list_path)$size),
    source_sha256 = digest::digest(
      file = raw_list_path, algo = "sha256", serialize = FALSE
    ),
    samples = raw_samples
  )
  partial <- paste0(private_cache_path, ".partial")
  saveRDS(private_cache, partial, compress = FALSE, version = 3)
  if (!file.rename(partial, private_cache_path)) {
    stop("Failed to atomically publish the private MV5-C raw cache.")
  }
}
read_seconds <- proc.time()[["elapsed"]] - started
for (sample_id in names(raw_samples)) {
  value <- raw_samples[[sample_id]]
  if ((!is.matrix(value) && !inherits(value, "Matrix")) ||
      is.null(rownames(value)) || is.null(colnames(value)) ||
      anyDuplicated(rownames(value)) || anyDuplicated(colnames(value)) ||
      ncol(value) != pilot$filtered_cells[match(sample_id, pilot$sample_id)]) {
    stop("Raw-count sample failed identity/shape validation: ", sample_id)
  }
}

private_cache_sha256 <- digest::digest(
  file = private_cache_path, algo = "sha256", serialize = FALSE
)

seeds <- 20260805:20260809
cell_rows <- list()
for (sample_id in pilot$sample_id) {
  for (seed in seeds) {
    selected <- select_matched_cells(
      colnames(raw_samples[[sample_id]]), n = 384L, seed = seed
    )
    cell_rows[[length(cell_rows) + 1L]] <- data.frame(
      contract_id = "mv05c_matched_cell_manifest_v1",
      pilot_id = "mv05c_one_tissue_feasibility_v1",
      sample_id = sample_id,
      seed = as.integer(seed),
      eligible_cells = ncol(raw_samples[[sample_id]]),
      selected_cells = 384L,
      selected_cell_sha256 = attr(selected, "selected_cell_sha256"),
      cell_id = as.character(selected),
      outcome_label_state = "closed",
      biological_outcomes_computed = FALSE,
      stringsAsFactors = FALSE
    )
  }
}
cell_manifest <- do.call(rbind, cell_rows)

sample_manifest <- pilot[c("sample_id", "study", "filtered_cells")]
sample_manifest$pilot_id <- "mv05c_one_tissue_feasibility_v1"
sample_manifest$outcome_label_state <- "closed"
sample_manifest$biological_outcomes_computed <- FALSE
sample_manifest <- sample_manifest[c(
  "pilot_id", "sample_id", "study", "filtered_cells",
  "outcome_label_state", "biological_outcomes_computed"
)]

folds <- sort(unique(pilot$study), method = "radix")
execution_rows <- list()
for (held_out in folds) {
  fold_id <- paste0("mv05c_loso_v1:", held_out)
  fit_scope_id <- paste0(fold_id, ":training_reference")
  for (seed in seeds) {
    for (representation in c("sct_fold", "inductive_integrated")) {
      for (view_id in c("cell_topology_v1", "gene_topology_v1")) {
        executable <- representation == "sct_fold" || view_id == "cell_topology_v1"
        execution_rows[[length(execution_rows) + 1L]] <- data.frame(
          contract_id = "mv05c_execution_manifest_v1",
          pilot_id = "mv05c_one_tissue_feasibility_v1",
          fold_id = fold_id,
          fit_scope_id = fit_scope_id,
          held_out_study = held_out,
          seed = as.integer(seed),
          representation = representation,
          view_id = view_id,
          sample_count = nrow(pilot),
          training_sample_count = sum(pilot$study != held_out),
          held_out_sample_count = sum(pilot$study == held_out),
          genes_target = 500L,
          cells_per_sample = 384L,
          pcs = 30L,
          homology_dimensions = "H0;H1",
          execution_policy = if (executable) "execute" else
            "structured_unavailable_no_inductive_gene_correction",
          outcome_label_state = "closed",
          biological_outcomes_computed = FALSE,
          stringsAsFactors = FALSE
        )
      }
    }
  }
}
execution_manifest <- do.call(rbind, execution_rows)

selection_path <- file.path(
  audit_dir, "mv05c-prospective-tissue-selection-2026-08-06.csv"
)
sample_path <- file.path(
  audit_dir, "mv05c-label-closed-sample-manifest-2026-08-06.csv"
)
cell_path <- file.path(
  audit_dir, "mv05c-matched-cell-manifest-2026-08-06.csv"
)
execution_path <- file.path(
  audit_dir, "mv05c-execution-manifest-2026-08-06.csv"
)
source_path <- file.path(
  audit_dir, "mv05c-source-preflight-2026-08-06.csv"
)
utils::write.csv(selection, selection_path, row.names = FALSE)
utils::write.csv(sample_manifest, sample_path, row.names = FALSE)
utils::write.csv(cell_manifest, cell_path, row.names = FALSE)
utils::write.csv(execution_manifest, execution_path, row.names = FALSE)
utils::write.csv(data.frame(
  contract_id = private_cache$contract_id,
  source_file = basename(raw_list_path),
  source_size_bytes = private_cache$source_size_bytes,
  source_sha256 = private_cache$source_sha256,
  source_read_mode = read_mode,
  source_read_seconds = read_seconds,
  private_cache_file = basename(private_cache_path),
  private_cache_size_bytes = unname(file.info(private_cache_path)$size),
  private_cache_sha256 = private_cache_sha256,
  selected_samples = nrow(pilot),
  selected_studies = length(unique(pilot$study)),
  minimum_filtered_cells = min(pilot$filtered_cells),
  all_samples_at_least_384_cells = all(pilot$filtered_cells >= 384L),
  outcome_label_state = "closed",
  biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
), source_path, row.names = FALSE)

stopifnot(
  identical(selected_tissue, "liver"),
  nrow(sample_manifest) == 6L,
  nrow(cell_manifest) == 6L * 5L * 384L,
  nrow(execution_manifest) == 2L * 5L * 2L * 2L,
  !any(c("tissue", "approach") %in% names(sample_manifest)),
  !any(c("tissue", "approach") %in% names(cell_manifest)),
  !any(c("tissue", "approach") %in% names(execution_manifest)),
  all(sample_manifest$outcome_label_state == "closed"),
  all(cell_manifest$outcome_label_state == "closed"),
  all(execution_manifest$outcome_label_state == "closed"),
  all(!sample_manifest$biological_outcomes_computed),
  all(!cell_manifest$biological_outcomes_computed),
  all(!execution_manifest$biological_outcomes_computed)
)
message("Frozen MV5-C one-tissue execution manifest without endpoint evaluation.")
