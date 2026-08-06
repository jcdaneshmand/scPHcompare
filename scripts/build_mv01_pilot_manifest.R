#!/usr/bin/env Rscript

options(warn = 2)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop(
    "Usage: build_mv01_pilot_manifest.R <historical-dir> <output-csv>",
    call. = FALSE
  )
}

historical_dir <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
output_csv <- args[[2]]

large_path <- file.path(historical_dir, "joined_metadata_cellcounts.csv")
bone_path <- file.path(
  historical_dir,
  "joined_metadata_cellcounts_bonemarrow.csv"
)

if (!file.exists(large_path) || !file.exists(bone_path)) {
  stop("Required historical metadata files are missing.", call. = FALSE)
}

required_large <- c(
  "orig.ident", "SRA", "Tissue.x", "Approach.x",
  "Number_of_Cells_After_Filtering"
)
required_bone <- c(
  "orig.ident", "SRA", "Tissue", "Approach",
  "Number_of_Cells_After_Filtering"
)

large_raw <- utils::read.csv(
  large_path,
  stringsAsFactors = FALSE,
  check.names = TRUE
)
bone_raw <- utils::read.csv(
  bone_path,
  stringsAsFactors = FALSE,
  check.names = TRUE
)

if (!all(required_large %in% names(large_raw))) {
  stop("Large-cohort metadata columns do not match the MV-01 contract.")
}
if (!all(required_bone %in% names(bone_raw))) {
  stop("Bone-marrow metadata columns do not match the MV-01 contract.")
}

normalize_rows <- function(values, cohort, tissue_col, approach_col) {
  result <- data.frame(
    cohort = cohort,
    sample_id = trimws(as.character(values[["orig.ident"]])),
    study = trimws(as.character(values[["SRA"]])),
    tissue = tolower(trimws(as.character(values[[tissue_col]]))),
    approach = trimws(as.character(values[[approach_col]])),
    filtered_cells = suppressWarnings(as.integer(
      values[["Number_of_Cells_After_Filtering"]]
    )),
    stringsAsFactors = FALSE
  )
  if (
    any(!nzchar(result$sample_id)) ||
      anyDuplicated(result$sample_id) ||
      any(!is.finite(result$filtered_cells)) ||
      any(result$filtered_cells <= 0L)
  ) {
    stop("Invalid sample IDs or filtered-cell counts in metadata.")
  }
  result
}

large <- normalize_rows(large_raw, "large", "Tissue.x", "Approach.x")
bone <- normalize_rows(bone_raw, "bone", "Tissue", "Approach")

nearest_to <- function(values, target) {
  values[order(abs(values$filtered_cells - target), values$sample_id), ][1L, ]
}

large_by_tissue <- split(large, large$tissue)
core <- do.call(rbind, lapply(large_by_tissue, function(values) {
  selected <- nearest_to(values, stats::median(values$filtered_cells))
  selected$pilot_role <- "core_tissue_median"
  selected$selection_rule <- paste0(
    "closest_to_within_tissue_median_cells_tie_sample_id"
  )
  selected
}))

large_order <- order(large$filtered_cells, large$sample_id)
minimum <- large[large_order[1L], ]
minimum$pilot_role <- "stress_cell_minimum"
minimum$selection_rule <- "global_minimum_filtered_cells_tie_sample_id"
maximum <- large[large_order[length(large_order)], ]
maximum$pilot_role <- "stress_cell_maximum"
maximum$selection_rule <- "global_maximum_filtered_cells_tie_sample_id"

bone_by_approach <- split(bone, bone$approach)
technical <- do.call(rbind, lapply(bone_by_approach, function(values) {
  ordered <- values[order(values$filtered_cells, values$sample_id), ]
  indices <- pmax(
    1L,
    floor((nrow(ordered) - 1L) * c(0.25, 0.75)) + 1L
  )
  selected <- ordered[indices, , drop = FALSE]
  selected$pilot_role <- c(
    "technical_approach_q25",
    "technical_approach_q75"
  )
  selected$selection_rule <- paste0(
    "within_approach_order_index_floor_nminus1_p_plus1_tie_sample_id"
  )
  selected
}))

manifest <- rbind(core, minimum, maximum, technical)
manifest <- manifest[order(
  match(
    manifest$pilot_role,
    c(
      "core_tissue_median", "stress_cell_minimum", "stress_cell_maximum",
      "technical_approach_q25", "technical_approach_q75"
    )
  ),
  manifest$cohort,
  manifest$tissue,
  manifest$approach,
  manifest$filtered_cells,
  manifest$sample_id
), ]

manifest$subsample_cells <- 384L
manifest$subsample_seeds <- "20260805;20260806;20260807;20260808;20260809"
manifest$primary_representations <- ifelse(
  manifest$cohort == "large",
  "sct_whole;seurat_integration",
  "sct_whole;integrated"
)
manifest$scientific_representative <- FALSE
manifest$manifest_purpose <- "feasibility_and_contract_validation_only"
manifest$source_metadata_file <- ifelse(
  manifest$cohort == "large",
  basename(large_path),
  basename(bone_path)
)
manifest$source_metadata_sha256 <- ifelse(
  manifest$cohort == "large",
  digest::digest(file = large_path, algo = "sha256"),
  digest::digest(file = bone_path, algo = "sha256")
)

if (nrow(manifest) != 14L || anyDuplicated(manifest$sample_id)) {
  stop("MV-01 pilot manifest must contain 14 unique samples.")
}
if (any(manifest$filtered_cells < manifest$subsample_cells)) {
  stop("A selected sample cannot supply the frozen 384-cell subsample.")
}

output_dir <- dirname(output_csv)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(manifest, output_csv, row.names = FALSE, na = "")
