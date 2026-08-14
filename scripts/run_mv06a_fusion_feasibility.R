args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(
    "Usage: run_mv06a_fusion_feasibility.R <repo-root> <input-dir> <output-dir>",
    call. = FALSE
  )
}

repo_root <- normalizePath(args[[1L]], mustWork = TRUE)
input_dir <- normalizePath(args[[2L]], mustWork = TRUE)
output_dir <- args[[3L]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, mustWork = TRUE)

source(file.path(repo_root, "R", "topological_distance_contract.R"))
source(file.path(repo_root, "R", "mv06_fusion.R"))

contract_id <- "mv06a_label_closed_fusion_feasibility_v1"
fit_scope_id <- "mv06a_descriptive_all_pilot_samples"

expected_sources <- data.frame(
  file = c(
    "bone__integrated__cell_topology_v1.rds",
    "bone__integrated__gene_topology_v1.rds",
    "bone__sct_whole__cell_topology_v1.rds",
    "bone__sct_whole__gene_topology_v1.rds",
    "large__sct_whole__cell_topology_v1.rds",
    "large__sct_whole__gene_topology_v1.rds",
    "large__seurat_integration__cell_topology_v1.rds",
    "large__seurat_integration__gene_topology_v1.rds"
  ),
  bytes = c(1374, 1381, 1375, 1374, 3798, 3794, 3795, 3810),
  sha256 = c(
    "b2988fbe2b69d350f1c219c5bf575b4883acaea563f097df34d30ac28bd31789",
    "cc79c500b1d2d2d5303913aae9ed68a2c3ccf1ac063d75ad2c4901e9ed2df7c8",
    "82991bd488d7dc33a383ede6de755e7d7a30e3b2f0b87ac977c2eb0961ca189a",
    "9616a928b9fbf6497c1d74a7e481ec9fc02064391f5e0e378134ca17af4588a2",
    "0a4d57baef8c6fb7498d44ca00efa97ded9173522cffafc2aa238147f372db8e",
    "fc3e4675e9b09287626e5174fd8f0c5f02eda078fe5167c364a0519bc2af0ac5",
    "7346d8b9571d5c997d5f673a3b38c47a0f9f0b64bf59a8493bca961246a37e30",
    "846a5efe764f4f444711af3e4bdaa7000a7530cfa2a82b0770c0bb2e5ae7bc1d"
  ),
  stringsAsFactors = FALSE
)
expected_sources$path <- file.path(input_dir, expected_sources$file)
if (any(!file.exists(expected_sources$path))) {
  stop("One or more frozen MV6-A source bundles are missing.", call. = FALSE)
}
actual_bytes <- as.numeric(file.info(expected_sources$path)$size)
actual_hashes <- vapply(
  expected_sources$path, digest::digest, character(1L),
  algo = "sha256", file = TRUE, serialize = FALSE
)
if (!identical(actual_bytes, expected_sources$bytes) ||
    !identical(unname(actual_hashes), expected_sources$sha256)) {
  stop("A frozen MV6-A source identity does not match the prefreeze.",
       call. = FALSE)
}
source_inventory <- expected_sources[c("file", "bytes", "sha256")]
source_inventory$contract_id <- contract_id
source_inventory <- source_inventory[
  c("contract_id", "file", "bytes", "sha256")
]

strata <- c(
  "bone__integrated", "bone__sct_whole",
  "large__sct_whole", "large__seurat_integration"
)
fits <- lapply(strata, function(stratum_id) {
  cell <- readRDS(file.path(
    input_dir, paste0(stratum_id, "__cell_topology_v1.rds")
  ))
  gene <- readRDS(file.path(
    input_dir, paste0(stratum_id, "__gene_topology_v1.rds")
  ))
  mv06_fit_fusion_components_v1(
    cell, gene,
    fit_scope_id = fit_scope_id,
    gene_weights = c(0, 0.25, 0.5, 0.75, 1)
  )
})
names(fits) <- strata

scales <- do.call(rbind, lapply(fits, mv06_scale_diagnostics_v1))
pairs <- do.call(rbind, lapply(fits, mv06_pair_diagnostics_v1))
weights <- do.call(rbind, lapply(fits, mv06_weight_diagnostics_v1))
correlations <- do.call(
  rbind, lapply(fits, mv06_correlation_diagnostics_v1)
)
neighbors <- do.call(rbind, lapply(fits, mv06_neighbor_diagnostics_v1))
matrix_hashes <- do.call(rbind, lapply(fits, mv06_matrix_hashes_v1))

neighbor_summary <- stats::aggregate(
  cbind(overlap_count, jaccard, exact_neighbor_set) ~
    contract_id + stratum_id + axis_id + k,
  data = transform(neighbors, exact_neighbor_set = as.numeric(exact_neighbor_set)),
  FUN = mean
)
names(neighbor_summary)[names(neighbor_summary) == "overlap_count"] <-
  "mean_overlap_count"
names(neighbor_summary)[names(neighbor_summary) == "jaccard"] <-
  "mean_jaccard"
names(neighbor_summary)[names(neighbor_summary) == "exact_neighbor_set"] <-
  "exact_neighbor_set_fraction"

weight_summary <- stats::aggregate(
  distance ~ contract_id + stratum_id + weight_id + gene_weight + cell_weight,
  data = weights,
  FUN = function(value) c(
    minimum = min(value), median = stats::median(value), maximum = max(value)
  )
)
weight_values <- weight_summary$distance
if (!is.matrix(weight_values) || ncol(weight_values) != 3L) {
  stop("Weight summary aggregation did not produce three frozen statistics.",
       call. = FALSE)
}
weight_summary$minimum <- weight_values[, "minimum"]
weight_summary$median <- weight_values[, "median"]
weight_summary$maximum <- weight_values[, "maximum"]
weight_summary$distance <- NULL

contract <- data.frame(
  contract_id = contract_id,
  accepted_base = "d0192d35a4ab52006aa83b0ad3b0ad6a19f066cb",
  fit_scope_id = fit_scope_id,
  stratum_count = length(fits),
  source_count = nrow(source_inventory),
  scale_count = nrow(scales),
  unordered_pair_count = nrow(pairs),
  weight_row_count = nrow(weights),
  correlation_count = nrow(correlations),
  neighbor_row_count = nrow(neighbors),
  matrix_hash_count = nrow(matrix_hashes),
  biological_values_read = 0L,
  endpoint_rows_computed = 0L,
  clusterings_computed = 0L,
  advanced_fusion_computed = 0L,
  stringsAsFactors = FALSE
)

write_stable_csv <- function(value, name) {
  utils::write.csv(
    value, file.path(output_dir, name),
    row.names = FALSE, na = "", quote = TRUE
  )
}

write_stable_csv(contract, "mv06a-contract.csv")
write_stable_csv(source_inventory, "mv06a-source-inventory.csv")
write_stable_csv(scales, "mv06a-scales.csv")
write_stable_csv(pairs, "mv06a-pairs.csv")
write_stable_csv(weights, "mv06a-weight-grid.csv")
write_stable_csv(weight_summary, "mv06a-weight-summary.csv")
write_stable_csv(correlations, "mv06a-correlations.csv")
write_stable_csv(neighbors, "mv06a-neighbors.csv")
write_stable_csv(neighbor_summary, "mv06a-neighbor-summary.csv")
write_stable_csv(matrix_hashes, "mv06a-matrix-hashes.csv")

artifact_files <- sort(list.files(
  output_dir, pattern = "^mv06a-.*[.]csv$", full.names = TRUE
))
artifact_manifest <- data.frame(
  contract_id = contract_id,
  file = basename(artifact_files),
  bytes = as.numeric(file.info(artifact_files)$size),
  sha256 = vapply(
    artifact_files, digest::digest, character(1L),
    algo = "sha256", file = TRUE, serialize = FALSE
  ),
  stringsAsFactors = FALSE
)
write_stable_csv(artifact_manifest, "mv06a-artifact-manifest.csv")

cat(
  "MV6-A completed:", length(fits), "strata,", nrow(pairs),
  "pairs,", nrow(weights), "weight rows,", nrow(neighbors),
  "neighbor rows\n"
)
