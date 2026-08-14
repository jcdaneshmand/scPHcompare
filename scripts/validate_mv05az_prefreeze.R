#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: validate_mv05az_prefreeze.R OUTPUT_DIR")
out <- normalizePath(args[[1L]], mustWork = TRUE)
read_out <- function(name) utils::read.csv(file.path(out, name), stringsAsFactors = FALSE)
inventory <- read_out("mv05az-resampling-asset-inventory-2026-08-12.csv")
target <- read_out("mv05az-matched-axis-target-2026-08-12.csv")
contract <- read_out("mv05az-stability-contract-2026-08-12.csv")
speed <- read_out("mv05az-acceleration-gate-2026-08-12.csv")
decision <- read_out("mv05az-continuation-decision-2026-08-12.csv")

checks <- data.frame(
  contract_id = "mv05az_validation_v1",
  validation_id = c("asset_inventory", "identifiability_boundary", "target_scope",
                    "matched_axes", "stability_rule", "landscape_contract",
                    "acceleration_equivalence", "rust_gate", "resources",
                    "closed_partitions", "closed_labels", "next_sprint"),
  passed = c(
    nrow(inventory) == 3L && identical(inventory$diagrams, c(80L, 150L, 56L)),
    all(!inventory$identifies_current_primary_k),
    nrow(target) == 4L && sum(target$diagrams_additional) == 160L &&
      sum(target$pairs_additional) == 720L &&
      sum(target$pairs_additional[target$view_id == "gene_topology_v1"]) == 360L,
    all(target$samples == 10L & target$seeds == 5L &
          target$additional_seed_count == 4L),
    all(c("seed_partitions", "stability", "uncertainty", "selection") %in%
          contract$rule_id) && grepl("ten pairwise adjusted Rand", contract$frozen_value[
            contract$rule_id == "stability"], fixed = TRUE),
    grepl("all finite intervals", contract$frozen_value[
      contract$rule_id == "landscape"], fixed = TRUE),
    all(c("fixtures", "exact_corpus", "adaptive_corpus", "determinism",
          "throughput") %in% speed$gate_id) && all(speed$abort_on_failure),
    grepl("begin Rust only if", speed$requirement[speed$gate_id == "rust_trigger"],
          fixed = TRUE),
    grepl("40 worker-hours", speed$requirement[speed$gate_id == "production_cap"],
          fixed = TRUE),
    isTRUE(decision$prefreeze_accepted) && !decision$partitions_authorized &&
      !decision$additional_seed_production_authorized,
    grepl("labels", contract$frozen_value[contract$rule_id == "labels"],
          ignore.case = TRUE) && grepl("unavailable", contract$frozen_value[
            contract$rule_id == "labels"], fixed = TRUE),
    isTRUE(decision$acceleration_benchmark_authorized) &&
      identical(decision$next_sprint,
                "MV5-BA_corrected_persim_equivalence_and_speed_benchmark")
  ), stringsAsFactors = FALSE
)
if (any(!checks$passed)) stop("MV5-AZ validation failed: ",
                              paste(checks$validation_id[!checks$passed], collapse = ", "))
utils::write.csv(checks, file.path(out, "mv05az-independent-validation-2026-08-12.csv"),
                 row.names = FALSE, quote = TRUE)
cat("MV5-AZ validation passed:", nrow(checks), "categories\n")
