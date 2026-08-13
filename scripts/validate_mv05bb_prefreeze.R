#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: validate_mv05bb_prefreeze.R OUTPUT_DIR")
out <- normalizePath(args[[1L]], mustWork = TRUE)
read_out <- function(name) utils::read.csv(file.path(out, name), stringsAsFactors = FALSE)
kernel <- read_out("mv05bb-rust-kernel-contract-2026-08-13.csv")
ffi <- read_out("mv05bb-ffi-contract-2026-08-13.csv")
corpus <- read_out("mv05bb-equivalence-corpus-2026-08-13.csv")
gate <- read_out("mv05bb-prototype-gate-2026-08-13.csv")
decision <- read_out("mv05bb-continuation-decision-2026-08-13.csv")
checks <- data.frame(
  contract_id = "mv05bb_validation_v1",
  validation_id = c("kernel_scope", "all_levels", "separate_dimensions",
                    "ffi", "equivalence", "speed", "memory", "fallback",
                    "package_boundary", "toolchain_boundary", "closed_scopes"),
  passed = c(
    nrow(kernel) == 13L && "forbidden" %in% kernel$boundary_id,
    grepl("all consecutive active", kernel$frozen_requirement[
      kernel$boundary_id == "levels"], fixed = TRUE),
    grepl("one H0 or H1", kernel$frozen_requirement[
      kernel$boundary_id == "dimensions"], fixed = TRUE),
    nrow(ffi) == 10L && identical(ffi$field_order, 1:10),
    nrow(corpus) == 5L && sum(corpus$required_results) == 443L &&
      all(c("MV5AY_exact_only", "MV5AY_adaptive_H1") %in% corpus$scope),
    grepl("3x median", gate$requirement[gate$gate_id == "speed"], fixed = TRUE),
    grepl("1 GiB", gate$requirement[gate$gate_id == "memory"], fixed = TRUE),
    grepl("fallback", gate$requirement[gate$gate_id == "fallback"], fixed = TRUE),
    grepl("Rust unavailable", gate$requirement[gate$gate_id == "package"], fixed = TRUE),
    decision$rust_prototype_authorized && !decision$rust_toolchain_install_authorized &&
      !decision$rust_production_adoption_authorized,
    decision$r_engine_retained && !decision$additional_seed_production_authorized &&
      !decision$partitions_authorized
  ), stringsAsFactors = FALSE
)
if (any(!checks$passed)) stop("MV5-BB validation failed: ",
                              paste(checks$validation_id[!checks$passed], collapse = ", "))
utils::write.csv(checks, file.path(out, "mv05bb-independent-validation-2026-08-13.csv"),
                 row.names = FALSE, quote = TRUE)
cat("MV5-BB validation passed:", nrow(checks), "categories\n")
