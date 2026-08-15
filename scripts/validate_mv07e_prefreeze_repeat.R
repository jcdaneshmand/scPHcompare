#!/usr/bin/env Rscript

options(warn = 2)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required.")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("usage: validate_mv07e_prefreeze_repeat.R PRIMARY_DIR REPEAT_DIR OUTPUT",
       call. = FALSE)
}
primary <- args[[1L]]; repeat <- args[[2L]]; output <- args[[3L]]
files <- c("mv07e-source-freeze.csv", "mv07e-accession-resolution.csv",
  "mv07e-approach-field-lineage.csv", "mv07e-canonical-approach.csv",
  "mv07e-approach-correction.csv", "mv07e-panel-availability.csv",
  "mv07e-panel-decision.csv", "mv07e-sample-seed-axis.csv",
  "mv07e-fit-scopes.csv", "mv07e-typed-topology.csv", "mv07e-pair-scope.csv",
  "mv07e-landscape-contract.csv", "mv07e-resource-resume-contract.csv",
  "mv07e-label-firewall.csv", "mv07e-stage-authorization.csv",
  "mv07e-acceptance-criteria.csv", "mv07e-decision.csv")
paths_a <- file.path(primary, files); paths_b <- file.path(repeat, files)
if (any(!file.exists(c(paths_a, paths_b))) || file.exists(output)) {
  stop("Repeat inputs are incomplete or output exists.", call. = FALSE)
}
sha <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
result <- data.frame(
  contract_id = "mv07e_byte_repeat_v1", artifact = files,
  primary_sha256 = vapply(paths_a, sha, character(1L)),
  repeat_sha256 = vapply(paths_b, sha, character(1L)),
  primary_bytes = as.numeric(file.info(paths_a)$size),
  repeat_bytes = as.numeric(file.info(paths_b)$size), stringsAsFactors = FALSE)
result$byte_identical <- result$primary_sha256 == result$repeat_sha256 &
  result$primary_bytes == result$repeat_bytes
if (!all(result$byte_identical)) stop("MV7-E repeat is not byte-identical.")
write.table(result, output, sep = ",", row.names = FALSE, col.names = TRUE,
            quote = TRUE, na = "")
message("MV7-E prefreeze repeat: 17/17 artifacts byte-identical")
