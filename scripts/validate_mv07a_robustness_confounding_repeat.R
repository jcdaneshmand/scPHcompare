#!/usr/bin/env Rscript

options(warn = 2)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required.")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("usage: validate_mv07a_robustness_confounding_repeat.R PRIMARY REPEAT OUTPUT")
primary <- args[[1L]]; repeat_dir <- args[[2L]]; output <- args[[3L]]
files <- c("mv07a-source-freeze.csv", "mv07a-landscape-contract.csv",
  "mv07a-robustness-coverage.csv", "mv07a-confounding-coverage.csv",
  "mv07a-computation-budget.csv", "mv07a-evidence-gap-registry.csv",
  "mv07a-claim-boundaries.csv", "mv07a-selection-firewall.csv",
  "mv07a-acceptance-criteria.csv", "mv07a-decision.csv")
sha <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
rows <- lapply(files, function(name) {
  a <- file.path(primary, name); b <- file.path(repeat_dir, name)
  data.frame(contract_id = "mv07a_byte_repeat_v1", artifact = name,
    primary_sha256 = if (file.exists(a)) sha(a) else NA_character_,
    repeat_sha256 = if (file.exists(b)) sha(b) else NA_character_,
    primary_bytes = if (file.exists(a)) file.info(a)$size else NA_real_,
    repeat_bytes = if (file.exists(b)) file.info(b)$size else NA_real_,
    identical = file.exists(a) && file.exists(b) && identical(sha(a), sha(b)) &&
      identical(as.numeric(file.info(a)$size), as.numeric(file.info(b)$size)),
    stringsAsFactors = FALSE)
})
rows <- do.call(rbind, rows)
if (!all(rows$identical)) stop("MV7-A repeat mismatch: ",
  paste(rows$artifact[!rows$identical], collapse = ", "), call. = FALSE)
if (file.exists(output)) stop("Refusing to overwrite: ", output, call. = FALSE)
write.table(rows, output, sep = ",", row.names = FALSE, col.names = TRUE,
  quote = TRUE, na = "")
message("MV7-A repeat: 10/10 artifacts byte-identical")
