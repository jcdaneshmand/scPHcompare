#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: validate_mv05au_corrected_consumer_prefreeze.R EVIDENCE_DIR OUTPUT_CSV")
e <- normalizePath(args[[1]], mustWork = TRUE); output <- normalizePath(args[[2]], mustWork = FALSE)
setwd("/mnt/e/Repositories/Jonah/PH Pipeline Repo/scPHcompare")
read_e <- function(id) utils::read.csv(file.path(e, paste0("mv05au-", id,
  "-2026-08-12.csv")), stringsAsFactors = FALSE, check.names = FALSE)
checks <- list(); record <- function(id, pass, evidence) {
  checks[[length(checks)+1L]] <<- data.frame(contract_id="mv05au_independent_validation_v1",
    validation_id=id, passed=isTRUE(pass), evidence=evidence, stringsAsFactors=FALSE)
  if (!isTRUE(pass)) stop("MV5-AU validation failed: ", id)
}
c <- read_e("consumer-inventory")
record("unique_first_consumer", sum(c$disposition == "implement_first") == 2L &&
  setequal(c$consumer[c$disposition == "implement_first"],
    c("verified_corrected_bundle_loader", "average_linkage_dendrogram")),
  "only verified loader and no-k average tree eligible")
record("blocked_choices", all(c$disposition[c$additional_choice != "none"] == "closed"),
  "all consumers requiring k, kernels, fusion, curves, or labels closed")
api <- read_e("api-contract")
record("api_boundary", nrow(api)==2L && all(api$default_off) &&
  !any(api$mutates_source) && !any(api$writes_legacy) &&
  all(api$combined_status=="descriptive_not_consumed"), "two additive read-only APIs")
v <- read_e("view-dimension-policy")
record("view_dimension", nrow(v)==6L && sum(v$may_build_tree)==4L &&
  !any(v$may_fuse_across_views) && !any(v$may_build_tree[v$dimension=="combined"]),
  "H0/H1 trees separate in both explicit views; combined excluded")
record("validation_plan", nrow(read_e("validation-plan"))==15L &&
  all(read_e("validation-plan")$required), "15 mandatory implementation validations")
record("abort_rules", nrow(read_e("abort-rules"))==14L, "14 stop rules")
record("migration", identical(read_e("migration-stages")$stage, 0:5) &&
  sum(read_e("migration-stages")$authorized_now)==1L, "six staged gates; current stage only")
freeze <- read_e("source-freeze")
sha <- vapply(freeze$path, function(p) sub(" .*", "", system2("sha256sum", p,
  stdout=TRUE)), character(1))
record("source_freeze", identical(unname(freeze$sha256), unname(sha)),
  "12 source, workflow, and contract hashes reproduce")
p <- read_e("prohibited-change-counters")
record("no_changes", nrow(p)==15L && all(p$value==0L), "all prohibited counters zero")
d <- read_e("continuation-decision")
record("bounded_decision", d$prefreeze_accepted && !d$partitions_authorized &&
  !d$evaluation_authorized && !d$default_change_authorized &&
  !d$legacy_change_authorized && d$next_sprint=="MV5-AV",
  "only loader and separate average-tree implementation authorized")
utils::write.csv(do.call(rbind, checks), output, row.names=FALSE, na="", quote=TRUE)
cat("MV5-AU independent validation passed:", length(checks), "categories\n")
