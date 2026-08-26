#!/usr/bin/env Rscript
options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("usage: build_mv17b_null_qualification_prefreeze.R <output>")
output <- args[[1L]]; if (dir.exists(output)) stop("MV17-B output exists")
dir.create(output, recursive = TRUE)
source("R/mv08z_landscape_production.R")
source("R/mv17_null_models.R")
sha <- .mv08z_sha256_file; atomic <- .mv08z_atomic_csv
head <- tolower(trimws(Sys.getenv("MV17B_PREFREEZE_HEAD", unset = "")))
if (!grepl("^[0-9a-f]{40}$", head)) stop("MV17B_PREFREEZE_HEAD required")
mv17a <- "docs/audits/mv17a-source-inventory-estimand-prefreeze-v4"
m <- .mv08z_read_csv(file.path(mv17a, "mv17a-artifact-manifest.csv"))
if (!all(vapply(file.path(mv17a, m$artifact), sha, character(1L)) == m$sha256))
  stop("MV17-B MV17-A binding drift")
families <- mv17b_null_registry_v1(); fixtures <- mv17b_fixture_registry_v1()
implementations <- c("R/mv17_null_models.R",
  "scripts/build_mv17b_null_qualification_prefreeze.R",
  "scripts/run_mv17b_null_qualification.R",
  "scripts/build_mv17b_null_qualification_closure.R",
  "tests/testthat/test-mv17b-null-qualification.R",
  "tests/testthat/test-mv17b-null-qualification-prefreeze.R")
binding <- data.frame(contract_id="mv17b_implementation_binding_v1",
  file=implementations, bytes=as.numeric(file.info(implementations)$size),
  sha256=vapply(implementations,sha,character(1L)))
contract <- data.frame(contract_id="mv17b_null_qualification_prefreeze_v1",
  execution_head=head, source_MV17A_manifest_sha256=sha(file.path(mv17a,"mv17a-artifact-manifest.csv")),
  null_families=4L, fixtures=3L, seeds="17001;17002;17003", points=128L,
  coordinates=3L, invariant_tolerance=1e-10, neighbor_k=8L,
  maximum_circle_neighbor_jaccard=0.80,
  determinism="byte_exact_RDS_digest", independent_repeat_required=TRUE,
  real_corpus_authorized=FALSE, PH_authorized=FALSE, labels_authorized=FALSE,
  outcomes_authorized=FALSE, clustering_authorized=FALSE, fusion_authorized=FALSE,
  biology_authorized=FALSE, manuscript_claims_authorized=FALSE)
validation <- data.frame(contract_id="mv17b_validation_v1",
  check_id=c("MV17A_closed","four_families","family_semantics","fixture_grid",
    "fixed_seeds","fixed_tolerances","implementation_bound","repeat_required",
    "real_corpus_closed","PH_closed","downstream_firewall","zero_null_execution"),
  passed=c(nrow(m)>0,nrow(families)==4L,all(nzchar(families$preserves)&nzchar(families$destroys)),
    nrow(fixtures)==9L,identical(sort(unique(fixtures$seed)),c(17001L,17002L,17003L)),
    contract$invariant_tolerance==1e-10 && contract$neighbor_k==8L,
    all(file.exists(implementations)),contract$independent_repeat_required,
    !contract$real_corpus_authorized,!contract$PH_authorized,
    !any(contract[c("labels_authorized","outcomes_authorized","clustering_authorized",
      "fusion_authorized","biology_authorized","manuscript_claims_authorized")]),TRUE))
if(!all(validation$passed)) stop("MV17-B prefreeze failed")
decision <- data.frame(contract_id="mv17b_decision_v1",
  decision="authorize_synthetic_null_qualification_after_commit",
  synthetic_execution_authorized_after_commit=TRUE, real_corpus_authorized=FALSE,
  next_after_closure="separate_MV17C_calibration_sentinel_prefreeze")
items <- list("mv17b-contract.csv"=contract,"mv17b-null-families.csv"=families,
  "mv17b-fixtures.csv"=fixtures,"mv17b-implementation-binding.csv"=binding,
  "mv17b-validation.csv"=validation,"mv17b-decision.csv"=decision)
for(n in names(items)) atomic(items[[n]],file.path(output,n))
writeLines(c("# MV17-B null-model qualification prefreeze","",
  "Four null families, three synthetic fixtures, three seeds, invariants,",
  "tolerances, and repeat rules are frozen. No null was generated.","",
  "Real-corpus calibration, PH, labels, outcomes, clustering, fusion, biology,",
  "and manuscript claims remain closed."),file.path(output,"MV17B_NULL_QUALIFICATION_PREFREEZE_2026-08-26.md"))
files<-sort(list.files(output)); manifest<-data.frame(contract_id="mv17b_artifact_manifest_v1",
  artifact=files,bytes=as.numeric(file.info(file.path(output,files))$size),
  sha256=vapply(file.path(output,files),sha,character(1L)))
atomic(manifest,file.path(output,"mv17b-artifact-manifest.csv"))
message("Built MV17-B prefreeze; checks=",nrow(validation))
