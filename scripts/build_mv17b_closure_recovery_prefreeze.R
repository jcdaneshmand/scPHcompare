#!/usr/bin/env Rscript
options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop(paste("usage: build_mv17b_closure_recovery_prefreeze.R",
  "<primary> <repeat> <rejected-closure> <output>"), call. = FALSE)
primary <- normalizePath(args[[1L]], mustWork = TRUE)
repeat_root <- normalizePath(args[[2L]], mustWork = TRUE)
rejected <- normalizePath(args[[3L]], mustWork = TRUE); output <- args[[4L]]
if (dir.exists(output)) stop("MV17-B recovery prefreeze exists", call. = FALSE)
dir.create(output, recursive = TRUE)
source("R/mv08z_landscape_production.R"); sha <- .mv08z_sha256_file; atomic <- .mv08z_atomic_csv
head <- tolower(trimws(Sys.getenv("MV17B_RECOVERY_HEAD", unset = "")))
if (!grepl("^[0-9a-f]{40}$", head)) stop("MV17B_RECOVERY_HEAD required")
files <- c(file.path(primary,c("mv17b-qualification.csv","mv17b-status.csv")),
  file.path(repeat_root,c("mv17b-qualification.csv","mv17b-status.csv")))
binding <- data.frame(contract_id="mv17b_recovery_production_binding_v1",
  role=c("primary_metrics","primary_status","repeat_metrics","repeat_status"),
  bytes=as.numeric(file.info(files)$size),sha256=vapply(files,sha,character(1L)))
implementation <- c("R/mv17_null_models.R","scripts/run_mv17b_null_qualification.R",
  "scripts/build_mv17b_null_qualification_prefreeze.R",
  "scripts/build_mv17b_null_qualification_closure.R")
impl <- data.frame(contract_id="mv17b_recovery_implementation_binding_v1",file=implementation,
  bytes=as.numeric(file.info(implementation)$size),sha256=vapply(implementation,sha,character(1L)))
contract <- data.frame(contract_id="mv17b_closure_recovery_prefreeze_v1",
  execution_head=head,primary_repeat_byte_identical=binding$sha256[[1L]]==binding$sha256[[3L]] && binding$sha256[[2L]]==binding$sha256[[4L]],
  production_rerun_authorized=FALSE,closure_v2_authorized_after_commit=TRUE,
  complete_metrics_required=TRUE,production_hash_binding_required=TRUE,
  real_corpus_authorized=FALSE,PH_authorized=FALSE,clustering_authorized=FALSE,
  fusion_authorized=FALSE,biology_authorized=FALSE,manuscript_claims_authorized=FALSE)
rejected_files <- sort(list.files(rejected)); rejected_binding <- data.frame(
  contract_id="mv17b_rejected_closure_binding_v1",artifact=rejected_files,
  bytes=as.numeric(file.info(file.path(rejected,rejected_files))$size),
  sha256=vapply(file.path(rejected,rejected_files),sha,character(1L)))
metrics <- read.csv(files[[1L]],stringsAsFactors=FALSE,check.names=FALSE)
validation <- data.frame(contract_id="mv17b_recovery_validation_v1",
  check_id=c("four_production_artifacts","primary_repeat_identical","36_jobs",
  "all_invariants","all_deterministic","rejected_closure_bound","implementation_bound",
  "production_frozen","closure_only","firewall"),passed=c(nrow(binding)==4L,
  contract$primary_repeat_byte_identical,nrow(metrics)==36L,all(metrics$passed_invariant),
  all(metrics$deterministic),nrow(rejected_binding)==3L,nrow(impl)==4L,
  !contract$production_rerun_authorized,contract$closure_v2_authorized_after_commit,
  !contract$real_corpus_authorized && !contract$PH_authorized && !contract$clustering_authorized &&
  !contract$fusion_authorized && !contract$biology_authorized && !contract$manuscript_claims_authorized))
if(!all(validation$passed)) stop("MV17-B recovery prefreeze failed")
items<-list("mv17b-recovery-contract.csv"=contract,"mv17b-recovery-production-binding.csv"=binding,
  "mv17b-recovery-implementation-binding.csv"=impl,"mv17b-rejected-closure-binding.csv"=rejected_binding,
  "mv17b-recovery-validation.csv"=validation)
for(name in names(items)) atomic(items[[name]],file.path(output,name))
writeLines(c("# MV17-B closure recovery prefreeze","",
  "The successful primary and independent repeat are frozen byte-for-byte.",
  "Production rerun is forbidden; only the provenance-complete v2 closure is authorized."),
  file.path(output,"MV17B_CLOSURE_RECOVERY_PREFREEZE_2026-08-26.md"))
out<-sort(list.files(output)); manifest<-data.frame(contract_id="mv17b_recovery_manifest_v1",
  artifact=out,bytes=as.numeric(file.info(file.path(output,out))$size),
  sha256=vapply(file.path(output,out),sha,character(1L)))
atomic(manifest,file.path(output,"mv17b-recovery-artifact-manifest.csv"))
message("Built MV17-B closure recovery prefreeze; checks=",nrow(validation))
