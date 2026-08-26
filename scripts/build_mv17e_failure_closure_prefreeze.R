#!/usr/bin/env Rscript
options(warn=2); args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=5L) stop("usage: build_mv17e_failure_closure_prefreeze.R <primary> <repeat> <primary-time> <repeat-time> <output>")
a<-normalizePath(args[[1]],mustWork=TRUE); b<-normalizePath(args[[2]],mustWork=TRUE); ta<-normalizePath(args[[3]],mustWork=TRUE); tb<-normalizePath(args[[4]],mustWork=TRUE); o<-args[[5]]
if(dir.exists(o)) stop("MV17-E failure prefreeze exists"); dir.create(o,recursive=TRUE)
source("R/mv08z_landscape_production.R"); h<-.mv08z_sha256_file; w<-.mv08z_atomic_csv
head<-tolower(trimws(Sys.getenv("MV17E_FAILURE_HEAD",unset=""))); if(!grepl("^[0-9a-f]{40}$",head)) stop("MV17E_FAILURE_HEAD required")
names<-c("sphere.rds","torus.rds","circle.rds","gaussian_cloud.rds","shuffled_sphere.rds","shuffled_torus.rds","shuffled_circle.rds","mv17e-summary.csv","mv17e-status.csv")
paths<-c(file.path(a,names),file.path(b,names),ta,tb); roles<-c(paste0("primary_",names),paste0("repeat_",names),"primary_GNU_time","repeat_GNU_time")
binding<-data.frame(contract_id="mv17e_failure_source_binding_v1",role=roles,bytes=as.numeric(file.info(paths)$size),sha256=vapply(paths,h,character(1L)))
payload_repeat<-all(binding$sha256[1:7]==binding$sha256[10:16]); status_repeat<-binding$sha256[[9]]==binding$sha256[[18]]
impls<-c("R/mv17_h2_fixtures.R","scripts/run_mv17e_h2_qualification.R","scripts/build_mv17e_failure_closure.R")
impl<-data.frame(contract_id="mv17e_failure_implementation_binding_v1",file=impls,bytes=as.numeric(file.info(impls)$size),sha256=vapply(impls,h,character(1L)))
contract<-data.frame(contract_id="mv17e_failure_closure_prefreeze_v1",execution_head=head,
 payload_repeat_exact=payload_repeat,status_repeat_exact=status_repeat,production_rerun_authorized=FALSE,
 thresholds_unchanged=TRUE,failure_closure_authorized_after_commit=TRUE,real_H2_authorized=FALSE,
 MV17F_authorized=FALSE,labels_authorized=FALSE,outcomes_authorized=FALSE,
 clustering_authorized=FALSE,fusion_authorized=FALSE,biology_authorized=FALSE,manuscript_claims_authorized=FALSE)
v<-data.frame(contract_id="mv17e_failure_prefreeze_validation_v1",check_id=c("20_sources_bound","payload_repeat_exact","status_repeat_exact","implementation_bound","production_frozen","thresholds_unchanged","failure_only","real_H2_closed","MV17F_closed","downstream_firewall"),passed=c(nrow(binding)==20L,payload_repeat,status_repeat,nrow(impl)==3L,!contract$production_rerun_authorized,contract$thresholds_unchanged,contract$failure_closure_authorized_after_commit,!contract$real_H2_authorized,!contract$MV17F_authorized,!any(contract[c("labels_authorized","outcomes_authorized","clustering_authorized","fusion_authorized","biology_authorized","manuscript_claims_authorized")])))
if(!all(v$passed)) stop("MV17-E failure prefreeze failed")
items<-list("mv17e-failure-contract.csv"=contract,"mv17e-failure-source-binding.csv"=binding,"mv17e-failure-implementation-binding.csv"=impl,"mv17e-failure-validation.csv"=v); for(n in names(items))w(items[[n]],file.path(o,n))
writeLines(c("# MV17-E failure-closure prefreeze","","All successful synthetic artifacts are frozen. No rerun or threshold change is authorized.","The closure must retain failed circle controls and close real-data H2 plus MV17-F."),file.path(o,"MV17E_FAILURE_CLOSURE_PREFREEZE_2026-08-26.md")); f<-sort(list.files(o)); w(data.frame(contract_id="mv17e_failure_prefreeze_manifest_v1",artifact=f,bytes=as.numeric(file.info(file.path(o,f))$size),sha256=vapply(file.path(o,f),h,character(1L))),file.path(o,"mv17e-failure-artifact-manifest.csv"))
