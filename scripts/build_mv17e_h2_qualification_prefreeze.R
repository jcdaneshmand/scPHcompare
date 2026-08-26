#!/usr/bin/env Rscript
options(warn=2); args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=1L) stop("usage: build_mv17e_h2_qualification_prefreeze.R <output>")
o<-args[[1L]]; if(dir.exists(o)) stop("MV17-E prefreeze exists"); dir.create(o,recursive=TRUE)
source("R/mv08z_landscape_production.R"); source("R/mv17_h2_fixtures.R")
h<-.mv08z_sha256_file; w<-.mv08z_atomic_csv; head<-tolower(trimws(Sys.getenv("MV17E_PREFREEZE_HEAD",unset="")))
if(!grepl("^[0-9a-f]{40}$",head)) stop("MV17E_PREFREEZE_HEAD required")
fixtures<-mv17e_fixture_registry_v1(); impl_files<-c("R/mv17_h2_fixtures.R",
 "scripts/build_mv17e_h2_qualification_prefreeze.R","scripts/run_mv17e_h2_qualification.R",
 "scripts/build_mv17e_h2_qualification_closure.R","tests/testthat/test-mv17e-h2-qualification.R",
 "tests/testthat/test-mv17e-h2-prefreeze.R")
impl<-data.frame(contract_id="mv17e_implementation_binding_v1",file=impl_files,
 bytes=as.numeric(file.info(impl_files)$size),sha256=vapply(impl_files,h,character(1L)))
contract<-data.frame(contract_id="mv17e_h2_qualification_prefreeze_v1",execution_head=head,
 primary_engine="ripserr",independent_engine="TDA_GUDHI",coefficient_field=2L,
 maximum_dimension=2L,complete_VR=TRUE,engine_tolerance=1e-5,
 sphere_minimum_H2_persistence=0.10,torus_minimum_H1_features=2L,
 torus_minimum_H2_persistence=0.05,circle_minimum_H1_persistence=0.50,
 circle_maximum_H2_persistence=1e-8,shuffled_parent_strict_attenuation=TRUE,
 repetitions=2L,one_worker=TRUE,timeout_seconds=1800L,RSS_cap_bytes=4294967296,
 primary_repeat_required=TRUE,independent_engine_required=TRUE,real_H2_authorized=FALSE,
 labels_authorized=FALSE,outcomes_authorized=FALSE,clustering_authorized=FALSE,
 fusion_authorized=FALSE,biology_authorized=FALSE,manuscript_claims_authorized=FALSE)
validation<-data.frame(contract_id="mv17e_validation_v1",check_id=c("seven_fixtures",
 "positive_negative_shuffle_classes","bounded_points","dual_engine","field_and_dimension",
 "fixed_tolerances","resource_caps","implementation_bound","repeat_required",
 "real_H2_closed","downstream_firewall","zero_H2_execution"),passed=c(nrow(fixtures)==7L,
 all(c("positive","negative","null_control","attenuated")%in%c(fixtures$expected_H2,"attenuated")),
 max(fixtures$points)<=48L,contract$independent_engine_required,
 contract$coefficient_field==2L&&contract$maximum_dimension==2L,
 contract$engine_tolerance==1e-5&&contract$circle_maximum_H2_persistence==1e-8,
 contract$one_worker&&contract$RSS_cap_bytes==4294967296,all(file.exists(impl_files)),
 contract$primary_repeat_required,!contract$real_H2_authorized,
 !any(contract[c("labels_authorized","outcomes_authorized","clustering_authorized",
 "fusion_authorized","biology_authorized","manuscript_claims_authorized")]),TRUE))
if(!all(validation$passed)) stop("MV17-E prefreeze failed")
decision<-data.frame(contract_id="mv17e_decision_v1",synthetic_execution_authorized_after_commit=TRUE,
 real_H2_authorized=FALSE,next_after_closure="MV17F_synthetic_resource_grid_prefreeze_only")
items<-list("mv17e-contract.csv"=contract,"mv17e-fixtures.csv"=fixtures,
 "mv17e-implementation-binding.csv"=impl,"mv17e-validation.csv"=validation,"mv17e-decision.csv"=decision)
for(n in names(items)) w(items[[n]],file.path(o,n)); writeLines(c("# MV17-E H2 qualification prefreeze","",
 "Sphere, torus, circle, Gaussian, and shuffled controls are frozen before H2 execution.",
 "Ripserr and TDA/GUDHI must agree under fixed qualitative, numerical, and resource gates.","",
 "Real-data H2 and every downstream scientific surface remain closed."),file.path(o,"MV17E_H2_QUALIFICATION_PREFREEZE_2026-08-26.md"))
f<-sort(list.files(o)); w(data.frame(contract_id="mv17e_artifact_manifest_v1",artifact=f,
 bytes=as.numeric(file.info(file.path(o,f))$size),sha256=vapply(file.path(o,f),h,character(1L))),file.path(o,"mv17e-artifact-manifest.csv"))
message("Built MV17-E prefreeze; checks=",nrow(validation))
