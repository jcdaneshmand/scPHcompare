#!/usr/bin/env Rscript
options(warn=2)
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=4L)stop("usage: build_mv17c_time_parser_failure_record.R <closure-v1> <primary-time> <repeat-time> <output>",call.=FALSE)
closure<-normalizePath(args[[1]],mustWork=TRUE);primary_time<-normalizePath(args[[2]],mustWork=TRUE);repeat_time<-normalizePath(args[[3]],mustWork=TRUE);output<-args[[4]]
if(dir.exists(output))stop("MV17-C time-parser failure record exists",call.=FALSE)
dir.create(output,recursive=TRUE)
source("R/mv08z_landscape_production.R");source("R/mv17_calibration.R")
r<-.mv08z_read_csv;w<-.mv08z_atomic_csv;h<-.mv08z_sha256_file
manifest<-r(file.path(closure,"mv17c-artifact-manifest.csv"));paths<-file.path(closure,manifest$artifact)
if(!all(as.numeric(file.info(paths)$size)==manifest$bytes)||!all(vapply(paths,h,character(1L))==tolower(manifest$sha256)))stop("MV17-C rejected closure manifest drift",call.=FALSE)
reported<-r(file.path(closure,"mv17c-resource-summary.csv"));correct<-rbind(cbind(mode="primary",mv17c_parse_gnu_time_v1(primary_time)),cbind(mode="repeat",mv17c_parse_gnu_time_v1(repeat_time)))
reported<-reported[match(correct$mode,reported$mode),,drop=FALSE]
failure<-data.frame(contract_id="mv17c_time_parser_failure_v1",mode=correct$mode,receipt_bytes=as.numeric(file.info(c(primary_time,repeat_time))$size),receipt_sha256=vapply(c(primary_time,repeat_time),h,character(1L)),closure_v1_reported_wall_seconds=reported$outer_wall_seconds,independent_wall_seconds=correct$wall_seconds,omitted_seconds=correct$wall_seconds-reported$outer_wall_seconds,scientific_artifact_changed=FALSE,production_rerun_authorized=FALSE,closure_v1_committable=FALSE,stringsAsFactors=FALSE)
whole_minutes<-failure$omitted_seconds/60
if(!all(failure$omitted_seconds>0)||!all(abs(whole_minutes-round(whole_minutes))<=1e-9)||any(failure$independent_wall_seconds<reported$aggregate_child_seconds))stop("MV17-C time-parser failure signature not reproduced",call.=FALSE)
binding<-data.frame(contract_id="mv17c_rejected_closure_binding_v1",artifact=manifest$artifact,bytes=as.numeric(file.info(paths)$size),sha256=vapply(paths,h,character(1L)),stringsAsFactors=FALSE)
w(failure,file.path(output,"mv17c-time-parser-failure.csv"));w(binding,file.path(output,"mv17c-rejected-closure-binding.csv"))
writeLines(c("# MV17-C closure v1 resource-accounting rejection","","The 195-job primary and 65-job repeat completed successfully and remain immutable.","Closure v1 is rejected because its GNU-time parser discarded minute fields from the two outer elapsed-time receipts.","The defect changes only public resource accounting: scientific metrics, empirical tails, queues, and source bindings are unchanged.","No production rerun is authorized. Recovery is limited to a committed parser correction, exact-head recovery prefreeze, and a fresh closure-only v2 root."),file.path(output,"MV17C_TIME_PARSER_FAILURE_2026-08-26.md"))
files<-sort(list.files(output));w(data.frame(contract_id="mv17c_time_parser_failure_manifest_v1",artifact=files,bytes=as.numeric(file.info(file.path(output,files))$size),sha256=vapply(file.path(output,files),h,character(1L))),file.path(output,"mv17c-time-parser-failure-artifact-manifest.csv"))
message("Built MV17-C time-parser failure record; production remains frozen")
