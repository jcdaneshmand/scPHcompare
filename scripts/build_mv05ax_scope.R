#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=TRUE);if(length(args)!=1L)stop("Usage: build_mv05ax_scope.R OUTPUT_DIR")
out<-normalizePath(args[[1]],mustWork=FALSE);dir.create(out,recursive=TRUE,showWarnings=FALSE)
script_arg<-grep("^--file=",commandArgs(trailingOnly=FALSE),value=TRUE)
if(length(script_arg)!=1L)stop("Unable to resolve the repository root.")
setwd(normalizePath(file.path(dirname(gsub("~+~"," ",sub("^--file=","",script_arg[[1L]]),fixed=TRUE)),".."),mustWork=TRUE))
m<-utils::read.csv("docs/audits/mv04-input-diagram-manifest-2026-08-05.csv",stringsAsFactors=FALSE)
verified<-vapply(seq_len(nrow(m)),function(i){x<-readRDS(m$result_file[[i]]);identical(
 digest::digest(m$result_file[[i]],algo="sha256",file=TRUE),m$result_file_sha256[[i]])&&
 identical(digest::digest(x$diagram,algo="sha256"),m$diagram_sha256[[i]])&&
 identical(x$provenance$diagram_sha256,m$diagram_sha256[[i]])&&isTRUE(x$provenance$scientific_eligible)},logical(1))
if(!all(verified))stop("MV5-AX source verification failed")
ids<-sort(unique(m$stratum_id),method="radix")
scope<-do.call(rbind,lapply(ids,function(id){x<-m[m$stratum_id==id,];n<-nrow(x);adaptive<-grepl("gene_topology",id)&&grepl("^large__",id)
 data.frame(contract_id="mv05ax_scope_v1",stratum_id=id,diagrams=n,pairs=choose(n,2),
 h0_min=min(x$h0_finite_intervals),h0_max=max(x$h0_finite_intervals),h1_min=min(x$h1_finite_intervals),h1_max=max(x$h1_finite_intervals),
 candidate_k=paste(2:min(10,n-1),collapse=";"),planned_route=if(adaptive)"exact_H0_adaptive_H1" else "exact_H0_exact_H1",
 planned_wall_seconds=30+choose(n,2)*if(adaptive)240 else 30, max_pairs=choose(n,2), max_rss_bytes=2*1024^3,workers=1L,stringsAsFactors=FALSE)}))
scope$execution_lane<-ifelse(scope$planned_route=="exact_H0_adaptive_H1","adaptive_concurrent_2","exact_bounded")
utils::write.csv(scope,file.path(out,"mv05ax-complete-scope-2026-08-12.csv"),row.names=FALSE,quote=TRUE)
summary<-data.frame(contract_id="mv05ax_scope_summary_v1",strata=nrow(scope),diagrams=sum(scope$diagrams),pairs=sum(scope$pairs),
 exact_only_pairs=sum(scope$pairs[scope$planned_route=="exact_H0_exact_H1"]),adaptive_h1_pairs=sum(scope$pairs[scope$planned_route=="exact_H0_adaptive_H1"]),
 all_sources_verified=all(verified),max_concurrent_processes=2L,projected_adaptive_sequential_hours=sum(scope$planned_wall_seconds[scope$planned_route=="exact_H0_adaptive_H1"])/3600,
 projected_adaptive_concurrent_wall_hours=max(scope$planned_wall_seconds[scope$planned_route=="exact_H0_adaptive_H1"])/3600,stringsAsFactors=FALSE)
utils::write.csv(summary,file.path(out,"mv05ax-scope-summary-2026-08-12.csv"),row.names=FALSE,quote=TRUE)
decision<-data.frame(contract_id="mv05ax_speed_decision_v1",production_authorized=TRUE,scheduling="two_adaptive_strata_concurrent_internal_workers_one",
 rust_required_before_production=FALSE,rust_candidate_after_production=TRUE,partition_authorized=FALSE,next_sprint="MV5-AY",stringsAsFactors=FALSE)
utils::write.csv(decision,file.path(out,"mv05ax-speed-decision-2026-08-12.csv"),row.names=FALSE,quote=TRUE)
cat("MV5-AX scope:",sum(scope$diagrams),"diagrams",sum(scope$pairs),"pairs\n")
