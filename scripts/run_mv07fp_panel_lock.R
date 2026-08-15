#!/usr/bin/env Rscript

Sys.setenv(OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1")
options(warn=2)
for(p in c("digest","Matrix","ps"))if(!requireNamespace(p,quietly=TRUE))stop(p," required.")
args<-commandArgs(trailingOnly=TRUE);if(length(args)!=6L)stop("usage: run_mv07fp_panel_lock.R PREFREEZE PRIMARY_CACHE ADDED_CACHE OLD_PANEL OUTPUT EXPECTED_HEAD")
prefreeze<-args[[1]];primary_cache<-args[[2]];added_cache<-args[[3]];old_panel_path<-args[[4]];output<-args[[5]];expected_head<-tolower(trimws(args[[6]]))
head<-tolower(trimws(system2("git",c("rev-parse","HEAD"),stdout=TRUE)));if(!identical(head,expected_head))stop("MV7-FP exact HEAD mismatch.")
if(dir.exists(output)&&length(list.files(output,all.files=TRUE,no..=TRUE)))stop("Output must be empty.")
source("R/provenance_utils.R");source("R/toy_baseline.R");source("R/dual_view_topology.R");source("R/mv05_resource_safe_execution.R");source("R/mv06_global_core_gate.R")
sha<-function(p)digest::digest(file=p,algo="sha256",serialize=FALSE);rss<-function()as.numeric(ps::ps_memory_info(ps::ps_handle())[["rss"]])
manifest<-read.csv(file.path(prefreeze,"mv07fp-cache-manifest.csv"),stringsAsFactors=FALSE,check.names=FALSE);caps<-read.csv(file.path(prefreeze,"mv07fp-resource-caps.csv"),stringsAsFactors=FALSE)
freeze<-read.csv(file.path(prefreeze,"mv07fp-source-freeze.csv"),stringsAsFactors=FALSE,check.names=FALSE)
public<-freeze[!as.logical(freeze$private_source),,drop=FALSE];if(any(vapply(seq_len(nrow(public)),function(i)!file.exists(public$artifact_locator[i])||sha(public$artifact_locator[i])!=public$sha256[i],logical(1))))stop("MV7-FP public source drift.")
accepted_panel<-freeze[freeze$source_id=="accepted_panel",,drop=FALSE];if(nrow(accepted_panel)!=1L||sha(old_panel_path)!=accepted_panel$sha256)stop("MV7-FP accepted panel drift.")
paths<-ifelse(manifest$source_tier=="primary90",file.path(primary_cache,manifest$private_cache_file),file.path(added_cache,manifest$private_cache_file))
started<-proc.time()[["elapsed"]];peak<-rss();base<-NULL;presence<-NULL;variances<-NULL;union_features<-character()
for(i in seq_len(nrow(manifest))){
 if(!file.exists(paths[i])||sha(paths[i])!=manifest$private_cache_sha256[i])stop("MV7-FP cache drift at row ",i)
 record<-readRDS(paths[i]);mv05d0_validate_normalization_cache_record_v2(record);value<-mv05d0_sct_matrix_from_cache_v1(record)
 if(record$identity$sample_id!=manifest$sample_id[i]||record$identity$seed!=manifest$seed[i]||record$cache_key!=manifest$normalization_cache_key[i]||record$payload_sha256!=manifest$payload_sha256[i]||ncol(value)!=384L||any(!is.finite(value@x)))stop("MV7-FP cache content drift at row ",i)
 features<-rownames(value);union_features<-union(union_features,features)
 if(is.null(base)){base<-features;presence<-integer(length(base));variances<-matrix(NA_real_,nrow=length(base),ncol=nrow(manifest),dimnames=list(base,NULL))}
 matched<-match(base,features);present<-which(!is.na(matched));presence[present]<-presence[present]+1L;variances[present,i]<-.mv03_row_variance(value[matched[present],,drop=FALSE])
 peak<-max(peak,rss());rm(record,value);if(i%%25L==0L)invisible(gc())
}
common<-sort(base[presence==nrow(manifest)],method="radix");variances<-variances[match(common,base),,drop=FALSE]
result<-mv06c_build_global_core_panel_v1(common,variances,as.integer(manifest$seed),panel_size=500L,samples_per_seed=124L)
elapsed<-proc.time()[["elapsed"]]-started
panel<-result$panel[c("panel_order","feature_id","gene","category","median_variance_rank","minimum_variance","finite_cache_count","positive_cache_count","selected_all_cache_core")];panel$contract_id<-"mv07fp_global_core_panel_v1";panel$panel_sha256<-result$panel_sha256;panel<-panel[c("contract_id","panel_sha256",setdiff(names(panel),c("contract_id","panel_sha256")))]
old<-read.csv(old_panel_path,stringsAsFactors=FALSE,check.names=FALSE);overlap<-length(intersect(panel$feature_id,old$feature_id));comparison<-data.frame(contract_id="mv07fp_panel_comparison_v1",old_panel_sha256=unique(old$panel_sha256),new_panel_sha256=result$panel_sha256,old_panel_size=nrow(old),new_panel_size=nrow(panel),feature_overlap=overlap,feature_jaccard=overlap/length(union(panel$feature_id,old$feature_id)),same_order=identical(panel$feature_id,old$feature_id),selection_basis="availability_and_within_cache_variance_only",stringsAsFactors=FALSE)
summary<-data.frame(contract_id="mv07fp_source_summary_v1",records_verified=620L,biological_samples=124L,seeds=5L,cells_per_record=384L,union_features=length(union_features),common_features=length(common),source_bytes=sum(manifest$private_cache_bytes),outcome_label_state="closed",biological_outcomes_computed=FALSE,stringsAsFactors=FALSE)
resource<-data.frame(contract_id="mv07fp_resource_v1",elapsed_seconds=elapsed,peak_process_rss_bytes=peak,source_bytes=sum(manifest$private_cache_bytes),elapsed_cap_seconds=caps$elapsed_cap_seconds,rss_cap_bytes=caps$rss_cap_bytes,elapsed_cap_pass=elapsed<=caps$elapsed_cap_seconds,rss_cap_pass=peak<=caps$rss_cap_bytes,stringsAsFactors=FALSE)
decision<-data.frame(contract_id="mv07fp_panel_decision_v1",decision=if(result$decision=="go_bounded_matched_sct_profile"&&nrow(panel)==500L&&all(panel$finite_cache_count==620L)&all(panel$positive_cache_count==620L)&&resource$elapsed_cap_pass&&resource$rss_cap_pass)"lock_exact_124_panel_and_authorize_MV7G_prefreeze" else "stop_panel_lock_failure",panel_sha256=result$panel_sha256,selected_panel_size=nrow(panel),eligible_unique_canonical_genes=result$eligibility$eligible_unique_canonical_genes,eligibility_margin=result$eligibility$eligibility_margin,pca_jobs=0L,ph_jobs=0L,landscape_jobs=0L,label_jobs=0L,outcome_jobs=0L,stringsAsFactors=FALSE)
dir.create(output,recursive=TRUE,showWarnings=FALSE);writep<-function(x,n){p<-file.path(output,n);write_provenance_csv(x,p);p}
det<-c(writep(summary,"mv07fp-source-summary.csv"),writep(result$eligibility,"mv07fp-eligibility.csv"),writep(panel,"mv07fp-panel.csv"),writep(result$seed_stability,"mv07fp-seed-stability.csv"),writep(comparison,"mv07fp-panel-comparison.csv"),writep(decision,"mv07fp-decision.csv"))
artifacts<-data.frame(contract_id="mv07fp_artifact_manifest_v1",file=basename(det),bytes=as.numeric(file.info(det)$size),sha256=vapply(det,sha,character(1)),stringsAsFactors=FALSE);writep(artifacts,"mv07fp-artifact-manifest.csv");writep(resource,"mv07fp-resource.csv")
if(decision$decision!="lock_exact_124_panel_and_authorize_MV7G_prefreeze")stop("MV7-FP panel lock failed.")
message("MV7-FP panel locked: ",result$panel_sha256)
