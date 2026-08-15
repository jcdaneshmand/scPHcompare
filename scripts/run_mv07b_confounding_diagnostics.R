#!/usr/bin/env Rscript
options(warn=2)
for(p in c("digest","pkgload"))if(!requireNamespace(p,quietly=TRUE))stop(p," required")
pkgload::load_all(".",quiet=TRUE,export_all=TRUE)
a<-commandArgs(TRUE);if(length(a)!=4L)stop("usage: run_mv07b_confounding_diagnostics.R PREFREEZE OUT HEAD METADATA")
pre<-a[[1]];out<-a[[2]];expected<-tolower(a[[3]]);metadata<-a[[4]]
head<-tolower(trimws(system2("git",c("rev-parse","HEAD"),stdout=TRUE)));if(head!=expected)stop("HEAD mismatch")
dir.create(out,recursive=TRUE,showWarnings=FALSE);if(length(list.files(out)))stop("output not empty")
sha<-function(x)digest::digest(file=x,algo="sha256",serialize=FALSE)
wc<-function(x,n)write_provenance_csv(x,file.path(out,n))
source_lock<-read.csv(file.path(pre,"mv07b-source-freeze.csv"),stringsAsFactors=FALSE)
if(any(vapply(seq_len(nrow(source_lock)),function(i)sha(source_lock$locator[i])!=source_lock$sha256[i],TRUE)))stop("source drift")
ml<-read.csv(file.path(pre,"mv07b-metadata-lock.csv"),stringsAsFactors=FALSE)
if(sha(metadata)!=ml$expected_sha256)stop("metadata hash drift")
wc(data.frame(contract_id="mv07b_label_access_receipt_v1",execution_head=head,
  prefreeze_source_sha256=sha(file.path(pre,"mv07b-source-freeze.csv")),
  metadata_expected_sha256=ml$expected_sha256,receipt_written_before_metadata_parse=TRUE,
  new_ph=FALSE,landscape=FALSE,distance=FALSE,ranking=FALSE),"mv07b-label-access-receipt.csv")

x<-read.csv("docs/audits/mv06h-outcome-evidence/mv06h-sample-method-summaries.csv",stringsAsFactors=FALSE)
methods<-mv07b_methods_v1()$method_id; endpoints<-mv07b_endpoints_v1()$endpoint_id
x<-x[x$method_id%in%methods,,drop=FALSE]
raw<-read.csv(metadata,stringsAsFactors=FALSE,check.names=FALSE)
idx<-match(x$query_sample_id,raw$Sample);if(anyNA(idx))stop("metadata join missing")
x$retained_cells<-as.numeric(raw$Number_of_Cells_After_Filtering[idx]);x$approach<-trimws(raw$Approach.x[idx])
if(any(!is.finite(x$retained_cells)|x$retained_cells<=0)|any(x$held_out_study!=raw$SRA[idx])|
   any(tolower(x$query_tissue)!=tolower(raw$Tissue.x[idx])))stop("metadata identity drift")
one<-x[x$method_id==methods[1],]
mixed<-names(which(tapply(one$approach,one$held_out_study,function(z)length(unique(z)))==2L))
design<-data.frame(contract_id="mv07b_design_inventory_v1",samples=nrow(one),studies=length(unique(one$held_out_study)),
  tissues=length(unique(one$query_tissue)),scrna=sum(one$approach=="scRNA-seq"),snrna=sum(one$approach=="snRNA-seq"),
  mixed_approach_studies=length(mixed),min_retained_cells=min(one$retained_cells),
  median_retained_cells=median(one$retained_cells),max_retained_cells=max(one$retained_cells),
  library_size_available=FALSE,cell_type_composition_available=FALSE)
wc(design,"mv07b-design-inventory.csv")

full<-do.call(rbind,lapply(methods,function(m)do.call(rbind,lapply(endpoints,function(e){z<-x[x$method_id==m,];data.frame(method_id=m,endpoint_id=e,estimate=mv07b_macro(z,e))}))))
study_ids<-sort(unique(one$held_out_study)); tissue_ids<-sort(unique(one$query_tissue))
delstudy<-do.call(rbind,lapply(methods,function(m)do.call(rbind,lapply(endpoints,function(e)do.call(rbind,lapply(study_ids,function(id){z<-x[x$method_id==m & x$held_out_study!=id,];est<-mv07b_macro(z,e);base<-full$estimate[full$method_id==m&full$endpoint_id==e];data.frame(method_id=m,endpoint_id=e,deleted_study=id,estimate=est,change_from_full=est-base)}))))))
deltissue<-do.call(rbind,lapply(methods,function(m)do.call(rbind,lapply(endpoints,function(e)do.call(rbind,lapply(tissue_ids,function(id){z<-x[x$method_id==m & x$query_tissue!=id,];means<-tapply(z[[e]],z$query_tissue,mean);est<-mean(means);base<-full$estimate[full$method_id==m&full$endpoint_id==e];data.frame(method_id=m,endpoint_id=e,deleted_tissue=id,estimate=est,change_from_full=est-base)}))))))
wc(full,"mv07b-full-summaries.csv");wc(delstudy,"mv07b-delete-study.csv");wc(deltissue,"mv07b-delete-tissue.csv")

contr<-mv07b_contrasts_v1(); ci<-list()
for(i in seq_len(nrow(contr)))for(e in endpoints){l<-contr$left_method[i];r<-contr$right_method[i];base<-full$estimate[full$method_id==l&full$endpoint_id==e]-full$estimate[full$method_id==r&full$endpoint_id==e]
 for(id in study_ids){v<-(delstudy$estimate[delstudy$method_id==l&delstudy$endpoint_id==e&delstudy$deleted_study==id]-delstudy$estimate[delstudy$method_id==r&delstudy$endpoint_id==e&delstudy$deleted_study==id]);ci[[length(ci)+1L]]<-data.frame(contrast_id=contr$contrast_id[i],endpoint_id=e,deletion_axis="study",deleted_id=id,full_estimate=base,deletion_estimate=v,sign_changed=(sign(v)!=sign(base)))}
 for(id in tissue_ids){v<-(deltissue$estimate[deltissue$method_id==l&deltissue$endpoint_id==e&deltissue$deleted_tissue==id]-deltissue$estimate[deltissue$method_id==r&deltissue$endpoint_id==e&deltissue$deleted_tissue==id]);ci[[length(ci)+1L]]<-data.frame(contrast_id=contr$contrast_id[i],endpoint_id=e,deletion_axis="tissue",deleted_id=id,full_estimate=base,deletion_estimate=v,sign_changed=(sign(v)!=sign(base)))}}
ci<-do.call(rbind,ci);wc(ci,"mv07b-contrast-influence.csv")

set.seed(20260817L);counts<-matrix(sample.int(length(study_ids),2000L*length(study_ids),replace=TRUE),nrow=2000L)
countrows<-list()
for(m in methods)for(e in endpoints){z<-x[x$method_id==m,];point<-mv07b_within_study_rank_correlation(z,e);boot<-apply(counts,1,function(ii)mv07b_within_study_rank_correlation(do.call(rbind,lapply(study_ids[ii],function(id)z[z$held_out_study==id,])),e));q<-mv07b_percentile(boot);countrows[[length(countrows)+1L]]<-data.frame(method_id=m,endpoint_id=e,estimate=point,ci_lower=q[1],ci_upper=q[2],replicates=2000L,seed=20260817L)}
countrows<-do.call(rbind,countrows);wc(countrows,"mv07b-retained-cell-association.csv")

arows<-list();for(m in methods)for(e in endpoints)for(id in mixed){z<-x[x$method_id==m&x$held_out_study==id,];arows[[length(arows)+1L]]<-data.frame(method_id=m,endpoint_id=e,study_id=id,snrna_minus_scrna=mean(z[[e]][z$approach=="snRNA-seq"])-mean(z[[e]][z$approach=="scRNA-seq"]))}
arows<-do.call(rbind,arows);set.seed(20260818L);ai<-matrix(sample.int(length(mixed),2000L*length(mixed),replace=TRUE),nrow=2000L)
asum<-do.call(rbind,lapply(methods,function(m)do.call(rbind,lapply(endpoints,function(e){v<-arows$snrna_minus_scrna[arows$method_id==m&arows$endpoint_id==e];boot<-apply(ai,1,function(ii)mean(v[ii]));q<-mv07b_percentile(boot);data.frame(method_id=m,endpoint_id=e,estimate=mean(v),ci_lower=q[1],ci_upper=q[2],mixed_studies=3L,replicates=2000L,seed=20260818L)}))))
wc(arows,"mv07b-mixed-approach-study-effects.csv");wc(asum,"mv07b-mixed-approach-summaries.csv")

studyflag<-aggregate(abs(delstudy$change_from_full),list(delstudy$method_id,delstudy$endpoint_id),max);names(studyflag)<-c("method_id","endpoint_id","max_abs_change")
tissueflag<-aggregate(abs(deltissue$change_from_full),list(deltissue$method_id,deltissue$endpoint_id),max);names(tissueflag)<-c("method_id","endpoint_id","max_abs_change")
flags<-data.frame(flag_id=c("large_study_influence","large_tissue_influence","contrast_sign_instability","retained_cell_association","approach_association"),
 triggered=c(any(studyflag$max_abs_change>=.05),any(tissueflag$max_abs_change>=.05),any(ci$sign_changed),
 any(abs(countrows$estimate)>=.3 & countrows$ci_lower*countrows$ci_upper>0),any(abs(asum$estimate)>=.1 & asum$ci_lower*asum$ci_upper>0)),
 consequence="narrow_claims_not_select_methods")
wc(flags,"mv07b-diagnostic-flags.csv")
disp<-if(flags$triggered[1]||flags$triggered[3]||flags$triggered[4])"narrow_claims_before_gene_rerun_or_new_data" else if(any(flags$triggered))"retain_secondary_claims_with_confounding_caveats" else "diagnostics_stable_pending_synthesis"
wc(data.frame(disposition=disp,triggered_flags=sum(flags$triggered),new_ph_authorized=FALSE,new_data_authorized=FALSE,method_selection_authorized=FALSE,fusion_authorized=FALSE,causal_technology_claim_authorized=FALSE,next_sprint="MV7-C_existing_data_synthesis"),"mv07b-decision.csv")
wc(data.frame(sample_rows=nrow(x),methods=6L,endpoints=2L,delete_study_rows=nrow(delstudy),delete_tissue_rows=nrow(deltissue),contrast_rows=nrow(ci),count_rows=nrow(countrows),approach_study_rows=nrow(arows),approach_summary_rows=nrow(asum),p_values=0L,new_ph_jobs=0L),"mv07b-production-summary.csv")
message("MV7-B diagnostics complete")
