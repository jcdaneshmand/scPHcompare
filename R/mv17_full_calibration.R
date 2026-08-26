# MV17-G full, label-blind H0/H1 calibration production-design contracts.

mv17g_family_registry_v1 <- function() {
  data.frame(view=c(rep("cell",3L),rep("gene",4L)),family_order=c(1:3,1:4),null_family=c("coordinate_permutation","covariance_gaussian","radial_density_cloud","coordinate_permutation","covariance_gaussian","radial_density_cloud","within_row_axis_shuffle"),stringsAsFactors=FALSE)
}

mv17g_seed_registry_v1 <- function(units_per_view=132L,replicates=99L) {
  units_per_view<-as.integer(units_per_view);replicates<-as.integer(replicates)
  if(length(units_per_view)!=1L||is.na(units_per_view)||units_per_view<1L||length(replicates)!=1L||is.na(replicates)||replicates<1L)stop("positive MV17-G units and replicates required",call.=FALSE)
  families<-mv17g_family_registry_v1();rows<-lapply(seq_len(nrow(families)),function(i){f<-families[i,];x<-expand.grid(unit_order=seq_len(units_per_view),replicate=seq_len(replicates),KEEP.OUT.ATTRS=FALSE,stringsAsFactors=FALSE);x$view<-f$view;x$family_order<-f$family_order;x$null_family<-f$null_family;x})
  out<-do.call(rbind,rows);out$view_order<-match(out$view,c("cell","gene"));out<-out[order(out$view_order,out$unit_order,out$family_order,out$replicate,method="radix"),,drop=FALSE];out$seed<-174000L+seq_len(nrow(out));out$contract_id<-"mv17g_null_seed_v1";rownames(out)<-NULL;out[c("contract_id","view","unit_order","null_family","replicate","seed")]
}

mv17g_grouped_queue_v1 <- function(mode=c("primary","repeat"),repeat_orders=NULL,units_per_view=132L,replicates=99L) {
  mode<-match.arg(mode);seeds<-mv17g_seed_registry_v1(units_per_view,replicates)
  if(mode=="repeat"){
    required<-c("view","unit_order","repeat_role");if(is.null(repeat_orders)||!all(required%in%names(repeat_orders))||nrow(repeat_orders)!=6L||anyDuplicated(repeat_orders[c("view","unit_order")])||!setequal(repeat_orders$view,c("cell","gene"))||!setequal(repeat_orders$repeat_role,c("minimum","median","maximum")))stop("six canonical MV17-G repeat orders required",call.=FALSE)
    keep<-paste(seeds$view,seeds$unit_order)%in%paste(repeat_orders$view,repeat_orders$unit_order);seeds<-seeds[keep,,drop=FALSE];units<-repeat_orders
  }else units<-expand.grid(view=c("cell","gene"),unit_order=seq_len(as.integer(units_per_view)),KEEP.OUT.ATTRS=FALSE,stringsAsFactors=FALSE)
  key<-unique(seeds[c("view","unit_order","null_family")]);null_groups<-do.call(rbind,lapply(seq_len(nrow(key)),function(i){z<-key[i,];take<-seeds$view==z$view&seeds$unit_order==z$unit_order&seeds$null_family==z$null_family;s<-seeds[take,,drop=FALSE];data.frame(view=z$view,unit_order=z$unit_order,null_family=z$null_family,replicate_count=nrow(s),seed_first=min(s$seed),seed_last=max(s$seed),stringsAsFactors=FALSE)}))
  observed<-data.frame(view=units$view,unit_order=units$unit_order,null_family="observed",replicate_count=1L,seed_first=0L,seed_last=0L,stringsAsFactors=FALSE)
  queue<-rbind(observed,null_groups);queue$view_order<-match(queue$view,c("cell","gene"));queue$family_order<-match(queue$null_family,c("observed","coordinate_permutation","covariance_gaussian","radial_density_cloud","within_row_axis_shuffle"));queue<-queue[order(queue$view_order,queue$unit_order,queue$family_order,method="radix"),,drop=FALSE];queue$job_order<-seq_len(nrow(queue));queue$mode<-mode;queue$scientific_runs<-queue$replicate_count;queue$contract_id<-"mv17g_grouped_queue_v1";rownames(queue)<-NULL;queue[c("contract_id","mode","job_order","view","unit_order","null_family","replicate_count","seed_first","seed_last","scientific_runs")]
}

mv17g_resource_projection_v1 <- function(mv17c_resource_summary,repeat_orders=data.frame(view=rep(c("cell","gene"),each=3L),unit_order=rep(c(1L,66L,132L),2L),repeat_role=rep(c("minimum","median","maximum"),2L))) {
  required<-c("mode","jobs","outer_wall_seconds","private_bytes");if(!all(required%in%names(mv17c_resource_summary))||!setequal(mv17c_resource_summary$mode,c("primary","repeat")))stop("canonical MV17-C resource summary required",call.=FALSE)
  primary<-mv17g_grouped_queue_v1("primary");repeat_queue<-mv17g_grouped_queue_v1("repeat",repeat_orders=repeat_orders);seconds_per_run<-max(mv17c_resource_summary$outer_wall_seconds/mv17c_resource_summary$jobs);private_bytes_per_run<-max(mv17c_resource_summary$private_bytes/mv17c_resource_summary$jobs);runs<-c(sum(primary$scientific_runs),sum(repeat_queue$scientific_runs));children<-c(nrow(primary),nrow(repeat_queue));projected_seconds<-runs*seconds_per_run;data.frame(contract_id="mv17g_resource_projection_v1",mode=c("primary","repeat"),grouped_children=children,scientific_runs=runs,conservative_seconds_per_run=seconds_per_run,projected_wall_seconds=projected_seconds,projected_wall_seconds_with_25pct_margin=projected_seconds*1.25,projected_private_bytes_with_25pct_margin=runs*private_bytes_per_run*1.25,child_timeout_seconds=1800L,aggregate_timeout_seconds=604800L,child_RSS_cap_bytes=8589934592,private_cap_bytes=12884901888,public_cap_bytes=67108864,workers=1L,retries=0L,stringsAsFactors=FALSE)
}

mv17g_design_contract_v1 <- function() {
  data.frame(contract_id="mv17g_design_v1",panel_id="exact500",seed=20260805L,units_per_view=132L,replicates=99L,minimum_attainable_probability=.01,worst_case_MCSE=.05,primary_grouped_children=1188L,primary_scientific_runs=91740L,repeat_grouped_children=27L,repeat_scientific_runs=2085L,H0_H1_separate=TRUE,summaries=4L,workers=1L,retries=0L,full_calibration_authorized=FALSE,real_localization_authorized=FALSE,labels_authorized=FALSE,outcomes_authorized=FALSE,clustering_authorized=FALSE,fusion_authorized=FALSE,biology_authorized=FALSE,manuscript_claims_authorized=FALSE,stringsAsFactors=FALSE)
}
