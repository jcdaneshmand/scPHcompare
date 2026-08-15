#!/usr/bin/env Rscript
options(warn = 2)
for (p in c("digest", "pkgload")) if (!requireNamespace(p, quietly=TRUE)) stop(p," required")
pkgload::load_all(".", quiet=TRUE, export_all=TRUE)
a <- commandArgs(TRUE); if(length(a)!=3L) stop("usage: build_mv07b_prefreeze.R OUT HEAD METADATA")
out <- a[[1]]; expected <- tolower(a[[2]]); metadata <- a[[3]]
head <- tolower(trimws(system2("git",c("rev-parse","HEAD"),stdout=TRUE)))
if(!identical(head,expected)) stop("HEAD mismatch")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
sha <- function(x) digest::digest(file=x,algo="sha256",serialize=FALSE)
sources <- c(sample_summaries="docs/audits/mv06h-outcome-evidence/mv06h-sample-method-summaries.csv",
  mv07a_decision="docs/audits/mv07a-synthesis-evidence/mv07a-decision.csv",
  helper="R/mv07b_confounding_diagnostics.R", runner="scripts/run_mv07b_confounding_diagnostics.R",
  builder="scripts/build_mv07b_prefreeze.R", validator="scripts/validate_mv07b_prefreeze.R",
  outcome_validator="scripts/validate_mv07b_confounding_diagnostics.R",
  spec="docs/specifications/MV07B_NO_NEW_PH_CONFOUNDING_DIAGNOSTICS_PREFREEZE_V1.md",
  tests="tests/testthat/test-mv07b-confounding-diagnostics.R")
if(any(!file.exists(sources))||!file.exists(metadata)) stop("source missing")
writep <- function(x,n) write_provenance_csv(x,file.path(out,n))
writep(data.frame(source_id=names(sources),locator=unname(sources),sha256=vapply(sources,sha,""),
  bytes=as.numeric(file.info(sources)$size),implementation_head=expected,stringsAsFactors=FALSE),"mv07b-source-freeze.csv")
writep(data.frame(locator=metadata,sha256=sha(metadata),bytes=as.numeric(file.info(metadata)$size),
  expected_sha256="e2b6f1e7d7dc05ba8e80627aad4ea1677c3b718a844397acd06596d7e22a18d0",
  fields="Sample;SRA;Tissue.x;Approach.x;Number_of_Cells_After_Filtering",
  contents_parsed=FALSE,stringsAsFactors=FALSE),"mv07b-metadata-lock.csv")
writep(mv07b_methods_v1(),"mv07b-methods.csv"); writep(mv07b_endpoints_v1(),"mv07b-endpoints.csv")
writep(mv07b_contrasts_v1(),"mv07b-contrasts.csv")
writep(data.frame(replicates=c(2000L,2000L),seed=c(20260817L,20260818L),
  diagnostic=c("retained_cell_study_block","mixed_approach_study_block"),
  p_values=FALSE,stringsAsFactors=FALSE),"mv07b-resampling.csv")
writep(data.frame(flag=c("delete_change","contrast_sign","cell_count","approach"),
  threshold=c(0.05,0,0.30,0.10),requires_interval_excludes_zero=c(FALSE,FALSE,TRUE,TRUE),
  consequence="narrow_claims_not_select_methods",stringsAsFactors=FALSE),"mv07b-flags.csv")
writep(data.frame(new_ph=FALSE,landscape=FALSE,distance=FALSE,ranking=FALSE,clustering=FALSE,
  fusion=FALSE,method_selection=FALSE,new_data=FALSE,default_change=FALSE,
  claim_promotion=FALSE,metadata_read=FALSE,outcomes_computed=FALSE),"mv07b-prefreeze-decision.csv")
message("MV7-B prefreeze built")
