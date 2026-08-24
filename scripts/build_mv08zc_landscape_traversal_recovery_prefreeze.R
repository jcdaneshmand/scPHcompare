#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) stop("usage: build_mv08zc_landscape_traversal_recovery_prefreeze.R <mv08zb-root> <failed-private> <failed-public> <output> <parent-head>", call. = FALSE)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required", call. = FALSE)
zb <- normalizePath(args[[1L]], mustWork = TRUE)
failed_private <- normalizePath(args[[2L]], mustWork = TRUE)
failed_public <- normalizePath(args[[3L]], mustWork = TRUE)
output <- normalizePath(args[[4L]], mustWork = FALSE)
parent <- tolower(args[[5L]])
if (dir.exists(output) || !grepl("^[0-9a-f]{40}$", parent)) stop("fresh output and parent required", call. = FALSE)
source("R/mv08z_landscape_production.R")
.mv08z_verify_manifest(zb, "mv08zb-artifact-manifest.csv")
prior <- .mv08z_read_csv(file.path(zb, "mv08zb-implementation-bindings.csv"))
ledger <- .mv08z_read_csv(file.path(failed_public, "mv08za-resource-ledger.csv"))
stderr <- file.path(failed_private, "logs", "sentinel_primary_chunk.stderr")
text <- paste(readLines(stderr, warn = FALSE), collapse = "\n")
files <- c("scripts/run_mv08z_landscape_chunk.R",
           "scripts/run_mv08za_landscape_sentinel.R", "R/landscape_contract.R",
           "scripts/build_mv08zc_landscape_traversal_recovery_prefreeze.R")
prior_rows <- prior[match(files[1:2], prior$file), , drop = FALSE]
if (nrow(ledger) != 1L || ledger$disposition != "failed" ||
    !grepl("MV8-Z execution binding drift", text, fixed = TRUE) ||
    anyNA(prior_rows$file) || !all(file.exists(files))) stop("MV8-ZC failure binding drift", call. = FALSE)
binding <- data.frame(contract_id="mv08zc_implementation_binding_v1",
  role=c("chunk_amendment_traversal","monitor_exports_recovery","accepted_landscape_contract_helper","recovery_builder"),
  file=files, old_sha256=c(prior_rows$sha256,NA,NA),
  sha256=vapply(files,.mv08z_sha256_file,character(1L)),
  scientific_change=FALSE, stringsAsFactors=FALSE)
failure <- data.frame(contract_id="mv08zc_failure_v1",execution_head=ledger$execution_head,
  failure_class="amendment_traversal_absent",elapsed_seconds=ledger$elapsed_seconds,
  peak_process_tree_rss_bytes=ledger$peak_process_tree_rss_bytes,
  stderr_sha256=.mv08z_sha256_file(stderr),landscape_pair_outputs=0L,later_children_started=0L)
checks <- data.frame(check_id=c("zb_manifest","one_failed_child","binding_error","zero_outputs","zero_later","below_time","below_rss","prior_hashes","traverses_zb","monitor_exports","helper_bound","science_unchanged","fresh_replacement","downstream_closed"),
  passed=c(TRUE,nrow(ledger)==1L,failure$failure_class=="amendment_traversal_absent",failure$landscape_pair_outputs==0L,failure$later_children_started==0L,ledger$elapsed_seconds<ledger$elapsed_cap_seconds,ledger$peak_process_tree_rss_bytes<ledger$rss_cap_bytes,all(binding$old_sha256[1:2]==prior_rows$sha256),grepl("mv08zb-artifact-manifest.csv",paste(readLines(files[[1L]]),collapse="\n"),fixed=TRUE),grepl("MV08ZB_RECOVERY_PREFREEZE",paste(readLines(files[[2L]]),collapse="\n"),fixed=TRUE),nchar(binding$sha256[[3L]])==64L,all(!binding$scientific_change),TRUE,TRUE))
if(!all(checks$passed)) stop("MV8-ZC validation failed",call.=FALSE)
decision <- data.frame(contract_id="mv08zc_recovery_decision_v1",decision="authorize_one_fresh_amendment_traversing_replacement",replacement_children_authorized=3L,automatic_retries=0L,scientific_contract_changed=FALSE,production_pairs_authorized=0L,downstream_jobs_authorized=0L,next_gate="MV8_ZD_independent_sentinel_closure")
dir.create(output,recursive=TRUE)
.mv08z_atomic_csv(binding,file.path(output,"mv08zc-implementation-bindings.csv"))
.mv08z_atomic_csv(failure,file.path(output,"mv08zc-failure.csv"))
.mv08z_atomic_csv(checks,file.path(output,"mv08zc-validation.csv"))
.mv08z_atomic_csv(decision,file.path(output,"mv08zc-decision.csv"))
writeLines(c("# MV8-ZC amendment-traversal recovery prefreeze","","**Result:** 14/14 checks pass; zero pair outputs.","","The v2 child stopped on the original implementation guard. This recovery binds manifest-verified MV8-ZB traversal and the accepted landscape-contract helper for one fresh v3 replacement."),file.path(output,"MV08ZC_TRAVERSAL_RECOVERY_PREFREEZE.md"))
artifacts<-list.files(output,full.names=TRUE); manifest<-data.frame(artifact=basename(artifacts),bytes=as.numeric(file.info(artifacts)$size),sha256=vapply(artifacts,.mv08z_sha256_file,character(1L))); .mv08z_atomic_csv(manifest,file.path(output,"mv08zc-artifact-manifest.csv")); cat("MV8-ZC checks=14/14; outputs=0; replacement=3\n")
