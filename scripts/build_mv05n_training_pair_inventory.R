#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 9L) {
  stop(paste(
    "usage: build_mv05n_training_pair_inventory.R",
    "SCT_PH_MANIFEST SCT_VIEW_METRICS INTEGRATED_PH_MANIFEST",
    "INTEGRATED_VIEW_METRICS IDENTITY_SUMMARY_OUTPUT GROUP_INVENTORY_OUTPUT",
    "CHUNK_INVENTORY_OUTPUT ADMISSION_REQUEST_OUTPUT ADMISSION_PAIRS_PER_DIM"
  ), call. = FALSE)
}

implementation_path <- Sys.getenv(
  "SCPH_MV05N_SOURCE", unset = "R/mv05n_clustering_gate.R"
)
source(implementation_path)
source("R/provenance_utils.R")

read_public <- function(path) {
  utils::read.csv(normalizePath(path, mustWork = TRUE), stringsAsFactors = FALSE,
                  check.names = FALSE)
}
file_sha <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}
write_immutable <- function(value, path) {
  if (file.exists(path)) stop("Refusing to overwrite ", path, call. = FALSE)
  write_provenance_csv(value, path)
}

sct_ph <- read_public(args[[1L]])
sct_metrics <- read_public(args[[2L]])
integrated_ph <- read_public(args[[3L]])
integrated_metrics <- read_public(args[[4L]])
identity_output <- args[[5L]]
group_output <- args[[6L]]
chunk_output <- args[[7L]]
admission_output <- args[[8L]]
admission_n <- as.integer(args[[9L]])
if (is.na(admission_n) || admission_n < 1L || admission_n > 64L) {
  stop("ADMISSION_PAIRS_PER_DIM must be in 1:64.", call. = FALSE)
}

inputs <- list(
  sct_whole = list(ph = sct_ph, metrics = sct_metrics),
  inductive_integrated = list(ph = integrated_ph, metrics = integrated_metrics)
)
group_rows <- list()
chunk_rows <- list()
admission_candidates <- list()
group_cursor <- chunk_cursor <- candidate_cursor <- 0L

for (representation in names(inputs)) {
  current <- inputs[[representation]]
  if (!setequal(current$ph$job_id, current$metrics$job_id) ||
      any(current$ph$representation != representation)) {
    stop("MV5-N source identities do not align for ", representation, ".",
         call. = FALSE)
  }
  source_groups <- unique(current$ph$group_id[order(current$ph$group_order)])
  for (source_group in source_groups) {
    ph_group <- current$ph[current$ph$group_id == source_group, , drop = FALSE]
    metrics <- current$metrics[current$metrics$job_id %in% ph_group$job_id,
                               , drop = FALSE]
    manifest <- mv05n_build_group_pair_manifest_v1(ph_group, metrics)
    group_cursor <- group_cursor + 1L
    group_rows[[group_cursor]] <- mv05n_group_inventory_v1(manifest)
    chunks <- split(manifest, manifest$chunk_id)
    for (chunk in chunks) {
      chunk_cursor <- chunk_cursor + 1L
      chunk_rows[[chunk_cursor]] <- data.frame(
        contract_id = "mv05n_training_pair_chunk_inventory_v1",
        chunk_id = chunk$chunk_id[[1L]], source_group_id = source_group,
        group_order = chunk$group_order[[1L]], fold_id = chunk$fold_id[[1L]],
        seed = chunk$seed[[1L]], representation = representation,
        homology_dimension = chunk$homology_dimension[[1L]],
        request_rows = nrow(chunk),
        request_identity_set_sha256 = .mv05n_digest(
          sort(chunk$pair_request_id, method = "radix")
        ),
        first_pair_request_id = min(chunk$pair_request_id),
        last_pair_request_id = max(chunk$pair_request_id),
        execution_disposition = "not_authorized_full_production",
        outcome_label_state = "closed", biological_outcomes_computed = FALSE,
        stringsAsFactors = FALSE
      )
    }
    ordered <- manifest[order(manifest$homology_dimension,
                              manifest$pair_request_id, method = "radix"), ]
    for (dimension in c("H0", "H1")) {
      eligible <- ordered[ordered$homology_dimension == dimension, , drop = FALSE]
      positions <- unique(as.integer(round(seq(1L, nrow(eligible),
                                               length.out = admission_n))))
      candidate_cursor <- candidate_cursor + 1L
      admission_candidates[[candidate_cursor]] <- eligible[positions, , drop = FALSE]
    }
  }
}

groups <- do.call(rbind, group_rows)
chunks <- do.call(rbind, chunk_rows)
if (nrow(groups) != 150L ||
    any(table(groups$representation) != 75L) ||
    any(tapply(groups$unordered_training_pairs, groups$representation, sum) !=
        262675L) ||
    any(tapply(groups$h0_h1_request_rows, groups$representation, sum) !=
        525350L)) {
  stop("MV5-N full request inventory does not reproduce the exact pair scope.",
       call. = FALSE)
}

profiles <- mv05n_select_admission_profiles_v1(groups)
admission <- do.call(rbind, admission_candidates)
admission <- merge(admission, profiles[c("profile", "fold_id")], by = "fold_id")
admission <- admission[admission$seed == min(.mv05n_required_seeds), , drop = FALSE]
admission <- admission[order(admission$profile, admission$representation,
                             admission$homology_dimension,
                             admission$pair_request_id, method = "radix"), ]
admission$admission_contract_id <- "mv05n_label_closed_admission_request_v1"
admission$admission_group_id <- paste(
  "mv05n_admission_v1", admission$profile, admission$representation,
  admission$homology_dimension, sep = ":"
)
admission$admission_seed_rule <- "canonical_first_frozen_seed_20260805"
admission$full_production_authorized <- FALSE
admission$clustering_jobs_executed <- 0L
rownames(admission) <- NULL

groups <- groups[order(groups$representation, groups$group_order,
                       method = "radix"), ]
chunks <- chunks[order(chunks$representation, chunks$group_order,
                       chunks$homology_dimension, chunks$chunk_id,
                       method = "radix"), ]
identity_summary <- do.call(rbind, lapply(names(inputs), function(representation) {
  selected_groups <- groups[groups$representation == representation, , drop = FALSE]
  selected_chunks <- chunks[chunks$representation == representation, , drop = FALSE]
  data.frame(
    contract_id = "mv05n_training_pair_identity_summary_v1",
    representation = representation, folds = 15L, seeds = 5L,
    groups = nrow(selected_groups), training_pairs = 262675L,
    h0_request_rows = 262675L, h1_request_rows = 262675L,
    h0_h1_request_rows = 525350L, chunks = nrow(selected_chunks),
    group_identity_set_sha256 = .mv05n_digest(
      sort(selected_groups$request_identity_set_sha256, method = "radix")
    ),
    chunk_identity_set_sha256 = .mv05n_digest(
      sort(selected_chunks$request_identity_set_sha256, method = "radix")
    ),
    source_ph_manifest_sha256 = if (representation == "sct_whole")
      file_sha(args[[1L]]) else file_sha(args[[3L]]),
    source_view_metrics_sha256 = if (representation == "sct_whole")
      file_sha(args[[2L]]) else file_sha(args[[4L]]),
    persistence_policy = "instantiate_validate_digest_discard_rows_until_authorized_v1",
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
}))
write_immutable(identity_summary, identity_output)
write_immutable(groups, group_output)
write_immutable(chunks, chunk_output)
write_immutable(admission, admission_output)
message("Built MV5-N identities: 525350 H0/H1 rows per representation; ",
        nrow(admission), " bounded admission rows.")
