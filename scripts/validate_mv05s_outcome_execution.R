#!/usr/bin/env Rscript

options(warn = 2)
for (package in c("aricode", "digest", "mclust")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(package, " is required for independent MV5-S validation.",
         call. = FALSE)
  }
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(paste(
    "usage: validate_mv05s_outcome_execution.R EXTERNAL_METADATA_PATH",
    "PRIVATE_ROOT AUDIT_DIR VALIDATION_OUTPUT"), call. = FALSE)
}
metadata_path <- normalizePath(args[[1L]], mustWork = TRUE)
private_root <- args[[2L]]
audit_dir <- args[[3L]]
validation_output <- args[[4L]]
read_public <- function(name) utils::read.csv(
  file.path(audit_dir, name), stringsAsFactors = FALSE, check.names = FALSE)
file_sha <- function(path) digest::digest(file = path, algo = "sha256",
                                          serialize = FALSE)
object_sha <- function(value) digest::digest(value, algo = "sha256",
                                            serialize = TRUE)
assert_close <- function(observed, expected, label, tolerance = 1e-12) {
  if (length(observed) != length(expected) ||
      any(is.na(observed) != is.na(expected)) ||
      any(abs(observed - expected) > tolerance, na.rm = TRUE)) {
    stop("Independent MV5-S mismatch: ", label, call. = FALSE)
  }
}

queue <- utils::read.csv(
  "docs/audits/mv05r-evaluation-queue-2026-08-10.csv",
  stringsAsFactors = FALSE, check.names = FALSE)
selected <- utils::read.csv(
  "docs/audits/mv05q-selected-training-partitions-2026-08-10.csv.gz",
  stringsAsFactors = FALSE, check.names = FALSE)
heldout <- utils::read.csv(
  "docs/audits/mv05q-heldout-assignments-2026-08-10.csv.gz",
  stringsAsFactors = FALSE, check.names = FALSE)
raw <- utils::read.csv(metadata_path, stringsAsFactors = FALSE, check.names = TRUE)
labels <- data.frame(
  sample_id = trimws(as.character(raw$orig.ident)),
  study = trimws(as.character(raw$SRA)),
  tissue = tolower(trimws(as.character(raw$Tissue.x))),
  approach = trimws(as.character(raw$Approach.x)), stringsAsFactors = FALSE)
labels <- labels[labels$tissue %in% c("bone marrow", "colon", "liver", "pbmc",
                                      "testis"), , drop = FALSE]
if (nrow(queue) != 2400L || nrow(selected) != 126000L ||
    nrow(heldout) != 9000L || nrow(labels) != 90L) {
  stop("Independent MV5-S source dimensions drifted.", call. = FALSE)
}

unit_dir <- file.path(private_root, "units")
results <- vector("list", nrow(queue))
status_rows <- vector("list", nrow(queue))
for (index in seq_len(nrow(queue))) {
  stem <- sub("^mv05r_eval_v1:", "", queue$evaluation_unit_id[[index]])
  artifact <- file.path(unit_dir, paste0(stem, ".rds"))
  status_path <- file.path(unit_dir, paste0(stem, ".status.rds"))
  if (!file.exists(artifact) || !file.exists(status_path)) {
    stop("Independent MV5-S unit pair is missing.", call. = FALSE)
  }
  status <- readRDS(status_path)
  if (!identical(status$evaluation_unit_id, queue$evaluation_unit_id[[index]]) ||
      !identical(status$state, "completed") ||
      !identical(status$artifact_sha256, file_sha(artifact))) {
    stop("Independent MV5-S unit hash/identity failed.", call. = FALSE)
  }
  results[[index]] <- readRDS(artifact)
  status_rows[[index]] <- status
}

is_training <- queue$evaluation_scope == "overlapping_training_partition_alignment"
training_seed <- do.call(rbind, lapply(results[is_training], `[[`, "seed"))
training_summary <- do.call(rbind, lapply(results[is_training], `[[`, "summary"))
public_training_seed <- read_public("mv05s-training-seed-outcomes-2026-08-10.csv")
public_training_summary <- read_public("mv05s-training-unit-summary-2026-08-10.csv")
order_seed <- order(training_seed$evaluation_unit_id, training_seed$seed,
                    method = "radix")
order_public_seed <- order(public_training_seed$evaluation_unit_id,
                           public_training_seed$seed, method = "radix")
if (!identical(training_seed[order_seed, setdiff(names(training_seed), "estimate")],
               public_training_seed[order_public_seed,
                                    setdiff(names(public_training_seed), "estimate")])) {
  stop("Independent MV5-S training public rows drifted.", call. = FALSE)
}
assert_close(training_seed$estimate[order_seed],
             public_training_seed$estimate[order_public_seed],
             "training public estimates")

expected_training <- numeric(nrow(training_seed))
for (index in seq_len(nrow(training_seed))) {
  row <- training_seed[index, , drop = FALSE]
  part <- selected[selected$analysis_group_id == row$analysis_group_id[[1L]] &
                     selected$algorithm_id == row$algorithm_id[[1L]] &
                     as.integer(selected$seed) == as.integer(row$seed[[1L]]), ,
                   drop = FALSE]
  label_index <- match(part$sample_id, labels$sample_id)
  truth <- labels[[row$label_axis[[1L]]]][label_index]
  expected_training[[index]] <- if (row$metric_id[[1L]] ==
                                    "adjusted_rand_index") {
    mclust::adjustedRandIndex(part$cluster, truth)
  } else {
    aricode::NMI(part$cluster, truth, variant = "max")
  }
}
assert_close(training_seed$estimate, expected_training,
             "all 9,000 training metrics")

expected_training_summary <- do.call(rbind, lapply(
  split(training_seed, training_seed$evaluation_unit_id), function(part) {
    values <- part$estimate
    loo <- vapply(seq_along(values), function(i) mean(values[-i]), numeric(1L))
    data.frame(evaluation_unit_id = part$evaluation_unit_id[[1L]],
               seed_mean = mean(values),
               seed_jackknife_se = sqrt((length(values) - 1) / length(values) *
                 sum((loo - mean(loo))^2)), stringsAsFactors = FALSE)
  }))
summary_index <- match(training_summary$evaluation_unit_id,
                       expected_training_summary$evaluation_unit_id)
assert_close(training_summary$seed_mean,
             expected_training_summary$seed_mean[summary_index],
             "training seed means")
assert_close(training_summary$seed_jackknife_se,
             expected_training_summary$seed_jackknife_se[summary_index],
             "training jackknife SE")
if (!setequal(training_summary$evaluation_unit_id,
              public_training_summary$evaluation_unit_id)) {
  stop("Independent MV5-S training summary identity drifted.", call. = FALSE)
}

heldout_results <- results[!is_training]
heldout_private <- do.call(rbind, lapply(heldout_results, `[[`, "private"))
heldout_seed <- do.call(rbind, lapply(heldout_results, `[[`, "seed"))
heldout_summary <- do.call(rbind, lapply(heldout_results, `[[`, "summary"))
for (index in which(!is_training)) {
  queue_row <- queue[index, , drop = FALSE]
  result <- results[[index]]
  for (seed_id in 20260805:20260809) {
    training <- selected[
      selected$analysis_group_id == queue_row$analysis_group_id[[1L]] &
      selected$algorithm_id == queue_row$algorithm_id[[1L]] &
      as.integer(selected$seed) == seed_id, , drop = FALSE]
    query <- heldout[
      heldout$analysis_group_id == queue_row$analysis_group_id[[1L]] &
      heldout$algorithm_id == queue_row$algorithm_id[[1L]] &
      as.integer(heldout$seed) == seed_id, , drop = FALSE]
    axis <- queue_row$label_axis[[1L]]
    train_truth <- labels[[axis]][match(training$sample_id, labels$sample_id)]
    counts <- as.data.frame(table(cluster = training$cluster,
                                  label = train_truth), stringsAsFactors = FALSE)
    counts <- counts[counts$Freq > 0, , drop = FALSE]
    maps <- lapply(split(counts, counts$cluster), function(part) {
      part <- part[order(-part$Freq, part$label, method = "radix"), , drop = FALSE]
      c(cluster = as.character(part$cluster[[1L]]),
        predicted = as.character(part$label[[1L]]),
        count = as.character(part$Freq[[1L]]),
        ties = as.character(sum(part$Freq == max(part$Freq))))
    })
    maps <- as.data.frame(do.call(rbind, maps), stringsAsFactors = FALSE)
    query_map <- match(as.character(query$cluster), maps$cluster)
    expected_prediction <- maps$predicted[query_map]
    expected_truth <- labels[[axis]][match(query$query_sample_id,
                                           labels$sample_id)]
    observed <- result$private[as.integer(result$private$seed) == seed_id, ,
                               drop = FALSE]
    observed <- observed[match(query$query_sample_id, observed$query_sample_id), ,
                               drop = FALSE]
    if (!identical(observed$predicted_label, expected_prediction) ||
        !identical(observed$true_label, expected_truth) ||
        !identical(observed$correct, expected_truth == expected_prediction) ||
        !identical(as.integer(observed$plurality_count),
                   as.integer(maps$count[query_map])) ||
        !identical(as.integer(observed$plurality_tie_size),
                   as.integer(maps$ties[query_map]))) {
      stop("Independent MV5-S held-out prediction reconstruction failed.",
           call. = FALSE)
    }
  }
}

public_heldout_seed <- read_public("mv05s-heldout-seed-outcomes-2026-08-10.csv")
public_heldout_summary <- read_public("mv05s-heldout-unit-summary-2026-08-10.csv")
expected_heldout_seed <- do.call(rbind, lapply(
  split(heldout_private, interaction(heldout_private$evaluation_unit_id,
                                     heldout_private$seed, drop = TRUE)),
  function(part) {
    class_values <- tapply(as.numeric(part$correct), part$true_label, mean)
    data.frame(evaluation_unit_id = part$evaluation_unit_id[[1L]],
               seed = part$seed[[1L]], estimate = mean(class_values),
               stringsAsFactors = FALSE)
  }))
public_key <- paste(public_heldout_seed$evaluation_unit_id,
                    public_heldout_seed$seed)
expected_key <- paste(expected_heldout_seed$evaluation_unit_id,
                      expected_heldout_seed$seed)
seed_index <- match(public_key, expected_key)
if (anyNA(seed_index)) stop("Independent MV5-S held-out seed keys drifted.",
                            call. = FALSE)
assert_close(public_heldout_seed$estimate,
             expected_heldout_seed$estimate[seed_index],
             "held-out seed balanced accuracy")
if (!setequal(heldout_seed$evaluation_unit_id,
              public_heldout_seed$evaluation_unit_id) ||
    !setequal(heldout_summary$evaluation_unit_id,
              public_heldout_summary$evaluation_unit_id)) {
  stop("Independent MV5-S held-out public identity drifted.", call. = FALSE)
}

private_sample <- readRDS(file.path(private_root,
                                    "heldout-sample-summary-private.rds"))
sample_columns <- c("representation", "distance_id", "algorithm_id",
                    "algorithm_role", "endpoint_id", "label_axis",
                    "query_sample_id", "heldout_study", "true_label")
expected_sample <- do.call(rbind, lapply(
  split(heldout_private, interaction(heldout_private[sample_columns],
                                     drop = TRUE, lex.order = TRUE)),
  function(part) data.frame(
    representation = part$representation[[1L]],
    distance_id = part$distance_id[[1L]], algorithm_id = part$algorithm_id[[1L]],
    algorithm_role = part$algorithm_role[[1L]], endpoint_id = part$endpoint_id[[1L]],
    label_axis = part$label_axis[[1L]], sample_id = part$query_sample_id[[1L]],
    study = part$heldout_study[[1L]], true_label = part$true_label[[1L]],
    seed_mean_correct = mean(part$correct), completed_seeds = nrow(part),
    stringsAsFactors = FALSE)))
sample_order <- order(private_sample$endpoint_id, private_sample$representation,
                      private_sample$distance_id, private_sample$algorithm_id,
                      private_sample$sample_id, method = "radix")
expected_order <- order(expected_sample$endpoint_id, expected_sample$representation,
                        expected_sample$distance_id, expected_sample$algorithm_id,
                        expected_sample$sample_id, method = "radix")
if (!identical(private_sample[sample_order, setdiff(names(private_sample),
                                                    "seed_mean_correct")],
               expected_sample[expected_order, setdiff(names(expected_sample),
                                                       "seed_mean_correct")])) {
  stop("Independent MV5-S private sample summary drifted.", call. = FALSE)
}
assert_close(private_sample$seed_mean_correct[sample_order],
             expected_sample$seed_mean_correct[expected_order],
             "held-out sample seed means")

independent_counts <- function(axis) {
  studies <- sort(unique(labels$study), method = "radix")
  counts <- matrix(0L, 2000L, length(studies), dimnames = list(NULL, studies))
  set.seed(20260810, kind = "Mersenne-Twister", normal.kind = "Inversion",
           sample.kind = "Rejection")
  if (axis == "tissue") {
    strata <- split(unique(labels[c("study", "tissue")])$study,
                    unique(labels[c("study", "tissue")])$tissue)
    for (i in seq_len(2000L)) {
      sampled <- unlist(lapply(strata, function(ids) sample(ids, length(ids),
                                                replace = TRUE)),
                        use.names = FALSE)
      counts[i, ] <- tabulate(match(sampled, studies), nbins = length(studies))
    }
  } else {
    for (i in seq_len(2000L)) {
      sampled <- sample(studies, length(studies), replace = TRUE)
      counts[i, ] <- tabulate(match(sampled, studies), nbins = length(studies))
    }
  }
  counts
}
tissue_counts <- independent_counts("tissue")
approach_counts <- independent_counts("approach")
macro_public <- read_public("mv05s-heldout-macro-outcomes-2026-08-10.csv")
bootstrap_public <- read_public("mv05s-bootstrap-audit-2026-08-10.csv")
context_columns <- c("representation", "distance_id", "algorithm_id",
                     "algorithm_role", "endpoint_id", "label_axis")
sample_contexts <- split(private_sample,
  interaction(private_sample[context_columns], drop = TRUE, lex.order = TRUE))
for (part in sample_contexts) {
  axis <- part$label_axis[[1L]]
  counts <- if (axis == "tissue") tissue_counts else approach_counts
  classes <- sort(unique(part$true_label), method = "radix")
  point <- mean(vapply(classes, function(class_id) {
    mean(part$seed_mean_correct[part$true_label == class_id])
  }, numeric(1L)))
  estimates <- vapply(seq_len(nrow(counts)), function(i) {
    weights <- counts[i, match(part$study, colnames(counts))]
    values <- vapply(classes, function(class_id) {
      keep <- part$true_label == class_id
      if (sum(weights[keep]) == 0) return(NA_real_)
      weighted.mean(part$seed_mean_correct[keep], weights[keep])
    }, numeric(1L))
    if (anyNA(values)) NA_real_ else mean(values)
  }, numeric(1L))
  keep <- macro_public$representation == part$representation[[1L]] &
    macro_public$distance_id == part$distance_id[[1L]] &
    macro_public$algorithm_id == part$algorithm_id[[1L]] &
    macro_public$endpoint_id == part$endpoint_id[[1L]]
  observed <- macro_public[keep, , drop = FALSE]
  estimable <- estimates[is.finite(estimates)]
  interval <- stats::quantile(estimable, c(0.025, 0.975), names = FALSE,
                              type = 7)
  if (nrow(observed) != 1L) stop("Independent MV5-S macro key drifted.",
                                 call. = FALSE)
  assert_close(observed$estimate, point, "held-out macro point")
  assert_close(c(observed$ci_lower, observed$ci_upper), interval,
               "held-out macro interval")
  if (observed$estimable_bootstrap_replicates != length(estimable)) {
    stop("Independent MV5-S estimable replicate count drifted.",
         call. = FALSE)
  }
  audit <- bootstrap_public[
    bootstrap_public$representation == part$representation[[1L]] &
    bootstrap_public$distance_id == part$distance_id[[1L]] &
    bootstrap_public$algorithm_id == part$algorithm_id[[1L]] &
    bootstrap_public$endpoint_id == part$endpoint_id[[1L]], , drop = FALSE]
  if (nrow(audit) != 1L ||
      !identical(audit$block_count_matrix_sha256[[1L]], object_sha(counts)) ||
      !identical(audit$replicate_estimates_sha256[[1L]], object_sha(estimates))) {
    stop("Independent MV5-S bootstrap hash reconstruction failed.",
         call. = FALSE)
  }
}

public_files <- list.files(audit_dir, pattern = "^mv05s-.*[.]csv$",
                           full.names = TRUE)
prohibited <- c("true_label", "predicted_label")
label_values <- tolower(unique(c(labels$tissue, labels$approach)))
for (path in public_files) {
  value <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (any(names(value) %in% prohibited) ||
      any(vapply(value, function(column) {
        any(tolower(as.character(column)) %in% label_values)
      }, logical(1L)))) {
    stop("Independent MV5-S public safety failed.", call. = FALSE)
  }
}

validation <- data.frame(
  contract_id = "mv05s_independent_validation_v1",
  category = c(
    "unit_identity_hashes", "all_training_ari_nmi",
    "training_seed_aggregation", "training_public_roundtrip",
    "heldout_training_only_plurality", "heldout_seed_balanced_accuracy",
    "heldout_sample_seed_aggregation", "tissue_block_bootstrap",
    "approach_mixed_study_bootstrap", "macro_points_intervals",
    "bootstrap_hashes", "public_label_safety"),
  checks = c(2400L, nrow(training_seed), nrow(training_summary),
             nrow(public_training_seed), nrow(heldout_private),
             nrow(public_heldout_seed), nrow(private_sample), 2000L, 2000L,
             nrow(macro_public), nrow(bootstrap_public), length(public_files)),
  passed = TRUE, stringsAsFactors = FALSE)
temporary <- paste0(validation_output, ".tmp-", Sys.getpid())
utils::write.csv(validation, temporary, row.names = FALSE, na = "")
if (file.exists(validation_output)) unlink(validation_output)
if (!file.rename(temporary, validation_output)) {
  stop("Independent MV5-S validation output rename failed.", call. = FALSE)
}
message("Independent MV5-S validation passed: categories=", nrow(validation),
        " training_metrics=", nrow(training_seed),
        " heldout_predictions=", nrow(heldout_private),
        " macro_contexts=", nrow(macro_public))
