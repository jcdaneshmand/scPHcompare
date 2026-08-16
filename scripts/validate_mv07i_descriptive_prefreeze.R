#!/usr/bin/env Rscript

options(warn = 2)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required.")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("usage: validate_mv07i_descriptive_prefreeze.R PREFREEZE REPEAT OUTPUT",
       call. = FALSE)
}
source("R/mv07h_full_topology.R")
read_csv <- function(path) utils::read.csv(
  path, stringsAsFactors = FALSE, check.names = FALSE)
prefreeze <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
repeat_root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
output <- args[[3L]]
if (dir.exists(output) && length(list.files(output, all.files = TRUE,
                                             no.. = TRUE))) {
  stop("MV7-I prefreeze validation output must be empty.", call. = FALSE)
}
files <- sort(list.files(prefreeze, pattern = "\\.csv$"), method = "radix")
repeat_files <- sort(list.files(repeat_root, pattern = "\\.csv$"),
                     method = "radix")
repeat_validation <- data.frame(
  contract_id = "mv07i_prefreeze_repeat_validation_v1", file = files,
  public_sha256 = vapply(file.path(prefreeze, files), .mv07h_sha256,
                         character(1L)),
  repeat_sha256 = vapply(file.path(repeat_root, files), .mv07h_sha256,
                         character(1L)),
  byte_identical = FALSE, stringsAsFactors = FALSE
)
repeat_validation$byte_identical <-
  repeat_validation$public_sha256 == repeat_validation$repeat_sha256
manifest <- read_csv(file.path(prefreeze, "mv07i-input-manifest.csv"))
representation <- read_csv(file.path(
  prefreeze, "mv07i-representation-registry.csv"))
aggregation <- read_csv(file.path(
  prefreeze, "mv07i-seed-aggregation-contract.csv"))
clustering <- read_csv(file.path(prefreeze, "mv07i-clustering-contract.csv"))
populations <- read_csv(file.path(prefreeze, "mv07i-populations.csv"))
metadata <- read_csv(file.path(prefreeze, "mv07i-metadata-registry.csv"))
boundaries <- read_csv(file.path(
  prefreeze, "mv07i-interpretation-boundaries.csv"))
stages <- read_csv(file.path(prefreeze, "mv07i-stage-authorization.csv"))
resources <- read_csv(file.path(
  prefreeze, "mv07i-resource-resume-contract.csv"))
builder_checks <- read_csv(file.path(prefreeze, "mv07i-prefreeze-checks.csv"))
decision <- read_csv(file.path(prefreeze, "mv07i-decision.csv"))
checks <- data.frame(
  contract_id = "mv07i_descriptive_prefreeze_independent_validation_v1",
  category = c("complete_output", "byte_repeat", "input_hashes",
               "representation_registry", "aggregation_policy",
               "label_free_clustering", "algorithm_boundary",
               "population_boundary", "metadata_lineage",
               "interpretation_boundary", "stage_firewall",
               "resource_resume", "builder_and_decision"),
  passed = c(
    identical(files, repeat_files) && length(files) == 11L,
    nrow(repeat_validation) == 11L &&
      all(repeat_validation$byte_identical),
    nrow(manifest) == 9L && all(file.exists(manifest$path)) &&
      identical(tolower(unname(vapply(
        manifest$path, .mv07h_sha256, character(1L)))),
        tolower(unname(manifest$sha256))),
    nrow(representation) == 6L &&
      sum(representation$scientific_role == "primary_component") == 4L &&
      sum(representation$scientific_role ==
            "secondary_descriptive_composite") == 2L,
    nrow(aggregation) == 5L &&
      all(aggregation$select_favorable_seed == FALSE) &&
      aggregation$uncertainty[aggregation$quantity == "distance"] ==
        "min;max;IQR;raw_MAD_constant_1",
    nrow(clustering) == 6L && all(clustering$candidate_k == "2:10") &&
      all(clustering$k_selection ==
            "smallest_k_within_one_SE_of_maximum_mean_stability") &&
      !any(as.logical(clustering$label_values_used)),
    all(clustering$primary_algorithm == "PAM_dissimilarity_v1") &&
      all(clustering$sensitivity_algorithm ==
            "average_linkage_at_PAM_selected_k") &&
      all(grepl("deferred", clustering$spectral_status, fixed = TRUE)),
    nrow(populations) == 3L &&
      populations$samples[populations$population_id ==
        "corrected_full_corpus_descriptive"] == 124L &&
      !populations$may_support_new_cross_study_tissue_claim[
        populations$population_id == "corrected_full_corpus_descriptive"],
    nrow(metadata) == 2L && all(metadata$expected_rows == 124L) &&
      !any(as.logical(metadata$may_select_method_k_view_dimension_seed)),
    nrow(boundaries) == 6L &&
      all(nzchar(boundaries$prohibited)),
    all(stages$authorized_now[1:2]) && !any(stages$authorized_now[3:4]) &&
      all(stages$label_access[1:2] == "FALSE"),
    resources$maximum_workers == 1L &&
      resources$candidate_pam_fits == 270L &&
      resources$selected_average_linkage_fits == 30L &&
      resources$deterministic_repeat_required &&
      resources$immutable_resume_required,
    nrow(builder_checks) == 13L && all(builder_checks$passed) &&
      decision$decision ==
        "authorize_label_closed_matrix_and_clustering_production_only" &&
      !as.logical(decision$labels_authorized) &&
      !as.logical(decision$outcomes_authorized)
  ),
  detail = c("11 expected public CSV artifacts", "11/11 byte-identical rebuild",
             "nine accepted source hashes", "four primary plus two secondary",
             "within-seed derivation and five-seed uncertainty",
             "PAM k=2:10 five-seed one-SE without labels",
             "average linkage only; spectral deferred",
             "90 cross-study context and 124 descriptive boundary",
             "canonical 124-sample tissue/study/approach sources",
             "six explicit claim prohibitions", "only label-closed stages open",
             "one worker, 270 PAM, 30 average, repeat and resume",
             "13/13 builder and label-closed production decision"),
  stringsAsFactors = FALSE
)
if (any(!checks$passed)) {
  stop("MV7-I independent prefreeze validation failed: ",
       paste(checks$category[!checks$passed], collapse = ", "),
       call. = FALSE)
}
validation_decision <- data.frame(
  contract_id = "mv07i_descriptive_prefreeze_validation_decision_v1",
  decision = "authorize_label_closed_matrix_and_clustering_production_only",
  validation_categories = 13L, byte_repeat_artifacts = 11L,
  labels_authorized = FALSE, outcomes_authorized = FALSE,
  claims_authorized = FALSE, external_data_authorized = FALSE,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
dir.create(output, recursive = TRUE, showWarnings = FALSE)
write.csv(repeat_validation, file.path(
  output, "mv07i-prefreeze-repeat-validation.csv"),
          row.names = FALSE, na = "")
write.csv(checks, file.path(output,
  "mv07i-prefreeze-independent-validation.csv"), row.names = FALSE, na = "")
write.csv(validation_decision, file.path(output,
  "mv07i-prefreeze-validation-decision.csv"), row.names = FALSE, na = "")
message("MV7-I descriptive prefreeze independent validation: 13/13 pass")
