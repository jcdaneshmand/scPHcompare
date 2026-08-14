#!/usr/bin/env Rscript

options(warn = 2)
args <- getOption("mv06f.args", commandArgs(trailingOnly = TRUE))
if (length(args) != 10L) {
  stop("usage: validate_mv06f_stage2_rebind.R QUEUE OLD_CONTRACT ",
       "NEW_CONTRACT OLD_GROUP NEW_GROUP NEW_REPEAT_GROUP NEW_R_ORACLES ",
       "NEW_PERSIM_ORACLES NEW_RESUME OUTPUT_DIR", call. = FALSE)
}
source("R/mv06f_production.R")
read_csv <- function(index) utils::read.csv(
  args[[index]], stringsAsFactors = FALSE, check.names = FALSE
)
queue <- read_csv(1L); old_contract <- read_csv(2L); new_contract <- read_csv(3L)
stage <- queue[queue$stage == "stage_1_maximum", , drop = FALSE]
if (nrow(stage) != 1L || nrow(old_contract) != 1L ||
    nrow(new_contract) != 1L) {
  stop("MV6-F rebind validation inputs are stale.", call. = FALSE)
}
mv06f_validate_group_directory_v1(
  args[[4L]], stage, old_contract$queue_root_sha256,
  old_contract$implementation_root_sha256, old_contract$rust_library_sha256
)
for (path in args[5:6]) mv06f_validate_group_directory_v1(
  path, stage, new_contract$queue_root_sha256,
  new_contract$implementation_root_sha256, new_contract$rust_library_sha256
)

old_manifest <- utils::read.csv(file.path(args[[4L]], "diagram-manifest.csv"),
                                stringsAsFactors = FALSE, check.names = FALSE)
new_manifest <- utils::read.csv(file.path(args[[5L]], "diagram-manifest.csv"),
                                stringsAsFactors = FALSE, check.names = FALSE)
manifest_fields <- c(
  "sample_id", "role", "view_id", "diagram_sha256", "finite_h0_intervals",
  "finite_h1_intervals", "essential_h0_intervals", "outcome_label_state",
  "biological_outcomes_computed"
)
manifest_equal <- identical(old_manifest[manifest_fields],
                            new_manifest[manifest_fields])
old_distances <- utils::read.csv(file.path(args[[4L]], "distances.csv"),
                                 stringsAsFactors = FALSE, check.names = FALSE)
new_distances <- utils::read.csv(file.path(args[[5L]], "distances.csv"),
                                 stringsAsFactors = FALSE, check.names = FALSE)
distance_equal <- identical(old_distances, new_distances)
old_records <- readRDS(file.path(args[[4L]], "diagrams.rds"))
new_records <- readRDS(file.path(args[[5L]], "diagrams.rds"))
diagram_equal <- identical(names(old_records), names(new_records)) &&
  all(vapply(names(old_records), function(name) {
    identical(old_records[[name]]$topology_result$diagram,
              new_records[[name]]$topology_result$diagram)
  }, logical(1L)))

scientific_files <- c("diagrams.rds", "diagram-manifest.csv", "distances.csv")
repeat_equal <- all(vapply(scientific_files, function(name) {
  .mv06f_sha256(file.path(args[[5L]], name)) ==
    .mv06f_sha256(file.path(args[[6L]], name))
}, logical(1L)))
r_rows <- read_csv(7L); persim_rows <- read_csv(8L); resume_rows <- read_csv(9L)
r_pass <- nrow(r_rows) == 12L && all(as.logical(r_rows$passed))
persim_pass <- nrow(persim_rows) == 12L && all(as.logical(persim_rows$passed))
resume_pass <- nrow(resume_rows) == 5L && all(as.logical(resume_rows$passed))
checks <- data.frame(
  contract_id = "mv06f_stage2_rebind_validation_v1",
  category = c(
    "old_and_new_group_contracts", "manifest_scientific_equivalence",
    "diagram_numeric_equivalence", "distance_exact_equivalence",
    "new_root_clean_repeat", "new_root_r_oracles",
    "new_root_persim_oracles", "new_root_immutable_resume"
  ),
  passed = c(TRUE, manifest_equal, diagram_equal, distance_equal,
             repeat_equal, r_pass, persim_pass, resume_pass),
  detail = c(
    "old and remediated stage-one directories validate under their roots",
    "sample/view roles, diagram hashes, and interval counts are identical",
    "all 180 H0/H1 topology diagram matrices are byte-value identical",
    "all 6,500 exact all-active-level component rows are identical",
    "three remediated scientific artifacts repeat byte-identically",
    "12 balanced cell/gene H0/H1 R oracles pass",
    "12 balanced cell/gene H0/H1 grouped-Persim oracles pass",
    "five remediated stage-one files survive zero-rebuild resume"
  ), outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
if (any(!checks$passed)) {
  print(checks)
  stop("MV6-F remediated stage-one rebind validation failed.", call. = FALSE)
}
admission <- data.frame(
  contract_id = "mv06f_stage2_admission_v1",
  queue_root_sha256 = new_contract$queue_root_sha256,
  implementation_root_sha256 = new_contract$implementation_root_sha256,
  rust_library_sha256 = new_contract$rust_library_sha256,
  stage1_scientific_equivalence_passed = TRUE,
  stage1_clean_repeat_passed = TRUE, stage1_r_oracles_passed = TRUE,
  stage1_persim_oracles_passed = TRUE, stage1_resume_passed = TRUE,
  stage2_authorized = TRUE, outcome_label_state = "closed",
  biological_outcomes_computed = FALSE, fusion_jobs = 0L,
  clustering_jobs = 0L, outcome_jobs = 0L, stringsAsFactors = FALSE
)
dir.create(args[[10L]], recursive = TRUE, showWarnings = FALSE)
utils::write.csv(checks, file.path(
  args[[10L]], "mv06f-stage2-rebind-validation.csv"
), row.names = FALSE, na = "")
utils::write.csv(admission, file.path(
  args[[10L]], "mv06f-stage2-admission.csv"
), row.names = FALSE, na = "")
message("Validated MV6-F remediated stage one: 8/8 categories pass.")
