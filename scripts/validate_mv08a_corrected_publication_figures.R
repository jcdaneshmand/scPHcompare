#!/usr/bin/env Rscript

options(warn = 2, digits = 17, scipen = 999)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required.")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop(paste(
    "usage: validate_mv08a_corrected_publication_figures.R",
    "PRODUCTION_DIR REPEAT_DIR PH_RDS H1_PRIVATE_CSV VALIDATION_DIR EXPECTED_HEAD"),
    call. = FALSE)
}

production_dir <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
repeat_dir <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
ph_path <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
h1_path <- normalizePath(args[[4L]], winslash = "/", mustWork = TRUE)
validation_dir <- args[[5L]]
expected_head <- tolower(trimws(args[[6L]]))
if (dir.exists(validation_dir) && length(list.files(validation_dir,
    all.files = TRUE, no.. = TRUE))) stop("MV8-A validation directory must be empty.")
dir.create(validation_dir, recursive = TRUE, showWarnings = FALSE)

readc <- function(path) utils::read.csv(path, stringsAsFactors = FALSE,
                                         check.names = FALSE)
sha <- function(path) digest::digest(file = path, algo = "sha256",
                                      serialize = FALSE)
truth <- function(x) if (is.logical(x)) !is.na(x) & x else
  tolower(trimws(x)) == "true"
write_csv <- function(x, path) utils::write.table(x, path, sep = ",",
  row.names = FALSE, col.names = TRUE, quote = TRUE, na = "",
  qmethod = "double")
png_dimensions <- function(path) {
  con <- file(path, open = "rb")
  on.exit(close(con))
  header <- readBin(con, what = "raw", n = 16L)
  if (length(header) != 16L ||
      !identical(rawToChar(header[13:16]), "IHDR")) return(c(NA, NA))
  readBin(con, what = integer(), n = 2L, size = 4L, endian = "big",
          signed = TRUE)
}

checks <- list()
add_check <- function(id, pass, detail) {
  checks[[length(checks) + 1L]] <<- data.frame(
    contract_id = "mv08a_validation_check_v1", check_id = id,
    pass = isTRUE(pass), detail = as.character(detail),
    stringsAsFactors = FALSE)
}

figure_ids <- paste0("figure-", 1:7, c("-corrected-architecture",
  "-cohort-confounding", "-diagram-to-landscape",
  "-h1-contribution-concordance", "-label-free-stability",
  "-descriptive-alignment", "-algorithm-sensitivity"))
root_files <- c(paste0(figure_ids, ".svg"), paste0(figure_ids, ".png"),
  "mv08a-source-manifest.csv", "mv08a-figure-manifest.csv",
  "mv08a-data-manifest.csv", "mv08a-renderer-provenance.csv")
data_rows <- c(
  "f1-nodes.csv" = 12L, "f1-edges.csv" = 12L,
  "f2-tissue-counts.csv" = 16L, "f2-cohort-flow.csv" = 5L,
  "f3-finite-diagram.csv" = 734L,
  "f3-display-landscape-levels.csv" = 7200L,
  "f4-h1-ecdf.csv" = 600L,
  "f4-h0-composite-concordance.csv" = 20L,
  "f5-complete-stability.csv" = 54L,
  "f6-complete-outcomes.csv" = 120L,
  "f7-algorithm-sensitivity.csv" = 30L)
complete <- all(file.exists(file.path(production_dir, root_files))) &&
  all(file.exists(file.path(repeat_dir, root_files))) &&
  all(file.exists(file.path(production_dir, "data", names(data_rows)))) &&
  all(file.exists(file.path(repeat_dir, "data", names(data_rows))))
add_check("complete_output_family", complete,
  "seven SVG/PNG pairs, four root manifests, and eleven data tables")
if (!complete) {
  report <- do.call(rbind, checks)
  write_csv(report, file.path(validation_dir, "mv08a-validation-checks.csv"))
  stop("MV8-A output family is incomplete.")
}

source <- readc(file.path(production_dir, "mv08a-source-manifest.csv"))
figures <- readc(file.path(production_dir, "mv08a-figure-manifest.csv"))
data_manifest <- readc(file.path(production_dir, "mv08a-data-manifest.csv"))
provenance <- readc(file.path(production_dir, "mv08a-renderer-provenance.csv"))

head_now <- tolower(trimws(system2("git", c("rev-parse", "HEAD"),
                                   stdout = TRUE)))
head_ok <- identical(head_now, expected_head) && nrow(source) == 16L &&
  all(tolower(source$expected_head) == expected_head)
add_check("prospective_head", head_ok,
  paste("expected", expected_head, "observed", head_now))

source_ok <- !anyDuplicated(source$source_id) &&
  all(!truth(source$new_ph)) && all(!truth(source$new_data))
for (i in seq_len(nrow(source))) {
  if (source$access_class[[i]] == "public_repository") {
    source_ok <- source_ok && file.exists(source$locator[[i]]) &&
      identical(tolower(sha(source$locator[[i]])),
                tolower(source$sha256[[i]]))
  }
}
private_paths <- c(fixed_ph_artifact = ph_path,
                   private_h1_pair_table = h1_path)
for (id in names(private_paths)) {
  row <- source[source$source_id == id, , drop = FALSE]
  source_ok <- source_ok && nrow(row) == 1L &&
    row$locator[[1L]] == "private_validated_artifact_not_published" &&
    identical(tolower(sha(private_paths[[id]])), tolower(row$sha256[[1L]]))
}
add_check("source_hashes_and_private_redaction", source_ok,
  "all public sources rehashed; two private sources supplied and redacted")

figure_ok <- nrow(figures) == 7L &&
  identical(figures$figure_id, figure_ids) &&
  all(figures$dpi == 300L) &&
  all(figures$png_width_pixels >= 2400L |
      figures$png_height_pixels >= 2400L) &&
  !any(truth(figures$manuscript_claim_authorized))
for (i in seq_len(nrow(figures))) {
  svg <- file.path(production_dir, figures$svg_filename[[i]])
  png <- file.path(production_dir, figures$png_filename[[i]])
  observed <- png_dimensions(png)
  figure_ok <- figure_ok && identical(tolower(sha(svg)),
    tolower(figures$svg_sha256[[i]])) && identical(tolower(sha(png)),
    tolower(figures$png_sha256[[i]])) &&
    identical(as.integer(observed), as.integer(c(
      figures$png_width_pixels[[i]], figures$png_height_pixels[[i]])))
}
add_check("figure_hashes_dimensions_and_dpi", figure_ok,
  "seven figure hashes and PNG IHDR dimensions independently checked")

svg_text <- vapply(file.path(production_dir, paste0(figure_ids, ".svg")),
  function(path) paste(readLines(path, warn = FALSE, encoding = "UTF-8"),
                       collapse = "\n"), character(1L))
svg_ok <- all(grepl("<svg", svg_text, fixed = TRUE)) &&
  all(grepl("<text", svg_text, fixed = TRUE)) &&
  !any(grepl("<image", svg_text, fixed = TRUE)) &&
  !any(grepl("100-point|100 point|first-level|first level", svg_text,
             ignore.case = TRUE))
add_check("editable_svg_and_no_legacy_claim", svg_ok,
  "live SVG text; no raster image tag or legacy first-level/100-grid language")

architecture <- svg_text[[1L]]
landscape <- svg_text[[3L]]
landscape_ok <- grepl("samples remain the independent comparison units",
  architecture, fixed = TRUE) && grepl("No fixed scientific grid",
  architecture, fixed = TRUE) && grepl("no level cap", architecture,
  fixed = TRUE) && grepl("every active level", landscape, fixed = TRUE) &&
  grepl("H0 and H1 remain separate", landscape, fixed = TRUE) &&
  grepl("Essential H0 class excluded", landscape, fixed = TRUE) &&
  grepl("visualization-only", landscape, fixed = TRUE)
add_check("corrected_landscape_contract_visible", landscape_ok,
  "all-active-level, separate H0/H1, essential-class and display-grid rules visible")

row_ok <- TRUE
for (name in names(data_rows)) {
  x <- readc(file.path(production_dir, "data", name))
  row_ok <- row_ok && nrow(x) == data_rows[[name]]
}
add_check("data_cardinality", row_ok,
  "12/12, 16/5, 734/7200, 600/20, 54, 120, and 30 rows")

stab <- readc(file.path(production_dir, "data", "f5-complete-stability.csv"))
outcome <- readc(file.path(production_dir, "data", "f6-complete-outcomes.csv"))
algorithm <- readc(file.path(production_dir, "data", "f7-algorithm-sensitivity.csv"))
families_ok <- length(unique(stab$representation_id)) == 6L &&
  all(table(stab$representation_id) == 9L) &&
  identical(sort(unique(stab$k)), 2:10) && all(stab$selected_k == 2L) &&
  length(unique(outcome$evaluation_unit_id)) == 120L &&
  all(outcome$status == "completed") && !any(truth(outcome$p_value_computed)) &&
  nrow(algorithm) == 30L &&
  !any(truth(algorithm$favorable_algorithm_selected))
add_check("complete_unselected_families", families_ok,
  "all stability, outcome, and algorithm rows included without winner selection")

data_files <- file.path(production_dir, "data", names(data_rows))
public_text <- paste(vapply(data_files, function(path) paste(readLines(path,
  warn = FALSE, encoding = "UTF-8"), collapse = "\n"), character(1L)),
  collapse = "\n")
privacy_ok <- !grepl("SRA701877|SRS3279688|Dear Dr|reviewer comment",
  public_text, ignore.case = TRUE) && nrow(data_manifest) == 11L &&
  !any(truth(data_manifest$sample_identifiers_published)) &&
  !any(truth(data_manifest$confidential_review_text)) &&
  !any(grepl("\\.pdf$", c(root_files, names(data_rows)), ignore.case = TRUE))
add_check("privacy_and_pdf_exclusion", privacy_ok,
  "fixed sample identity, confidential comments, and PDFs absent from public bundle")

prov_ok <- nrow(provenance) == 1L && provenance$figures_rendered == 7L &&
  provenance$figure8_status == "deferred_cross_stage_estimand" &&
  provenance$scientific_distance_grid == "none_exact_or_error_controlled" &&
  provenance$display_grid_only_figure3 == 600L &&
  !truth(provenance$new_ph) && !truth(provenance$new_data) &&
  !truth(provenance$p_values) && !truth(provenance$rankings) &&
  !truth(provenance$confidential_review_published)
add_check("scope_and_figure8_deferral", prov_ok,
  "no new PH/data/ranking; figure 8 explicitly deferred")

production_files <- sort(list.files(production_dir, recursive = TRUE,
  full.names = FALSE))
repeat_files <- sort(list.files(repeat_dir, recursive = TRUE,
  full.names = FALSE))
repeat_ok <- identical(production_files, repeat_files)
if (repeat_ok) {
  repeat_ok <- all(vapply(production_files, function(path) identical(
    tolower(sha(file.path(production_dir, path))),
    tolower(sha(file.path(repeat_dir, path)))), logical(1L)))
}
add_check("byte_identical_repeat", repeat_ok,
  paste(length(production_files),
        "production files independently rerendered byte-identically"))

report <- do.call(rbind, checks)
write_csv(report, file.path(validation_dir, "mv08a-validation-checks.csv"))
decision <- data.frame(
  contract_id = "mv08a_validation_decision_v1",
  checks_passed = sum(report$pass), checks_total = nrow(report),
  all_checks_pass = all(report$pass),
  decision = if (all(report$pass)) "authorize_author_review_of_figures" else
    "block_author_review",
  publication_submission_authorized = FALSE,
  external_data_download_authorized = FALSE,
  new_ph_authorized = FALSE, expected_head = expected_head,
  stringsAsFactors = FALSE)
write_csv(decision, file.path(validation_dir, "mv08a-validation-decision.csv"))
manifest_files <- sort(list.files(validation_dir, full.names = TRUE))
validation_manifest <- data.frame(
  contract_id = "mv08a_validation_manifest_v1",
  filename = basename(manifest_files),
  sha256 = vapply(manifest_files, sha, character(1L)),
  bytes = as.numeric(file.info(manifest_files)$size),
  stringsAsFactors = FALSE)
write_csv(validation_manifest,
  file.path(validation_dir, "mv08a-validation-manifest.csv"))
if (!all(report$pass)) {
  stop("MV8-A independent validation failed: ",
       paste(report$check_id[!report$pass], collapse = ", "))
}
message(sum(report$pass), "/", nrow(report),
        " MV8-A checks pass; author figure review is authorized")
