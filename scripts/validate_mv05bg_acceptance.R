summary <- utils::read.csv(
  "docs/audits/mv05bg-hardening-summary-2026-08-13.csv",
  stringsAsFactors = FALSE, check.names = FALSE
)
actions <- utils::read.csv(
  "docs/audits/mv05bg-action-pins-2026-08-13.csv",
  stringsAsFactors = FALSE, check.names = FALSE
)
candidate <- utils::read.csv(
  "docs/audits/mv05bg-local-linux-candidate-manifest-2026-08-13.csv",
  stringsAsFactors = FALSE, check.names = FALSE
)
diagnostics <- utils::read.csv(
  "docs/audits/mv05bg-diagnostic-revisions-2026-08-13.csv",
  stringsAsFactors = FALSE, check.names = FALSE
)
decision <- utils::read.csv(
  "docs/audits/mv05bg-continuation-decision-2026-08-13.csv",
  stringsAsFactors = FALSE, check.names = FALSE
)
report <- paste(readLines(
  "docs/audits/MV05BG_PREHOSTED_RUST_CI_HARDENING_2026-08-13.md",
  warn = FALSE
), collapse = "\n")

value <- function(measure) summary$value[summary$measure == measure]
checks <- c(
  nrow(summary) == 16L,
  all(summary$accepted),
  identical(value("native_matrix_rows"), "4"),
  identical(value("job_timeout_minutes"), "180"),
  identical(value("locked_R_records_matched"), "265/265"),
  identical(value("commit_pinned_R_records_matched"), "29/29"),
  identical(value("sparse_PDF_or_RDS_files"), "0"),
  identical(value("normalized_R_package_files"), "219"),
  identical(value("normalized_R_package_tree"),
            "b4079f4ba7526ca711199eef71b2479a90581dad"),
  identical(value("Rust_unit_tests"), "4/4"),
  identical(value("R_candidate_fixtures"), "5/5"),
  identical(value("static_validation"), "31/31"),
  identical(value("hosted_runs"), "0"),
  identical(value("pushes"), "0"),
  nrow(actions) == 4L,
  all(nchar(actions$commit_sha) == 40L),
  all(actions$verified_before_edit),
  nrow(candidate) == 1L,
  identical(candidate$source_hash_basis, "git:HEAD"),
  identical(candidate$library_sha256,
            "51d3fca4808c6620c46c7b1100ddfad4753af02d37f87ec3eb278008846b160d"),
  identical(candidate$dependency_inventory_sha256,
            "b0441ee1cd02d495c77b4a6c43ea557da3065dec34cc0b433cea0f824d20ff48"),
  all(nchar(candidate[c(
    "cargo_lock_sha256", "cargo_manifest_sha256", "rust_source_sha256",
    "public_abi_header_sha256", "r_shim_sha256"
  )]) == 64L),
  identical(candidate$private_data_used, "false"),
  identical(candidate$published_release, "false"),
  nrow(diagnostics) == 5L,
  !any(diagnostics$scientific_contract_changed),
  isTRUE(decision$local_hardening_accepted),
  !isTRUE(decision$push_performed),
  !isTRUE(decision$hosted_certification_complete),
  !isTRUE(decision$release_authorized),
  !isTRUE(decision$R_default_changed),
  !isTRUE(decision$private_data_used),
  grepl("branch should\\s+remain local", report, perl = TRUE),
  grepl("Windows and macOS remain implementation\\s+only", report,
        perl = TRUE),
  grepl("not a release or default\\s+change", report, perl = TRUE)
)
if (!all(checks)) {
  stop("MV5-BG acceptance validation failed at checks: ",
       paste(which(!checks), collapse = ", "), call. = FALSE)
}
cat("MV5-BG acceptance validation: ", sum(checks), "/",
    length(checks), " passed\n", sep = "")
