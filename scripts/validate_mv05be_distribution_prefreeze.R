inventory <- utils::read.csv(
  "docs/audits/mv05be-current-boundary-inventory-2026-08-13.csv",
  stringsAsFactors = FALSE, check.names = FALSE
)
targets <- utils::read.csv(
  "docs/audits/mv05be-target-matrix-2026-08-13.csv",
  stringsAsFactors = FALSE, check.names = FALSE
)
decision <- utils::read.csv(
  "docs/audits/mv05be-distribution-decision-2026-08-13.csv",
  stringsAsFactors = FALSE, check.names = FALSE
)
specification <- paste(readLines(
  "docs/specifications/MV05BE_RUST_CROSS_PLATFORM_DISTRIBUTION_PREFREEZE_V1.md",
  warn = FALSE
), collapse = "\n")

checks <- c(
  nrow(inventory) == 5L,
  identical(inventory$item, c(
    "r_package_ci", "Rbuildignore", "R_prototype_wrapper", "Cargo_lock",
    "Rust_kernel_source"
  )),
  nrow(targets) == 6L,
  sum(as.logical(targets$first_release_required)) == 4L,
  sum(as.logical(targets$current_build_certified)) == 1L,
  sum(as.logical(targets$current_R_runtime_certified)) == 1L,
  all(targets$fallback_if_absent == "canonical_R"),
  all(c("linux", "windows", "macos") %in% targets$platform),
  all(c("x86_64-unknown-linux-gnu", "x86_64-pc-windows-msvc",
        "aarch64-apple-darwin", "x86_64-apple-darwin") %in%
      targets$rust_target[as.logical(targets$first_release_required)]),
  isTRUE(decision$two_track_distribution_accepted),
  !isTRUE(decision$R_source_package_contains_Rust_binary),
  !isTRUE(decision$ordinary_R_install_requires_Cargo),
  !isTRUE(decision$automatic_download_authorized),
  isTRUE(decision$explicit_opt_in_artifacts_specified),
  isTRUE(decision$checksum_and_provenance_required),
  !isTRUE(decision$R_default_changed),
  !isTRUE(decision$Rust_production_adoption_authorized),
  identical(decision$next_sprint,
            "MV5-BF_cross_platform_CI_certification_implementation"),
  grepl("tools::R_user_dir", specification, fixed = TRUE),
  grepl("never run", specification, fixed = TRUE),
  grepl("private data is never uploaded", specification, fixed = TRUE),
  grepl("R source tarball contains no executable Rust binary", specification,
        fixed = TRUE),
  grepl("default remains `engine = \"r\"`", specification, fixed = TRUE),
  grepl("MV5-BF may implement CI certification without publishing artifacts",
        specification, fixed = TRUE)
)

if (!all(checks)) {
  stop("MV5-BE validation failed at checks: ",
       paste(which(!checks), collapse = ", "), call. = FALSE)
}
cat("MV5-BE independent validation: ", sum(checks), "/",
    length(checks), " passed\n", sep = "")
