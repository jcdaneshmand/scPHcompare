workflow <- paste(readLines(
  ".github/workflows/rust-accelerator-certification.yml", warn = FALSE
), collapse = "\n")
specification <- paste(readLines(
  "docs/specifications/MV05BF_RUST_CANDIDATE_CERTIFICATION_WORKFLOW_V1.md",
  warn = FALSE
), collapse = "\n")
manifest_script <- paste(readLines(
  "scripts/mv05bf_build_candidate_manifest.py", warn = FALSE
), collapse = "\n")
fixture_script <- paste(readLines(
  "scripts/mv05bf_rust_ci_fixture.R", warn = FALSE
), collapse = "\n")
windows_driver <- paste(readLines(
  "scripts/mv05bf_run_windows_harness.ps1", warn = FALSE
), collapse = "\n")
rbuildignore <- readLines(".Rbuildignore", warn = FALSE)

required_workflow <- c(
  'permissions:\n  contents: read',
  "ubuntu-22.04", "windows-2022", "macos-14", "macos-15-intel",
  "x86_64-unknown-linux-gnu", "x86_64-pc-windows-msvc",
  "aarch64-apple-darwin", "x86_64-apple-darwin",
  'RUST_TOOLCHAIN: "1.97.1"', "--locked", "-D warnings",
  "clean-build-a.sha256", "Clean release builds are not byte-identical",
  "mv05bc_ffi_sanitizer.c", "mv05bf_run_windows_harness.ps1",
  "mv05bf_rust_ci_fixture.R",
  "Check exact R source package with accelerator absent and present",
  "r-check-absent.log", "r-check-present.log",
  "retention-days: 7", "no-release-guard"
)
forbidden_workflow <- c(
  "contents: write", "gh release", "create-release", "release-action",
  "tmp/mv05", "docs/private", "\\.rds", "\\.pdf"
)
candidate_workflow <- sub("(?s)\n  no-release-guard:.*$", "", workflow,
                          perl = TRUE)

checks <- c(
  all(vapply(required_workflow, grepl, logical(1L), x = workflow,
             fixed = TRUE)),
  !any(vapply(forbidden_workflow, grepl, logical(1L),
              x = candidate_workflow, fixed = TRUE)),
  !grepl("uses:[^\n]*release", workflow, perl = TRUE),
  !grepl("^  release:", workflow, perl = TRUE),
  lengths(regmatches(workflow, gregexpr("runner: ", workflow,
                                        fixed = TRUE))) == 4L,
  grepl('certification_class": "candidate-only"', manifest_script,
        fixed = TRUE),
  grepl('private_data_used": "false"', manifest_script, fixed = TRUE),
  grepl('published_release": "false"', manifest_script, fixed = TRUE),
  grepl("r_package_checks_absent_and_present_passed", manifest_script,
        fixed = TRUE),
  grepl("error <= 1e-12", fixture_script, fixed = TRUE),
  grepl("missing_library_fallback", fixture_script, fixed = TRUE),
  grepl("corrupt_library_fallback", fixture_script, fixed = TRUE),
  grepl("scripts\\mv05bf_windows_ffi_harness.c", windows_driver,
        fixed = TRUE),
  !grepl("tests\\native\\mv05bf_windows_ffi_harness.c", windows_driver,
         fixed = TRUE),
  "^rust-candidate-evidence$" %in% rbuildignore,
  grepl("A local Linux run cannot certify Windows", specification,
        fixed = TRUE),
  grepl("does not meet the MV5-BE glibc 2.17-compatible", specification,
        fixed = TRUE),
  grepl("R package remains the canonical implementation", specification,
        fixed = TRUE),
  grepl("no hosted run or publication authorized", specification,
        fixed = TRUE)
)

if (!all(checks)) {
  stop("MV5-BF static validation failed at checks: ",
       paste(which(!checks), collapse = ", "), call. = FALSE)
}
cat("MV5-BF static contract validation: ", sum(checks), "/",
    length(checks), " passed\n", sep = "")
