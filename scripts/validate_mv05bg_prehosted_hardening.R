workflow <- paste(readLines(
  ".github/workflows/rust-accelerator-certification.yml", warn = FALSE
), collapse = "\n")
manifest <- paste(readLines(
  "scripts/mv05bf_build_candidate_manifest.py", warn = FALSE
), collapse = "\n")
specification <- paste(readLines(
  "docs/specifications/MV05BG_PREHOSTED_RUST_CI_HARDENING_V1.md",
  warn = FALSE
), collapse = "\n")

pins <- c(
  "actions/checkout@d23441a48e516b6c34aea4fa41551a30e30af803",
  "actions/cache@0057852bfaa89a56745cba8c7296529d2fc39830",
  "actions/upload-artifact@ea165f8d65b6e75b540449e92b4886f43607fa02",
  "r-lib/actions/setup-r@d3c5be51b12e724e68f33216ca3c148b66d5f0b6"
)
paths <- c(
  '- ".Rbuildignore"', '- "DESCRIPTION"', '- "renv.lock"',
  '- "renv/**"', '- "scripts/restore_locked_dependencies.R"'
)
boundary_hashes <- c(
  "cargo_manifest_sha256", "public_abi_header_sha256", "r_shim_sha256"
)
tracked <- system2("git", "ls-files", stdout = TRUE)
forbidden_tracked <- grepl(
  "(^docs/private/|\\.pdf$|^example_run\\.r$|^tmp/|^rust-candidate-evidence/)",
  tracked
)

checks <- c(
  all(vapply(pins, grepl, logical(1L), x = workflow, fixed = TRUE)),
  !grepl("uses: actions/checkout@v", workflow, fixed = TRUE),
  !grepl("uses: actions/cache@v", workflow, fixed = TRUE),
  !grepl("uses: actions/upload-artifact@v", workflow, fixed = TRUE),
  !grepl("uses: r-lib/actions/setup-r@v", workflow, fixed = TRUE),
  grepl("needs: no-release-guard", workflow, fixed = TRUE),
  grepl("timeout-minutes: 180", workflow, fixed = TRUE),
  grepl("Restore and verify exact locked R dependencies", workflow,
        fixed = TRUE),
  grepl("scripts/restore_locked_dependencies.R", workflow, fixed = TRUE),
  grepl("dependency-bootstrap-manifest.csv", workflow, fixed = TRUE),
  grepl("Normalize native dependency inventory", workflow, fixed = TRUE),
  grepl("mv05bg_normalize_dependencies.py", workflow, fixed = TRUE),
  grepl("native-dependencies-raw.txt", workflow, fixed = TRUE),
  grepl("RENV_PATHS_LIBRARY", workflow, fixed = TRUE),
  grepl("hashFiles('renv.lock', 'renv/activate.R',", workflow, fixed = TRUE),
  !grepl("setup-r-dependencies", workflow, fixed = TRUE),
  all(vapply(paths, grepl, logical(1L), x = workflow, fixed = TRUE)),
  grepl("if: always()", workflow, fixed = TRUE),
  grepl("if-no-files-found: warn", workflow, fixed = TRUE),
  grepl("Forbidden tracked candidate input", workflow, fixed = TRUE),
  grepl("sparse-checkout-cone-mode: false", workflow, fixed = TRUE),
  grepl("Unexpected binary document/data in sparse worktree", workflow,
        fixed = TRUE),
  grepl("persist-credentials: false", workflow, fixed = TRUE),
  !any(forbidden_tracked),
  all(vapply(boundary_hashes, grepl, logical(1L), x = manifest,
             fixed = TRUE)),
  grepl("git_source_sha256", manifest, fixed = TRUE),
  grepl('"source_hash_basis"', manifest, fixed = TRUE),
  grepl("contents: read", workflow, fixed = TRUE),
  !grepl("contents: write", sub("  no-release-guard:.*", "", workflow),
         fixed = TRUE),
  grepl("local-only", specification, fixed = TRUE),
  grepl("No push, PR, workflow run", specification, fixed = TRUE)
)

if (!all(checks)) {
  stop("MV5-BG static validation failed at checks: ",
       paste(which(!checks), collapse = ", "), call. = FALSE)
}
cat("MV5-BG pre-hosted static validation: ", sum(checks), "/",
    length(checks), " passed\n", sep = "")
