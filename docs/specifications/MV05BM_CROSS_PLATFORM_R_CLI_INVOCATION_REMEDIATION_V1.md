# MV5-BM cross-platform R CLI invocation remediation v1

## Purpose

MV5-BM corrects a P08-only hosted certification defect discovered after the
MV5-BL expression-context repair allowed the Rust matrix to execute. It changes
only how PowerShell invokes the R command-line application; it does not change
R, Rust, landscape, topology, data, analysis, or package source.

## Hosted evidence and root cause

P08 Rust run `31742331945`, macOS ARM64 job `94588557542`, successfully
completed exact dependency restoration, pinned Rust installation, unit tests,
documentation tests, formatting, clippy, two byte-identical release builds,
the native C ABI harness, dependency normalization, and all five R/Rust
fixtures. The package-check step then failed before `R CMD build` executed:

`Invoke-History: A positional parameter cannot be found that accepts argument 'build'.`

PowerShell's built-in `R` alias resolves to `Invoke-History` before external
applications. Therefore bare `R CMD ...` is not a portable PowerShell command,
even when the R application is installed and `Rscript` works.

## Exact correction

Within the existing `Check exact R source package with accelerator absent and
present` PowerShell step:

1. resolve the external R application with
   `Get-Command R -CommandType Application -ErrorAction Stop`;
2. select its first `Source` path into `$rExecutable`;
3. replace `R CMD build .` with `& $rExecutable CMD build .`;
4. replace `R CMD check ...` with `& $rExecutable CMD check ...`.

No other workflow or repository source change is accepted.

## Forward-only publication rule

Add one ordinary commit after the current P08 head and push it normally only if
the current remote head is the expected MV5-BL head. Never amend, reset, rebase,
replace, or force-push. P01 through P07 remain unchanged. Update the P08 draft
body and remote-publication ledger with the new head/tree and replacement run.

## Validation and stop conditions

- Require official `actionlint` v1.7.12 to report zero diagnostics.
- Require the MV5-BL P08 head and the original P08 scientific candidate to
  remain ancestors of the new head.
- Require the new commit to change only the Rust certification workflow and
  only the portable R CLI invocation lines described above.
- Require the P08 scientific slice patch to remain identical.
- Monitor the replacement P08 R run, no-release guard, and four Rust candidates.
- Stop on any additional failure for root-cause diagnosis; do not merge.

MV5-BM does not authorize merge, branch deletion, force-push, tag, release,
binary publication, Rust adoption/default, DOI/Zenodo action, new calculation,
or manuscript claim.
