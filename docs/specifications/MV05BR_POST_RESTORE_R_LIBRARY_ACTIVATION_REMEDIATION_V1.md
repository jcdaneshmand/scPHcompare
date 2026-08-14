# MV5-BR post-restore R library activation remediation contract v1

Date frozen: 2026-08-13
Status at freeze: prospective, forward-only CI correction; all PRs remain draft and unmerged

## Triggering evidence

P02 hosted R run `31745847174` restored and independently verified all 265
records in `renv.lock` at the explicit library
`/home/runner/work/scPHcompare/scphcompare-renv-library`. The immediately
following `R CMD build` and `R CMD check` step received that path as
`R_LIBS_USER`, yet `renv` warned that locked packages were not installed and
the check reported every declared dependency unavailable.

The committed `renv/activate.R` establishes the cause: `RENV_PATHS_LIBRARY` is
treated as a library *root* and receives a platform/profile suffix during
project activation, whereas the restore script intentionally passes the exact
unsuffixed directory as `renv::restore(library = ...)`. The autoloader also
documents `R_LIBS` as the subprocess signal that prevents redundant project
activation. Exporting `R_LIBS_USER` alone therefore cannot bind later R
processes to the exact directory that was verified.

## Prospective correction

Only after `scripts/restore_locked_dependencies.R` exits successfully:

1. export the exact verified `RENV_PATHS_LIBRARY` as both `R_LIBS` and
   `R_LIBS_USER` through `GITHUB_ENV`; and
2. export `RENV_CONFIG_AUTOLOADER_ENABLED=FALSE` through `GITHUB_ENV` so later
   `R`, `Rscript`, and `R CMD` processes use the already restored library
   instead of recomputing an unused project-library path.

Apply this correction to the general R-package workflow at P01 and propagate
it by ordinary parent merges through P08. Apply the same post-restore exports
to P08's separate Rust candidate workflow on every native matrix host. Existing
step-local `R_LIBS_USER` assignment may remain as an explicit nested-check
guard, but it is not the primary binding.

## Frozen boundaries

- Do not alter the restore script, lockfile, dependency versions/sources,
  platform matrix, R 4.4.1, Rust 1.97.1, package source, fixtures, landscapes,
  tolerances, inputs, results, evidence, APIs, defaults, or scientific slices.
- Keep the exact restore and 265-record identity check ahead of every export.
- Do not disable renv during restoration; disable only redundant activation in
  later steps after the exact library has passed verification.
- Preserve read-only workflow permissions, the private-data/no-release guard,
  all R/Rust/ABI/fallback/package gates, and short artifact retention.
- Advance every affected branch by normal fast-forward push only. No force,
  rewrite, merge, tag, release, binary publication, default change, DOI action,
  or manuscript claim is authorized.

## Acceptance gate

1. Both workflow files parse and pass pinned `actionlint` v1.7.12.
2. Static inspection proves the exports occur only after successful exact
   restore and before all build/check/fixture processes.
3. Each updated parent remains an ancestor of its child; every scientific slice
   stays patch-identical to its retained pre-publication ref.
4. Final general R-package runs for P01 through P08 all pass.
5. Final P08 Rust no-release guard and all four native jobs pass, including
   Windows raw DLL byte identity and both accelerator-absent/present R checks.
6. Final heads, trees, run IDs, conclusions, PR metadata, and unchanged `main`
   are recorded and independently rechecked before this goal closes.

## Nonclaims

This remediation makes hosted processes consume the dependency library they
already restored and verified. It does not change the reproducibility target,
certify a release, authorize Rust adoption, add biological evidence, or make a
scientific/manuscript claim.
