# MV5-BP R package CI library-visibility remediation v1

## Purpose

MV5-BP makes the general R package workflow's post-restore steps see the exact
library that `restore_locked_dependencies.R` just restored and verified. It
changes only CI environment propagation and does not change the lockfile,
dependencies, package source, Rust source, landscapes, data, or analyses.

## Hosted evidence and root cause

P07 replacement R run `31744074169` completed the full 365-record restore after
MV5-BN fixed its Ubuntu system libraries. `R CMD build` succeeded, but
`R CMD check` then reported the restored Imports and Suggests as unavailable.
The job log showed `R_LIBS_USER` still pointed at setup-r's temporary library,
while the verified packages reside at `RENV_PATHS_LIBRARY`.

The restore script sets its selected library inside its own R process. That
does not mutate the environment of later shell steps. Those later R processes
therefore cannot see the exact verified library unless the workflow exports it.

## Exact correction

After the restore command succeeds, append exactly:

`printf 'R_LIBS_USER=%s\n' "${RENV_PATHS_LIBRARY}" >> "${GITHUB_ENV}"`

GitHub applies that binding only to subsequent steps, so the build/check and
realistic-fixture steps inherit the already-verified library. No other source
line or repository file may change.

## Forward-only stack and validation contract

Add one ordinary P01 commit, propagate it through P02 to P08 with ordinary
merge commits, and push normal fast-forwards only after remote-head checks.
Require official `actionlint` v1.7.12 to report zero diagnostics, preserve all
original and prior remediation commits in ancestry, retain patch-identical
scientific slices, and constrain all new candidate deltas to the R workflow.
Update all draft bodies and the ledger, then monitor replacement R checks on
all eight PRs and the P08 Rust matrix.

MV5-BP does not authorize merge, branch deletion, force-push, tag, release,
binary publication, Rust adoption/default, DOI/Zenodo action, new calculation,
or manuscript claim.
