# MV5-BN hosted R system-dependency remediation v1

## Purpose

MV5-BN provisions the Ubuntu system libraries that the exact committed R
lockfile requires before `renv::restore()`. It is CI infrastructure only and
does not change the lockfile, R packages, R source, Rust source, landscape
definitions, data, analyses, results, or scientific review boundaries.

## Hosted evidence and root cause

The retained R package runs for P02, P05, P06, and P07 failed identically while
loading restored `igraph`: `libglpk.so.40` was absent. Before restoration,
`renv` explicitly reported these missing Ubuntu prerequisites:

- `libcurl4-openssl-dev` for `curl`;
- `libglpk-dev` for `igraph` and `leidenAlg`;
- `libgmp3-dev` for `TDA`;
- `libhdf5-dev` for `hdf5r`;
- `pandoc` for `knitr` and `rmarkdown`.

The missing GLPK runtime caused the observed failure. The remaining packages
are provisioned at the same time because they are reported requirements of the
same frozen lockfile and would otherwise become later cold-run failures.

## Exact correction

1. In the Ubuntu-only R package workflow, add one pre-restore step that runs
   `sudo apt-get update` and installs the five reported packages with
   `--no-install-recommends`.
2. In P08's cross-platform Rust certification workflow, add the equivalent
   step guarded by `if: runner.os == 'Linux'` so only its Ubuntu candidate runs
   `apt-get`; macOS and Windows candidates remain unchanged.
3. Make no other workflow or repository source changes.

## Forward-only stack construction

Add the R-workflow correction as one ordinary P01 commit. Propagate the updated
parent through P02 to P07 with ordinary merge commits. Merge updated P07 into
the current P08 head, preserving all MV5-BL and MV5-BM commits, then add the
P08-only Linux prerequisite step as a separate ordinary commit. Push only
normal fast-forwards after confirming each expected remote head. Never amend,
reset, rebase, replace, or force-push.

## Validation and stop conditions

- Require official `actionlint` v1.7.12 to report zero diagnostics.
- Require all original candidate commits and all prior remediation heads to
  remain ancestors of their current branches.
- Require P02 through P07 scientific slice patches to remain identical.
- Require P08's scientific slice patch to remain identical.
- Require remediation commits to touch workflow files only.
- Monitor the replacement R checks on all eight PRs and the replacement P08
  no-release guard plus four-platform Rust matrix.
- Diagnose any additional hosted failure before another source change.

MV5-BN does not authorize merge, branch deletion, force-push, tag, release,
binary publication, Rust adoption/default, DOI/Zenodo action, new calculation,
or manuscript claim.
