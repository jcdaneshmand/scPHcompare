# MV5-BO Rust CI check-library visibility remediation v1

## Purpose

MV5-BO makes P08's nested `R CMD check` processes see the exact library already
restored and verified by the Rust candidate workflow. It changes only a
PowerShell environment binding; it does not add, remove, or update any package
and does not change R, Rust, landscape, topology, data, or analysis source.

## Hosted evidence and root cause

Replacement Rust run `31744074509`, macOS ARM64 job `94594349625`, passed exact
dependency restoration, pinned Rust installation, tests, formatting, clippy,
two byte-identical release builds, the native C ABI harness, dependency
normalization, and all five R/Rust fixtures. The external R application also
executed successfully after MV5-BM. `R CMD check` then reported the restored
package dependencies as unavailable.

The job restores the lockfile to `RENV_PATHS_LIBRARY`. The package checks run
from nested evidence directories, outside the project directory that activates
renv. Their child R processes therefore retain the setup action's temporary
`R_LIBS_USER` and do not include the verified renv library.

## Exact correction

In the existing P08 PowerShell package-check step, immediately after resolving
the external R application, set:

`$env:R_LIBS_USER = $env:RENV_PATHS_LIBRARY`

The subsequent `R CMD build` and both absent/present `R CMD check` invocations
then inherit the same exact library verified earlier in the job. No other line
or repository file may change.

## Forward-only and validation contract

Add one ordinary commit after the current P08 head. Require official
`actionlint` v1.7.12 to report zero diagnostics, all previous P08 remediation
heads and the original P08 scientific candidate to remain ancestors, and the
new commit to change only the Rust certification workflow. Confirm the P08
non-workflow scientific slice remains patch-identical, verify the remote P08
head before a normal non-force push, update the draft body and ledger, and
monitor the replacement P08 R and Rust runs.

MV5-BO does not authorize merge, branch deletion, force-push, tag, release,
binary publication, Rust adoption/default, DOI/Zenodo action, new calculation,
or manuscript claim.
