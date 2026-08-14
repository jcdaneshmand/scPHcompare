# MV5-BG pre-hosted Rust CI hardening

Date: 2026-08-13
Status: local hardening accepted; push and hosted execution remain paused

## Outcome

MV5-BG continued locally after the owner paused pushing. It found and corrected
five CI/provenance defects that could be diagnosed without a hosted run. No
scientific implementation, landscape definition, R API, package default,
kernel source, private input, GitHub branch, PR, workflow run, artifact upload,
or release changed.

The workflow now restores the exact 265-record `renv.lock` baseline rather than
asking a generic DESCRIPTION installer to reconstruct GitHub and historical
Bioconductor dependencies. The active Ubuntu R 4.4.1 library matched 265/265
versions and all 29 commit-pinned records. The cache identity includes OS,
architecture, R version, lockfile, renv bootstrap, and restoration script. The
native timeout is 180 minutes, matching the established package CI risk scale.

All four external actions are pinned to immutable commits read directly from
their official Git repositories. The nonrelease/private-data job must pass
before native jobs start, checkout credentials are not persisted, and failed
jobs may retain seven-day partial diagnostics without converting missing
evidence into another failure.

## Source and data boundary

The first hardened guard revealed that the repository intentionally tracks
public generated RDS results; treating every tracked RDS as private was
incorrect. Native jobs now sparse-check out only the package, Rust, test, and CI
sources they need. That selection contains 550 tracked files and zero PDF/RDS
files. The guard still rejects `docs/private`, PDFs, `example_run.r`, `tmp`, and
candidate-evidence paths if tracked, and separately proves the sparse worktree
contains no PDF/RDS files.

The sparse selection produced the same normalized 219-file R source package
tree as MV5-BF: `b4079f4ba7526ca711199eef71b2479a90581dad`.
Because workflow and helper scripts remain build-excluded and that content tree
is identical, MV5-BF's two authoritative absent/present `Status: OK` logs remain
applicable; they were hash-checked and reused rather than rerunning an unchanged
eleven-minute double package check.

GitHub Checkout's official v6 documentation confirms both single-file and
multi-path sparse checkout inputs and recommends `contents: read`. It also says
the partial-clone filter overrides sparse checkout, so this workflow uses
sparse checkout alone. Reference checked 2026-08-13:
<https://github.com/actions/checkout/blob/v6/README.md>.

## Reproducible provenance

Candidate source hashes now use canonical bytes from the exact Git revision,
not OS-transformed working-tree bytes. The manifest binds Cargo.lock,
Cargo.toml, Rust source, public ABI header, R shim, toolchain, target, command,
binary, dependency inventory, fixtures, and R-check logs.

Raw `ldd`, `otool`, or `dumpbin` output remains available for diagnosis, while
a platform-specific normalizer produces the inventory used for hashing. Three
analytical normalizer fixtures passed. Two independent Linux `ldd` executions
with address-varying raw output produced the same normalized SHA-256:
`b0441ee1cd02d495c77b4a6c43ea557da3065dec34cc0b433cea0f824d20ff48`.

## Local parity results

- Hardened static validation: 31/31.
- Independent YAML four-row matrix parse: pass.
- Windows PowerShell parse: pass.
- Exact active R library: 265/265 records, including 29/29 commit pins.
- Rust format, strict Clippy, and unit tests: pass; 4/4 unit tests.
- Two clean release builds: byte-identical at accepted SHA-256 `51d3fca...`.
- Linux native C ABI harness: pass.
- Analytical/missing/corrupt R fixtures: 5/5.
- Sparse R-package invariance: 219/219 files, identical tree.
- Sparse CI source selection: 550 files, zero PDF/RDS.
- Git push, PR, hosted run, upload, and release count: zero.

## Decision and next gate

The local MV5-BG hardening is accepted. Windows and macOS remain implementation
only, and the Linux glibc 2.17 release baseline remains open. The branch should
remain local until the owner again authorizes pushing. At that point, the next
gate is the first hosted four-platform candidate run—not a release or default
change.
