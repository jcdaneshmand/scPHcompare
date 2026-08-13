# MV5-BG pre-hosted Rust CI hardening contract v1

Date frozen: 2026-08-13
Status at freeze: local-only; push, PR, hosted workflow, and release remain closed

## Purpose

Harden the accepted MV5-BF candidate-certification workflow against failure
modes discoverable before its first hosted run. This sprint may change only CI,
audit, and build-excluded helper files. It may not change the Rust numerical
kernel, R landscape implementation, scientific inputs, package API, defaults,
or published state.

## Prospective corrections

1. Replace generic DESCRIPTION-based dependency installation with the exact
   `renv.lock` restoration already used by package CI. The package includes
   GitHub and historical Bioconductor records that a generic installer cannot
   reconstruct from DESCRIPTION alone.
2. Cache the locked library per OS, architecture, R version, and lock/bootstrap
   hash; retain a machine-readable dependency verification report.
3. Increase the job timeout from 45 to 180 minutes because four native hosts
   may compile the locked R dependency graph before running two package checks.
4. Pin every third-party GitHub Action to the exact official commit observed at
   freeze, retaining a human-readable major-version comment.
5. Run the nonrelease/private-data gate before all native jobs and reject
   forbidden tracked paths in addition to forbidden workflow capabilities.
6. Expand pull-request path filters to cover DESCRIPTION, `.Rbuildignore`, the
   lock/bootstrap inputs, and all directly executed helper sources.
7. Preserve partial evidence after a failure without treating a missing
   evidence directory as an additional error.
8. Extend candidate manifests with `Cargo.toml`, public ABI header, and R shim
   hashes so the manifest describes the complete executable boundary in
   addition to the Git commit.

## Frozen external action identities

- `actions/checkout@v6`:
  `d23441a48e516b6c34aea4fa41551a30e30af803`
- `actions/cache@v4`:
  `0057852bfaa89a56745cba8c7296529d2fc39830`
- `actions/upload-artifact@v4`:
  `ea165f8d65b6e75b540449e92b4886f43607fa02`
- `r-lib/actions@v2`:
  `d3c5be51b12e724e68f33216ca3c148b66d5f0b6`

These identities were read directly from the official Git repositories with
`git ls-remote` before workflow edits.

## Local acceptance gate

- YAML and PowerShell parse successfully.
- Static validation proves the exact four hosts, immutable action pins,
  dependency restore/cache, 180-minute timeout, guard dependency, expanded
  paths, complete manifest hashes, and nonpublishing behavior.
- The Linux command-parity path passes using the isolated Rust 1.97.1 toolchain
  and already restored exact R 4.4.1 library.
- Candidate manifest generation passes with all required evidence and the new
  boundary hashes.
- The R source-package content tree remains identical and excludes candidate
  evidence.
- Public acceptance artifacts contain no private input, binary, PDF, RDS, or
  reviewer content.

## Explicit nonclaims

Local hardening cannot certify Windows or macOS, validate hosted runner images,
or satisfy the Linux glibc 2.17 release baseline. No push, PR, workflow run,
artifact upload, release, installer, or production Rust adoption is authorized.
