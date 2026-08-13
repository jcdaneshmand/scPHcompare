# MV5-BF Rust candidate-certification workflow contract v1

Date frozen: 2026-08-13
Status at freeze: local implementation only; no hosted run or publication authorized

## Purpose

Implement the MV5-BE distribution decision as a nonpublishing, auditable CI
candidate gate. The workflow may prove that the accepted Rust landscape kernel
builds and behaves correctly on the four required native host targets, but it
must not publish a release, install an accelerator for users, or change the
canonical R execution path.

## Immutable boundaries

- The R package remains the canonical implementation and default path.
- `.github`, `scripts`, `rust`, `docs`, and all accelerator evidence remain
  excluded from the R source tarball under the existing `.Rbuildignore` rules.
- No Rust binary is committed to the package or repository.
- No Cargo or Rust toolchain is required to build, install, load, or use the R
  package without the accelerator.
- The workflow has read-only repository permissions and no release event,
  release action, mutable `latest` URL, external upload, package downloader, or
  installation step.
- Only public analytical/generated fixtures enter hosted CI. Private diagrams,
  result corpora, PDFs, RDS files, and `tmp/mv05*` evidence are forbidden.
- Candidate artifacts expire after seven days and are explicitly marked
  unpublished and candidate-only.

## Required native matrix

| Host runner | Rust target | Candidate library |
|---|---|---|
| `ubuntu-22.04` | `x86_64-unknown-linux-gnu` | `.so` |
| `windows-2022` | `x86_64-pc-windows-msvc` | `.dll` |
| `macos-14` | `aarch64-apple-darwin` | `.dylib` |
| `macos-15-intel` | `x86_64-apple-darwin` | `.dylib` |

Each row must use R 4.4.1 and Rust 1.97.1, verify that the native rustc host
matches the declared target, and use Cargo's locked dependency graph.

## Per-host acceptance sequence

1. Check Rust formatting, run Rust unit tests, and run Clippy with warnings
   denied.
2. Produce two clean release builds and require byte-identical candidate
   libraries.
3. Compile and run a native C ABI harness against the dynamic library.
4. Record native dynamic dependencies with `ldd`, `dumpbin`, or `otool`.
5. Call the library from R on three exact analytical fixtures, including a
   narrow feature, at absolute error at most `1e-12`.
6. Prove canonical-R fallback for both a missing and a corrupt library.
7. Build the exact R source package and require `R CMD check` success with the
   accelerator absent and present.
8. Emit a normalized candidate manifest containing the commit, runner, target,
   Rust identities, source and lock hashes, candidate hash and size, dependency
   inventory hash, fixture hash, R-check log hashes, and explicit
   `candidate-only`, `private_data_used=false`, and
   `published_release=false` declarations.
9. Upload only the short-lived candidate and evidence directory.

## Local implementation gate

Before any hosted run, the repository must pass all of the following locally:

- YAML parsing and exact matrix/static-policy validation;
- PowerShell parsing for the Windows harness driver;
- Python byte-code compilation for the manifest builder;
- the existing Rust format, unit, strict-Clippy, and clean-build checks;
- Linux native C ABI execution and dependency inventory;
- all five R analytical/fallback fixtures;
- candidate-manifest generation and content validation using local stand-in
  R-check logs;
- exact R-package source-tree invariance against the last accepted package
  input tree; and
- an independent validator covering the frozen contract.

## Explicit nonclaims

A local Linux run cannot certify Windows or either macOS host. Even a passing
Ubuntu 22.04 job does not meet the MV5-BE glibc 2.17-compatible *release*
baseline. This workflow therefore creates candidate evidence only. Hosted
four-platform execution, Linux release-baseline construction, signed
attestation, GitHub release publication, installer work, and production Rust
adoption remain separate future gates requiring authorization and evidence.
