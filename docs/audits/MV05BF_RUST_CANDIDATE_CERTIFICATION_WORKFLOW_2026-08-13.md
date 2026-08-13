# MV5-BF Rust candidate-certification workflow

Date: 2026-08-13
Status: local implementation accepted; hosted cross-platform certification and release remain closed

## Outcome

MV5-BF implemented the nonpublishing four-host candidate workflow frozen in
MV5-BE. It covers native Linux x86-64, Windows x86-64 MSVC, macOS ARM64, and
macOS x86-64. Every row pins Rust 1.97.1 and R 4.4.1, verifies the native Rust
host, uses the locked Cargo graph, runs format/unit/strict-Clippy checks,
requires two clean byte-identical release builds, executes a native C ABI
harness, inventories dynamic dependencies, exercises R analytical and failure
fixtures, and checks the exact R package with the accelerator absent and
present. Evidence artifacts expire after seven days. The workflow has
read-only repository permissions and no push or release trigger.

The workflow was parsed independently as YAML with four matrix rows. Its
Windows PowerShell driver and Python manifest builder also parsed cleanly. The
nonpublishing/private-data guard was executed against the exact workflow after
correcting an initial self-match defect. The final static contract validator
passed 19/19 checks.

## Local native result

The isolated official Rust 1.97.1 toolchain reported host
`x86_64-unknown-linux-gnu`. Rust unit tests passed 4/4; formatting and strict
Clippy passed; the native C ABI harness passed; and two fresh release builds
were byte-identical. The library SHA-256 remained
`51d3fca4808c6620c46c7b1100ddfad4753af02d37f87ec3eb278008846b160d`,
the same accepted MV5-BC binary. Its local dynamic dependencies were only
`libgcc_s`, `libc`, the loader, and Linux's virtual DSO.

R passed all five candidate fixtures. The three analytical results included a
sign-changing pair and a narrow feature; maximum absolute error was
`1.76e-24`. Missing and corrupt libraries both produced status 9001 and used
the canonical R fallback. The normalized manifest records the toolchain,
source, Cargo.lock, binary, dependency, fixture, and package-check hashes and
declares `candidate-only`, `private_data_used=false`, and
`published_release=false`.

## R package boundary

Both authoritative `R CMD check --no-manual --no-build-vignettes` executions
returned `Status: OK`, first with Rust absent and then with the candidate path
present. The 219-file source content tree matched the accepted package exactly
after removing only R's generated `Packaged:` timestamp and normalizing regular
file modes. The old WSL-on-Windows tarball had spuriously marked all files 755;
clean Git staging correctly produced 644 files and passed R's permission check.

A final review found that the root `rust-candidate-evidence` directory needed
an explicit `.Rbuildignore` entry. A clean build containing a dummy candidate
proved that the directory was excluded and that the normalized content tree
remained `b4079f4ba7526ca711199eef71b2479a90581dad`. Thus neither the Rust crate,
candidate binary, workflow, scripts, nor audit material enters the R source
package, and ordinary R installation remains Rust- and network-free.

## Diagnostic corrections

The source-freeze trail records six non-scientific implementation revisions:
three validator assertions were narrowed to test their intended files, the
safety guard was made non-self-matching, the CI evidence directory was added
to `.Rbuildignore`, and the validator was extended to assert that exclusion.
A first direct package build was also stopped after it began traversing the
large ignored `tmp` history; its exact process group was terminated and its
826 MB temporary copy removed. Authoritative package checks used a clean Git
archive staging directory. No private data entered any build or tracked output.

## Decision and nonclaims

The MV5-BF local implementation is accepted. This does **not** certify Windows
or either macOS target, because those native jobs have not run. It also does
not make Ubuntu 22.04 a portable Linux release: MV5-BE's glibc 2.17-compatible
release baseline remains open. No artifact was pushed, uploaded, signed,
attested, or published; no installer exists; and the R engine remains the
default and scientific reference.

The next gate is an explicitly authorized push followed by the first hosted
four-platform candidate run. Only after all four native rows pass should the
project design the Linux glibc 2.17 build baseline and release attestation.
