# MV5-BE Rust cross-platform distribution prefreeze v1

Date: 2026-08-13
Authorization: safe prefreeze after MV5-BD `b486779`

## Decision

Use a two-track distribution model. The `scPHcompare` R source package remains
self-contained, R-canonical, and free of compiled Rust binaries and toolchain
requirements. The Rust accelerator is a separate, optional set of immutable
GitHub release artifacts installed only by explicit user action and activated
only by an explicit engine choice. Missing, incompatible, corrupt, or failing
artifacts always return to the canonical R engine.

This prefreeze does not edit CI, build or publish any platform artifact, add a
network installer, change package contents/defaults, or authorize production
Rust use.

## Required first-release targets

| Platform | Rust target | Library | Build constraint |
|---|---|---|---|
| Linux x86-64 | `x86_64-unknown-linux-gnu` | `libscph_landscape_kernel.so` | Build/test against glibc 2.17-compatible baseline, not the Ubuntu 22.04 host baseline |
| Windows x86-64 | `x86_64-pc-windows-msvc` | `scph_landscape_kernel.dll` | Native Windows host build; test from current x86-64 R/UCRT runtime |
| macOS ARM64 | `aarch64-apple-darwin` | `libscph_landscape_kernel.dylib` | Native Apple Silicon build with explicit deployment target |
| macOS x86-64 | `x86_64-apple-darwin` | `libscph_landscape_kernel.dylib` | Native Intel build with explicit deployment target |

Linux ARM64 and Windows ARM64 are later optional targets. Their absence must
select R rather than selecting a mismatched artifact.

## Immutable artifact contract

Each artifact is bound to engine version, Git tag and commit, Rust 1.97.1 host
toolchain identity, target triple, Cargo.lock SHA-256, Rust source SHA-256,
exact locked build command, library byte count and SHA-256, native dependency
inventory, unit/lint/ABI/R-fixture results, and build-run identity. Release tags
and files are immutable. A signed or GitHub-attested release manifest is
required before adoption; checksum alone detects corruption but does not prove
publisher identity.

Artifact names are
`scph-landscape-kernel-v1_<target-triple>.<so|dll|dylib>`. The manifest maps
exact supported OS/architecture combinations to one artifact and never guesses
across targets.

## CI and certification gates

Create a new Rust workflow later; do not overload current R-package CI. It will
run on pull requests for source tests and only publish on explicit version tags
after approval. Every required target must pass:

1. checksum-verified pinned minimal toolchain and `cargo --locked`;
2. format, 4/4 unit tests, strict Clippy, and release build;
3. exported-symbol/header ABI inspection and native C success/error harness;
4. R analytical/sign-crossing fixtures through the built library;
5. unavailable, wrong-target, corrupt, nonzero-status, and invalid-output R
   fallback tests with no partial installation;
6. two clean same-runner builds with source/result determinism recorded;
7. native dependency inventory and platform load test;
8. R package check unchanged with accelerator absent and optionally supplied.

The private 408-result MV5-BD scientific corpus remains bound to the exact Rust
source hash. Cross-platform certification uses the same source and public
analytical/generated fixtures; private data is never uploaded to CI.

## Explicit installation and runtime policy

A later installer may accept either a local artifact/manifest pair (air-gapped
path) or an explicit release version to download over HTTPS. It must never run
on package install, load, startup, examples, tests, or ordinary analysis. It
downloads to a temporary file, verifies manifest provenance and SHA-256 before
atomic placement under `tools::R_user_dir("scPHcompare", "cache")`, validates
the exact target and ABI, and preserves prior valid versions. No background
update or mutable `latest` URL is allowed.

The R package default remains `engine = "r"`. A later user must explicitly
select `engine = "rust"` or an equally explicit option after installation.
Rust failure emits certified diagnostics and recomputes in R. Cache identities
must include engine and ABI version while scientific identities remain tied to
the same landscape definition.

## R/CRAN boundary

The R source tarball contains no executable Rust binary. Rust source remains in
Git but build-excluded unless a later source-build design is separately
approved. Ordinary R installation cannot require Cargo or network access.
Examples and checks remain offline and R-canonical. This follows current R
guidance that source packages should avoid binary executables and current CRAN
policy that maintainer-supplied binary packages are not accepted.

## Adoption boundary

MV5-BF may implement CI certification without publishing artifacts. Publication
and any installer/runtime integration require a later explicit release gate.
Production adoption remains false until all four required targets and manifest,
provenance, fallback, and package gates pass. R removal or Rust-by-default is
not authorized by any success in this sequence.
