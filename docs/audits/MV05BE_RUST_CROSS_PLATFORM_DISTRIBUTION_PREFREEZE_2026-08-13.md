# MV5-BE Rust cross-platform distribution prefreeze

Date: 2026-08-13
Status: accepted prefreeze; implementation and production adoption remain closed

## Repository finding

The existing package CI has one Ubuntu 22.04/R 4.4.1 job. It contains no Rust
toolchain, cross-platform matrix, artifact build, or release path. The Rust
crate is deliberately excluded from the R source tarball, and the prototype R
wrapper requires an explicit library path. Thus MV5-BD's complete numerical
acceptance is not yet a distribution mechanism.

## Frozen design

Adopt a two-track model. `scPHcompare` remains an R-canonical, portable source
package with no embedded Rust executable, no Cargo requirement, and no network
activity during install/load/startup. Separately, immutable opt-in accelerator
artifacts may later be published for Linux x86-64, Windows x86-64, macOS ARM64,
and macOS x86-64. Unsupported targets use R.

Each artifact must be tied to the exact source, Cargo.lock, toolchain, target,
build command, native dependencies, byte count, SHA-256, Git commit/tag, test
results, and signed or attested manifest. A future installer must require an
explicit version or local file, verify provenance and checksum before atomic
placement in the R user cache, perform an ABI/load self-test, and never use a
mutable `latest` selector. No background update is allowed.

Runtime remains explicitly selected and default-off. Absence, mismatch,
corruption, load failure, nonzero status, or invalid output returns to canonical
R. The private 408-result corpus remains local; platform CI receives only the
same source and public analytical/generated fixtures.

## Platform boundary

Rust officially provides host tools for the required Linux x86-64, Windows
x86-64 MSVC, and macOS ARM64 targets; macOS x86-64 remains an official supported
host target at a lower tier. Native builds and R runtime calls are nevertheless
required because target availability alone does not certify this ABI or R
integration. Linux release builds must use a glibc 2.17-compatible baseline;
the Ubuntu 22.04 prototype binary is evidence, not a portable release artifact.

## R distribution rationale

Current *Writing R Extensions* guidance says source packages should avoid
binary executables because they are nonportable and a security risk, and states
that CRAN will not accept source submissions containing them. Current CRAN
policy also says maintainer-built binary R packages are not accepted and that
platform binaries are produced by CRAN. Therefore a compiled Rust library
cannot simply be committed inside the R source tarball. An opt-in external
accelerator preserves a clean R source package while allowing carefully
verified GitHub releases for users who want the speedup.

Official references checked 2026-08-13:

- <https://doc.rust-lang.org/rustc/platform-support.html>
- <https://doc.rust-lang.org/cargo/reference/build-scripts.html>
- <https://cran.r-project.org/doc/manuals/r-release/R-exts.html>
- <https://cran.r-project.org/web/packages/policies.html>

## Next gate

MV5-BF may implement a nonpublishing cross-platform certification workflow:
four required target builds, unit/lint/ABI/R-fixture/fallback tests, native
dependency inventories, repeat-build evidence, and candidate checksums. It may
upload short-lived CI artifacts for inspection but may not create a GitHub
release, add a package downloader, change the R default, or declare production
adoption.
