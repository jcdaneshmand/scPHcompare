# MV5-BC bounded Rust landscape-kernel prototype

Date: 2026-08-13
Status: accepted prototype; production adoption remains closed

## Outcome

MV5-BC implemented the frozen pair-bounded Rust prototype for the exact
persistence-landscape squared-L2 kernel. It preserves the dissertation-aligned
definition: all finite intervals, all consecutive active levels, zero padding,
separate H0/H1 calls, exact piecewise-linear integration, no grid, and no level
cap. R still owns scientific validation, provenance, certification, fallback,
and every public API. No package default changed.

The Rust implementation constructs the same exact critical-pair landscapes as
the accepted independent Persim oracle, but indexes births with a deterministic
segment tree instead of repeatedly scanning and removing Python list entries.
It then merges corresponding level functions and integrates their squared
difference with compensated summation. The crate has no external dependencies
and uses one thread.

## Frozen toolchain and build

The owner authorized an isolated Ubuntu toolchain installation. Official
`rustup-init` 1.29.0 was checksum-verified before execution; its SHA-256 is
`4acc9acc76d5079515b46346a485974457b5a79893cfb01112423c89aeb5aa10`.
Rust `1.97.1-x86_64-unknown-linux-gnu` was installed under ignored MV5-BC
staging without changing Ubuntu packages, the home directory, or shell files.
The locked release build used no external crates.

Two clean build/test/lint cycles passed 4/4 Rust tests and strict Clippy with
warnings denied. Both produced the identical shared-library SHA-256
`51d3fca4808c6620c46c7b1100ddfad4753af02d37f87ec3eb278008846b160d`.
The versioned C ABI success and error paths also passed an AddressSanitizer and
UndefinedBehaviorSanitizer harness.

## Scientific equivalence

All prototype equivalence tiers passed:

- Tier A: 3/3 analytical and sign-crossing fixtures.
- Tier B: 20/20 prospectively selected tractable exact R references; maximum
  squared-distance absolute error `2.22826201934367e-11`.
- Tier C: 12/12 frozen MV5-BA worst-depth H0/H1 certificates; maximum absolute
  error `3.9249430528604e-13`.

Normalized fixture, Tier B, Tier C, and environment outputs were byte-identical
across independent clean builds/runs. Private diagrams, pair identities,
measurements, and toolchain/build outputs remain under ignored `tmp/` storage.

## Performance and memory

On the same six worst-depth pairs used in MV5-BA, every Rust run was faster than
the accepted R timing. Speedups ranged from `599.287619047619x` to
`689.425531914891x`, with median `643.190881516968x`, far above the frozen `3x`
gate. Whole-run peak RSS was `240254976` bytes (about 229 MiB), below the 1 GiB
pair-bounded limit. This large gain is plausible because the R exact reference
enumerates pairwise tent-segment crossings and the Persim implementation spends
most time in repeated Python list scans; the Rust prototype uses indexed exact
critical-pair events.

## R integration and package health

The R shim is internal and optional. It canonicalizes finite interval arrays,
loads an explicitly supplied shared library, resolves a `NativeSymbolInfo` from
that exact DLL, and falls back to the accepted R calculation on absence, status
failure, lookup error, or invalid output. An unavailable library produced the
stable status `9001` and executed canonical fallback. Package-namespace testing
found and corrected an initial name-based `.C` resolution restriction; no Rust
result was accepted during the failure.

The complete source-loaded suite passed with Rust absent (1,411 passes, one
expected optional-library skip) and present (1,423 passes, no skips). The final
checked package-input tree `ef0d8d4edca662765407170ab5f72129940f6a42` produced
source tarball SHA-256
`6f28a83231a79935696f562b98020e8c911a3e197d9cdfdbe2ba418407a495bb`.
`R CMD check --no-manual --no-build-vignettes` returned `Status: OK` with Rust
absent and present. A first diagnostic check produced one NOTE for the obsolete,
ignored `.C(DUP=FALSE)` argument; it was removed before both authoritative
checks. This report and final log hashes are under build-excluded `docs/` and do
not alter that checked package input or tarball.

## Decision and next gate

The bounded prototype is accepted. This does **not** authorize production Rust
use: the R engine remains canonical, the Rust crate is excluded from R source
builds, and no public default or API changed. Added seeds, partitions, labels,
and outcomes remain closed.

The next safe sprint is MV5-BD: run the complete Tier D/E equivalence corpus
(318 exact and 90 adaptive-certified results), specify build/distribution and
fallback certification for supported platforms, and make a separate auditable
production-adoption decision. That sprint may evaluate adoption but must not
change defaults unless every full gate passes.
