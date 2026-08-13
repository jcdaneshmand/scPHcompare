# MV5-BD Rust full equivalence and adoption gate

Date: 2026-08-13
Status: full numerical equivalence accepted; production adoption deferred

## Outcome

MV5-BD evaluated the accepted MV5-BC Rust kernel against the complete frozen
MV5-AY corpus: all 408 H0/H1 dimension references from 56 existing,
raw-hash-verified diagrams. No diagrams, seeds, partitions, labels, outcomes,
or biological analyses were added.

The full numerical gate passes. Tier D passes 318/318 exact R references, with
maximum squared-distance absolute error `7.27595761418343e-11`. Tier E passes
90/90 adaptive-certified H1 references, with maximum error
`3.68900873390365e-12`. The tightest Tier E result remains
`4.46097915185136e-12` inside its accepted R error certificate.

## Invariants and reproducibility

All 408 pair directions were reversed. Every reversed squared distance is
bit-identical, finite-interval counts swap correctly, and level/segment
diagnostics match. All 112 diagram/dimension self-comparisons are exactly zero.
The normalized 408-result table, self table, and environment record repeat
identically across two clean executions.

The two runs completed in 69.01 and 68.91 seconds. Forward calculation for all
408 accepted results consumed 28.068 measured kernel-call seconds in run A;
the remaining time includes 408 reverse calculations, 112 self-calculations,
R orchestration, CSV parsing, hashing, and atomic output. Peak whole-run RSS was
239,700 KiB (`245452800` bytes), well below the 1 GiB gate.

## Adoption decision

The Rust kernel is now scientifically equivalent across the entire accepted
existing-data corpus. Production adoption is nevertheless deferred because
the engineering/distribution gate is incomplete: the pinned Linux build and
runtime are certified, but Windows and macOS builds, platform artifact
provenance, artifact verification/selection, installed-package fallback, and
release distribution are not yet certified.

R therefore remains canonical. The Rust crate remains excluded from R source
builds, no binary is embedded, no public API or default is changed, and every
cache identity remains unchanged. This is an engineering deferral, not a
scientific rejection.

## Next sprint

MV5-BE should be a cross-platform build and distribution prefreeze. It should
inspect current CI and package layout, specify source-versus-prebuilt artifact
policy, supported targets, checksums and provenance, failure/fallback behavior,
and release/CRAN constraints before any workflow or packaging implementation.
It must not enable Rust by default or remove the R path.
