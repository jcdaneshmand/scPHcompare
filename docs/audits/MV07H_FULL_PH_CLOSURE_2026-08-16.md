# MV7-H full-PH closure

Date: 2026-08-16

Status: accepted; one prospectively selected landscape stress group authorized

## Scope and scientific boundary

The complete label-closed 124-sample, five-seed, dual-view PH stage is closed.
It contains 620 cell-topology and 620 gene-topology records. Every record
retains H0 and H1 separately. No landscape, clustering, label, outcome,
dimension-combination, or manuscript-claim job was run during PH production.

The persistence-landscape estimand is unchanged: finite positive intervals;
the essential H0 class excluded; all active consecutive landscape levels; H0
and H1 calculated separately; exact squared-L2 integration; no uniform grid;
no universal level cap; and streamed pair groups rather than saved dense
landscape-function matrices.

## Production result

The production decision is
`full_PH_complete_await_independent_validation`. Five source bundles yielded
1,240 valid PH records: 1,238 selected from Ripserr and two selected from the
approved exact TDA/GUDHI resource fallback. Both fallbacks occurred for sample
`SRA628554_SRS2664364`, gene topology, in seeds 20260805 and 20260806 after
their primary Ripserr attempts exceeded the frozen RSS cap. The other three
seeds completed that sample under Ripserr, demonstrating why the engine gate
was applied independently per seed.

The accepted aggregate elapsed time is 18,980.127 seconds against a
172,800-second cap. Private state is 1,082,955,320 bytes against a
4,294,967,296-byte cap. All 60 accepted MV7-G sentinel views reproduce exactly.

## Audited recovery history

The first seed-20260806 fallback repeat failed in 0.544 seconds before
scientific computation because the dynamic repeat queue had not regenerated
that seed's support source. Its stderr records the exact missing-source path,
and its failed ledger row remains preserved. Commit `9e00c38` added a narrowly
gated, named recovery attempt and dynamic support-source derivation. Commit
`49578a7` added strict immutable-resume ownership for shared output paths: an
RSS-failed primary can accept an existing output only when exactly one
completed fallback receipt matches both SHA-256 and byte size.

The independently regenerated seed-20260806 support source is 146,569,565
bytes with SHA-256
`f820a59ff62a4fc5c5109492db234aea7eb39b9f00060b03f19a62c0621ae757`.
The named recovery completed in 359.534 seconds at 6,654,717,952 bytes peak
RSS. Its 92,527-byte artifact has SHA-256
`24ab5e77566feb4f0167a6736db1b67b491c5bd077dd337a112238a9fcf820aa`,
exactly matching the accepted production artifact. No failed evidence was
deleted or rewritten.

## Independent validation

- All 12 independent validation categories pass.
- All 1,240 record identities, orientations, hashes, sizes, H0 MST oracles, and
  finite-positive H1 requirements pass.
- All 20 fresh Ripserr/GUDHI subset comparisons pass with maximum absolute
  H0/H1 interval error 0.
- All 15 prospectively required repeat artifacts are byte- and SHA-identical.
- Whole-stage immutable resume passes 7/7 across 3,802 files: 3,793 private
  artifacts/ledgers/logs and nine public production/validation files retained
  identical file axes, SHA-256 values, byte sizes, and mtimes.
- The canonical repository suite passes with only the two established optional
  MV5-BC skips.

The accepted closure implementation HEAD is
`c8d658e3153d3be16cc7d39168e91c7deaf852d4`. Evidence is in
`docs/audits/mv07h-full-ph-evidence/`.

## Decision

The independent decision is
`authorize_one_MV7H_landscape_stress_group`. It authorizes only the frozen
seed-20260807, gene-topology, H1 group containing all 7,626 unordered pairs.
The other 19 landscape groups, H0/H1 combinations, clustering, labels,
outcomes, and result-dependent claims remain closed until that group passes
its Rust/reference-oracle, deterministic-repeat, resource, completeness, and
immutable-resume gates.
