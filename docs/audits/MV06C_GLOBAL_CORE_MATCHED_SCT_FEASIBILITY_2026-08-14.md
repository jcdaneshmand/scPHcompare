# MV6-C global-core matched SCT feasibility

| Field | Value |
|---|---|
| Date | 2026-08-14 |
| Prefreeze commit | `b52c857` |
| Implementation commit | `2414024` |
| Contract | `mv06c_global_core_matched_sct_feasibility_v1` |
| Inputs | 450 accepted MV5-D0 v2 SCT caches |
| Biological samples / seeds | 90 / 5 |
| Cells per cache | 384 |
| Labels/outcomes | Closed / zero |
| Decision | `go_bounded_matched_sct_profile` |

## Outcome

Option A is technically viable. A single fixed 500-gene SCT panel can be used
for every accepted sample-seed instance without zero filling, sample removal,
or sample-specific gene sets.

All 450 accepted private cache hashes, cache keys, payload identities, selected
cell axes, and label-closed states verified. The two-pass production inventory
found:

| Quantity | Result |
|---|---:|
| Union of cache feature IDs | 33,914 |
| Exact feature intersection across all 450 caches | 2,730 |
| Retained-category common features | 2,596 |
| Excluded mitochondrial/ribosomal/hemoglobin features | 134 |
| Common features nonfinite in any cache | 0 |
| Finite but nonpositive-variance in any cache | 62 |
| Eligible features before canonical deduplication | 2,536 |
| Duplicate canonical features removed | 0 |
| Eligible unique canonical genes | 2,536 |
| Frozen panel size | 500 |
| Eligibility margin above 500 | 2,036 |

The ordered panel SHA-256 is
`7be22cdb9056fed427c78d58be2b19e258c7c6807e6f7ac3900bd1bfa1d19eb8`.
Every selected gene is present, finite, and nonconstant across all 450 caches.
The compact ordered panel is retained as
`docs/audits/mv06c-global-core-evidence/mv06c-panel.csv`; no per-sample
expression values, variances, ranks, or cell identities are published.

## Seed stability

The five seed-specific top-500 panels are highly stable against the global
panel:

| Seed | Eligible genes | Shared top-500 | Jaccard | Rank correlation |
|---:|---:|---:|---:|---:|
| 20260805 | 2,572 | 491 | 0.96464 | 0.99908 |
| 20260806 | 2,572 | 492 | 0.96850 | 0.99905 |
| 20260807 | 2,581 | 493 | 0.97239 | 0.99909 |
| 20260808 | 2,579 | 492 | 0.96850 | 0.99904 |
| 20260809 | 2,576 | 493 | 0.97239 | 0.99910 |

This supports a fixed all-seed panel rather than seed-specific panel tuning. It
does not show biological usefulness or fusion improvement.

## Scientific interpretation

The global core resolves the MV6-B structural missing-feature problem for SCT.
It is deliberately transductive technical harmonization: feature presence and
within-sample variance from all existing samples/seeds define the fixed panel.
No biological label or result selects a gene, but the panel is not a
future-study-inductive feature-selection rule.

Accordingly, a later blocked analysis can estimate downstream study-held-out
behavior conditional on this fixed existing-data panel. It cannot claim that
the panel itself was learned without seeing held-out technical availability.
For external validation, this exact panel must be frozen before the new data
are opened, and unavailable genes require a prespecified disposition.

The selected panel also merits later biological/technical interpretation.
The existing technical-feature rules remove mitochondrial, ribosomal-protein,
and hemoglobin genes, but globally variable housekeeping, activation, stress,
or cell-composition genes can still dominate. That is a robustness/confounding
question for MV-07, not a reason to alter the prospectively selected panel now.

## Validation, determinism, and resources

The independent validator uses a separate one-pass aggregation and does not
call the production panel selector. It independently rehashes all 450 caches
and reconstructs the union/intersection axes, all eligibility counts, exact
ordered panel and aggregate values, panel hash, all five seed panels and
stability statistics, workload, decision, resource gate, public schema, and
zero-execution boundary. All 12 categories pass.

The production inventory used 840.012 seconds, 913,883,136 bytes peak process
RSS, and 2,991,811,724 unique source bytes. It performed two deserialization and
two hash passes, each totaling 5,983,623,448 bytes. The run passes the frozen
1,800-second and 8-GiB caps. Nondeterministic timing/RSS remain in a separate
resource sidecar and are excluded from scientific identities.

The complete source-loaded R test suite passes, including all 17 focused MV6-C
expectations. Its two skips are the established optional Rust-library and
public-audit-in-build exclusions, not MV6-C failures. The installed-package
PH-child-process boundary found during MV6-B remains a separate
publication-readiness issue; MV6-C changes neither that launcher nor its toy
baseline.

A clean repeat reproduced all eight immutable scientific/decision files with
identical filename axes, byte counts, and SHA-256 hashes. The repeat used
857.449 seconds and 923,742,208 bytes peak process RSS, independently passing
the same 1,800-second and 8-GiB caps. MV6-C therefore satisfies its prospective
determinism and resource acceptance rules.

## What can and cannot be reused

Reusable without scientific change:

- the 450 accepted MV5-D0 SCT normalization caches;
- the exact 384 selected cells per sample and seed;
- the 90-sample, 15-study, five-seed blocked axis;
- the corrected typed view, PH, landscape, provenance, and fusion
  implementations; and
- the dissertation-aligned landscape definition.

Must be rebuilt or re-evaluated because the panel changed:

- fold-local training standardization and cell PCA;
- all global-core cell and gene typed views;
- cell and gene H0/H1 diagrams;
- separate all-active-level landscape distances;
- cell-energy and panel-matched expression baselines;
- training-fitted four-component scales; and
- any later fusion, clustering, retrieval, robustness, or outcome result.

The accepted old-panel cell analysis remains valid for its original estimand,
but it is not the matched cell component of the new global-core fusion
estimand.

## Future workload and next gate

The complete candidate would contain:

| Unit | Count |
|---|---:|
| Global-core cell views | 6,750 |
| Global-core gene views | 6,750 |
| H0/H1 diagram components | 27,000 |
| Directed held-out-query/training-reference pairs | 35,350 |
| Four component-specific landscape distances | 141,400 |
| Fold-seed component scales | 300 |
| Five-weight fusion pair rows | 176,750 |

These counts are not authorization. The next sprint is MV6-D: prospectively
select bounded label-closed fold-seed/source/PH sentinels, reconstruct matched
global-core SCT cell and gene sources, profile diagram and landscape cost, and
stop before full production or outcomes. Gene PH construction—not the already
validated landscape formula—must be treated as the primary unknown bottleneck.
The optional Rust landscape prototype may be profiled only after canonical R
correctness and cannot change the landscape definition.

## Gate disposition

MV6-C passes and rescues a coherent route toward fusion. It does not pass
G-MV6. Full source construction, PH, landscapes, fusion, clustering, outcomes,
new data, integrated gene topology, defaults, claims, release, DOI, binaries,
PDFs, reviewer material, and `example_run.r` remain closed.

The accepted landscape definition remains separate H0/H1, finite
positive-persistence intervals, every active consecutive level, and exact or
error-controlled L2 integration without a universal grid or level cap.
