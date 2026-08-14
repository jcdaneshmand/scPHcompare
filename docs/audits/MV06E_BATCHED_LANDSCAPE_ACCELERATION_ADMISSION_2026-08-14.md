# MV6-E batched landscape acceleration admission

| Field | Result |
|---|---|
| Date | 2026-08-14 |
| Prefreeze commits | `c201c8b`; workload correction `c458a59` |
| Frozen inputs | 20 accepted MV6-D cell/gene H0/H1 diagrams |
| Reference / throughput rows | 20 canonical R / 180 cross-engine primary |
| Invariant rows per engine | 20 reverse; 40 exact self-zero |
| Clean scientific repeats | 2/2 grouped Persim; 2/2 Rust, byte-identical |
| Independent validation | 8/8 categories pass |
| Resume validation | 4/4 scientific/metrics artifacts unchanged |
| Labels/outcomes | Closed / zero |
| Decision | `admit_both_with_rust_preferred_private` |

## Outcome

MV6-E removes the landscape-computation barrier without changing the scientific
landscape definition. Both the portable grouped Persim engine and the accepted
Rust kernel reproduce the frozen canonical R references, agree across the full
throughput corpus, satisfy reverse and self invariants, repeat deterministically,
and remain below the prefrozen resource caps.

The Rust kernel is admitted only as an explicit, hash-verified private WSL
production accelerator. Grouped Persim remains the canonical portable exact
fallback. This decision does not embed or distribute a binary, change an R
package default, replace the R oracle, or authorize biological evaluation.

## Frozen corpus and checks

The private staging file contains 18,733 finite positive-persistence intervals
from ten cell-view and ten gene-view dimension records. Public evidence exposes
only technical diagram identities, interval counts, frozen pair requests,
references, scientific results, resource metrics, projections, and the decision;
the private intervals are not tracked.

For each view, the throughput corpus contains all 45 unordered pairs among ten
diagrams in H0 and all 45 in H1. This gives 180 primary rows. Each engine also
computed the reverse of all 20 canonical reference requests and 40
diagram/dimension self-comparisons, for 240 normalized scientific rows per run.

The validator established:

- 20/20 R exact or adaptive-certified references for both engines;
- 180/180 cross-engine throughput rows within the frozen exact tolerance;
- 20/20 reverse rows per engine with swapped interval counts;
- 40/40 self rows per engine exactly zero;
- exact H0/H1 separation, finite nonnegative distances, and all-active-level
  provenance; and
- byte-identical A/B scientific outputs for each engine.

The maximum absolute squared-distance discrepancies were
`4.093e-12` for grouped Persim versus R, `1.933e-12` for Rust versus R, and
`1.228e-11` between the two exact candidates across all 180 throughput rows.
All are within the prospectively frozen absolute-plus-relative certificates.

## Measured resources

| Engine/run | Total engine time | Peak process RSS | Scientific SHA-256 |
|---|---:|---:|---|
| Grouped Persim A | 260.729 s | 981,184,512 B | `b1ff51f566b4fdb667d6ff70f9733c0c9c4bf82da001475abae508c80a0c226f` |
| Grouped Persim B | 256.743 s | 982,118,400 B | `b1ff51f566b4fdb667d6ff70f9733c0c9c4bf82da001475abae508c80a0c226f` |
| Rust A | 5.387 s | 199,958,528 B | `1157be5b7c1248e6d96fe75d3a18f66ed6116719c964f5852029ee37168ce8f5` |
| Rust B | 5.903 s | 199,917,568 B | `1157be5b7c1248e6d96fe75d3a18f66ed6116719c964f5852029ee37168ce8f5` |

The Rust engine times are measured after the accepted R environment is active;
observed end-to-end launches were approximately 13.5--14.1 seconds. A production
worker may amortize that startup, and even charging roughly nine seconds to each
of 75 groups adds under 0.2 worker-hours. The admission decision therefore does
not depend on omitting startup overhead.

Both candidates remain below the 3,600-second and 8-GiB per-process guards.
Grouped Persim constructs each of the 40 diagram/dimension landscapes once and
reuses it. Rust integrates interval representations directly; its zero landscape-
build count is an engine-design distinction, not a zero-cost or altered-landscape
claim.

## Full-workload reconstruction

The corrected full workload contains 27,000 diagram/dimension preparations
(13,500 diagrams times H0/H1) and 141,400 component pair rows (35,350 biological
pairs times two views times H0/H1), assigned to 75 groups.

| Engine | Maximum projected landscape worker-hours | Maximum total with MV6-D source and PH | Landscape cap |
|---|---:|---:|---:|
| Grouped Persim | 42.814 h | 51.290 h | Pass (<60 h) |
| Rust | 0.873 h | 9.349 h | Pass (<60 h) |

The reconstructed maximum projection gives Rust a 49.015-fold landscape
advantage over grouped Persim. Both routes pass the 60-hour landscape and
72-hour combined planning caps. The prior 837.828-hour estimate measured one
public-R process per pair and is superseded for production planning, not for
the R oracle's correctness role.

## Resume and failure behavior

Both runners publish scientific and metrics files atomically and refuse to
overwrite stale output. Resume now requires both companion files, the exact
row count, and the expected contract IDs. A deliberately incomplete artifact
pair fails closed.

After the clean runs, each A runner was invoked again. Both printed their
validated-reuse message. The four accepted A artifacts retained identical
SHA-256 hashes, byte sizes, and modification times; the public
`mv06e-resume-validation.csv` records those observations.

## Landscape contract preserved

The dissertation-aligned primary definition remains:

- H0 and H1 are calculated and retained separately;
- only finite positive-persistence intervals enter the landscape;
- every active consecutive landscape level is retained;
- absent depth is zero-padded between diagrams;
- squared L2 is integrated exactly or under an explicit error certificate;
- no universal landscape-level cap or uniform grid is introduced; and
- any combined H0/H1 or cell/gene quantity remains derived and fully
  decomposable into recorded components.

Rust and grouped Persim accelerate this same quantity. Neither substitutes a
persistence image, fixed grid, top-k approximation, bottleneck distance, or
Wasserstein distance.

## Gate disposition

MV6-E completes as `admit_both_with_rust_preferred_private`. The next action is
a separate immutable MV6-F full matched-production prefreeze. That contract
should stream source/view/PH artifacts, invoke the accepted Rust library only
after verifying SHA-256
`51d3fca4808c6620c46c7b1100ddfad4753af02d37f87ec3eb278008846b160d`,
retain grouped Persim as the portable fallback, and sample canonical R oracle
checks before any blocked fusion evaluation.

MV6-E alone does not authorize full production, outcome access, clustering,
fusion fitting, G-MV6, new data, a public accelerator default, binary release,
manuscript claims, PDFs, confidential reviewer material, or `example_run.r`.
