# MV5-AY complete corrected landscape matrix production

Date: 2026-08-12  
Status: accepted for label-closed stability prefreeze; partitions remain closed

## Outcome

MV5-AY completed the owner-authorized existing-data expansion frozen by MV5-AX.
All eight strata, 56 eligible persistence diagrams, 204 unordered sample pairs,
and 408 separate H0/H1 distance results completed. Every result uses all finite
intervals and all consecutive active landscape levels. H0 and H1 remain the
primary separate matrices. The Euclidean H0/H1 combination is retained only as
a secondary descriptive matrix. No universal landscape-level cap or uniform
grid is used.

All 204 pair artifacts are atomic, versioned, hash-addressed, resumable shards.
An independent validator reconstructed every H0, H1, and combined matrix entry
from those shards, reverified all 56 source file and diagram hashes, recomputed
the canonical public input identities, checked strict numerical certificates,
and verified completion-last and immutable-resume behavior. Sixteen validation
categories pass. A second clean validation produces four byte-identical public
evidence files.

## Numerical contract

- Scientific contract: `full_l2_error_controlled_v1`.
- Intervals: finite deaths only; essential intervals are excluded and counted.
- Levels: all consecutive active levels, with zero padding at missing depth.
- H0/H1: calculated, certified, stored, and interpreted separately.
- Exact route: streamed exact breakpoint integration.
- Adaptive route: partitioned QUADPACK with independent refinement and frozen
  absolute/relative tolerances of `1e-8`.
- Routing: H0 exact in every pair; H1 adaptive only for the two large gene-view
  strata whose interval depths exceed the exact guard of 500.
- Combined distance: exactly `sqrt(H0^2 + H1^2)`, secondary and descriptive.
- Legacy grid/K1 behavior, clustering, labels, outcomes, visualization, Betti
  curves, and workflow defaults remain unchanged and closed.

## Recovery event

The first concurrent adaptive attempt retained 20 SCT-gene and 18
integrated-gene shards before SCT pair 21 encountered QUADPACK's specific
“extremely bad integrand behaviour” condition. The difficult H1 pair contained
1,625 and 1,471 finite intervals. This exposed a numerical partitioning edge
case, not a change in the scientific definition.

Commit `c82a2d6` added a narrow deterministic remediation: the original
integration call is always attempted first; only that named condition triggers
midpoint bisection; child absolute-error budgets sum to the parent budget; and
unrelated errors remain fatal. Thirty-seven focused expectations pass. The real
failed pair certified at an achieved squared-integral error estimate of
`1.728983e-10`, below the frozen `1e-8` threshold, and repeated with identical
normalized result and cache identity. Both production units then resumed from
their retained shards and completed all 90 adaptive-H1 pairs.

## Resource result

All eight processes stayed below the frozen 2 GiB RSS cap. The conservative
peak recorded by the public validator is 945,364,992 bytes (about 0.88 GiB).
The two adaptive units are charged for the entire 3,718-second failed concurrent
attempt plus their successful recovery runs: 7,614 seconds cumulative for
SCT-gene and 9,715 seconds for integrated-gene. Both remain below their frozen
10,830-second budgets. The six exact-only units completed between 41.92 and
286.04 seconds of process wall time.

## Interpretation and next boundary

These matrices establish corrected landscape geometry for the larger existing
sample panels. They do not identify a partition. MV5-AZ is authorized only as a
prospective, label-closed stability/resampling prefreeze. Its first obligation
is to distinguish technical distance reproducibility from matched-axis
partition stability, audit reusable five-seed assets, and freeze the minimum
additional calculation and acceleration equivalence gates before any `k` is
selected.

The measured R adaptive path is correct and bounded, but high-depth gene H1 is
the dominant scaling cost. The repository already contains a corrected Persim
critical-pair implementation whose project-controlled segment integral avoids
the rejected built-in Persim norm. A prospective acceleration sprint should
benchmark that mature exact path against these accepted R matrices first. Rust
is justified only if numerical equivalence passes and the verified mature path
still misses the throughput, memory, packaging, or maintenance target.

## Public evidence

- `mv05ay-production-summary-2026-08-12.csv`
- `mv05ay-resource-summary-2026-08-12.csv`
- `mv05ay-independent-validation-2026-08-12.csv`
- `mv05ay-adaptive-recovery-2026-08-12.csv`
- `mv05ay-continuation-decision-2026-08-12.csv`

Private diagrams, pair shards, matrices, process logs, and source paths remain
under ignored `tmp/` storage and are not publication artifacts.
