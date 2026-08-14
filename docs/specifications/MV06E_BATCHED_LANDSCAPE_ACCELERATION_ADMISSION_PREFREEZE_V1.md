# MV6-E batched landscape acceleration admission prefreeze v1

| Field | Frozen value |
|---|---|
| Date | 2026-08-14 |
| Parent decision | MV6-D `revise_for_targeted_acceleration` |
| Scientific inputs | 20 accepted MV6-D cell/gene H0/H1 diagrams |
| Accepted reference pairs | 10 held-out/training pairs, 20 dimension results |
| Throughput corpus | All 45 unordered pairs within each 10-diagram view, H0/H1 separate |
| Candidate engines | Existing grouped Persim exact engine; existing certified Rust exact kernel |
| Canonical oracle | R `full_l2_error_controlled_v1` results |
| Labels/outcomes | Closed / zero |
| Required stop | Before 75-fold source/PH production, full distances, fusion, clustering, or outcomes |

## Purpose

MV6-E determines whether the accepted all-active-level landscape definition can
be executed in production style on the new matched cell and gene diagrams. It
does not reconsider the landscape mathematics, PH orientation, global-core
panel, fusion estimand, or biological endpoints.

MV6-D showed that fold source/PCA and cell/gene PH are feasible, while starting
one R reference process per diagram pair projects to 837.828 worker-hours. The
accepted MV5-D4 grouped exact calculation previously completed all 35,350 cell
pairs in 1.165 worker-hours. MV6-E therefore tests reusable grouped/compiled
execution before any new algorithm is invented.

## Frozen inputs and implementations

| Artifact | SHA-256 |
|---|---|
| MV6-D PH metrics | `69a665b8fa484da603acf92ba5877a96938083f00e5fddee3017410d5c51bd16` |
| MV6-D landscape metrics | `44d3ba89d7441ab03328907cbd6ab75842129560402726d7aaeb02d92a7004bd` |
| MV6-D decision | `b534e0d6a582fffd7a5d2985a593c23939675d4c9179d33de8787d3115178089` |
| Rust kernel source | `35979730d533005aee5bf9fa79083f17c531123329b89e78e87c1968d1428381` |
| Rust lockfile | `b5656db9c642426754e11124cedb85c32e79be871a5c84da44c7d1dd74baf7f9` |
| R Rust shim | `4771ec396759dacaf36e860819b49978b37d8cc8d4224e54cc72b0077f8a5eb1` |
| MV5-D4 grouped Persim runner | `870cef3f3622faaa2774e54bb9a41931ee42a6010e77067ad62b9866eecf98f3` |
| MV5-W grouped exact runner | `590642095d18143db0b0ce9f6725bed31b3e4e50cb0bd1df5d70430803f77a59` |

The private diagram files must match every MV6-D PH metric hash, cache key,
view, role, sample, seed, point count, and interval count. Only finite positive-
persistence H0/H1 intervals are staged. Essential H0 is audited and excluded.
No expression matrix, cell identity, or topology-view distance object enters
the acceleration corpus.

## Immutable landscape definition

Every engine must compute the same dimension-specific squared L2 distance:

- H0 and H1 are separate calls and separate output rows;
- every finite positive-persistence interval is included;
- every consecutive active landscape level is included;
- the shallower diagram is zero-padded by level;
- integration is exact over piecewise-linear critical segments, or compared
  against the accepted R error-controlled certificate;
- there is no level cap, interval cap, or uniform grid; and
- any square-rooted or H0/H1-combined value is derived only after preserving
  the squared dimension components.

An engine that is fast only after truncation, gridding, or changing H1 policy
fails automatically.

## Frozen corpora

### Reference corpus

The ten MV6-D held-out/training comparisons produce 20 dimension rows: cell
H0/H1 and gene H0/H1 in each of five fold-seeds. These rows retain the accepted
R distance, exact/adaptive route, achieved error estimate, finite-interval
counts, diagram identities, and view/dimension identity.

For R-exact rows, a candidate squared distance must satisfy
`abs(candidate - reference) <= 1e-10 + 1e-10 * abs(reference)`. For an
R-adaptive-certified row, it must fall within the recorded achieved absolute
error estimate plus `100 * .Machine$double.eps * max(1, abs(reference))`.
Rust and grouped Persim must also agree with each other under the exact-row
tolerance on all 20 rows.

### Throughput corpus

Within each view, the ten accepted diagrams are canonically ordered by diagram
hash. All 45 unordered pairs are evaluated for H0 and H1, yielding 90 dimension
rows per view and 180 total. The pair axis is label closed and is used only for
correctness/resource profiling.

Each grouped candidate must build or canonicalize a diagram/dimension at most
once per execution and reuse that representation across its pair requests when
the engine supports reuse. Rust may instead reconstruct per call only if its
measured throughput remains competitive. Timing must separate staging,
landscape construction, pair integration, serialization, and total process
cost wherever the implementation exposes those phases.

## Candidate rules

### Grouped Persim exact

Reuse the accepted critical-pair construction and exact linear-segment
integration; do not use Persim's rejected built-in norm. Cache each diagram's
H0/H1 critical-pair landscape once within the group. Record Persim/Python
versions and implementation hashes.

### Rust exact

Use only the existing dependency-free locked crate and explicit R shim. Build
from the frozen source/lock with the already installed isolated toolchain when
available; no network, system-package, shell-profile, or package-default change
is authorized. The library hash, compiler/toolchain, status code, engine
version, interval counts, active levels, and event segments are mandatory.

MV5-BD already established full-corpus numerical equivalence, while public
production adoption was deferred for cross-platform distribution. MV6-E may
admit Rust only as an explicit private WSL production engine with verified
library identity and canonical fallback. It cannot embed a binary, enable a
default, or claim public cross-platform distribution.

## Correctness, determinism, and execution guards

Both engines must pass:

1. 20/20 accepted reference rows;
2. 180/180 finite nonnegative throughput rows with identical interval counts;
3. 20/20 reversed reference pairs, with bit-identical or tolerance-equivalent
   squared distance and swapped interval counts;
4. 40/40 diagram/dimension self-comparisons exactly zero;
5. exact H0/H1 separation and all-active-level provenance;
6. two clean normalized executions with byte-identical scientific output; and
7. atomic output, stale-identity refusal, and zero-rebuild resume.

| Guard | Frozen value |
|---|---:|
| Heavy workers | 1 |
| One engine execution elapsed cap | 3,600 seconds |
| One engine process-tree RSS cap | 8 GiB |
| Combined bounded private storage cap | 2 GiB |
| Labels/outcomes/fusion/clustering jobs | 0 |

## Production projection and decision

The benchmark projects:

- 27,000 diagram/dimension constructions if no reuse exists, or 27,000 staged
  diagram/dimensions amortized within 75 groups when caching is supported; and
- 141,400 component-distance rows (35,350 pairs × cell/gene × H0/H1).

Projection must include process/group startup and serialization rather than
kernel time alone. A production candidate passes the landscape resource gate
only if its observed-maximum projection is at most 60 worker-hours and its
peak RSS is at most 8 GiB. Together with MV6-D's 8.48-hour maximum source/PCA
plus PH projection, this retains margin under the 72-hour overall gate.

The prospective disposition is:

- `admit_grouped_persim`: grouped Persim passes correctness/determinism/resume
  and the 60-hour landscape cap; Rust adds no required material advantage;
- `admit_explicit_rust_wsl`: Rust passes every gate and has a material
  throughput or memory advantage needed for the cap; use remains explicit,
  hash-verified, private-production only with canonical fallback;
- `admit_both_with_rust_preferred_private`: both pass and Rust provides at
  least 3× total-throughput speedup while grouped Persim remains the portable
  canonical fallback;
- `revise_batch_design`: correctness passes but neither production projection
  meets the cap; or
- `stop_acceleration_candidate`: identity or numerical correctness fails.

Passing MV6-E authorizes only a separate immutable source/PH/landscape
production prefreeze. It does not launch full production or pass G-MV6.
Fusion, clustering, outcomes, new data, integrated gene topology, public
defaults, binary distribution, manuscript claims, release actions, PDFs,
reviewer material, and `example_run.r` remain closed.
