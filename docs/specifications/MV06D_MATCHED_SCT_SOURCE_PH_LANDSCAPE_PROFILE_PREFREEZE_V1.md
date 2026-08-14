# MV6-D matched SCT source/PH/landscape profile prefreeze v1

| Field | Frozen value |
|---|---|
| Date | 2026-08-14 |
| Parent decision | MV6-C `go_bounded_matched_sct_profile` |
| Global panel | 500 genes; scientific panel SHA-256 `7be22cdb9056fed427c78d58be2b19e258c7c6807e6f7ac3900bd1bfa1d19eb8` |
| Candidate axis | 90 samples, 15 studies, five seeds, LOSO folds |
| Maximum fold sentinels | Five, one per seed |
| Maximum profiled samples | Ten, one held-out and one training sample per fold sentinel |
| Maximum PH jobs | 20: 10 cell and 10 gene |
| Maximum landscape pair calls | 10: five cell and five gene; H0/H1 retained separately |
| Labels/outcomes | Closed / zero |
| Required stop | Before full fold/view production, fusion, clustering, or outcomes |

## Purpose

MV6-D measures the real corrected matched-SCT pipeline before any 6,750-view
cell/gene rerun is authorized. It reconstructs a bounded number of fold-local
sources from the accepted MV5-D0 v2 SCT caches, executes typed H0/H1 PH for
matched cell and gene views, evaluates the dissertation-aligned landscape
calculation on bounded pairs, and projects the full workload.

This is a resource and correctness profile, not a fusion evaluation. It cannot
select a weight, compare a biological label, cluster samples, or support a
biological claim.

## Frozen inputs

| Input | File SHA-256 |
|---|---|
| MV5-D0 v2 resource ledger | `73f757a91c202a8e38dfa746fc8816f8e272caf912c33b4b959ef55102c68308` |
| Label-closed candidate manifest | `842c047ba821f8eca317da52504910733509fb4fddd11d6f54f7e79d9f29d0b7` |
| LOSO fold-seed plan | `50379f98cd4927c5c8cb19dbd9ca8ecc7b7b3a9af2e04eb9c8358ecb0b722c6d` |
| Public MV6-C panel CSV | `b3a5aff1a0bc01e871751fb9db0b3babfaf18835e68c5699346d8476d903d0ab` |

Hashes are compared case-insensitively. The scientific panel hash is computed
from the ordered feature/gene identity and is distinct from the public CSV file
hash.

Every selected cache must exist, match its accepted file hash and cache key,
contain exactly the accepted selected-cell identity, and remain label closed.
No source reconstruction may proceed after an input mismatch.

## Deterministic sentinel selection

The 15 studies are sorted by held-out sample count and then study ID using
radix order. Ranks 1, 4, 8, 12, and 15 are paired in order with seeds 20260805
through 20260809. This gives five distinct fold-seed sentinels spanning the
observed fold-size range without using tissue, assay approach, endpoints,
outcomes, PH, landscapes, or fusion results.

Within each selected fold-seed:

1. the held-out sentinel is the held-out sample with the largest accepted
   private-cache byte count, breaking ties by sample ID;
2. the training sentinel is selected by the same rule among training samples;
3. the rank-8/seed-20260807 sentinel executes first as stage 1; and
4. the other four sentinels may execute only if stage 1 passes every correctness
   guard and both gene-PH jobs remain within the continuation limits below.

The resulting public manifest may contain study, sample, seed, role, cache
hashes/keys, and technical counts. It may not contain expression values, cell
IDs, biological labels, or outcomes.

## Matched source contract

For each selected fold-seed, all 90 caches for that seed are read and validated.
The exact ordered MV6-C panel is extracted with no zero filling, feature
substitution, panel shrinking, or sample removal.

The global panel is an explicitly transductive technical harmonization learned
in MV6-C. Conditional on that fixed panel, all downstream transformations here
remain training-only within the LOSO fold:

- gene centers and scales are fitted from pooled training-sample cells only;
- every scale must be finite and greater than `sqrt(.Machine$double.eps)`;
- the same positive training-fitted affine transform is applied to training
  and held-out sources;
- a deterministic 30-component `stats::prcomp` model is fitted only on the
  training sources; and
- all 90 standardized matrices must satisfy the scientific 500-gene by
  384-cell typed-source contract, although only the two sentinel cell and gene
  views are retained in the bounded private bundle.

The cell view is 384 cells in the shared 30-PC coordinates. The gene view is
500 genes in the exact same matched 384 cells, centered and unit-normalized
within gene, with explicit Pearson-correlation chord distance. These are
separate observation spaces and may not be pooled into one diagram.

## PH and diagram contract

All jobs call `run_topology_view_ph()` on a validated typed view:

- complete Vietoris–Rips filtration;
- cell metric: Euclidean distance in shared training-fitted PCs;
- gene metric: explicit Pearson-correlation chord distance;
- `max_dim = 1`, `threshold = -1`, coefficient field 2; and
- exactly one essential H0 interval, retained in the diagram but excluded from
  landscapes.

For a view with `n` observations, H0 must contain `n - 1` finite positive merge
intervals plus one essential interval. Finite H0 deaths must match a separately
implemented Prim minimum-spanning-tree oracle on the view's actual metric.
Every H1 interval must be finite and have positive persistence. Stage-1 cell
and gene results are independently rerun and must be byte-identical after
excluding timing/RSS sidecars.

## Landscape contract

Each fold compares its held-out sentinel with its training sentinel separately
for the cell view and gene view. The canonical public API is used in
`full_l2_error_controlled_v1` mode with automatic exact/error-controlled
routing. The following dissertation-aligned definition is immutable:

- H0 and H1 are calculated and reported separately;
- only finite positive-persistence intervals enter landscapes;
- every active consecutive landscape level is retained;
- missing levels are zero-padded between diagrams;
- integration is exact or error-controlled on dimension-specific support;
- there is no universal level cap or fixed uniform grid; and
- any combined distance is secondary to its recorded H0 and H1 components.

The stage-1 landscape calculations are repeated and must reproduce scientific
distances and cache identities exactly. Runtime fields are sidecars and do not
enter immutable scientific identity.

## Resource guards and staged stop rules

Only one heavy worker may run at a time; BLAS/OpenMP thread counts are one.

| Unit | Hard cap |
|---|---:|
| One fold source/PCA bundle | 1,800 seconds and 8 GiB process-tree RSS |
| One cell PH job | 600 seconds and 4 GiB process-tree RSS |
| One gene PH job | 1,800 seconds and 8 GiB process-tree RSS |
| One cell or gene landscape pair call | 1,800 seconds and 8 GiB process-tree RSS |
| Entire bounded stage | 14,400 worker-seconds and 10 GiB private storage |

Stage 2 is prohibited if stage 1 has any correctness, identity, serialization,
or resource failure, or if either stage-1 gene-PH job exceeds 900 seconds or
6 GiB. A stopped stage remains evidence; caps may not be relaxed after seeing
the result within MV6-D.

## Full-workload projection and decision

Observed median, type-8 p90, and observed-maximum scenarios are projected for:

- 75 fold source/PCA builds;
- 6,750 cell PH jobs;
- 6,750 gene PH jobs;
- 35,350 cell-diagram pair calls; and
- 35,350 gene-diagram pair calls.

Private source/diagram storage is projected separately. The public record may
contain aggregate counts, times, RSS, sizes, interval counts, distances, and
hashes, but not private matrices, diagrams, expression values, or cell IDs.

The prospective disposition is:

- `go_prefreeze_full_matched_sct`: all correctness/determinism guards pass,
  stage 2 completes, observed-maximum total is at most 72 worker-hours, and
  projected private storage is at most 10 GiB;
- `revise_for_targeted_acceleration`: correctness passes but the maximum
  projection exceeds 72 worker-hours or one component is an obvious dominant
  bottleneck; any Rust work must reproduce the canonical R scientific result;
- `stop_matched_sct_scaleup`: a scientific/identity guard fails or stage 1
  breaches a hard cap; or
- `insufficient_profile`: stage 2 cannot complete for a reason that does not
  invalidate the estimand but leaves the projection unsupported.

No disposition passes G-MV6 or authorizes fusion. A `go` authorizes only a new
immutable full-run prefreeze. An acceleration disposition authorizes only a
correctness-first accelerator sprint. Outcomes, clustering, new data,
integrated gene topology, default changes, manuscript claims, and release
actions remain closed.
