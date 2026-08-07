# MV5-D0 stage-1 completion audit

| Field | Value |
|---|---|
| Date | 2026-08-07 |
| Contract | `mv05d0_production_cache_validation_v1` |
| Scope | 90 existing-data samples × 5 frozen seeds |
| Outcome-label state | Closed |
| SCT caches | 450/450 complete and independently validated |
| Downstream folds, PH, landscapes, distances, clustering, integration, outcomes | 0 |
| Decision | MV5-D0 stage 1 complete; stop and review before any fold execution |

## Outcome

MV5-D0 stage 1 is complete. The previously blocked monolithic-source route was
replaced with the existing per-sample Seurat files already stored under the
historical project. No new dataset was added and no outcome label was opened.

Exactly one source file was found for each of the 90 frozen candidates. Each
job extracted only the RNA count matrix, removed the deterministic sample-name
prefix from cell IDs, converted counts to a canonical sparse representation,
and published an atomic private shard. All 53 samples recovered during the
failed monolithic attempt have exact canonical count identity with their
individual-source versions. The remaining 37 samples were recovered directly
from their existing individual files.

## Raw source migration

| Measurement | Result |
|---|---:|
| Raw shards | 90/90 |
| Exact recovered-monolithic comparisons | 53/53 |
| Newly recovered individual sources | 37 |
| Failures | 0 |
| Worker-hours | 0.7052 h |
| Maximum process-tree RSS | 3,096,002,560 B |
| Total compressed raw-shard storage | 1,048,827,846 B |

The source migration used at most two workers and stayed below the same
1,800-second/8-GiB/40-GB envelope. Source paths and private objects remain
outside Git; the public audit contains sample IDs, basenames, dimensions,
hashes, resources, and closed-label state only.

## Frozen selections

All 450 sample-seed selections were regenerated from the v2 raw shards using
the existing deterministic 384-cell rule and seeds `20260805`–`20260809`.
The selection-set SHA-256 is
`b1210a68f063d1d64aae8109cc5f4dbe1e34067d66d28695e8af9a905327a6f5`.
No tissue or approach column is present.

## SCT production cache

| Measurement | Result |
|---|---:|
| Entries | 450/450 |
| Samples × seeds | 90 × 5 |
| Build failures | 0 |
| Normalization worker-hours | 2.5623 h |
| SCTransform operation-hours | 1.7642 h |
| Median entry elapsed | 20.113 s |
| Maximum entry elapsed | 37.516 s |
| Maximum process-tree RSS | 1,942,548,480 B |
| Total cache storage | 2,991,811,724 B |
| Largest cache file | 22,570,179 B |
| Runtime SHA-256 | `480da9fc132f94a5c901c9fa3fac52d77f602f61615e226981d1490a8d4f3948` |
| Cache-set SHA-256 | `a7e1809930b3b9b6314de9325c554808fe16c57dbdaded5ef0aa7a6470fedfb0` |

An independent post-build pass reopened every cache and checked the v2 record,
runtime identity, sample, seed, selected-cell hash, payload hash, file size,
and file SHA-256 against the public manifests. All 450 are valid resumable
entries. The cache occupies 7.48% of the 40-GB cap; the largest observed RSS is
22.61% of the 8-GiB cap; and the longest entry is 2.08% of the 30-minute cap.

## Mandatory feasibility reprojection

The old normalization estimate was 10.717 worker-hours. Replacing it with the
measured 2.562 worker-hours gives:

| Scenario | Updated lower bound | 21.6-h planning cap | Disposition |
|---|---:|---:|---|
| SCT cell-primary | 10.525 h | Pass | Feasible for a future label-closed fold stage after review |
| SCT cell + gene | 14.097 h | Pass numerically | Deferred: gene eligibility/scope gate remains open |
| All planned views | 17.668 h before integrated mapping | Incomplete | Prohibited: integrated mapping is unmeasured and not authorized |
| Naive full MV5-D | 242.516 h | Fail | Prohibited |

This is a resource decision, not a biological result. It does not authorize a
gene-view expansion, integration, complete-matrix clustering, fusion, or any
outcome calculation.

## Landscape boundary

The dissertation-aligned landscape definition is unchanged: H0 and H1 remain
separate; every active consecutive landscape level is retained; and integration
is exact or error-controlled. No persistence diagram or landscape was computed
in this stage.

## Privacy and scope

- Outcome labels remained closed.
- No tissue/approach endpoint join or method comparison was performed.
- No fold, PH, landscape, distance, clustering, integration, or outcome job ran.
- Raw shards, SCT caches, detailed logs, PDFs, reviewer correspondence, and
  `example_run.r` remain outside Git.
- Nothing was pushed.

## Decision table

| Question | Disposition |
|---|---|
| Scientific contract coherent? | Pass for stage-1 SCT inputs; landscape definition unchanged |
| Correctness demonstrated? | Pass: 90 raw shards, 450 selections, 450 independently validated v2 caches |
| Computation feasible? | Yes for future SCT cell-primary fold execution under the measured projection |
| Biological interpretation permitted? | Prohibited |
| Next action | Owner review, then specify MV5-D1 label-closed SCT cell-primary fold execution |
