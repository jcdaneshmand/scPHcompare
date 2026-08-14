# MV5-D1 label-closed SCT cell-fold specification v1

| Field | Frozen value |
|---|---|
| Date | 2026-08-07 |
| Stage | MV5-D1 |
| Input set | 90 existing-data samples × 5 frozen seeds |
| Outer split | 15 leave-one-study-out folds |
| Representation | SCT |
| View | Cell topology only |
| Outcome-label state | Closed |
| Authorized output | Training-fitted cell coordinates and provenance |
| Required stop | Before PH, landscapes, distances, clustering, integration, gene views, or outcomes |

## Purpose

MV5-D1 converts the 450 independently validated MV5-D0 normalization caches
into 75 immutable fold-seed coordinate bundles. It establishes the exact
cells-as-observations input that a later persistence-homology stage may consume.
It does not compute topology or make a biological comparison.

## Frozen inputs

1. The public 90-sample candidate manifest and 75-row fold-seed plan from
   MV5-C2.
2. The 450 runtime-complete MV5-D0 SCT matrix caches and their public resource
   manifest.
3. Seeds `20260805`–`20260809`.
4. Exactly 384 deterministically selected cells per sample and seed.
5. Exactly 15 held-out studies and 90 samples in every fold-seed job.

Only `sample_id`, `study`, the outer-fold role, and immutable technical
provenance may enter this stage. Tissue, approach, other endpoints, and derived
outcomes are prohibited.

## Fold-local transformation

For each held-out study and seed:

1. Validate and load the exact 90 SCT caches for that seed.
2. Define training and held-out sample IDs only from the frozen study split.
3. Form the feature universe from the intersection of **training matrices
   only**. Held-out feature availability must not affect selection.
4. Exclude non-biological feature categories under the existing MV-03 rule.
5. Rank each retained gene within each training sample by variance; order genes
   by median training rank with deterministic gene/feature tie breaks; retain
   500 genes.
6. Estimate gene-wise centers and scales from pooled training cells only.
7. Apply those frozen parameters to both training and held-out matrices. If a
   held-out matrix lacks a training-selected gene, map that feature to the
   training mean (zero after z-scoring); never refit the panel from held-out
   availability. This constant coordinate contributes zero to every within-
   sample pairwise cell distance and therefore cannot invent Rips topology.
8. Fit deterministic `stats::prcomp` on training cells only, without additional
   centering or scaling, and retain 30 components.
9. Project training and held-out cells through the one training-fitted rotation.
10. Store typed 384-cell × 30-PC cell-coordinate views. Do not invoke Ripser or
    any other PH implementation.

The observation definition is therefore explicit: each sample contributes 384
cell points in the same 30-dimensional coordinate system fitted without the
held-out study.

## Identity and resumability

Every private fold cache must bind:

- fold, fit scope, held-out study, and seed;
- sorted training and query sample IDs;
- all 90 MV5-D0 normalization cache keys;
- candidate-manifest and fold-plan file hashes;
- implementation hash and numerical runtime identity;
- panel size, component count, and transformation contract IDs; and
- payload, PCA-model, standardization, and coordinate-view hashes.

Publication is atomic. A missing cache may be built; a valid identity match may
be reused; an unreadable, invalid, or stale cache must never be overwritten
silently.

## Resource and execution gates

| Gate | Limit |
|---|---:|
| Heavy workers | At most 2 |
| Per-entry elapsed time | 1,800 seconds |
| Per-entry process-tree RSS | 8 GiB |
| Total MV5-D1 private cache | 40 GB |
| Expected entries | 75 |

A bounded pilot must pass functional validation, determinism, held-out
perturbation leakage tests, and resource gates before the full queue is
admitted. The queue stops admitting work after the first failed job or measured
cap violation, and completed valid entries remain resumable.

## Public evidence

Public manifests may contain sample/study split identifiers, dimensions,
counts, hashes, runtime versions, resource measurements, and explicit zero
counters. They must not contain expression matrices, cell coordinates, feature
panels, tissue or approach labels, source paths, or biological outcomes.

An independent validation pass must reopen every private cache and verify its
identity, payload, PCA fit scope, all 90 coordinate views, file hash/size,
label-closed status, and these zero counters:

- PH jobs;
- landscape jobs;
- distance jobs;
- clustering jobs;
- integration jobs;
- gene-view jobs; and
- biological-outcome jobs.

## Landscape boundary

The dissertation-aligned landscape definition is unchanged and is not executed
here: H0 and H1 remain separate, every active consecutive level is retained,
and integration is exact or error-controlled. A later stage must consume the
MV5-D1 coordinates without changing that definition.

## Completion decision

MV5-D1 completes only if all 75 fold-seed caches pass independent validation,
the resource envelope holds, downstream counters remain zero, the full test
suite and package check pass, and a post-fold feasibility projection is
recorded. Completion authorizes review of a later cell-primary PH stage; it
does not authorize that stage automatically.
