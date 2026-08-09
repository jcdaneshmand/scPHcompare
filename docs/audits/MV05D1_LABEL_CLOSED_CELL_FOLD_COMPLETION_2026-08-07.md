# MV5-D1 label-closed SCT cell-fold completion audit

| Field | Value |
|---|---|
| Date | 2026-08-07 |
| Contract | `mv05d1_production_cache_validation_v1` |
| Scope | 15 LOSO studies × 5 seeds × 90 samples |
| Outcome-label state | Closed |
| Fold-seed caches | 75/75 complete and independently validated |
| Typed cell-coordinate views | 6,750 |
| PH, landscapes, distances, clustering, integration, gene views, outcomes | 0 |
| Decision | MV5-D1 complete; stop before production cell PH |

## Outcome

MV5-D1 is complete. The 450 frozen MV5-D0 SCT caches were transformed into 75
immutable, label-closed fold-seed bundles. Each bundle contains 90 typed cell
views with exactly 384 cell observations in a shared 30-PC coordinate system
fitted from the training studies only.

No persistence diagram, landscape, between-sample distance, clustering,
integration, gene view, or biological result was computed.

## Leakage-safe transform

Within every fold and seed:

1. the feature universe and 500-gene median-variance-rank panel were derived
   only from training samples;
2. gene centers and scales were estimated from pooled training cells only;
3. PCA was fitted only to training cells; and
4. the one frozen transform was applied to both training and held-out cells.

The code audit corrected an earlier helper that intersected feature
availability across training and held-out matrices before training-only
ranking. The production contract now derives the feature universe strictly
from training matrices.

## Held-out feature compatibility

An initial strict pilot admitted four jobs, completed two, and stopped on two
structured incompatibilities when held-out matrices lacked a small number of
training-selected genes. It opened no labels or downstream artifacts.

The replacement policy maps an absent held-out feature to its training mean,
which is zero after the frozen training z-score. Held-out availability never
changes the selected panel or fitted transform. Because the mapped feature is
constant over all cells in that sample, it contributes zero to every
within-sample pairwise cell distance and cannot create Rips topology.

Across the completed production run:

| Measurement | Result |
|---|---:|
| Held-out sample-seed views | 450 |
| Views with at least one mapped feature | 71 |
| Missing feature-view instances | 111 |
| Maximum absent features in one 500-gene view | 4 (0.8%) |

This policy is authorized for the cell-topology coordinate/PH route only. It
does not automatically authorize a cross-sample coordinate baseline, gene
topology, or integrated representation.

## Determinism and production resources

The corrected real-data pilot repeated one fold independently. All 90 cell
views, 1,036,800 coordinate values, the complete payload, and the serialized
file were exactly identical; maximum absolute coordinate difference was zero.

| Measurement | Result |
|---|---:|
| Fold-seed entries | 75/75 |
| Build failures | 0 |
| Independently reopened/validated | 75/75 |
| Coordinate views validated | 6,750 |
| Worker-hours | 2.3756 h |
| Transformation operation-hours | 1.1658 h |
| Median entry elapsed | 114.119 s |
| Maximum entry elapsed | 120.846 s |
| Maximum process-tree RSS | 2,098,720,768 B (1.95 GiB) |
| Total private cache storage | 894,548,038 B (0.895 GB) |
| Largest private cache | 11,928,875 B |
| Implementation SHA-256 | `0252c7da173976337d0f643026098a956080579af8fdf2546ccee00747c7045a` |
| Cache-set SHA-256 | `829a92ef27da19a522e1d6237414b6f24ed3cf27173624bd4158afb67af7ffc2` |
| Coordinate-bundle SHA-256 | `4452739f23970b6a314f7aa402fbcf8a0e1880dccdb21917f3daeaa01a70e620` |

The queue used at most two workers. The maximum job consumed 24.43% of the
8-GiB cap, 6.71% of the 30-minute cap, and the complete cache used 2.24% of the
40-GB cap.

## Independent validation

A separate pass reconstructed each expected identity from the candidate,
fold, SCT-resource, implementation, and runtime manifests. It reopened and
validated every record, PCA fit scope, payload hash, all 90 views, coordinate
set hash, file hash, and file size. All 75 records passed. Every downstream
execution counter remained zero.

## Post-fold feasibility projection

The earlier 4.3906-hour “cached fold” estimate bundled coordinate fitting, PH,
gene work, and baselines. It is not a valid estimate of the remaining cell-PH
stage and was not subtracted as if PH had already run.

For SCT cell-primary, the measured/previously projected known components are:

| Component | Worker-hours |
|---|---:|
| MV5-D0 normalization, measured | 2.5623 |
| MV5-D1 cell coordinates, measured | 2.3756 |
| Exact all-level landscape distances, prior projection | 3.5718 |
| Known-components lower bound | 8.5097 |
| Production cell PH | Unmeasured |

Therefore the total feasibility decision is incomplete even though the known
components are well below the 21.6-hour planning cap. The next stage must be a
bounded, label-closed production cell-PH profiling sprint before authorizing
all 6,750 diagram calculations.

## Landscape boundary

The dissertation-aligned definition is unchanged: H0 and H1 remain separate;
all active consecutive landscape levels are retained; and integration is exact
or error-controlled. MV5-D1 produced coordinates only, so no landscape choice
was silently approximated or activated.

## Privacy and scope

- Outcome labels remained closed.
- Tissue and approach were absent from the fold and resource artifacts.
- Private SCT caches, coordinate bundles, logs, PDFs, reviewer correspondence,
  and `example_run.r` remain outside Git.
- No dependency or lockfile was changed.
- Nothing was pushed.
- The final source-loaded suite passed 529/529 expectations, and the source
  package check reported `Status: OK`.

## Decision table

| Question | Disposition |
|---|---|
| Scientific contract coherent? | Pass for SCT cell coordinates; missing-feature policy restricted to cell topology |
| Correctness demonstrated? | Pass: leakage tests, exact repetition, 75/75 independent validation |
| Computation feasible? | Yes for coordinates; total cell-primary feasibility incomplete until cell PH is measured |
| Biological interpretation permitted? | Prohibited |
| Next action | Specify and run bounded MV5-D2 label-closed cell-PH profiling; do not launch full PH automatically |
