# MV8-N Pearson-residual migration sensitivity prefreeze

**Date:** 2026-08-21

**Contract:** `mv08n_pearson_residual_migration_prefreeze_v1`

**Accepted source commit:** `547f5b7c721a223f042fe4c1c473d0875801f2db`

**Result:** prefreeze passed; only the bounded source/view sentinel is authorized

## Decision

The owner-approved candidate gene representation is SCT Pearson residuals from
one standard SCTransform model fitted to every cell passing the already-frozen
QC rule for each sample or external unit. Gene topology is then evaluated on
the unchanged deterministic 384-cell selections. Pearson residuals are a
candidate, not yet the project default.

This prefreeze does not execute SCTransform, recover residual values, calculate
persistence, calculate landscapes, compare results, open labels or outcomes,
cluster, fuse views, or replace historical evidence.

## Why the paired design is necessary

MV8-M showed that exact-500 feasibility depends on the transformed
representation. It also exposed two changes that must not be conflated:

1. fitting SCTransform on all frozen-QC cells instead of the selected 384;
2. using Pearson residuals instead of SCT `data` values.

The internal sensitivity therefore compares the existing selected-384-fit SCT
`data` baseline with all-QC-fit SCT `data` and all-QC-fit Pearson residuals on
the same 384 cells, exact 500 genes, seeds, and topology contract. This
separates fit-scope, representation-layer, and net migration effects.

The external sensitivity adds common-475 bridges and an exact-500 residual
view. All-QC-fit SCT `data` on exact500 is retained only as a source diagnostic
because MV8-M found that geometry invalid in HCA_BM_002; no PH job is queued for
that view.

## Frozen inventories

- Internal: 124 raw-count shards (90 primary plus 34 added), 620 existing
  sample-seed selections, five seeds (`20260805`--`20260809`), and 384 cells per
  selection. Every live raw file was independently rehashed and matched its
  accepted audit.
- External: eight current exact-reference HCA filtered matrices. Every matrix
  was independently rehashed. Frozen-QC cell counts range from 2,571 to 4,614,
  so every unit clears the 384-cell observation requirement.
- Panels: ordered exact500 SHA-256
  `48e881ee753893bfaecd7101fc16fbcb552dd30a2a791fedfc7204d7b322a32e`
  and ordered common475 SHA-256
  `b7b802ca862a63d7a4bbcaeab5af1192577663992a5ebde831371b6efafbc0ba`.
- New all-QC model fits: 132 (124 internal plus eight external).
- Proposed new gene views: 1,280 (1,240 internal plus 40 external).
- Planned PH records: 1,272. The eight external all-QC-fit SCT-data exact500
  diagnostic views are excluded.
- Planned exact all-active-level landscape groups: 28 (20 internal plus eight
  external).
- Planned descriptive paired comparison strata: 40 (30 internal plus ten
  external).

Cell topology and its 30-PC coordinates are immutable comparators and will not
be recomputed.

## Landscape contract

The migration preserves the dissertation-aligned definition: finite positive
intervals only, the essential H0 class excluded, all active consecutive
landscape levels retained, and exact critical-breakpoint squared-L2
integration. There is no uniform grid and no universal level cap. H0 and H1
remain separate throughout.

## Stage firewall and resource policy

Execution is serial with one worker and no retry. Each source fit/view sentinel
has an 1,800-second and 12-GiB cap. Future PH jobs retain an 1,800-second and
12-GiB cap; future landscape work retains a 3,600-second and 12-GiB cap. The
prospective private-artifact ceiling is 40 GiB.

Only these three source/view sentinel units become executable after this
prefreeze is committed:

| Scope | Unit | QC cells | Frozen selection axes | Repeat |
|---|---|---:|---:|---|
| Internal minimum | `SRA701877_SRS3279688` | 396 | 5 | no |
| Internal maximum | `SRA742961_SRS3565197` | 11,475 | 5 | yes |
| External bridge | `HCA_BM_002` | 4,614 | 1 | yes |

The remaining 129 source fits, all PH, all landscapes, and all comparisons stay
closed. Passing a stage never automatically opens the next stage. A failure
preserves partial evidence and stops without retry or scientific-contract
substitution.

## Validation and privacy

The metadata-only builder passed 18/18 independent checks. Public tables contain
source identities, hashes, queue contracts, and aggregate dimensions only. Raw
counts, expression values, residuals, barcodes, labels, outcomes, and private
cache paths are not published. The builder itself does not open matrix payloads
or call SCTransform, `GetResidual`, PH, landscape, clustering, or fusion code.

## Next gate

MV8-O is the bounded source/view sentinel. It must freeze the seven remaining
external selected-cell axes, run the three authorized units under monitoring,
repeat the maximum internal and HCA_BM_002 units, and demonstrate exact500
finite/nonzero-variance residual geometry plus common475 bridges before any
full-source-production decision.
