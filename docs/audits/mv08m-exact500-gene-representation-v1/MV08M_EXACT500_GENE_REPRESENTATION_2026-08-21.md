# MV8-M exact-500 gene-representation audit

**Date:** 2026-08-21
**Status:** complete; exact-500 Pearson-residual geometry is feasible, but adoption is owner-gated

## Decision

The 500-gene failure is not caused by absent genes, too few QC-eligible cells,
or an inherent inability to form an exact-500 gene topology.  It is specific
to the later production representation: the `SCT` `data` layer (log1p
corrected UMI counts) makes one exact-panel gene constant across all 4,614
frozen-QC cells.  Pearson correlation, and therefore correlation-chord
distance, is undefined for that gene.

The original corrected contract's SCT Pearson residuals retain all 500 genes,
are finite, have no effectively constant gene, and produce a valid exact-500
correlation-chord matrix.  A fixed ordinary log-normalized-count diagnostic
also passes, but it is not an adoption candidate.  Pearson residuals are the
recommended scientific direction because they restore the explicit
`DUAL_VIEW_TOPOLOGY_SPECIFICATION_V1` and MV03 method rather than inventing a
new repair.

Do **not** silently replace the existing MV05-D0/MV07 or HCA common-475
results.  Adopting Pearson residuals requires an explicit owner decision and a
prospective re-prefreeze covering the internal 124-sample reference plus all
eight external units.  PH, landscapes, clustering, fusion, labels, outcomes,
and manuscript claims remain closed here.

## What drifted

The corrected dual-view specification says that `sct_whole` uses SCT Pearson
residuals and that missing residual rows must be recovered from the recorded
SCT model.  MV03 follows that rule and explicitly rejects `SCT@data` as a
substitute.  MV05-D0 later froze complete `SCT data` matrices; MV07's exact-500
reference and the MV8 external runners inherited that later cache contract.

Official Seurat semantics distinguish these layers:

- `SCT data`: log1p corrected UMI counts;
- `SCT scale.data`: Pearson residuals.

MV8-M therefore diagnoses a representation-contract drift in the evolving
code, not a failure of the dissertation's updated landscape definition and
not evidence that a 475-gene biological panel is intrinsically preferable.

## Frozen execution

All three representations use the same HCA_BM_002 filtered H5, same ordered
500-gene panel, same ordered common-475 subset, and all 4,614 cells passing the
already-frozen QC rule.  Barcodes are deterministically ordered and QC occurs
before panel inspection.  One standard SCTransform fit supplies both the
later-production comparator and residual candidate; `Seurat::GetResidual`
materializes the exact residual panel from that fitted model.  The ordinary
log-normalized comparator uses scale factor 10,000 on the same RNA counts.

The frozen MV07 reference is used for identity and lineage only.  Its
data-layer center/scale vectors are not transferred to a different
representation.  Correlation-chord geometry is evaluated within each
representation on its unambiguously defined values.

## Exact-500 result

| Representation | Constant genes | Exact-500 chord geometry | Role |
|---|---:|---|---|
| SCT `data` / log1p corrected UMIs | 1 | invalid | later production comparator |
| SCT Pearson residuals | 0 | valid | corrected-contract candidate |
| RNA LogNormalize(10,000) | 0 | valid | diagnostic only |

The residual exact-500 distance hash is
`ae149db84c760c156d7f0f1ca5d3f38b947aba010bce03fe126ae8accbbc87f5`.
No matrix or failed-gene identity is published.

## Common-475 structural sensitivity

All three representations form valid common-475 geometries over 112,575
unique gene pairs.  The later SCT-data versus Pearson-residual comparison has:

- Pearson distance correlation: 0.9503;
- Spearman distance correlation: 0.9268;
- median absolute chord-distance difference: 0.00974;
- 95th-percentile absolute difference: 0.03062; and
- mean per-gene top-10-neighbor overlap: 0.6813.

Thus the representations are globally similar but not locally
interchangeable.  These measurements are descriptive.  MV08-G's thresholds
were frozen for 500-versus-475 panel sensitivity after PH and are not reused
as representation-selection thresholds.

## Determinism, resources, and failed attempt

The corrected primary and fresh repeat reproduce the representation summary
and stability table byte-for-byte.  Primary and repeat elapsed times are
559.74 and 567.16 seconds; peak RSS is 3,879,828 and 4,289,784 KiB.  Both are
within the one-worker, 1,800-second, 12-GiB envelope.

The first execution is retained as a non-scientific implementation failure.
All transforms completed, but scalar-first R `pmax`/`pmin` dropped the
correlation matrix dimensions before diagonal assignment.  A
dimension-preserving clipping correction passed a targeted regression.  The
failed attempt produced no admitted scientific output and stayed within caps.

## Recommendation and next gate

Recommend restoring SCT Pearson residuals as the candidate gene-side
representation.  Before adoption, prospectively freeze a migration/sensitivity
gate that:

1. rebuilds or recovers exact-panel Pearson residuals without `SCT@data`
   fallback;
2. verifies exact-500 viability and resources across the internal reference
   scope and all eight HCA units;
3. quantifies paired changes in gene H0/H1 persistence and all-active-level
   landscapes while keeping cell topology unchanged;
4. treats common-475 SCT-data results as preserved historical sensitivity
   evidence, not deleted or relabeled evidence; and
5. stops for owner review before clustering, fusion, outcomes, or claims.

## Firewalls

No labels, outcomes, biological comparisons, persistence, landscapes,
clustering, fusion, manuscript claims, other-unit matrices, feature dropping,
imputation, or deletion were opened.  Published artifacts contain aggregate
counts, metrics, and hashes only.
