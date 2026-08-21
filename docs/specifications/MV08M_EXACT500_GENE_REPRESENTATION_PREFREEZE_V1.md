# MV8-M exact-500 gene-representation prefreeze v1

**Status:** prospectively frozen for one label-closed HCA_BM_002 sentinel only

**Date:** 2026-08-21

**Parent decisions:** D-103 through D-109

## Purpose

MV8-L showed that all 500 ordered genes are present and retained when the
gene-side source uses all 4,614 cells passing frozen QC, but one gene is
exactly constant in the `SCT` `data` layer.  The unchanged Pearson
correlation-chord distance is therefore undefined.  MV8-M tests whether this
failure belongs to the 500-gene panel or to the expression representation.

This is a representation feasibility audit, not a topology rerun or an
automatic method change.  It compares three named representations on exactly
the same cells and genes and stops before persistence, landscapes, clustering,
fusion, labels, outcomes, or manuscript claims.

## Reconciled representation history

The project has two incompatible recorded SCT contracts:

1. `DUAL_VIEW_TOPOLOGY_SPECIFICATION_V1` and the corrected MV03 pilot define
   `sct_whole` as SCT Pearson residuals.  MV03 explicitly rejects historical
   `SCT@data` values as substitutes and uses `Seurat::GetResidual` when a
   requested feature is absent from `SCT@scale.data`.
2. MV05-D0 later freezes complete `SCT` `data` matrices, whose Seurat semantics
   are log1p corrected UMI counts.  MV07's exact-500 reference and the MV8-H
   through MV8-L external runners consume that later cache contract.

MV8-M records this as method-contract drift.  It does not retroactively erase
either evidence chain.  A passing Pearson-residual candidate would support the
original corrected definition, but adoption would require a separately
approved re-prefreeze and paired internal/external sensitivity analysis.

## Frozen inputs and common scope

- Unit: `HCA_BM_002` only.
- Filtered Cell Ranger 8.0.1 H5 SHA-256:
  `7b0a9f90b05d68c14af29cea39926f74f3c52d07c8825bdc4f01303bb6ad0e2f`.
- Ordered exact-500 panel SHA-256:
  `48e881ee753893bfaecd7101fc16fbcb552dd30a2a791fedfc7204d7b322a32e`.
- Ordered common-475 subset SHA-256:
  `b7b802ca862a63d7a4bbcaeab5af1192577663992a5ebde831371b6efafbc0ba`.
- Frozen unit-level QC: 500--9,000 detected features, mitochondrial fraction
  at most 20%, ribosomal fraction above 5%, and positive total counts.
- Expected eligible scope: all 4,614 cells, deterministically ordered by
  barcode after QC and before panel inspection.
- The MV07 reference record is read only to bind the exact ordered panel and
  production lineage.  Its center/scale vectors are not applied across the
  three representations because they were fitted under the later `SCT data`
  contract.

No gene-aware cell selection, feature dropping, padding, imputation,
replacement, reordering, or identifier substitution is permitted.

## Frozen representations

One standard `Seurat::SCTransform` fit is run on all eligible cells with
`return.only.var.genes=FALSE`, `variable.features.n=3000`, `min_cells=5`, one
worker, and otherwise documented defaults.  It supplies two candidates from
the same fitted model:

1. **`sct_data_log1p_corrected_umi`** -- the ordered panel extracted from the
   `SCT` assay `data` layer.  This is the later MV05-D0/MV07 production
   representation and the direct MV8-L comparator.
2. **`sct_pearson_residual`** -- the ordered panel materialized with
   `Seurat::GetResidual` from the fitted SCT model and RNA counts, then
   extracted from `SCT@scale.data`.  No `SCT@data` fallback is allowed.
3. **`rna_lognormalize_10000`** -- ordinary library-size normalization of the
   same RNA counts with `Seurat::NormalizeData(method="LogNormalize",
   scale.factor=10000)`, followed by the RNA `data` layer.  This is a fixed
   diagnostic comparator, not an adoption candidate.

## Frozen geometry and outputs

Within each representation independently:

1. require exactly 500 ordered rows and 4,614 ordered columns;
2. require all values finite;
3. center each gene across the same cells for the variance test;
4. reject any gene with standard deviation at or below
   `sqrt(.Machine$double.eps)`;
5. compute the 500-by-500 Pearson correlation matrix across cells and the
   unchanged chord distance `sqrt(2 * (1 - r))`;
6. require finite, symmetric, nonnegative distances with a zero diagonal; and
7. extract the ordered common-475 principal submatrix.

Only aggregate summaries and hashes may be published.  No expression matrix,
gene-level failure identity, barcode, correlation matrix, or distance matrix
may leave ignored execution storage.

For each pair of valid representations, MV8-M descriptively records Pearson
and Spearman correlation of the 112,575 unique common-475 off-diagonal
distances, median and 95th-percentile absolute distance differences, and mean
per-gene overlap of the ten nearest neighbors.  These are descriptive and
have no adoption threshold.  MV08-G's 0.95 median-Spearman and 0.80 median
top-10 thresholds were frozen for a different panel-sensitivity estimand and
are not repurposed here.

## Admission, repeat, and resource gates

A representation is exact-500 viable only if it retains the exact ordered
panel, is finite, has zero effective constant genes, and forms a valid
correlation-chord matrix and valid common-475 submatrix.  The sentinel must
also:

- complete a fresh deterministic repeat with identical scientific aggregate
  fields and hashes, excluding elapsed time and peak RSS;
- remain within one worker, 12 GiB peak RSS, and 1,800 seconds per run; and
- preserve all firewalls below.

A viable Pearson-residual result supports advancing to an owner decision on a
new residual-based representation contract.  It does not itself authorize
that contract.  Log-normalized viability is diagnostic only.  If no
representation is viable, exact-500 remains blocked without repair.

## Stop rules and firewalls

- Stop rather than drop or alter a failed panel feature or change the topology
  metric.
- Do not inspect another HCA unit at matrix level.
- Do not compute persistence, landscapes, fusion, clustering, labels,
  outcomes, biological comparisons, or manuscript claims.
- Do not modify, delete, or relabel existing common-475 or exact-500 evidence.
- A representation adoption, internal/external rerun, or reinterpretation of
  existing results requires explicit project-owner approval after MV8-M.

## Authoritative method sources

- Hafemeister C, Satija R. *Normalization and variance stabilization of
  single-cell RNA-seq data using regularized negative binomial regression*.
  Genome Biology 20, 296 (2019). <https://doi.org/10.1186/s13059-019-1874-1>
- Seurat SCTransform reference, including `data` versus `scale.data` semantics
  and residual recovery controls.
  <https://satijalab.org/seurat/reference/sctransform>
- Seurat GetResidual reference.
  <https://satijalab.org/seurat/reference/getresidual>
- Seurat NormalizeData reference.
  <https://satijalab.org/seurat/reference/normalizedata>
- Skinnider MA, Squair JW, Foster LJ. *Evaluating measures of association for
  single-cell transcriptomics*. Nature Methods 16, 1082--1090 (2019).
  <https://doi.org/10.1038/s41592-019-0612-7>
