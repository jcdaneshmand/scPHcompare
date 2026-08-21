# MV8-L exact-500 view-specific observation-scope prefreeze v1

**Status:** prospectively frozen for one HCA_BM_002 feasibility sentinel only

**Date:** 2026-08-21

**Parent decisions:** D-103 through D-108

## Purpose

MV8-K established that the recovered HCA_BM_002 raw-read matrix contains all
500 ordered panel features, but the fixed, panel-agnostic 384-cell selection
does not supply nonzero transformed variance for three features.  This
prefreeze asks a narrower question: can the existing exact-500 *gene-side*
correlation-chord representation be formed from all cells that pass the
already-frozen QC rule, while the existing exact-500 384-cell/PCA cell-side
diagnostic remains unchanged?

This is a feasibility test, not an adoption of a new topology contract.  The
current typed dual-view source remains 500 genes by 384 cells; it will not be
modified or bypassed in MV8-L.

## Scientific rationale and limits

SCTransform fits regularized negative-binomial models to UMI data to control
technical depth variation; its documented default drops features expressed in
fewer than five cells.  Using all QC-eligible cells is therefore a direct test
of whether the MV8-K loss was induced by a fixed computational subsample rather
than an absent source feature.  It does not impute expression, rescue a gene by
symbol substitution, or select cells because of a gene's expression.

Gene–gene associations are estimated across cells, whereas the existing cell
view uses cells as observations and an immutable 30-PC projection.  A large
benchmark of single-cell association measures found that larger cell counts
were associated with more functionally coherent gene networks, while also
showing that association-metric choice matters.  Therefore MV8-L keeps the
frozen gene correlation-chord metric unchanged and tests only its observation
scope.

This rationale supports feasibility only.  It does **not** establish that
cell-side and gene-side views with different observation counts may later be
combined, fused, or described as the same estimand.

## Frozen sentinel inputs and operations

- Unit: `HCA_BM_002`, using the existing filtered Cell Ranger 8.0.1 H5 whose
  SHA-256 is `7b0a9f90b05d68c14af29cea39926f74f3c52d07c8825bdc4f01303bb6ad0e2f`.
- Exact panel: the existing ordered 500 stable IDs, SHA-256
  `48e881ee753893bfaecd7101fc16fbcb552dd30a2a791fedfc7204d7b322a32e`.
- Frozen unit-level QC rule: 500–9,000 detected features, mitochondrial
  fraction at most 20%, ribosomal fraction above 5%, and positive total counts.
  The QC rule must be evaluated before panel expression is inspected.  The
  expected sentinel aggregate is 4,614 eligible cells.
- Cell-side comparator: MV8-K's existing panel-agnostic deterministic 384-cell
  selection (seed `20260805`) and immutable 30-PC model.  It is read as an
  identity-bound comparator only; MV8-L does not rerun cell topology.
- Gene-side candidate source: all QC-eligible cells, in deterministic barcode
  order.  No gene-aware cell filtering, stratification, or re-sampling is
  permitted.
- Transform comparison: SCTransform with `return.only.var.genes=FALSE`, first
  at documented `min_cells=5`, then at `min_cells=1` only if needed to diagnose
  a retention threshold.  Both retain the same QC cells and all other options.
- Standardization: use the frozen exact-500 center and positive scale vectors
  only to verify finite values and an immutable identity binding.  No PCA is
  refit or used for the gene-side candidate.
- Gene-distance diagnostic: calculate the existing correlation-chord distance
  only if all 500 standardized gene vectors are finite and have effective
  nonzero variance.  Also calculate the corresponding common-475 submatrix
  comparison.  Do not compute persistence, landscapes, clustering, or fusion.

## Admission criteria

The sentinel is feasible only if one declared transform configuration has all
of the following:

1. all 500 ordered stable IDs retained exactly once;
2. 4,614 QC-eligible cells (or a documented deterministic count from the
   frozen QC rule), with no gene-aware filtering;
3. finite standardized values and no effectively zero-variance gene;
4. a finite symmetric 500-by-500 correlation-chord distance matrix with zero
   diagonal and nonnegative off-diagonal values;
5. finite corresponding common-475 distances and a documented structural
   comparison; and
6. resource use within one worker, 12 GiB peak RSS, and 1,800 seconds.

Passing MV8-L only authorizes a cross-unit **aggregate feasibility** prefreeze;
it does not authorize a new typed topology source, PH, landscapes, or adoption
of mixed-scope views.

## Stop rules and firewalls

- A failed retention, variance, distance, resource, determinism, input-binding,
  or firewall check closes the candidate configuration.  It may not be repaired
  by gene-aware sampling, padding, imputation, dropping genes, or changing the
  correlation-chord metric.
- No labels, outcomes, cell identities, expression values, matrices, diagrams,
  persistence, landscapes, clustering, fusion, manuscript claims, or deletions
  are opened or published.
- No other HCA unit may be inspected at matrix level unless the single sentinel
  passes all admission criteria and an auditable cross-unit prefreeze is first
  written.
- Any subsequent implementation of a source with different cell counts for the
  two views requires explicit project-owner authorization.

## Authoritative method sources

- Hafemeister C, Satija R. *Normalization and variance stabilization of
  single-cell RNA-seq data using regularized negative binomial regression*.
  Genome Biology 20, 296 (2019).
  <https://doi.org/10.1186/s13059-019-1874-1>
- Seurat SCTransform reference, including the documented feature-retention and
  fitted-model controls. <https://satijalab.org/seurat/reference/sctransform>
- Skinnider MA, Squair JW, Foster LJ. *Evaluating measures of association for
  single-cell transcriptomics*. Nature Methods 16, 1082–1090 (2019).
  <https://doi.org/10.1038/s41592-019-0612-7>
