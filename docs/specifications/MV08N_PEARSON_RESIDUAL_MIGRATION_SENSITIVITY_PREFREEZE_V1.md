# MV8-N Pearson-residual migration and sensitivity prefreeze v1

**Status:** prospectively frozen; only the bounded source/view sentinel is
authorized after this prefreeze is committed

**Date:** 2026-08-21

**Parent decisions:** D-084 through D-110

## Purpose

MV8-M established that the exact-500 failure is representation-specific:
later-production SCT log1p corrected UMI values leave one constant gene in
HCA_BM_002, while SCT Pearson residuals form a valid 500-gene
correlation-chord geometry.  The project owner approved Pearson residuals as
the candidate gene-side representation and authorized a prospective
internal-124/external-eight migration sensitivity sprint.

MV8-N freezes that sprint before any new persistence or landscape execution.
It does not make residuals the package default, reinterpret or replace prior
evidence, open outcomes, or authorize the complete workload.

## Scientific estimand and unchanged axes

The migration is gene-side only.  Cell topology, the immutable 30-PC cell
projection, selected-cell identities, panel order, filtration, homology
dimensions, and landscape mathematics remain unchanged.

- Internal scope: 124 retained samples, five deterministic seeds
  (`20260805`--`20260809`), and exactly 384 selected cells per sample and seed.
- External scope: eight exact-reference HCA bone-marrow units, seed
  `20260805`, and exactly 384 selected cells per unit after the frozen HCA QC
  rule.
- Exact panel: ordered 500-gene axis SHA-256
  `48e881ee753893bfaecd7101fc16fbcb552dd30a2a791fedfc7204d7b322a32e`.
- Common subset: ordered 475-gene axis SHA-256
  `b7b802ca862a63d7a4bbcaeab5af1192577663992a5ebde831371b6efafbc0ba`.

Gene topology continues to treat genes as points and uses Pearson
correlation-chord distance `sqrt(2 * (1 - r))` across the frozen 384 cells.
No gene may be dropped, padded, imputed, substituted, or rescued by changing
the topology metric.

## Normalization fit scope and representations

One standard `Seurat::SCTransform` model is fit per internal sample or external
unit using all cells in that unit's already-frozen QC-eligible raw-count
source.  The model uses `min_cells=5`, `variable.features.n=3000`,
`return.only.var.genes=FALSE`, one worker, and the bound Seurat/sctransform
runtime.  The exact panel's Pearson residuals are materialized with
`Seurat::GetResidual` from that model and RNA counts.  Each seed-specific view
then subsets the model output to its pre-existing 384 selected cells.

This all-QC model fit is used uniformly across all 132 units.  It avoids a
unit-specific low-count fallback and makes model estimation invariant to the
five internal cell-selection seeds.  The topology observation axis remains
384 cells; cells outside that axis contribute only to normalization-model
estimation.

The paired sensitivity distinguishes three effects:

1. `sct_data_selected384_fit`: the existing selected-384 SCT-data contract;
2. `sct_data_all_qc_fit_selected384`: log1p corrected UMIs from the new
   all-QC model, evaluated on the same 384 cells; and
3. `sct_pearson_residual_all_qc_fit_selected384`: Pearson residuals from that
   same all-QC model and same 384 cells.

Thus existing versus all-QC SCT data isolates fit-scope effects, all-QC SCT
data versus all-QC residuals isolates layer effects, and existing versus
residuals measures the net migration effect.

## Internal-124 workload

The existing exact-500 selected-fit SCT-data gene PH and landscapes are
immutable comparators and must rehash before reuse.  MV8-N proposes two new
exact-500 gene views for each of 620 sample–seed axes: all-QC-fit SCT data and
all-QC-fit Pearson residuals.

After source admission, the internal candidate workload is:

- 124 all-QC SCT model fits reused across five seeds;
- 1,240 new gene views (620 per new representation);
- 1,240 new complete VR PH records, each containing separate H0/H1;
- 20 new landscape-distance groups (two representations by five seeds by
  H0/H1), each containing all 7,626 unordered sample pairs; and
- 30 paired comparison strata (fit-scope, layer, and net migration effects by
  five seeds by H0/H1).

The exact existing 620 gene records and ten existing gene landscape groups
remain read-only comparators.  Cell PH and cell landscapes are not rerun.

## External-eight workload

The existing MV8-I common-475 result came from the historical reference H5
axis and is retained as historical replication evidence, not used as the
primary paired comparator for the newly processed exact-reference matrices.
For each exact-reference HCA unit, MV8-N freezes:

- selected-384-fit SCT data on common475;
- all-QC-fit SCT data on common475;
- all-QC-fit Pearson residuals on common475;
- all-QC-fit Pearson residuals on exact500; and
- all-QC-fit SCT data on exact500 as a feasibility diagnostic only.

The first four stacks form 32 new PH records.  The exact500 SCT-data
diagnostic is not advanced to PH unless all eight views are valid, and it is
never required for residual migration.  New residual/common475/exact500
landscapes retain every active consecutive level and exact integration.

External paired comparisons are fit-scope, layer, common475 net migration,
residual panel, and total exact500 migration effects, each evaluated
separately in H0 and H1 across the 28 unordered unit pairs.

## Persistence and landscape contract

Every admitted gene view receives complete Vietoris--Rips H0/H1 persistence
over field 2 with no filtration truncation.  Essential H0 is excluded before
landscape construction.  Finite positive-persistence intervals contribute to
all consecutive active landscape levels; paired diagrams are zero-padded only
above their active depths; squared-L2 distance is integrated exactly at
critical breakpoints.  There is no fixed grid, universal level cap, or dense
landscape-function matrix.

Ripserr remains primary.  The previously approved GUDHI exact fallback may be
used only after a recorded Ripserr resource failure and within a 12-GiB cap.
The R exact implementation remains the numerical oracle; the hash-bound Rust
kernel is execution-only.

## Comparison metrics and interpretation

For every paired landscape-distance stratum, publish:

- Pearson and Spearman correlation across the complete unordered-pair axis;
- median and 95th-percentile absolute distance change;
- mean per-unit top-10-neighbor overlap (or all other seven units externally);
- distance-matrix and ranking hashes; and
- per-stratum determinism and completeness checks.

These metrics are descriptive, not equivalence tests.  MV08-G's panel
thresholds are not repurposed.  Technical migration admission requires valid
exact500 residual geometry, complete pair axes, deterministic repeat, and
resource compliance—not a favorable structural result.  Final default
adoption and any biological interpretation remain owner-gated after results.

## Staged execution and authorization

1. **Input closure:** rehash 124 internal raw shards, 620 selected-cell axes,
   eight exact-reference H5 files, both panels, and immutable comparators.
2. **Source/view sentinel (authorized next):** run the minimum-cell and
   maximum-cell internal units across all five seeds plus HCA_BM_002; repeat
   the maximum internal unit and HCA_BM_002.  No PH is computed.
3. **Full source production (closed):** build 132 resumable residual caches and
   validate all proposed selected views only after a separate sentinel
   closure.
4. **PH sentinel and full PH (closed):** require separate prospective engine,
   cross-engine, repeat, and resource admission.
5. **Landscape stress and completion (closed):** select the deterministic
   maximum interval-burden group after PH, then require R-oracle, repeat,
   resource, and immutable-resume admission.
6. **Paired structural comparison (closed):** run only after all source, PH,
   and landscape artifacts are immutable and complete.

No later stage is implicitly authorized by an earlier pass.

## Resource and storage policy

- Source/model sentinel and production: one worker, no retry, 1,800 seconds
  and 12 GiB per unit.
- PH: one worker, no retry, 1,800 seconds and 12 GiB per gene view, with
  tighter measured caps frozen after the PH sentinel.
- Landscape group: one worker, no retry, 3,600 seconds and 12 GiB.
- Aggregate heavy execution is serial.  Existing unrelated processes are not
  killed, paused, or modified.
- Private residual caches, matrices, diagrams, and landscapes remain ignored
  and resumable by exact identity.  Public artifacts contain aggregate
  metrics and hashes only.

The source sentinel must measure cache size and project total retained state
before full production.  A projected total above 40 GiB, a per-process cap
breach, a source identity mismatch, or any non-deterministic repeat stops the
sprint for a revised resource decision.

## Firewalls and stop rules

- Labels, outcomes, tissues, studies, approaches, benchmark results, and prior
  cluster assignments are prohibited from execution and stopping decisions.
- No clustering, fusion, method selection, manuscript claim, default change,
  or deletion is authorized.
- A failed exact-panel feature or correlation geometry stops the applicable
  global representation stack; no per-unit feature dropping or fallback is
  allowed.
- Existing SCT-data/common475 evidence is preserved and never relabeled as
  Pearson-residual evidence.
- Public outputs may not contain expression values, residual matrices,
  barcodes, private paths, diagrams, or donor-private metadata.
