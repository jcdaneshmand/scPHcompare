# MV-03 corrected pilot diagrams and feasibility boundary

| Field | Value |
|---|---|
| Date | 2026-08-05 |
| Scope | Existing-data technical feasibility only |
| Contract | `DUAL_VIEW_TOPOLOGY_SPECIFICATION_V1.md` |
| Stages | A, B, and C completed in frozen order |
| PH jobs | 132 completed; 0 failed |
| Biological interpretation | Prohibited |
| Landscape/distance/clustering work | Not performed; deferred to MV-04/MV-05 |
| Gate | G-MV3 passes for technical advancement to MV-04 |

## Outcome

MV-03 generated the first scientifically eligible corrected H0/H1 diagrams for
both `cell_topology_v1` and `gene_topology_v1`. Every job used the frozen
500-gene panel, 384 matched cells, representation-specific fit-scope
standardization, 30 shared PCs for the cell view, correlation-chord geometry
for the gene view, and complete Vietoris-Rips filtration (`threshold = -1`).

The staged execution completed 4 Stage A jobs, 48 Stage B jobs, and 80 Stage C
jobs. No sample, representation, view, seed, or failed result was removed or
replaced.

## Representation-native source resolution

The historical `expr_list_sctWhole*.Rds` files contain `SCT@data`, not Pearson
residuals, and are therefore ineligible for the corrected source contract.
Direct inspection of both saved merged Seurat objects established that:

- each object contains one recorded `SCTModel`;
- `SCT@scale.data` contains 3,000 Pearson-residual features for every cell;
- every frozen sample has the exact manifest cell count;
- all extracted values and canonical cell identifiers are valid.

For bone marrow, the retained residual features supplied more than 500 genes
common to the 6,000-feature integrated representation. For the large cohort,
only 364 retained residual features overlapped the independently retained 3,000
integrated features. The preflight stopped as required. The missing 2,636
features were then calculated with `Seurat::GetResidual` from the recorded SCT
model and RNA counts after subsetting the object to the ten frozen pilot
samples. No `SCT@data` values, zero padding, or alternate representation was
substituted.

Historical SCT and integrated cell prefixes differ. Pairing therefore uses the
canonical barcode suffix following `__`, proves uniqueness and exact set
equality per sample, and reindexes by identifier. Pairing by column position is
prohibited.

## Feature panel and fitted transformations

Panels were selected independently for the large and bone cohorts under
`descriptive_all_pilot_samples`:

1. require presence in every pilot sample and both primary representations;
2. require finite variance above `.Machine$double.eps` in every sample and
   representation;
3. exclude the versioned `^MT-`, `^RP[SL]`, and `^HB(?!P)` categories;
4. rank SCT Pearson-residual variance within each sample;
5. aggregate equally weighted sample ranks by median;
6. break ties by canonical gene symbol and retain the first 500.

For each cohort, representation, and seed, exactly 384 cells from every fit
sample contributed to fitted gene means, standard deviations, and shared PCA.
The same selected cell identifiers and gene order feed paired representations
and both views. Full panel, rank, exclusion, subsampling, scaling, and PCA audit
tables are retained under `docs/audits/`.

## Frozen stage results

| Stage | Scope | Jobs | Failures | Worker time |
|---|---|---:|---:|---:|
| A | Large global min/max; SCT whole; first seed; both views | 4 | 0 | 5.23 s |
| B | Eight large core plus four bone samples; both representations; first seed; both views | 48 | 0 | 73.22 s |
| C | Four sentinels; both representations; five seeds; both views | 80 | 0 | 83.16 s |
| Total | Frozen MV-03 execution | 132 | 0 | 161.61 s |

The slowest job was a large-cohort Seurat-integration gene view at 14.64
seconds. The maximum measured process-tree RSS was 1,121,476,608 bytes (about
1.04 GiB), also in that stratum. Both remain far below the frozen 30-minute and
8-GiB per-job caps. All diagrams contain one explicit essential H0 interval and
nontrivial finite H1 intervals.

The heavy cost is representation-native extraction, not PH at pilot scale. The
large 16.9-GB Seurat object required 496 seconds to read during model recovery,
and residual recovery required another 30 seconds. Ignored, hash-bound pilot
caches prevent repeating this cost during MV-03. Later production design should
materialize scientifically eligible representation-native inputs directly.

## Seed sensitivity boundary

Stage C treats the five seeds as repeated measurements, not independent
biological samples. Across 32 sample-by-representation-by-view-by-dimension
strata, the median coefficient of variation for total persistence was 0.023 and
the maximum was 0.120. The largest observed value occurred for large-cohort
SCT-whole cell-view H1 in the maximum-cell sentinel. Detailed interval-count,
total-persistence, maximum-persistence, and selected-cell-overlap ranges are in
`mv03-seed-stability-2026-08-05.csv`.

These descriptive summaries show bounded numerical/subsampling behavior within
the pilot. They do not establish biological reproducibility and do not replace
the all-level landscape and diagram-distance validation required by MV-04.

## Gate decision

- Scientifically eligible H0/H1 diagrams exist for both views: **pass**.
- Artifact point counts, axes, metrics, filtration, field, essential H0 policy,
  hashes, sample IDs, representations, and seeds are explicit: **pass**.
- Runtime and peak process-tree memory are measured under enforced caps:
  **pass**.
- Five-seed sensitivity is retained and bounded at the summary level:
  **pass for MV-04 technical advancement; not a biological claim**.
- Pairwise distances, landscapes, clustering, fusion, and biological analysis:
  **not authorized by G-MV3 and not performed**.

| Question | Disposition |
|---|---|
| Scientific contract coherent? | approve for corrected pilot artifacts |
| Correctness demonstrated? | pass at typed input/PH/provenance boundary |
| Computation feasible? | yes at frozen pilot scale |
| Biological interpretation permitted? | prohibited |
| Next action | advance to MV-04 distance and production-engine validation |

## Auditable evidence map

- Source/model audits: `mv03-sct-residual-*.csv` and
  `mv03-integrated-*.csv`.
- Feature eligibility and ranks: `mv03-gene-*.csv` and
  `mv03-sct-variance-ranks-*.csv`.
- Cell/PCA/scaling provenance: `mv03-stage-*-matched-cells-*.csv`,
  `mv03-stage-*-standardization-*.csv`, and
  `mv03-stage-*-pca-summary-*.csv`.
- Immutable job/result identities and resource measurements:
  `mv03-stage-*-job-manifest-*.csv` and `mv03-stage-*-ph-metrics-*.csv`.
- Artifact-by-artifact hash, contract, axis, filtration, interval, and identity
  validation: `mv03-bundle-validation-2026-08-05.csv` (132/132 pass).
- Diagram summaries and stability: `mv03-diagram-summary-2026-08-05.csv` and
  `mv03-seed-stability-2026-08-05.csv`.
- Corrected diagram RDS artifacts: `results/mv03/diagrams/`.
