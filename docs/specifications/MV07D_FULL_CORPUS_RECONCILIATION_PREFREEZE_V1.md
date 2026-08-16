# MV7-D full-corpus reconciliation and omitted-source feasibility prefreeze

Date: 2026-08-15  
Status: prospective prefreeze; no expanded PH or outcome calculation

## Purpose

Reconcile every historical sample against the corrected pipeline before a
manuscript claim map is treated as final. This sprint separates populations by
estimand instead of calling all of them “the full data.” It also freezes a
bounded source/SCT feasibility check for samples omitted from the primary
90-sample benchmark.

## Frozen sample flow

1. The public heterogeneous-tissue metadata contains 127 candidate samples.
2. The historical `MIN_CELLS = 250` post-QC rule excludes three substantia
   nigra samples before PH, leaving 124 retained samples.
3. Five tissues occur in at least two studies. Their 90 samples support the
   primary leave-one-study-out tissue-retrieval estimand.
4. Pancreatic islets (12), prostate (16), and retained substantia nigra (6)
   each occur in one study. These 34 samples are not primary failures or
   missing results: they are structurally inestimable for cross-study tissue
   retrieval and remain eligible for corrected descriptive topology.
5. The three below-250 samples are a separate threshold-sensitivity population.
   They cannot satisfy the fixed 384-cell contract and are not admitted here.

The separately stored 25-library GSE120221 bone-marrow processing path remains
a technical validation population. MV8-B later established that its 25
biological libraries are the same `SRA779509` libraries already inside the
127/124 flow; it is not an additional independent cohort.

Study and tissue fields agree between the public candidate metadata and the
historical retained ledger. Approach labels disagree for 16 retained samples
(14 public scRNA/historical snRNA and two public snRNA/historical scRNA). The
reconciliation preserves both values and flags every disagreement. It carries
the historical retained value forward solely to remain consistent with the
accepted MV7-B confounding diagnostic; no technology claim may use these rows
until an accession-level source audit resolves the provenance.

## Prospective six-sample feasibility panel

The panel contains the minimum- and maximum-post-QC-depth sample within each of
the three single-study tissues. It uses seed `20260805` and exactly 384 selected
cells. Selection uses only tissue stratum and retained-cell depth; no PH,
landscape, distance, clustering, retrieval, or method-ranking value is read.

The authorized execution reuses the accepted individual-Seurat-to-raw-shard,
deterministic cell-selection, and v2 SCT-cache code paths. It must run
sequentially with one heavy child, a 1,800-second per-job cap, and an 8-GiB
process-tree RSS cap. Private source paths, caches, and logs remain ignored.

This check does not authorize a persistence diagram. A corrected full-corpus
cell view first needs a frozen descriptive fit scope; the gene view needs a
frozen 124-sample panel-availability and transform scope. Those choices must be
defined before topology is calculated.

## Landscape invariant

Any later 124-sample expansion must use the dissertation-aligned corrected
definition without modification:

- all finite positive-persistence intervals, with the essential H0 interval
  excluded;
- all consecutive active landscape levels;
- H0 and H1 calculated and reported separately;
- exact or error-controlled squared-L2 integration on the dimension-specific
  support;
- no universal uniform grid and no universal landscape-level cap; and
- streamed or chunked computation rather than dense landscape materialization.

## Prohibited actions

- Recalculate or alter the accepted primary 90-sample result.
- Treat the 34 single-study samples as primary endpoint failures.
- Use legacy feature-as-point PH as corrected evidence.
- Fit a 124-sample PCA or gene panel, calculate PH, open expanded outcomes, or
  promote manuscript claims in this prefreeze.
- Admit the three below-250 samples by silently changing the 384-cell depth or
  the historical 250-cell rule.
- Add new biological data.

## Acceptance criteria

- A one-row-per-candidate 127-row reconciliation is produced.
- Counts reconcile exactly to 124 retained, 90 primary, 34 descriptive-only,
  and three threshold-sensitivity candidates.
- All 127 sparse sources and all 127 individual Seurat sources are uniquely
  present locally; the 90 accepted corrected sample identities are complete.
- The six-sample feasibility panel is deterministic and prospectively frozen.
- The revised landscape contract is carried forward unchanged.
- All 16 approach-label disagreements are explicit rather than silently
  resolved.
- Independent validation passes before the bounded feasibility run begins.

## Gate

On prefreeze acceptance, authorize only `MV7-D1`: the six-sample source/SCT
feasibility run. If it passes, the next artifact is a prospective 124-sample
descriptive cell/gene fit-scope and resource plan. Expanded PH and outcomes
remain separately gated.
