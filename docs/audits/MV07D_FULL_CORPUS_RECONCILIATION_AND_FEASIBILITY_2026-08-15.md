# MV7-D full-corpus reconciliation and omitted-source feasibility audit

Date: 2026-08-15
Status: complete; prospective full-corpus expansion plan only

## Outcome

The historical corpus and corrected benchmark now reconcile sample by sample:

`127 candidates -> 124 post-QC retained -> 90 primary cross-study samples`

The 34-sample difference between 124 and 90 is not a Ripser failure set.
Those samples belong to pancreatic islets (12), prostate (16), and retained
substantia nigra (6), each represented by only one study. They cannot estimate
the primary leave-one-study-out tissue-retrieval endpoint because no same-tissue
training study exists. They are eligible for a separately identified corrected
124-sample descriptive-topology analysis.

The other three candidates are substantia nigra samples with 166, 169, and 197
post-QC cells. The historical code excluded them explicitly before PH under
`MIN_CELLS = 250`. Existing evidence does not prove that these exact samples
threw a captured `ripserr` exception. They also cannot satisfy the corrected
384-cell depth contract, so they remain a separate, currently unauthorized
threshold-sensitivity population.

## Estimand-specific populations

| Population | N | Permitted role |
|---|---:|---|
| Historical candidates | 127 | Sample-flow denominator |
| Historical retained/reproduction | 124 | Legacy reproduction only; historical PH orientation is scientifically ineligible |
| Corrected primary cross-study | 90 | Accepted blocked tissue-retrieval estimand |
| Corrected full-corpus descriptive | 124 | Prospective topology, distances, and clustering without cross-study claims for single-study tissues |
| Below-threshold sensitivity | 3 | Separate reduced-depth feasibility/stability question only |
| GSE120221 validation processing | 25 | Separate historical processing path; MV8-B later proved exact overlap with 25 `SRA779509` libraries already inside the 124-sample flow, so it is technical rather than independent validation |

The complete one-row-per-candidate ledger is
`mv07d-prefreeze-evidence/mv07d-sample-reconciliation.csv`.

## Source and depth feasibility

All 127 candidate sparse matrices and all 127 individual Seurat sources exist
locally with exactly one individual source per sample. Corrected accepted
artifacts cover the primary 90; they do not cover the 34 descriptive-only
samples.

Every one of the 34 retained omitted samples supports 384 selected cells:

| Tissue | N | Minimum | Median | Maximum |
|---|---:|---:|---:|---:|
| Pancreatic islets | 12 | 396 | 1,885 | 2,992 |
| Prostate | 16 | 4,751 | 6,335 | 11,475 |
| Substantia nigra | 6 | 584 | 713 | 1,015 |

A prospectively frozen six-sample panel selected the minimum and maximum depth
within each tissue. All six passed the accepted individual-source extraction,
canonical sparse raw-cache, deterministic seed-`20260805` 384-cell selection,
and v2 SCT matrix-cache path. Observed source counts exactly matched the
historical ledger; every SCT payload had 384 cells and finite values.
Independent validation reopened and rehashed all private source/raw/SCT
artifacts, reconstructed the resource projection, and passed 9/9 categories.

| Measurement | Observed maximum or total |
|---|---:|
| Raw-source child time, six total | 138.162 s |
| SCT child time, six total | 119.562 s |
| Maximum raw process-tree RSS | 1,559,584,768 B |
| Maximum SCT process-tree RSS | 1,581,699,072 B |
| Six raw caches | 65,650,937 B |
| Six SCT caches | 36,360,528 B |

Simple linear planning projects 34 raw shards plus 170 five-seed SCT caches at
about 1.158 worker-hours and 1.402 GB. This is an upstream planning estimate,
not an authorization and not a PH/landscape projection.

## Approach-label discrepancy

Study and tissue fields agree between the public 127-row metadata and the
historical retained ledger. Approach labels disagree for 16 retained samples:
14 are public `scRNA-seq`/historical `snRNA-seq`, and two are public
`snRNA-seq`/historical `scRNA-seq`. MV7-B used the historical ledger. MV7-D
preserves both values and flags each row; no causal or stratified technology
claim should proceed until accession-level source provenance resolves them.

This discrepancy does not change the tissue/study eligibility counts or the
primary 90-sample identity.

## Landscape contract

The full-corpus plan carries forward the revised dissertation-aligned
definition unchanged: finite positive-persistence intervals, essential H0
excluded, every consecutive active level, H0/H1 separate, exact or
error-controlled squared-L2 integration on dimension-specific support, no
universal grid, no universal level cap, and streamed/chunked computation.

No PH, persistence landscape, distance, clustering, retrieval result, or
biological outcome was calculated in MV7-D.

## Operational deviation

The first Windows wrapper was given a ten-second client timeout. Its WSL parent
continued temporarily and atomically completed five raw entries, then ended
when the detached wrapper session closed. No scientific output was published.
The unchanged resume-safe runner was relaunched with a long-lived wrapper,
validated/reused existing atomic entries, completed all six raw and SCT jobs,
and published only after the complete gate passed. The final evidence and
independent validation are complete; no partial artifact entered Git.

## Decision

MV7-D closes the sample-count ambiguity and demonstrates that extending the
accepted upstream source/SCT path to the 34 retained single-study samples is
technically reasonable. It does not authorize the expansion itself.

The next sprint is MV7-E: resolve the 16 approach-label discrepancies as far as
public accession evidence permits and freeze the 124-sample descriptive cell
fit scope, gene-panel rule, pair scope, resource sentinels, and label firewall.
Only after that prospective contract passes should the 34 raw shards and 170
SCT caches be built.

## Verification

- MV7-D helper expectations: 9/9 pass.
- Prefreeze independent categories: 11/11 pass.
- Source/SCT feasibility independent categories: 9/9 pass.
- Complete package source-test suite: exit zero with two established optional
  skips and no failures or warnings.
