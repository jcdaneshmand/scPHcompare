# MV7-E metadata provenance and full-124 descriptive-estimand prefreeze audit

Date: 2026-08-15  
Exact prefreeze commit: `4e1feea0f6bea626a44ffc39911ed8deb9c844f4`  
Status: complete; MV7-F upstream production only authorized

## Outcome

MV7-E resolves every sequencing-approach discrepancy, corrects the scope of an
earlier approach-only diagnostic, and prospectively freezes the separate
124-sample descriptive analysis. It computes no PCA, PH, landscapes,
distances, clustering, or biological outcomes.

The resulting continuation is:

`MV7-F upstream caches -> exact 124-derived panel lock -> MV7-G PH sentinel`

Full PH remains closed.

## Why the approach fields disagreed

The historical `joined_metadata_cellcounts.csv` retained two joined columns:

| Field | Origin | Agreement with public metadata | Scientific disposition |
|---|---|---:|---|
| `Approach.x` | Seurat metadata assigned by an expression heuristic | 108/124 | Prohibited |
| `Approach.y` | Public candidate metadata | 124/124 | Authoritative |
| canonical | Public value plus official accession confirmation where disputed | 124/124 | Authoritative |

The historical heuristic classified technology from median detected genes,
mitochondrial fraction, and UMI count. Those expression properties cannot
establish whether cells or nuclei were experimentally isolated.

Official GEO method records resolve all 16 conflicts in favor of the public
value: 14 are explicitly single-cell and two are explicitly nuclei extraction
with DAPI singlet sorting. The checked-in accession table contains every SRS,
SRR, GSM, GSE, resolution, evidence summary, access date, and official URL.
Representative primary records include [PBMC GSE96583](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96583),
[PBMC GSE106543](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE106543),
[testis GSE109037](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109037),
[pancreatic islets GSE114297](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114297),
[prostate GSE117403](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117403),
[bone marrow GSE120221](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120221),
and [substantia nigra GSE126836](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126836).

## Correction to MV7-B

The corrected 90-sample primary axis contains 90 scRNA-seq and zero snRNA-seq
samples. It therefore contains zero mixed-approach studies. MV7-B's
mixed-study approach association is not estimable and its earlier numeric
result is superseded as a metadata-field error.

This does not require recalculating topology or landscapes. MV7-B's study,
tissue, retained-cell, deletion-influence, and cell-versus-gene results are
unaffected. No technology claim had been authorized. In the full 124 samples,
118 are scRNA-seq and six are snRNA-seq; the six nuclei samples are nested in
one substantia-nigra study/tissue, so a causal or separable technology effect
is still not identifiable.

## Panel-availability gate

The accepted 90-derived MV6-C panel contains 500 features with semantic SHA-256
`7be22cdb9056fed427c78d58be2b19e258c7c6807e6f7ac3900bd1bfa1d19eb8`.
The audit inspected only source feature axes after reproducing Seurat's
historical underscore-to-dash normalization; expression values were not used
or exported, and topology and labels were not read.

| Result | Count |
|---|---:|
| Added samples checked | 34 |
| Complete for all 500 features | 33 |
| Incomplete | 1 |
| Missing-feature occurrences | 1 |

`SRA701877_SRS3279688` contains 499/500 and lacks
`KLF2-ENSG00000127528.5`. This single availability result activates the
already specified fallback. After MV7-F, the unchanged MV6-C global-core
algorithm will be applied to all 620 SCT caches. The exact 124-derived panel
must then be committed, compared with the 90-derived panel, and independently
validated before any PCA or PH.

## Frozen descriptive estimand

- Exact population: 124 retained samples; the 90-sample LOSO benchmark remains
  primary and immutable.
- Seeds: `20260805` through `20260809`.
- Depth: 384 matched cells per sample and seed; 620 sample-seed states.
- Fit: one transductive global center/scale and 30-PC PCA per seed over all
  124 equal-depth samples (47,616 cells per fit).
- Cell view: 384 cells as points in 30 PCs with Euclidean geometry.
- Gene view: 500 genes as points across the same 384 cells with
  correlation-chord geometry.
- PH: complete Vietoris-Rips H0/H1; essential H0 excluded.
- Pairs: 7,626 unordered pairs per seed, 38,130 per view, 76,260 across views,
  and 152,520 separate H0/H1 component rows; no cross-seed pairs.

The global descriptive transform is not the same estimand as training-only
LOSO PCA and may not be reported as such.

## Landscape definition

The dissertation-aligned corrected definition remains unchanged:

- all finite positive-persistence intervals;
- essential H0 excluded;
- every consecutive active landscape level;
- H0 and H1 calculated and reported separately;
- exact or error-controlled squared-L2 integration on dimension-specific
  support;
- no universal grid and no universal level cap; and
- streamed/chunked computation without dense landscape matrices.

## Resource and stage gate

MV7-F is authorized for 34 atomic raw shards and 170 SCT caches only. Each
child has an 1,800-second and 8-GiB cap; the aggregate has a four-hour
worker-time and 4-GiB storage gate. Existing valid entries are immutable,
resume is identity-checked, and partial state is not publishable.

MV7-F's measured result must set the MV7-G cap. MV7-G's measured result must set
the MV7-H cap. Tissue, study, approach, and outcomes remain outside panel,
transform, topology, PH, and landscape selection.

## Validation and determinism

- Focused MV7-E helper expectations: 18/18 pass.
- Acceptance criteria: 16/16 pass.
- Independent categories: 17/17 pass, including reopening all 34 sources.
- Clean exact-commit repeat: 17/17 builder artifacts byte-identical.
- Complete R source test suite: exit zero with two established optional skips.
- PCA, PH, landscape, clustering, and outcome jobs: zero.

During development, the first repeat-validator version failed to parse because
it used R's reserved word `repeat` as a variable. That evidence root was
discarded. The validator was corrected, every MV7-E script was parse-checked,
and public source locators were normalized to repository-relative paths. The
complete primary/repeat/independent sequence was regenerated at the exact
commit above. No scientific value changed.

## Decision

MV7-E passes. Proceed to MV7-F upstream production, then stop calculation long
enough to freeze and independently validate the exact 124-derived panel. Do
not run PCA or PH until that panel gate passes.
