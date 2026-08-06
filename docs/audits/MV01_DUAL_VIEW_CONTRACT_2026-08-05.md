# MV-01 dual-view scientific-contract audit — 2026-08-05

## Outcome

MV-01 freezes a bounded technical contract for implementing and piloting two different topological views:

- `cell_topology_v1`: cells as points in a shared 30-PC coordinate system;
- `gene_topology_v1`: the same 500 genes as points under correlation-chord geometry across the same 384 selected cells.

The same gene panel and cell subsample feed both views. This makes the view orientation scientifically interpretable instead of confounding it with different input genes, cells, or sample sizes.

G-MV1 passes for MV-02 implementation and MV-03 feasibility work. Nothing in this sprint activates corrected production PH, authorizes a full biological rerun, or validates the views' biological performance.

## Evidence inspected

### Dissertation and existing audits

- The dissertation describes a cell-by-cell Euclidean distance matrix for each sample (local Git-ignored PDF page 54; printed page 40).
- The orientation audit demonstrated that current and historical feature-by-cell assay matrices were passed directly to a row-as-point PH API.
- The landscape audits approved all-level exact/error-controlled H0/H1 landscape L2 as the corrected target and invalidated fixed-grid/level-cap reuse.
- Historical filtered-cell inventory contains 149 samples across the two cohorts, ranging from 396 to 11,475 cells; large-cohort tissue/study/approach imbalance remains substantial.

### Actual current code

- `R/PH_Calculation.R` extracts RNA, SCT, integrated, and Harmony assay matrices in feature-by-cell orientation.
- `R/PH_Functions.R::process_and_monitor()` serializes the matrix and passes it without an axis contract to `ripserr::vietoris_rips()`.
- Current retry logic can derive a sample-specific threshold from a PCA/kNN calculation; that behavior is prohibited for corrected cross-sample primary comparisons.
- Current code constructs `SCT_Whole`, `SCT_Individual`, Seurat-integrated, and Harmony-derived matrices. The pilot limits itself to SCT Whole and Seurat integration.
- Current Harmony code materializes a corrected gene-by-cell assay and checks its orientation. This is not automatically equivalent to using the method's native integrated embedding and requires a later audit.
- Current sample clustering applies multiple hierarchical linkages and spectral clustering to PH matrices, while conventional k-means is applied to cell embeddings. The new contract freezes PAM for arbitrary sample dissimilarities and excludes direct k-means/unsupported Ward use.

### Locked environment

| Component | Version relevant to MV-01 |
|---|---:|
| R | 4.4.1 |
| Seurat | 5.3.0 |
| SeuratObject | 5.0.2 |
| sctransform | 0.4.1 |
| harmony | 1.2.4 |
| ripserr | 0.3.0 |
| TDA | 1.9.4 |
| cluster | 2.1.8.1 |
| mclust | 6.1.1 |

The tracked lockfile remains authoritative. MV-01 installed or upgraded no dependencies.

## Key decisions and rationale

### Same 500 genes and 384 cells in both views

Using the same source matrix isolates the scientific effect of treating rows versus columns as points. A 500-gene pilot is large enough to exercise multi-level gene landscapes but remains below the historical 3,000–tens-of-thousands feature-point workloads that produced extreme depth and memory. It is a feasibility setting, not a final biological optimum.

The matched cell count is 384 because every currently inventoried post-QC sample has at least 396 cells. Equal cells prevent larger samples from receiving more coordinate dimensions in the gene view or denser point sampling in the cell view. Five frozen seeds make cell-sampling instability observable.

### Cell geometry: shared PCA with Euclidean distance

The dissertation intended cells as observations, but direct Euclidean geometry in thousands of expression features is vulnerable to sparsity, scaling, and compute. The contract therefore fits a shared 30-PC basis within each cohort/representation/repetition and projects every sample into identical coordinates. Seurat's `RunPCA()` documentation confirms that its conventional orientation is cells by genes and that zero-variance features are removed. Twenty and fifty PCs are prespecified sensitivities.

The pilot uses equal cells per fit sample so large samples cannot dominate the basis. A descriptive all-pilot fit is allowed only for feasibility; blocked evaluation requires training-only fitting.

### Gene geometry: correlation-chord metric

The accidental gene orientation was not accompanied by an intentional gene metric. The new view centers and unit-normalizes each gene's 384-cell expression vector, then uses Euclidean chord distance. It equals `sqrt(2*(1-r))`, is bounded, and is a metric because it is ordinary Euclidean distance after a specified embedding.

This represents gene coexpression-pattern geometry. Single-cell association measures are sensitive to the data-generating and preprocessing properties, so normalized RMS Euclidean and technical-feature inclusion are mandatory sensitivities rather than a large unstructured metric search.

### Feature selection

The panel is derived from unintegrated SCT Whole data using median within-sample variance rank, with samples weighted equally. It is intersected across samples and the paired integrated representation, then reused unchanged. Selection in the unintegrated reference avoids allowing integration output to define which features are later judged to have changed.

Mitochondrial, ribosomal-protein, and hemoglobin genes are excluded from the primary panel because they are already explicit technical/QC categories in the project. Their inclusion is retained as a prespecified sensitivity so a technical gene-topology signal is not concealed.

### Representation strata

The first comparison is SCT Whole versus Seurat integration, separately in the large and bone-marrow cohorts. This directly addresses the dissertation's central preprocessing contrast while bounding computation. Raw counts, SCT Individual, and Harmony are deferred sensitivities rather than multiplying the initial design.

PCA/standardization models are representation-specific, so absolute point-cloud or diagram scales are not compared naively across representations. Sample-distance structures and prespecified within-stratum endpoints are compared instead.

### Complete filtration and separate H0/H1

Primary MV-03 PH uses complete Vietoris-Rips filtration (`threshold=-1`) with H0/H1 calculated together and reported separately. The prior sample-specific retry threshold is incompatible with a clean primary cross-sample interpretation. If complete H1 is infeasible, MV-03 must stop and make a separate common-censoring/approximation decision.

The dissertation-aligned all-level landscape L2 remains primary. Bottleneck and one Wasserstein-family method are sensitivities; generic Wasserstein calculations on Betti/Euler/landscape curve arrays are not persistence-diagram Wasserstein distances.

### PAM and label-free `k`

PAM is frozen as the primary sample clustering method because it accepts arbitrary sample dissimilarities and returns sample medoids. Candidate `k` is selected from repeated-subsample stability without biological labels. Average linkage and a separately specified spectral affinity are sensitivities. Oracle `k` is historical sensitivity only.

Direct k-means on a sample-distance matrix is excluded. Ward linkage is excluded unless a valid Euclidean embedding is demonstrated.

## Pilot manifest

`scripts/build_mv01_pilot_manifest.R` read the two historical metadata tables without modifying them and generated `docs/audits/mv01-pilot-sample-manifest-2026-08-05.csv`.

Two consecutive locked-environment builds produced the same manifest SHA-256: `435e9f9f2887bc2d91250555193ed79e2084afe62527282c39d60571dceb4495`.

Source identities:

| Metadata source | Rows | SHA-256 |
|---|---:|---|
| `joined_metadata_cellcounts.csv` | 124 | `e2b6f1e7d7dc05ba8e80627aad4ea1677c3b718a844397acd06596d7e22a18d0` |
| `joined_metadata_cellcounts_bonemarrow.csv` | 25 | `602fb2baee73ab7a1042f7d361e8a6748820f6fcfa48244b941269ff1caa8015` |

The 14-sample manifest contains:

- eight large-cohort samples closest to their tissue-specific median cell count;
- the global 396-cell and 11,475-cell stress samples;
- deterministic lower/upper cell-count order-statistic samples within each bone-marrow assay approach, using one-based index `floor((n - 1) * p) + 1` for `p = 0.25, 0.75` after ordering by cell count and sample ID.

Every selected sample can supply 384 historical post-QC cells. Fresh MV-03 QC must reconcile these counts; the manifest does not override new provenance.

## Staged feasibility boundary

MV-03 begins with the two large-cohort cell-count extremes in SCT Whole, one seed, and both views. It advances to the 12 core/technical samples and paired representations only when every Stage A job completes within 30 minutes and 8 GiB sampled process-tree RSS.

At most two PH workers may run concurrently. Stage caps are 24 worker-hours for the first multi-sample corpus, 40 worker-hours for five-seed stability, and 20 GiB retained artifacts. These are stop thresholds, not optimization targets.

## Prospective scientific questions

The contract supports later tests of:

- whether cell topology captures biological identity beyond sample-level non-topological baselines;
- whether gene topology contains independent or merely technical structure;
- whether integration reduces study/approach structure while retaining biology;
- whether H1 adds stable information beyond H0;
- whether cell/gene fusion is complementary, redundant, or leakage-prone.

The pilot is prohibited from answering those questions. Its role is eligibility, correctness, numerical stability, and feasibility.

## Failure and honesty boundary

Twelve failure hypotheses are frozen in the specification. The important scientific stop cases are:

- either view remains unstable after matched repeated subsampling;
- cell topology changes mainly with cell count or PCA choice;
- gene topology changes mainly with technical-feature inclusion;
- H1 is empty/unstable or requires incomparable thresholds;
- study/technology explains the apparent tissue structure;
- fusion merely averages redundant matrices or improves only through label tuning;
- a common 500-gene panel cannot be formed without zero padding;
- resource use exceeds the staged envelope.

Any of those outcomes narrows, revises, or rejects a view. None authorizes tuning around a negative result.

## MV-01 gate disposition

| Question | Disposition |
|---|---|
| Scientific contracts coherent? | `cell_topology_v1` and `gene_topology_v1` technically frozen |
| Observation axes explicit? | Pass |
| Primary point metrics valid? | Pass by construction; implementation tests still required in MV-02 |
| Cell/gene matching specified? | Pass: same 500 genes and 384 cells |
| Filtration/H0/H1 specified? | Pass |
| Pilot bounded and deterministic? | Pass: public 14-sample manifest and staged resource gate |
| Biological claims permitted? | No; feasibility only |
| Production change authorized? | No |
| Immediate next action | MV-02 orientation-safe constructors and fixtures |

## References

- Seurat `RunPCA`: <https://satijalab.org/seurat/reference/runpca>
- Seurat `SCTransform`: <https://satijalab.org/seurat/reference/sctransform>
- Chazal et al., *Subsampling Methods for Persistent Homology*: <https://proceedings.mlr.press/v37/chazal15.html>
- Bubenik, *Statistical Topological Data Analysis using Persistence Landscapes*: <https://www.jmlr.org/papers/v16/bubenik15a.html>
- Skinnider et al., *Evaluating measures of association for single-cell transcriptomics*: <https://www.nature.com/articles/s41592-019-0372-4>
- Luecken et al., *Benchmarking atlas-level data integration in single-cell genomics*: <https://www.nature.com/articles/s41592-021-01336-8>
