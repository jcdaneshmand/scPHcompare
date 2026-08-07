# Dual-View Topology Specification v1

## Document control

| Field | Value |
|---|---|
| Status | Frozen and enforced by MV-02 typed constructors/tests; ready for MV-03 feasibility pilot; not activated for confirmatory science |
| Date | 2026-08-05 |
| Contract IDs | `cell_topology_v1`, `gene_topology_v1` |
| Legacy ID | `legacy_gene_view_v0` — historical reproduction/stress use only |
| Landscape contract | `full_l2_error_controlled_v1` |
| Data scope | Existing large multi-tissue and bone-marrow cohorts |
| Approval boundary | Project-owner-directed technical freeze; broader author-team review required before a confirmatory full run |

## 1. Purpose and scientific hierarchy

This specification defines two different topological objects from the same sample expression matrix. It resolves the historical observation-orientation conflict without discarding the possibility that gene geometry is scientifically useful.

1. `cell_topology_v1` is the corrected confirmatory view. Cells are points, and topology describes the organization of a sample's cell population.
2. `gene_topology_v1` is a prespecified secondary view. A common gene panel supplies the points, and topology describes within-sample gene-expression relationships across matched cells.
3. `legacy_gene_view_v0` is not an implementation of `gene_topology_v1`. It passed feature-by-cell matrices directly to PH without the matching, scaling, metric, provenance, or interpretation controls below.

Cell and gene results must be interpretable independently before fusion. Fusion remains exploratory until both views pass MV-02 through MV-05.

## 2. Three distance layers

The term “distance” must always carry a layer identifier.

| Layer | Symbol | Object compared | Primary v1 choice |
|---|---|---|---|
| Within-sample geometry | `point_metric` | Cells with cells, or genes with genes, inside one sample | Cell: Euclidean in 30 shared PCs; gene: correlation-chord/Euclidean chord |
| Between-sample topology | `topology_metric` | Same-view, same-stratum persistence diagrams/landscapes | Separate all-level H0 and H1 landscape L2 |
| Sample analysis | `sample_method` | Samples represented by one distance matrix | PAM with label-free stability-selected `k` |

No code, artifact, table, or manuscript prose may call an object merely “the distance” when more than one layer is in scope.

## 3. Common source-matrix contract

For sample `s`, representation `q`, fit scope `f`, and subsample repetition `r`, let

`X[s,q,r]` be a numeric matrix with rows equal to the frozen gene panel and columns equal to the frozen cell subsample.

Required shape and identity:

- rows: exactly 500 unique human gene symbols in the frozen order;
- columns: exactly 384 unique cell IDs selected without replacement;
- the same cell IDs and gene order feed both views for a given `(s,q,r)`;
- all values finite;
- no duplicated or anonymous axis identifiers;
- source sample ID, cohort, assay/layer, representation, and input digest recorded;
- an axis-role field explicitly equal to `genes_by_cells`.

Representation-native continuous values are converted to a canonical standardized source matrix:

1. `sct_whole`: obtain SCT Pearson residuals for the frozen genes. If residuals are absent because `SCTransform()` retained only variable genes, calculate them from the recorded SCT model rather than substituting another slot silently.
2. `seurat_integration`/bone `integrated`: obtain the integrated assay's corrected feature values for the frozen genes and run the explicitly recorded centering/scaling operation required for PCA.
3. Within each cohort-by-representation-by-fit-scope stratum, estimate gene means and standard deviations from equal 384-cell contributions per fit sample. No sample contributes more cells to fitting than another.
4. Center and scale every retained gene with those fit-scope parameters. Reject nonfinite or effectively zero-variance features.

The pilot compares `sct_whole` with Seurat integration because that is the central dissertation contrast and both offer a cohort-shared representation. Raw, SCT Individual, and Harmony are deferred sensitivities. Harmony requires a separate native-output audit because current code materializes a corrected gene-by-cell assay rather than consuming a conventional low-dimensional Harmony embedding.

Large multi-tissue and bone-marrow cohorts are separate comparison strata. Representations are also separate strata. Diagrams or point distances are never pooled across cohorts or representations merely because they share a shape.

## 4. Frozen gene-panel contract

The pilot panel contains 500 genes and is selected once per cohort fit scope from the unintegrated `sct_whole` representation, then reused unchanged in the paired integrated stratum and both topology views.

Selection procedure:

1. Begin with genes present in every fit-scope sample and available in both primary representations.
2. Require finite, non-negligible variance in every fit-scope sample and representation.
3. For the primary panel, exclude the exact recorded set of mitochondrial, ribosomal-protein, and hemoglobin genes matching the implementation's versioned feature-category rules. Retain that excluded set for an explicit technical-feature sensitivity rather than discarding it from provenance.
4. Rank remaining genes by within-sample variance in `sct_whole`; aggregate ranks across equally weighted samples using median rank.
5. Select the first 500 by increasing median rank, breaking ties by canonical gene symbol.
6. Store the ordered panel, excluded candidates, per-sample ranks, fit-scope IDs, and SHA-256 digest.

Feature selection is label-free but still a fitted transformation. In descriptive MV-03 pilot work, the fit scope is all pilot samples and must be labeled `descriptive_all_pilot_samples`. In held-out or blocked evaluation, feature selection is fit on the training partition only and then applied unchanged to test samples.

If fewer than 500 eligible genes remain, the stratum fails preflight. It must not fill missing positions with zeros, sample-specific genes, or lower-ranked genes selected under a different rule.

Panel-size sensitivities are 250 and 1,000 genes, admitted only after the 500-gene Stage A/B pilot passes its resource gate.

## 5. Matched-cell subsampling contract

The primary pilot uses 384 post-QC cells per sample because the historical eligible range begins at 396 cells. This preserves all currently inventoried samples while leaving a 12-cell preflight margin in the smallest historical sample.

Rules:

- sample exactly 384 cells without replacement;
- use seeds `20260805`, `20260806`, `20260807`, `20260808`, and `20260809`;
- sort eligible cell IDs before applying the seeded sample operation;
- reuse exactly the same selected cells for cell and gene views and paired representations whenever the source object contains those cell IDs;
- store the ordered selected IDs and digest;
- never replace a failed sample with another sample after seeing topology or label outcomes;
- if a fresh rerun yields fewer than 384 eligible cells, record `ineligible_matched_cell_count` and stop for a design decision.

The five repetitions estimate sensitivity to which cells were observed. They are repeated measurements of one sample, not five independent biological samples.

Future cell-count sensitivities may use 256, 512, and 1,024 cells where eligibility permits. They are not part of the initial pilot core.

## 6. Cell topology: `cell_topology_v1`

Let `Z[s,q,r]` be the standardized 500-gene by 384-cell source matrix. Fit one PCA model per cohort-by-representation-by-repetition fit scope using the pooled, equally weighted cells from its fit samples. The PCA is computed in the conventional orientation: cells by genes (`rev.pca = FALSE`).

For the primary cell view:

- retain 30 variance-weighted cell PC scores;
- require exactly the same ordered loading coordinates for every sample in the stratum;
- represent each sample as a `384 x 30` cell-by-PC matrix;
- use ordinary Euclidean distance between cell PC scores;
- pass either the explicitly validated cell-by-PC matrix or its named `dist` object to PH;
- record the PCA seed, algorithm, center, scale, loadings, singular values, fit samples, genes, and model digest.

The primary geometry therefore describes cell-population structure in a shared linear coordinate system. It does not claim that PCA preserves every biological manifold or that a persistent H1 loop is automatically biological.

Prespecified sensitivities:

- 20 and 50 PCs using the same fitted decomposition;
- cosine-chord geometry after unit-normalizing nonzero 30-PC cell vectors, using Euclidean chord distance rather than unverified `1-cosine` dissimilarity.

Zero-norm cells, missing loadings, fewer than 30 usable PCs, or inconsistent coordinate names are hard failures.

## 7. Gene topology: `gene_topology_v1`

For each sample source matrix `Z[s,q,r]`, center every gene row over its 384 selected cells and divide it by its L2 norm. Call the resulting row vector `u_g`. Genes with a zero or nonfinite norm are preflight failures because the point set must remain the same across samples.

The primary gene metric is

`d_gene(g,h) = ||u_g - u_h||_2 = sqrt(2 * (1 - r[g,h]))`,

where `r[g,h]` is Pearson correlation over the matched cells. This is Euclidean chord distance between centered unit vectors, lies in `[0,2]`, and supplies a metric-space representation suitable for Vietoris-Rips PH.

For the primary gene view:

- represent every sample by the same 500 named gene points;
- coordinates are the 384 selected cells only for constructing within-sample gene relationships;
- calculate and pass an explicit named `dist` object to PH, avoiding matrix-orientation inference;
- require invariance to a joint permutation of cell columns;
- record row centering, norm, constant-gene disposition, cell IDs, and distance digest.

Cells are not aligned across different samples, nor do they need to be for this contract: each diagram summarizes the invariant within-sample geometry of the common gene set. Equal cell counts and repeated matched subsampling control the precision and sample-size dependence of that geometry.

The interpretation is gene coexpression-pattern topology, not a regulatory network, causal interaction graph, or cell-population topology.

Prespecified sensitivity:

`d_gene_rms(g,h) = sqrt(mean_c (Z[g,c] - Z[h,c])^2)`,

using the same standardized source matrix and selected cells. This is a scaled Euclidean metric that retains expression-magnitude differences suppressed by correlation-chord geometry.

The technical-feature-inclusion analysis and 250/1,000-gene panels remain labeled sensitivities.

## 8. Filtration and persistence-diagram contract

Primary pilot PH uses:

- Vietoris-Rips filtration;
- `ripserr` 0.3.0 in the locked R 4.4.1 environment;
- `max_dim = 1`, yielding H0 and H1 in one run;
- complete filtration with `threshold = -1`;
- field and engine defaults recorded explicitly in provenance;
- one essential H0 interval retained in diagram provenance but excluded from landscape evaluation under the landscape contract;
- invalid, zero-persistence, and nonfinite interval counts recorded.

Sample-specific finite thresholds are prohibited for corrected pairwise comparison. If complete H1 calculation exceeds the MV-03 resource envelope, the job is a failure, not an invitation to choose a favorable per-sample threshold. MV-03 may propose a separately versioned common-censoring sensitivity selected without labels and applied identically within a stratum, but it cannot silently replace `full_vr_v1`.

Absolute metric scale is primary within each stratum. Scale-normalized filtrations are separate sensitivities. A diagram from one representation is not directly pooled with a diagram from another representation.

## 9. Landscape and between-sample distances

Primary sample comparison follows `full_l2_error_controlled_v1`:

- all active consecutive landscape levels;
- exact or independently validated error-controlled integration;
- H0 distance reported separately;
- H1 distance reported separately;
- secondary combined distance `sqrt(d_H0^2 + d_H1^2)` accompanied by component values and H1 squared-distance contribution;
- no universal landscape-level cap or fixed uniform-grid primary calculation.

`landscape_reference_v1` remains the R correctness oracle. The corrected Persim critical-pair calculation is only a production prototype until it is batch-validated on eligible MV-03 diagrams.

Prespecified diagram-distance sensitivities:

1. bottleneck distance, with analytical fixtures and corrected-diagram provenance;
2. one Wasserstein-family implementation, preferably sliced Wasserstein if a mature isolated dependency trial passes correctness, licensing, packaging, and resource review.

The current generic curve “Wasserstein” functions are not substitutes for persistence-diagram Wasserstein distance.

## 10. Sample clustering and continuous tasks

The frozen primary clustering method is PAM/k-medoids because it accepts a general valid sample dissimilarity matrix without pretending the matrix itself is a Euclidean feature table.

Primary label-free cluster-count rule:

1. candidate `k` values are `2:min(10, n_samples - 1)`;
2. cluster every repeated-cell-subsample distance matrix at every candidate `k`;
3. calculate pairwise adjusted Rand agreement among repetitions;
4. select the largest mean stability, using the smallest `k` within one Monte Carlo standard error of the maximum;
5. if repeated inputs are unavailable or every solution is degenerate, report `no_stable_k` rather than importing the known class count.

Sensitivity clustering methods:

- average-linkage hierarchical clustering;
- spectral clustering after a frozen distance-to-affinity transformation and eigengap/stability audit.

`ward.D2` is excluded for arbitrary topological dissimilarities unless a valid Euclidean embedding is independently demonstrated. K-means is excluded from direct distance-matrix clustering.

An oracle `k` equal to a known tissue, study, or approach class count is retained only as an explicitly labeled historical-comparison sensitivity. Biological labels do not choose the primary `k`, metric, feature panel, PCA dimension, filtration, or fusion weight.

Continuous sample retrieval and blocked nearest-neighbor classification are required later so conclusions do not depend entirely on one clustering partition.

## 11. Pilot comparison strata and manifest

The frozen public candidate manifest is `docs/audits/mv01-pilot-sample-manifest-2026-08-05.csv`, generated by `scripts/build_mv01_pilot_manifest.R` from read-only historical metadata.

It contains 14 unique samples:

- eight large-cohort core samples: one per tissue, closest to the within-tissue median filtered-cell count;
- two large-cohort stress samples: global minimum and maximum filtered-cell count;
- four bone-marrow technical samples: deterministic lower/upper cell-count order-statistic selections within scRNA-seq and snRNA-seq, using one-based index `floor((n - 1) * p) + 1` for `p = 0.25, 0.75` after ordering by cell count and sample ID.

The manifest is a feasibility and contract-validation design, not a biologically representative sample or an analysis set selected to produce favorable clustering.

Primary representation pairs:

- large cohort: `sct_whole` versus `seurat_integration`;
- bone marrow: `sct_whole` versus `integrated`.

Stratum identity is:

`cohort × representation × view × homology_dimension × subsample_seed × fit_scope_id`.

Only samples sharing all required stratum fields may enter one sample-distance matrix.

## 12. Staged MV-03 execution and resource envelope

| Stage | Workload | Purpose | Advance condition |
|---|---|---|---|
| A | Large-cohort global min/max; `sct_whole`; one seed; both views | Orientation, H0/H1, and worst cell-count-ratio smoke test | Four jobs complete under per-job caps with valid provenance |
| B | Eight tissue-median plus four bone technical samples; both primary representations; first seed; both views | First eligible distance-matrix corpus | Complete artifact/identifier checks and estimated total cost inside stage cap |
| C | Four stability sentinels; both representations; five seeds; both views | Subsampling stability and production estimate | Stability report supports continuation or explicitly narrows a view |

Four Stage C sentinels are the two large-cohort cell-count stress samples plus the bone-marrow `technical_approach_q25` selections for scRNA-seq and snRNA-seq.

Resource caps:

- PH wall time: 30 minutes per sample-view-representation job;
- peak process-tree RSS: 8 GiB per worker;
- concurrency: at most two PH workers and 16 GiB aggregate sampled worker RSS;
- Stage A: at most two worker-hours;
- Stage B: at most 24 worker-hours;
- Stage C: at most 40 worker-hours;
- total retained pilot artifacts: at most 20 GiB before an explicit storage decision.

Exceeding a cap stops the stage. It does not authorize changing thresholds, removing failed samples, increasing parallelism, or weakening provenance. MV-03 must write a revised feasibility decision first.

## 13. Estimands and interpretation boundaries

Prospective estimands for later statistical planning:

1. **Biological distance contrast:** within-tissue versus between-tissue sample topology distance, with study blocking and uncertainty.
2. **Technical distance contrast:** same versus different study/approach topology distance, reported separately from biology.
3. **Integration change:** paired change in the biological and technical contrasts from unintegrated SCT Whole to integrated representation.
4. **View complementarity:** whether gene topology adds stable information beyond cell topology and matched non-topological sample baselines.
5. **H1 contribution:** whether H1 provides stable information beyond H0 rather than merely a nonzero numerical component.

The MV-03 pilot cannot estimate these biological quantities reliably; it only verifies eligibility, numerical behavior, stability, and cost.

Prohibited interpretations without additional evidence:

- an H1 feature is a biological cycle because it is present;
- integration caused a geometric mechanism because a distance changed;
- gene topology represents a regulatory or causal network;
- tissue retrieval proves removal of batch effects;
- fusion is superior because one weight performs best on the same labels used to choose it.

## 14. Failure hypotheses and prespecified dispositions

| ID | Failure hypothesis | Required diagnostic | Default disposition |
|---|---|---|---|
| F-01 | Cell topology is dominated by cell count | Matched-count versus count sensitivity; distance/count association | Keep matched design; reject uncontrolled result |
| F-02 | Gene topology is unstable to selected cells | Five-seed diagram/distance stability | Narrow gene claim, increase controlled sampling, or reject view |
| F-03 | Cell topology is unstable to PCA/subsampling | Seed and 20/30/50-PC sensitivity | Revise coordinate contract before full run |
| F-04 | H1 is empty, censored, or unstable | H1 counts, persistence, threshold, and repeat stability | Report H0 primary; retain H1 as null/unsupported |
| F-05 | Representation scales prevent paired interpretation | Fit/scale diagnostics and normalized sensitivity | Compare within strata only; prohibit raw cross-scale claims |
| F-06 | Tissue signal is study/technology signal | Study-blocked and approach-stratified evaluation | Narrow or withdraw biological claim |
| F-07 | Gene signal is driven by technical gene categories | Excluded/included technical-feature sensitivity | Label technical diagnostic or revise panel |
| F-08 | Complete H1 PH is computationally infeasible | Stage time/RSS/failure scaling | Evaluate approved approximation/common censoring; do not silently truncate |
| F-09 | Cell and gene views are redundant | Distance/neighborhood/component comparison | Report stronger single view; do not force fusion |
| F-10 | Fusion gains are tuning leakage or one-study effects | Nested/blocked evaluation and full weight grid | Reject learned fusion claim |
| F-11 | Common 500-gene panel cannot be formed | Per-sample/representation availability report | Stop affected stratum; no zero padding |
| F-12 | Historical metadata no longer match rerun QC | Manifest reconciliation | Regenerate pilot manifest before PH |

Null, negative, and failed results remain in the audit bundle.

## 15. Provenance and cache identity

Every corrected view object and downstream artifact must record at least:

- contract ID and version;
- cohort, representation, sample ID, and stratum ID;
- source object/file digest and assay/layer;
- axis roles, dimensions, ordered gene digest, ordered cell digest;
- QC and exclusion policy;
- subsample count, seed, and selected IDs;
- fit-scope ID and train/test role;
- feature-selection rule and excluded technical genes;
- normalization/standardization parameters and digest;
- PCA algorithm, seed, dimensions, loadings/model digest where applicable;
- point metric ID and explicit distance-object digest;
- PH engine/version, maximum homology dimension, field, filtration policy, threshold, runtime, peak process-tree RSS, and failure state;
- diagram digest and finite/infinite/invalid interval counts by H0/H1;
- topology-distance engine/version, specification, numerical tolerance/error, and matrix digest;
- clustering method/version, candidate `k`, selection rule, seeds, and assignments;
- creation revision and software/session manifest.

The cache key is a SHA-256 digest of a canonical serialization of scientific inputs and parameters. Paths, timestamps, and display labels are not sufficient cache identities. Any change to view, axis roles, cells, genes, fit scope, representation, metric, filtration, homology, engine, or numerical policy must invalidate the cache.

## 16. Leakage rules

- Tissue, study, approach, and manuscript outcome labels are never inputs to feature selection, PCA fitting, metric choice, `k` selection, or fusion weights.
- A descriptive all-sample fit is allowed only for MV-03 feasibility and must be labeled transductive/descriptive.
- In blocked evaluation, every fitted transformation uses training samples only. Test samples are projected or transformed with frozen training parameters.
- Hyperparameter selection is nested inside the training partition.
- Labels may be opened only after immutable primary distance matrices and cluster assignments are written for the evaluation split.
- The pilot sample-selection rules use tissue/approach only to ensure coverage for feasibility; the manifest is prohibited from confirmatory biological effect estimation.

## 17. MV-01 gate decision

G-MV1 passes as a technical contract freeze for MV-02/MV-03 because:

- both axes and scientific objects are explicit;
- primary within-sample metrics are mathematical metrics;
- equal cells and a common gene panel isolate the view difference;
- primary representations and strata are bounded;
- complete filtration and H0/H1 handling are explicit;
- the method shortlist is frozen before corrected biological results exist;
- provenance, leakage, resource, failure, and stop rules are prespecified.

This pass does not approve a confirmatory full run or biological claims. Broader author-team review remains required at G-MV7 and before manuscript activation.

## 18. References and evidence

Implementation evidence: `R/dual_view_topology.R`, `tests/testthat/test-dual-view-topology.R`, and `docs/audits/MV02_ORIENTATION_SAFE_CONSTRUCTORS_2026-08-05.md`. The implementation provides a strict `scientific` profile and a reduced `analytical_fixture` profile that is always marked scientifically ineligible. The corrected PH entry point refuses bare matrices and bare `dist` objects.

- Local orientation and engine audit: `docs/audits/LANDSCAPE_ORACLE_AND_DIAGRAM_ELIGIBILITY_2026-08-05.md`.
- Dissertation cell-by-cell intent: `docs/Dissertation_SubmissionReady_October.pdf`, PDF page 54 (printed page 40); local and Git-ignored.
- ripserr row-as-point contract: <https://cran.r-project.org/web/packages/ripserr/refman/ripserr.html>.
- Seurat PCA orientation and zero-variance behavior: <https://satijalab.org/seurat/reference/runpca>.
- Seurat SCTransform output contract: <https://satijalab.org/seurat/reference/sctransform>.
- Chazal et al., subsampling for PH: <https://proceedings.mlr.press/v37/chazal15.html>.
- Bubenik, persistence landscapes: <https://www.jmlr.org/papers/v16/bubenik15a.html>.
- Skinnider et al., association measures for single-cell transcriptomics: <https://www.nature.com/articles/s41592-019-0372-4>.
- Luecken et al., integration benchmarking and biology/batch separation: <https://www.nature.com/articles/s41592-021-01336-8>.
- Current-data/code evidence and the exact frozen method table: `docs/audits/MV01_DUAL_VIEW_CONTRACT_2026-08-05.md` and `docs/audits/mv01-method-freeze-2026-08-05.csv`.
