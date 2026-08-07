# MV-05 statistical and fair single-view benchmark plan v1

## Document control

| Field | Value |
|---|---|
| Contract ID | `mv05_statistical_benchmark_contract_v1` |
| Date frozen | 2026-08-06 |
| Status | Frozen for implementation review; biological execution not authorized |
| Existing-data scope | 124-sample large cohort plus 25-sample bone-marrow cohort |
| Confirmatory view | Cell topology, H0 and H1 separately |
| Secondary view | Gene topology, H0 and H1 separately |
| Fusion | Prohibited until G-MV5 passes |
| Pilot role | Technical smoke testing only |
| New-data policy | Existing data first; no trigger yet |

## 1. Purpose and decisive boundary

MV-05 asks whether a correctly defined sample-level topological representation
adds stable information beyond matched sample-level non-topological
representations. It does not ask whether copied sample labels score well when
repeated over thousands of cells.

The immutable MV-04 pilot bundle is numerically valid but cannot estimate the
confirmatory biological endpoints in this plan:

1. its 14 samples were selected for feasibility rather than biological
   representativeness;
2. no pilot tissue spans two studies;
3. its fitted transformations are explicitly
   `descriptive_all_pilot_samples`; and
4. the current Seurat integration code jointly fits all supplied samples and
   has no audited held-out query mapping path.

MV-04 artifacts may therefore verify API behavior, identifiers, matrix
handling, and resource estimates. They may not be presented as held-out tissue
validation or used to choose a winning method.

## 2. Scientific questions and hierarchy

The frozen questions are:

1. Does cell topology preserve cross-study tissue similarity better than a
   distributional baseline built from the same cells and shared coordinates?
2. Does integration change biological conservation and residual technical
   structure, reported as two separate axes rather than one composite score?
3. Does H1 add stable information beyond H0, rather than merely being nonzero?
4. Does the intentionally defined gene-topology view carry reproducible
   secondary information beyond a matched gene-correlation baseline?
5. Are label-free sample clusters stable across cell subsamples, and do their
   post hoc biological and technical alignments tell the same or conflicting
   stories?

Cell topology is confirmatory. Gene topology is secondary. No result from the
gene view can rescue a failed cell-view primary endpoint, and the stronger of
H0/H1 is not selected after outcomes are seen.

## 3. Experimental and inferential units

| Object | Unit | Permitted use |
|---|---|---|
| Cell | Subsample within a biological sample | Construct one sample representation |
| Gene | Common named feature | Construct the gene view or matched gene baseline |
| Sample/donor preparation | Primary observation | Distance, retrieval, prediction, clustering |
| Study/accession | Dependence block and outer fold | Generalization and uncertainty |
| Cell-subsample seed | Repeated technical realization | Stability/technical uncertainty, never independent n |
| Sample pair | Dependent derived observation | Descriptive distance contrasts; aggregate by study pair before inference |

Every inferential table contains one row per sample, study summary, or study
pair as declared. Sample assignments are never copied to cells for primary
metrics. Cell counts never weight sample-level accuracy or external-cluster
agreement.

This boundary follows evidence that cells from one individual are correlated
subsamples rather than independent experimental units and that sample-level
aggregation is a defensible comparison unit in multi-sample single-cell work.

## 4. Existing-data eligibility

The large cohort contains 124 samples, 18 studies, eight tissues, and two
approaches. Five tissues have at least two studies and form the current
cross-study candidate subset:

- bone marrow: 31 samples, 3 studies;
- PBMC: 12 samples, 4 studies;
- liver: 6 samples, 2 studies;
- testis: 28 samples, 4 studies;
- colon: 13 samples, 2 studies.

Together they provide 90 candidate samples across 15 studies. Pancreatic
islets, substantia nigra, and prostate are each represented by one study and
are descriptive-only for tissue endpoints. They remain in failure/confounding
tables and may support technical diagnostics; they are not silently discarded
from the project.

The bone-marrow cohort contains 25 samples from one study and one tissue (7
scRNA-seq and 18 snRNA-seq). It cannot validate biological tissue recovery or
cross-study generalization. It is a technical approach/distinctness task only.

Eligibility is frozen before outcome calculation. A tissue enters the primary
cross-study endpoint only when it occurs in at least two studies and every
evaluated held-out sample has at least one training-study sample of that tissue.

## 5. Representation and baseline panel

### 5.1 Topological methods

- `cell_landscape_h0_v1`: all-level exact/error-controlled landscape L2, H0.
- `cell_landscape_h1_v1`: all-level exact/error-controlled landscape L2, H1.
- `gene_landscape_h0_v1`: same summary contract, secondary gene view, H0.
- `gene_landscape_h1_v1`: same summary contract, secondary gene view, H1.

The raw H0/H1 Euclidean combination is descriptive only because MV-04 showed
strong H0 dominance. H0 and H1 are never collapsed before their separate
results and uncertainty are written.

### 5.2 Matched non-topological baselines

1. `cell_distribution_energy_shared_pca_v1` is the matched cell-view baseline.
   It uses the same 384 cells and training-fit shared 30-PC coordinates and
   computes the square root of the empirical V-statistic energy divergence,
   `sqrt(max(0, 2 mean d(X,Y) - mean d(X,X) - mean d(Y,Y)))`, including the
   zero self-pairs in both within-sample means. It preserves distributional
   information without persistent homology.
2. `pseudobulk_shared_panel_euclidean_v1` is a context baseline. It computes
   Euclidean distance between sample means over the same 500-gene panel after
   training-fit scaling. It is sample-level and equally weighted by sample.
3. `gene_correlation_frobenius_v1` is the matched gene-view baseline. It uses
   the same 500 genes and matched cells and computes RMS Frobenius distance
   between sample gene-correlation matrices over all `500^2` entries.
4. `composition_hellinger_conditional_v1` is a sensitivity only if a harmonized,
   externally frozen cell-label vocabulary exists for every included study.
   Cell labels may not be inferred using tissue or test-fold information.

Pseudobulk and cell-centroid distances are algebraically close under the same
linear PCA model; they are not counted as two independent wins. The centroid
may be retained as an equivalence diagnostic, not as another headline method.

Cell-level Seurat/Louvain results are a separate biological analysis and are
excluded from the primary same-unit benchmark. K-means is prohibited on a
distance matrix. Ward linkage is prohibited unless Euclidean embeddability is
demonstrated independently.

## 6. Outer validation and fit/transform rules

### 6.1 Primary outer split

The primary split is `large_leave_one_study_out_v1`. Each of the 18 large-cohort
studies is held out once. The five cross-study tissues determine endpoint
eligibility, but tissue labels do not determine fold assignment.

For every outer fold:

1. Freeze training and held-out sample IDs by study.
2. Fit feature eligibility, feature ranking, means, variances, PCA, distance
   scales, and any integration/reference model using training samples only.
3. Apply the frozen transformation to held-out samples without refitting.
4. Generate the same five cell subsamples (`20260805:20260809`) in every method.
5. Write and hash view objects, diagrams, distances, baseline matrices, and
   predictions before tissue/approach outcomes are opened.
6. Evaluate only held-out samples whose tissue is represented in another
   training study; record all other held-out samples as
   `single_study_tissue_not_estimable`.

No sample from one study may appear on both sides of a fold. If future methods
add tunable hyperparameters, tuning occurs in study-blocked inner folds inside
the outer training set. The current 1-NN retrieval task and method registry are
fixed and require no outcome-driven tuning.

### 6.2 Integration induction gate

The current `perform_integration()` route calls `FindIntegrationAnchors()` and
`IntegrateData()` jointly and remains transductive. MV5-B added a separate
`mv05_inductive_mapping_v1` route using reference-only feature selection/PCA,
`FindTransferAnchors(reduction="pcaproject")`, and label-free
`IntegrateEmbeddings()`. Its synthetic two-study fixture passed twice with
identical embeddings, but full-data fold feasibility is not yet established.
Consequently:

- current all-sample integrated objects support only explicitly transductive,
  descriptive integration diagnostics;
- a confirmatory integrated LOSO result must use the new
  training-reference/held-out-query contract and pass MV5-C on an eligible
  existing-data tissue; and
- if a held-out mapping cannot reproduce the approved feature/coordinate
  contract, the integrated confirmatory endpoint is `not_executable`, not
  replaced by the transductive result.

This integration induction gate must be solved before biological MV-05
execution. It does not require new data.

## 7. Label firewall

| Stage | Allowed labels | Prohibited use |
|---|---|---|
| Eligibility audit | Study, tissue, approach | Outcome calculation or method ranking |
| Outer split | Study only | Tissue-balanced outcome-informed folds |
| Feature/PCA/integration fitting | None | Any biological or technical label |
| Distance scaling | None | Label-dependent component weights |
| Method selection | Prespecified registry only | Choosing method from outcomes |
| Primary `k` selection | Repeated distance stability only | Known class count |
| Prediction writing | None | Opening labels before hashes are written |
| Biological evaluation | Tissue | Refitting any upstream object |
| Technical evaluation | Study and approach | Refitting any upstream object |
| Oracle-`k` sensitivity | Tissue/study/approach | Promotion to primary result |

The software label-use registry and validator enforce these boundaries.

## 8. Frozen endpoints

### 8.1 Biological conservation

The primary continuous endpoint is cross-study tissue retrieval macro mean
reciprocal rank (`cross_study_tissue_mrr_v1`). For each held-out sample, rank
training samples by distance after excluding its study. Record the reciprocal
rank of the first sample with the same tissue. Average first within tissue and
then equally across eligible tissues, so bone marrow and testis do not dominate
smaller classes.

Supportive endpoints are:

- fixed 1-nearest-neighbor tissue balanced accuracy;
- between-tissue minus within-tissue cross-study distance, after reducing
  sample pairs to study-pair summaries; and
- tissue-specific results with study and sample counts, never winner-only
  averages.

Ties use canonical sample ID order and are reported. A held-out tissue absent
from training is not scored as an error; it is a prespecified non-estimable
case.

### 8.2 Technical signal

Technical structure is reported separately:

- cross-study minus same-study distance within tissue (lower residual study
  separation is preferred after integration);
- different-approach minus same-approach distance within study and tissue where
  identifiable; and
- study/approach retrieval as explicitly technical diagnostics.

Biological conservation and technical removal are not combined into a weighted
overall score. Integration can improve one and harm the other.

### 8.3 Integration change

For each fold, seed, method, view, and dimension, compute paired integrated
minus SCT changes in biological and technical endpoints. The central
integration interpretation uses both coordinates:

| Biological change | Technical change | Disposition |
|---|---|---|
| improves | improves | favorable within stated endpoint limits |
| improves | worsens | biology/technical trade-off |
| worsens | improves | probable overcorrection or biological loss |
| worsens | worsens | unfavorable |
| uncertain | any | uncertain; do not rank |

No causal geometric mechanism is inferred from a distance change alone.

## 9. Clustering panel and cluster count

Clustering is secondary to continuous retrieval.

1. Primary: PAM/k-medoids on each valid sample-distance matrix.
2. Sensitivity: average-linkage hierarchical clustering.
3. Sensitivity after implementation gate: Gaussian-affinity normalized spectral
   clustering with bandwidth equal to the median positive off-diagonal distance,
   zero graph diagonal, normalized Laplacian, deterministic eigensolver and
   clustering seed, and eigengap diagnostics.

The current `perform_spectral_clustering()` helper is not eligible because its
seed, eigengap, and artifact provenance are not frozen.

The primary candidate set is `k=2:min(10,n-1)`. Every method/distance must have
all five subsample-seed matrices with identical sample axes. For each `k`,
compute all pairwise adjusted Rand agreements across seeds. Estimate Monte
Carlo SE by leaving out one seed at a time. Select the smallest `k` whose mean
stability is within one SE of the largest mean. Missing seeds, inconsistent
sample axes, or wholly degenerate partitions produce `no_stable_k`.

Oracle `k` equal to the known number of tissue, study, or approach classes is a
historical sensitivity only. External agreement uses one row per sample; ARI is
the primary external-cluster metric and NMI is supportive. Purity, Jaccard, VI,
and cell-replicated metrics do not form a confirmatory multiplicity family.

## 10. Method comparisons and multiplicity

Three families are frozen:

1. `F1_primary_retrieval`: cell H0/H1 versus the matched cell-energy baseline
   within SCT and integrated representations.
2. `F2_primary_technical`: the corresponding residual-study contrasts.
3. `F3_integration_change`: cell H0/H1 paired integrated-minus-SCT biological
   and technical changes.

Holm adjustment is applied within each family. Gene-view, pseudobulk-context,
clustering, per-tissue, composition, oracle-`k`, and technical bone-marrow
results are labeled secondary/sensitivity and reported completely without
promotion based on nominal significance.

A statement that topology outperforms matched conventional structure requires
the directionally favorable effect and interval against the matched baseline;
being better than one weak baseline while worse than the matched distributional
baseline is not success.

## 11. Uncertainty and randomization

- Use all five fixed cell-subsample seeds; they are repeated technical
  realizations, not five independent samples.
- Primary 95% intervals use 2,000 paired, tissue-stratified study-block
  bootstrap replicates. All samples from a resampled study travel together.
- Pairwise-distance endpoints are reduced to study-pair summaries before
  inferential resampling; raw sample pairs are not treated as independent.
- Optional randomization inference uses paired study-level sign flips with
  `B=9,999` and reports `(b+1)/(B+1)`. Exact zero p-values are impossible.
- Report point estimate, interval, number of tissues, studies, samples, folds,
  completed seeds, failures, and adjusted p-value where applicable.
- If fewer than four independent study blocks contribute to a contrast, report
  the effect descriptively with no p-value and mark uncertainty insufficient.

## 12. Failures, missingness, and null results

Every attempted fold/method/seed receives one of:

- `completed`;
- `single_study_tissue_not_estimable`;
- `training_feature_contract_failed`;
- `held_out_mapping_unavailable`;
- `distance_or_ph_resource_cap`;
- `baseline_resource_cap`;
- `incomplete_seed_set`;
- `no_stable_k`;
- `degenerate_distance_matrix`;
- `metric_not_identifiable`; or
- `software_failure`.

Failures are never replaced with another sample, study, method, feature panel,
or favorable seed after labels are opened. Denominators and reasons accompany
every summary. Null, negative, contradictory, and technically excluded results
remain in the immutable bundle.

## 13. Execution stages and resource gates

| Stage | Work | Advance condition |
|---|---|---|
| MV5-A | Contract fixtures and pilot technical smoke test | No biological outputs; IDs, folds, label firewall, baseline formulas pass |
| MV5-B | Two-study inductive integration/mapping and baseline feasibility | Held-out query contract passes; numerical provenance complete |
| MV5-C | One complete eligible tissue across all five seeds/methods | Matched axes, failures retained, total cost projected |
| MV5-D | Full 90-sample/15-study candidate execution | All frozen methods attempted under approved caps |
| MV5-E | Labels opened; endpoints/inference only | Predictions/matrices already immutable; complete reporting |

Initial caps for stages B/C are 8 GiB and 30 minutes per sample-view job, at
most two heavy workers, and 24 aggregate worker-hours. The MV-04 landscape
profile informs caching and scheduling; it does not authorize truncation. Stage
C must write a revised full-run estimate and storage budget before D.

## 14. Existing-data and new-data decision

Existing data are sufficient to attempt cross-study endpoints for five tissues
and to determine whether the inductive integration contract is executable. No
new-data trigger is opened now.

A later trigger may be proposed only if an executed existing-data result leaves
a named decision-relevant gap, such as too few independent studies for a tissue
or no valid held-out integration mapping. The trigger must record the estimand,
candidate data, confounding structure, access rights, expected value, and cost,
and requires owner/author-team approval.

## 15. Evidence and implementation references

- Luecken et al. (2022), integration benchmarking with biological conservation
  and batch-removal evaluated separately:
  <https://doi.org/10.1038/s41592-021-01336-8>.
- Crowell et al. (2020), multi-sample single-cell inference and pseudobulk
  aggregation: <https://doi.org/10.1038/s41467-020-19894-4>.
- Zimmerman et al. (2021), single-cell pseudoreplication and experimental-unit
  dependence: <https://doi.org/10.1038/s41467-021-21038-1>.
- Bernau et al. (2014), cross-study validation and the optimism of ordinary
  within-study cross-validation: <https://doi.org/10.1093/bioinformatics/btu279>.
- nestedcv (2023), training-only feature selection and nested validation:
  <https://doi.org/10.1093/bioadv/vbad048>.
- Tibshirani and Walther (2005), prediction-strength concepts for cluster
  validation: <https://doi.org/10.1198/106186005X59243>.
- Székely and Rizzo (2013), energy statistics based on distances:
  <https://doi.org/10.1016/j.csda.2013.03.018>.

These sources justify design principles; they do not establish biological
performance for scPHcompare.

## 16. Authorization boundary

Freezing this document authorizes implementation of fold-fit transformations,
matched baselines, deterministic clustering helpers, provenance, and staged
technical smoke tests. It does not authorize opening biological outcomes,
changing methods after results, cell/gene fusion, manuscript claims, new data,
or a public scientific default. G-MV5 passes only after the full frozen
single-view benchmark is executed and each view receives an independent
works/fails/uncertain disposition.
