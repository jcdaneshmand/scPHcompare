# scPHcompare Project Evidence and Scientific Audit Baseline

## Document control

| Field | Value |
|---|---|
| Status | Draft evidence baseline; findings require verification before manuscript use |
| Version | 0.3.14 |
| Created | 2026-08-03 |
| Canonical project | `jcdaneshmand/scPHcompare` |
| Evidence owner | Jonah Daneshmand |
| Purpose | Preserve what the dissertation and preprint establish, identify discrepancies and risks, and ground the modernization/resubmission roadmap |

This file separates four kinds of statements:

- **Documented:** explicitly stated in the dissertation or preprint.
- **Observed:** directly observed in the local repository or source code.
- **Inference:** a reasoned interpretation that must be tested.
- **Proposed:** a recommended next step, not a completed result.

No source PDF was edited. Page references below use PDF page numbers, followed by the printed dissertation page where useful.

## 1. Evidence sources

### Primary project documents

1. `docs/Dissertation_SubmissionReady_October.pdf`
   - Title metadata: *A Persistent Homology Framework for scRNA-seq: Assessing Clustering Performance and Biological Signal Across Integration Strategies*.
   - 159 PDF pages.
   - Relevant sections: research questions and objectives (PDF 18-20; printed 4-6), methods (PDF 46-65; printed 32-51), clustering interpretation (PDF 82-84; printed 68-70), curve interpretation (PDF 123-126; printed 109-112), future directions (PDF 127-133; printed 113-119), and limitations (PDF 139-140; printed 125-126).
2. `docs/Jonah-BioRxiv_v2.pdf`
   - Title: *Impact of integration on persistent homology clustering and biological signal detection in scRNA-seq data*.
   - 19 PDF pages.
   - BioRxiv preprint, also indexed by PubMed/PMC as a preprint (PMID 40766507).
   - Relevant sections: methods (PDF 2-6), results (PDF 6-17), discussion and conclusion (PDF 17-18).

### Local repositories inspected

1. `E:\Repositories\Jonah\PH Pipeline Repo\scPHcompare`
   - Local `main` at `64ec4fd` during inspection.
   - Cached `origin/main` was two commits ahead at `3910b15`.
   - `docs/` and `example_run.r` were untracked.
   - R package structure with `DESCRIPTION`, `NAMESPACE`, `renv.lock`, `R/`, `man/`, `tests/`, `README.md`, `LICENSE`, and `CITATION.cff`.
   - No GitHub Actions workflow was present; `.github` contained only `copilot-instructions.md`.
2. `E:\Repositories\Jonah\PH_ClusteringApp`
   - Earlier prototype/research workspace and original GitHub remote.
   - Contains the early Shiny application and historical pipeline generations.
   - Its relationship to the paper package must be documented, but it should not be assumed to be the current manuscript repository.

### Current literature sampled for roadmap updates

The following are initial high-value sources, not a complete systematic review:

- Luecken et al. (2022), [Benchmarking atlas-level data integration in single-cell genomics](https://doi.org/10.1038/s41592-021-01336-8). Establishes the scIB framework and separates batch removal from biological conservation across multiple integration methods and preprocessing choices.
- Dann et al. (2023), [Population-level integration of single-cell datasets enables multi-scale analysis across samples](https://doi.org/10.1038/s41592-023-02035-2). Provides a population-level integration perspective relevant to sample-aware evaluation.
- Zhang et al. (2025), [Recovery of biological signals lost in single-cell batch integration with CellANOVA](https://doi.org/10.1038/s41587-024-02463-1). Demonstrates that integration can remove meaningful biological variation and introduce distortion, directly relevant to the paper's central question.
- Huynh and Cang (2024), [Topological and geometric analysis of cell states in single-cell transcriptomic data](https://doi.org/10.1093/bib/bbae176). Introduces scGeom and combines curvature with persistent homology on cell and gene networks.
- Wang et al. (2023), [K-nearest-neighbors induced topological PCA for single-cell RNA-sequence data analysis](https://pubmed.ncbi.nlm.nih.gov/37961744/). Introduces persistent-Laplacian/topological PCA alternatives and evaluates them on 11 scRNA-seq datasets.
- Hernandez-Lemus (2025), [Topological data analysis in single cell biology](https://doi.org/10.3389/fimmu.2025.1615278). A recent review useful for updating the field context and terminology.
- Bauer (2021), [Ripser: efficient computation of Vietoris-Rips persistence barcodes](https://doi.org/10.1007/s41468-021-00071-5). Important for understanding the compiled backend and why a Rust rewrite must be evidence-gated.
- Zhang, Xiao, and Wang (2020), [GPU-accelerated computation of Vietoris-Rips persistence barcodes](https://doi.org/10.4230/LIPIcs.SoCG.2020.70). Reports the Ripser++ approach and substantial GPU speedups for appropriate workloads.
- Solomon et al. (2025), [The Flood Complex: Large-Scale Persistent Homology on Millions of Points](https://papers.nips.cc/paper_files/paper/2025/hash/aba03e6f25decb32bda9c5bf81c58305-Abstract-Conference.html). A recent scalable-complex alternative motivated by limitations of Rips, alpha, sparse Rips, and witness complexes.

## 2. Project identity and repository relationship

### Documented

- Both the dissertation and preprint identify `https://github.com/jcdaneshmand/scPHcompare` as the code repository.
- The project analyzes whole scRNA-seq samples as point clouds and compares their sample-level topological signatures.
- `PH_ClusteringApp` contains an earlier Shiny prototype and historical scripts.

### Observed

- `scPHcompare` has already undergone meaningful packaging and cleanup work that is absent from the historical `PH_ClusteringApp` repository.
- The existing modernization plan was initially written in `PH_ClusteringApp` before the documents revealed that `scPHcompare` is the repository cited by the dissertation and preprint.

### Confirmed repository lineage

The owner confirmed on 2026-08-03 that `scPHcompare` is the public-facing evolved version of the original scripts and prototype under `PH_ClusteringApp`. Therefore:

- `scPHcompare` is the canonical repository for scientific reproduction, modernization, and resubmission.
- `PH_ClusteringApp` is the historical source/prototype lineage.
- Any future Shiny application work should be scoped as a deliberate interface for `scPHcompare`, not treated as the authoritative scientific implementation.

### Contributor credit

The owner confirmed that both Dr. Eric Rouchka and Dr. Akshitkumar Mistry will receive credit on the project. Exact authorship order and CRediT roles remain decisions for the author team. Repository and manuscript work must preserve historical contributions while separately recording new contributions made during modernization.

## 3. Scientific intent recovered from the documents

### Central hypothesis

The dissertation hypothesizes that integration transforms scRNA-seq geometry in a way that may suppress fine-scale/local structure while making global topological structure more detectable by persistent-homology representations (dissertation PDF 18; printed 4). It further proposes that full pairwise sample relationships may retain biological information even when group-average topological curves converge.

### Research questions

The dissertation asks whether PH can:

- recover biologically meaningful sample groupings across heterogeneous scRNA-seq collections;
- remain useful in raw or minimally normalized data containing batch effects;
- quantify how Raw, SCT Individual, SCT Whole, and Seurat-integrated representations alter topology;
- compare favorably with K-means and Seurat/Louvain clustering;
- reveal information in full pairwise landscape distances that is lost in mean-curve summaries;
- identify structural relationships complementary to local-density and correlation-based methods.

### Original objectives

- Analyze 124 sample datasets across eight tissues and four preprocessing configurations.
- Compute H0/H1 persistence diagrams, Betti curves, Euler curves, and persistence landscapes.
- Build bottleneck, spectral, and landscape distance matrices.
- Apply hierarchical and spectral clustering to sample-level topological distances.
- Compare against conventional K-means and Seurat/Louvain methods.
- Assess biological alignment with ARI, NMI, Jaccard, purity, and VI relative to random baselines.
- Statistically compare topological curves across tissues and preprocessing iterations.

## 4. Data and analysis design recovered from the documents

### Primary heterogeneous collection

- Eight tissue labels: bone marrow, colon, liver, pancreatic islets, PBMC, prostate, substantia nigra, and testis.
- 124 sample-level point clouds according to the dissertation.
- Local metadata currently contains the following sample counts:

| Tissue | Samples |
|---|---:|
| Bone marrow | 31 |
| Colon | 13 |
| Liver | 6 |
| Pancreatic islets | 12 |
| PBMC | 12 |
| Prostate | 16 |
| Substantia nigra | 9 |
| Testis | 28 |

- Data were sourced through PanglaoDB with GEO/SRA identifiers.
- The collection is unbalanced across tissues and contains samples nested within SRA studies.
- The dissertation table marks the substantia nigra study as single-nucleus RNA-seq; the remaining listed studies appear to be scRNA-seq.

### Homogeneous collection

- GSE120221, 25 separately processed bone-marrow libraries from 20 donors.
- MV8-B accession closure proves these are the same 25 `SRA779509` libraries
  already present in the 124-sample corpus, not an independent external set.
- Retained as a cell-view technical/reprocessing check. Its primary frozen
  gene view is non-estimable because the historical axis retains 467/500
  panel genes.
- The topological clustering evaluation sets `k=25`, matching the number of samples.

### Quality control

- Retain cells with 500-9,000 detected genes.
- Require less than 20% mitochondrial reads and greater than 5% ribosomal reads.
- Remove genes detected in fewer than three cells.
- Technical covariates were intentionally not regressed during SCTransform.

### Four preprocessing configurations

1. **Raw:** filtered counts without normalization.
2. **SCT Individual:** SCTransform applied independently to each sample.
3. **SCT Whole:** samples merged and SCTransform applied jointly, without formal integration.
4. **Integrated:** Seurat SCT integration using RPCA anchors, with CCA fallback and `k.weight=50`.

### Topological calculation

- Each whole sample is treated as a cell-by-gene point cloud.
- Euclidean distance is used.
- Vietoris-Rips persistent homology is computed for H0 and H1 using `ripserr`.
- Derived representations include Betti curves, Euler curves, and persistence landscapes.
- Sample-to-sample matrices include bottleneck distance (BDM), a spectral transformation of BDM (SDM), and persistence-landscape distance (LDM).

### Reported principal finding

The paper and dissertation report a performance inversion: conventional methods perform better for tissue recovery on raw/minimally processed data, while landscape-based sample clustering improves relative to conventional methods after Seurat integration. The documents interpret this as evidence that integration weakens local geometry while clarifying global sample-level structure.

This finding is scientifically interesting, but the audit items below must be resolved before it is treated as confirmed.

## 5. High-priority scientific and implementation audit findings

The findings in this section are issues to verify, not accusations about how the work was produced.

### S-001 — Conventional and PH comparisons use different analytical units

**Priority:** Critical

**Evidence:** The dissertation explicitly describes K-means and Seurat/Louvain as cell-level clustering and PH methods as sample-level clustering (PDF 58-59; printed 44-45). In the current code, PH sample assignments are copied to every cell before ARI/NMI and other metrics are calculated (`assign_ph_clusters()` and `run_cluster_comparison()`).

**Risk:** Cell-level conventional labels and replicated sample-level PH labels are not directly comparable observations. Samples with more retained cells receive more weight in PH metrics, and the effective sample size differs between method families. This can change both effect sizes and null distributions.

**Required resolution:** Rebuild the core comparison at a common unit of analysis. The preferred primary comparison is sample-level PH versus sample-level non-topological baselines. Cell-level clustering can remain a separate secondary analysis with separate claims.

### S-002 — The bone-marrow `k=25` result is a technical distinctness check, not biological validation

**Priority:** Critical

**Evidence:** The paper states that hierarchical clustering with 25 samples and `k=25` is expected to score perfectly whenever distinct samples have nonzero pairwise distances (paper PDF 6-8). The resulting perfect score therefore follows from assigning one cluster per sample.

**Risk:** Calling this a validation dataset or evidence of clustering quality overstates what the result demonstrates.

**Required resolution:** Relabel it as a pipeline integrity/distinctness check or replace it with a validation task that has replicated biological classes, held-out samples, and nontrivial cluster structure.

### S-003 — Tissue, study, technology, and sample size are confounded

**Priority:** Critical

**Evidence:** Samples are nested in SRA studies; tissue representation is unbalanced; and the dissertation acknowledges that the scRNA-seq versus snRNA-seq result may reflect tissue rather than technology because additional snRNA-seq tissues are absent (dissertation PDF 123; printed 109). The shorter paper reports persistent technology differences but omits this caveat.

**Risk:** Tissue recovery may partly measure study identity, platform, cell count, or other dataset-specific structure rather than tissue biology.

**Required resolution:** Use study-blocked resampling, leave-one-study-out evaluation, balanced sensitivity analyses, explicit variance/confounding models, and datasets in which biological and technical variables are not perfectly nested.

### S-004 — Persistence-landscape definitions disagree across paper, dissertation, and code

**Priority:** Critical

**Evidence:**

- The preprint describes using the first landscape function (`k=1`) for H0 and H1, then combining dimensions by an L2 norm (paper PDF 5).
- The dissertation describes a pointwise L2 norm over all available landscape functions within each homology dimension, followed by an L2 combination across H0/H1 (dissertation PDF 55; printed 41).
- Current `compute_landscape_curve()` uses `rowMeans()` when a landscape is a matrix, then combines H0 and H1 by a pointwise Euclidean norm.
- The pre-audit public path used TDA's implicit `KK=1` default and a fixed `[0,1]` grid. The compatibility path now requests `KK=1` explicitly without changing its scientific behavior.
- Eight saved historical analyses contain five levels. The 124-sample artifacts are oriented `5x500`, while the 25-sample artifacts are `500x5`; their saved LDM entries exactly reproduce index-wise five-level Frobenius distances.
- Historical finite filtration endpoints differ across samples by factors from approximately 4.2 to 52,410, so index-wise comparisons on per-diagram grids do not compare functions at common filtration values.
- Legacy combined LDMs are generally H0-dominated; the median H1 energy fraction is below 0.004 in six of eight profiles and reaches 0.058 in the highest profile.
- Convergence profiling found exact active depths of 2,999–34,839 for H0 and 4–4,327 for H1 in the diagnostic strata. No universal level cap represented every preprocessing condition.
- Although 250-point distance ranks were nearly identical to 1,000-point ranks at a high cap, an analytic energy audit showed that even 1,000 uniform points could miss more than 1%—and sometimes all—of a narrow sample landscape. Rank agreement between fixed grids is therefore not a sufficient numerical-accuracy check.

**Risk:** The principal representation behind the strongest result is not defined consistently. Grid mismatch, orientation-dependent averaging, and fixed unit-grid truncation can materially change or erase the representation, making legacy landscape-derived clustering and claims unsuitable for reuse.

**Required resolution:** The project owner approved the dissertation-aligned all-level, exact/error-controlled L2 target on 2026-08-05. The R exact oracle now agrees with corrected Persim critical-pair integration and independent SciPy quadrature to floating precision. Persim 0.3.8's built-in norm fails sign-changing-difference fixtures and is rejected. Fresh eligible diagrams, production batch validation, and broader author-team confirmation remain required before activation or regeneration. H0 and H1 must be reported separately before any secondary combination.

### S-005 — Null models and p-values require statistical redesign

**Priority:** Critical

**Evidence:** The dissertation describes 100 curve bootstraps, 10 sets of 5 random groups for some curve nulls, and 10,000 label permutations for other comparisons. Current helper defaults include 50 or 1,000 replicates in different locations. Cluster p-values use the fraction of random metrics at least as extreme as observed and may return exact zero.

**Risks:**

- Exact zero is not an appropriate finite Monte Carlo p-value.
- Randomization may not preserve sample nesting, tissue imbalance, cluster sizes, or cell-count weights.
- Multiple families of pairwise and pointwise tests are not defined consistently.
- Pointwise curve values are strongly dependent across filtration thresholds.
- Some uses of Wasserstein distance appear to treat generic curves as probability distributions; assumptions and normalization need formal verification, especially for signed Euler curves.

**Required resolution:** Write a statistical analysis plan before rerunning results. Use `(b+1)/(B+1)` Monte Carlo p-values, blocked/stratified permutations, uncertainty intervals, explicit multiplicity families, and functional-data or global curve tests whose assumptions match the objects being compared.

### S-006 — Sample-specific filtration truncation may compromise distance comparability

**Priority:** High

**Evidence:** The dissertation reports adaptive SCT Whole thresholds based on the 90th percentile of `k=100` nearest-neighbor distances (PDF 64; printed 50). Current code can select thresholds per expression matrix and stores the selected threshold per job.

**Risk:** Persistence diagrams truncated at different filtration maxima may not be comparable without explicit censoring/truncation treatment. Preprocessing iterations may also be compared over normalized `[0,1]` grids despite having different physical distance scales.

**Required resolution:** Define a common-scale or scale-normalized filtration policy, quantify sensitivity to truncation, and record thresholds as first-class provenance in every result.

### S-007 — The strongest causal/mechanistic claims exceed the current evidence

**Priority:** High

**Evidence:** The paper attributes conventional-method decline to integration warping local geometry and PH improvement to clarified global topology. It uses phrases such as “significantly outperform” while presenting normalized scores, but does not show a direct paired inferential comparison between method families.

**Risk:** A plausible post hoc mechanism may be presented as demonstrated causality. PH stability to bounded perturbations does not by itself imply robustness to batch effects that alter topology.

**Required resolution:** Recast the mechanism as a hypothesis unless directly tested. Add local/global geometry diagnostics, direct method-difference uncertainty, controlled simulations, ablations, and integration methods with different geometrical behavior.

### S-008 — Dataset inventory contains apparent cell-count inconsistencies

**Priority:** High

**Evidence:** Several dissertation Table 1 rows report more cells after filtering than before filtering, including examples in PBMC, prostate, substantia nigra, and testis entries (PDF 49-51; printed 35-37).

**Risk:** The table may combine counts from different stages, contain transcription errors, or reveal metadata mismatches.

**Required resolution:** Regenerate dataset tables directly from an auditable manifest and pipeline outputs. Validate accession uniqueness, cell counts, tissue labels, sequencing approach, and exclusion reasons.

### S-009 — The paper lacks a genuinely independent biological validation task

**Priority:** High

**Evidence:** The homogeneous dataset validates technical distinctness; the heterogeneous collection is used to develop and interpret the primary finding. No held-out multi-study biological validation is described.

**Risk:** The principal result may be specific to the chosen studies, labels, preprocessing, or oracle `k`.

**Required resolution:** Add at least one external heterogeneous validation collection and a study-held-out design. Controlled simulation should complement, not replace, real external validation.

### S-010 — Cluster count uses ground-truth labels

**Priority:** High

**Evidence:** The dissertation and code set `k` to the number of known tissue, SRA, or approach groups.

**Risk:** This is an oracle-`k` benchmark, not fully unsupervised discovery. It is valid only if labeled accurately and accompanied by non-oracle sensitivity analyses.

**Required resolution:** Report oracle-`k` results explicitly, then add cluster-number selection/stability analyses and continuous sample-retrieval or classification tasks that do not require the known class count.

### S-011 — Raw-count Euclidean geometry needs justification and sensitivity analysis

**Priority:** High

**Evidence:** Raw filtered counts are used as a PH input with Euclidean distance. Alternative distances are deferred to future work.

**Risk:** Library size, sparsity, feature scale, and cell number can dominate sample topology. Comparing point clouds with different numbers of cells can alter PH summaries even without biological differences.

**Required resolution:** Include cell-count-matched subsampling, repeated landmark sampling, library-size controls, reduced-space inputs, and alternative metrics such as correlation/cosine or other biologically justified distances.

### S-012 — Current package execution and CI need verification

**Priority:** High

**Evidence:** The package has tests but no GitHub Actions workflow. `run_postprocessing_pipeline()` calls `assignRandomGroup()`, but no implementation was found in tracked `R/` files; the test suite mocks it.

**Risk:** A clean installation may not run the advertised full pipeline, and mocked tests may conceal missing production dependencies.

**Required resolution:** Run package checks in a clean environment, add an end-to-end toy-data smoke test without mocks, add CI, and resolve missing or implicitly sourced functions.

### S-013 — Persistence homology used features rather than cells as observations

**Priority:** Critical

**Evidence:** The dissertation states that each sample's expression matrix produced a cell-by-cell Euclidean distance matrix (PDF page 54; printed page 40). Current `R/PH_Calculation.R` and the historical Part 1 script extract Seurat assays as feature-by-cell matrices. Current and historical `PH_Functions.R` pass those matrices directly, without transposition, to `ripserr::vietoris_rips()`. The ripserr point-cloud contract interprets rows as points. All nine audited historical diagram artifacts carry this code-path conflict; integrated H0 counts of 2,999 and much larger raw/SCT counts are consistent with features acting as points.

**Risk:** The persistence diagrams describe gene/feature geometry, not the intended cell geometry. Every downstream bottleneck, spectral, landscape, clustering, statistical, figure, and manuscript result derived from those diagrams is scientifically ineligible for the corrected analysis.

**Resolution:** MV-01 froze the prospective definitions, MV-02 enforces them
through versioned typed objects, and MV-03 has now generated 132 scientifically
eligible corrected H0/H1 pilot diagrams for both views under monitored resource
caps. MV-03 recovered missing SCT Pearson residuals from the recorded SCT model,
proved cell pairing by canonical identifiers, fit label-free cohort panels and
equal-contribution transformations, and retained five-seed sensitivity. All
jobs completed; the slowest took 14.64 seconds and peak process-tree RSS was
about 1.04 GiB. This passes G-MV3 for MV-04 technical distance validation only.
Clustering, fusion, and biological interpretation remain prohibited.
Historical artifacts remain restricted to labeled reproduction or stress
testing.

## 6. Manuscript, reference, and figure audit findings

### M-001 — Preliminary citation mismatches

**Priority:** High

- Preprint reference 1 (Lum et al., 2013) is a general Mapper/TDA paper but is cited after a sentence about scRNA-seq revolutionizing transcriptomics.
- Preprint reference 3 (Wang et al., 2019) is a single-cell topology paper but is cited after a statement about the growth of public scRNA-seq datasets.
- The DOI printed for reference 3 appears malformed and unrelated (`10.1249/MSS.0000000000000901.Physical`); the authoritative article record should replace it.

These examples justify a complete sentence-to-source citation audit rather than a DOI-only cleanup.

### M-002 — Version and identifier inconsistencies

**Priority:** Medium

- The paper lists Seurat 5.3.0 in one location and describes imports using Seurat 5.2.1 elsewhere.
- Figure 2's caption uses `GSE12022` instead of `GSE120221`.
- Package versions in prose, the dissertation table, `DESCRIPTION`, and `renv.lock` need one generated source of truth.

### M-003 — Figure provenance and regeneration

**Priority:** High

The owner reports that some content was created with earlier LLMs and some figures were made manually. Therefore every figure should receive:

- source/provenance classification: generated by code, manually assembled, AI-assisted, or unknown;
- source data and generating script where applicable;
- exact parameter/configuration record;
- accessibility and legibility review;
- check that caption, labels, equations, and manuscript interpretation match the rendered figure;
- reconstruction in a reproducible workflow for every quantitative figure.

The synthetic PH illustration should be regenerated from deterministic source points and code. Manuscript figures should avoid displaying degenerate comparisons as if they were informative results.

### M-004 — LLM-assisted text requires claim-level review

**Priority:** High

Each substantive sentence should be classified as:

- supported directly by project results;
- supported by an external source;
- a hypothesis/inference;
- background/common knowledge;
- unsupported or too strong.

This review should focus especially on claims of deformation invariance, batch-effect robustness, biological interpretation of H1 features, causal explanations of integration geometry, and suggested clinical applications.

### M-005 — Full reference verification protocol

For every reference:

1. Match title, authors, venue, year, volume, pages/article number, DOI, and PubMed identifier against an authoritative record.
2. Check correction, expression-of-concern, and retraction status.
3. Confirm that the cited source actually supports the sentence.
4. Prefer the peer-reviewed version over a preprint when available.
5. Replace outdated background claims with newer reviews or primary work where appropriate.
6. Record verification date and verifier in a reference ledger.

## 7. Recommended scientific reframing

### Current broad claim

Integration creates a state in which PH landscape clustering outperforms conventional clustering for tissue recovery.

### More defensible provisional framing

**Proposed:** scPHcompare is a sample-level topological auditing framework for measuring how integration changes biological conservation and technical mixing in collections of scRNA-seq point clouds.

This framing preserves the genuinely interesting contribution while reducing dependence on a potentially non-comparable cell-level-versus-sample-level performance claim. It also aligns naturally with modern integration benchmarking, which evaluates both batch removal and biological conservation.

### Candidate primary questions for revision

1. Do sample-level PH representations add information beyond matched non-topological sample representations when evaluated on the same samples?
2. How do different integration methods change within-biological-group cohesion and between-study mixing in topological space?
3. Which PH representation, distance, filtration policy, and sampling strategy is most stable under repeated cell subsampling?
4. Can topology detect integration-induced distortion that established scIB metrics miss, and vice versa?
5. Does the result generalize across studies, tissues, sequencing technologies, and external datasets?

### Matched sample-level baselines

Potential baselines should operate on the same sample units as PH:

- pseudobulk expression distances;
- sample centroids in PCA or integrated latent space;
- cell-type composition distances;
- distributional distances between sample embeddings;
- kernel or optimal-transport sample distances;
- graph-derived sample descriptors;
- conventional clustering applied to each matched sample-distance representation.

Cell-level Louvain/K-means remains valuable, but should be presented as a separate analysis rather than a directly commensurate sample-level competitor.

## 8. Updated technical priorities

1. Establish a canonical sample-level data contract and manifest.
2. Correct and test the persistence-landscape definition.
3. Establish common filtration/scaling and repeated subsampling policies.
4. Build matched sample-level baselines and study-blocked evaluation.
5. Add external validation and controlled simulations.
6. Convert file-existence caching to dependency-aware workflow management (`targets` remains a strong candidate).
7. Profile the corrected pipeline before optimizing.
8. Benchmark existing Ripser, GPU, sparse-complex, landmark, and newer scalable-complex approaches before writing Rust.
9. Use Rust only for a measured, isolated bottleneck not already solved by an established implementation.

## 9. Literature-update workstreams

### Integration benchmarking

- Compare Seurat RPCA with a justified subset such as Harmony, LIGER, FastMNN, Scanorama, and scVI/scANVI where data and compute permit.
- Adopt the scIB distinction between biological conservation and batch removal.
- Include recent evidence that integration may remove biological signal or introduce distortion.
- Evaluate at both cell and sample levels without conflating the units.

### Topology and geometry in single-cell analysis

- Position scPHcompare relative to scTDA, scGeom, topological/persistent PCA, persistent community detection, and recent TDA reviews.
- Clarify whether the novelty is the PH representation, sample-level comparison, integration auditing, or the combined benchmark.
- Avoid claiming PH is intrinsically batch-proof; test sensitivity explicitly.

### Scalable PH

- Benchmark cell subsampling and landmark strategies with stability estimates.
- Evaluate reduced spaces and alternative filtrations/complexes.
- Compare current `ripserr`/Ripser with GPU or other established implementations.
- Consider recent alternatives such as Flood complexes only after checking implementation maturity, guarantees, and fit to the biological question.

## 10. Evidence needed from the project owner

- Exact submitted manuscript version if different from `Jonah-BioRxiv_v2.pdf`.
- Supplemental tables and figures referenced by the paper.
- Any response-to-reviewers draft.
- Source files for manually assembled figures.
- Notes identifying which text or figures were LLM-assisted.
- The exact analysis run/configuration used for the BioRxiv v2 figures and tables.
- Access/location information for large intermediate artifacts that cannot be regenerated cheaply.

The original editorial decision and reviewer reports have now been recovered. Because the letter states that the reports are confidential and must not be published without editor permission, the verbatim source and detailed response matrix are stored under the Git-ignored `docs/private/` directory. Public planning documents use only generalized workstreams derived from them.

## 11. MV5-P complete-training-distance evidence

**Observed:** The complete label-closed MV5-O distance scope has been
materialized and independently validated: 150 groups, 4,340 separate exact
all-active-level H0/H1 landscape chunks, 150 representation-specific energy
units, 75 shared-pseudobulk units, 4,565 total immutable units, 1,838,725
values, and 525 complete symmetric zero-diagonal finite/nonnegative training
matrix components.

**Observed:** Twelve frozen exact R oracles pass; both prespecified maximum
groups reproduce 32 landscape outputs plus energy byte-for-byte; the extra SCT
shared-pseudobulk output passes a separately reported supplemental repeat; and
all 4,565 output/status pairs survive completed-run resume with unchanged
hashes, sizes, and timestamps and zero rebuild.

**Observed:** Execution uses `12.044379` worker-hours, reaches
`492163072` bytes maximum process-tree RSS, and stores
`4570070656` bytes in the private production root, within all frozen caps.
The prior 1.278-GB estimate is now classified as output-focused and incomplete
because it omitted full repeated group-local interval staging; the discrepancy
was disclosed during execution and remained below a conservative
`6.011`-GB reserve forecast and the 10-GiB hard guard.

**Scientific boundary:** These are validated training-distance inputs, not a
biological result. No production clusters, held-out assignments, labels,
ARI/NMI, tissue or method selection, robustness result, gene topology, fusion,
new data, optimization, or manuscript claim were produced in MV5-P.

## 12. MV5-Q label-closed clustering-artifact evidence

**Observed:** All 150 frozen fold/representation/distance groups completed from
the accepted MV5-P training matrices and MV5-D5/MV5-J query distances. The run
fit 6,750 candidate PAM models and produced 567,000 private candidate-partition
rows, 1,350 public stability rows, 126,000 selected PAM/average-linkage
training-partition rows, and 9,000 held-out assignments across 750 complete
analysis-matrix contexts.

**Observed:** Independent validation reproduces all 150 five-seed one-SE k
selections, canonical sorted-member cluster identities, one-medoid-per-PAM-
cluster identities, and every PAM and average-linkage held-out assignment. The
canonical maximum-size fold (`large_loso_v1:SRA713577`, 89 training samples)
reproduces all 40 artifacts byte-for-byte. Full resume reuses all 150 groups and
leaves all 753 paths, hashes, sizes, and timestamps unchanged with zero rebuild.

**Observed:** Clustering itself uses `134.125` elapsed seconds across groups,
with `1.838` seconds maximum group time and `1219973120` bytes peak process RSS,
inside the frozen limits. Ten query aliases are explicitly audited across all
15 folds and five seeds. The stopped preliminary root is not accepted; it
exposed and led to correction of an SCT pseudobulk representation alias before
the accepted `production-v2` run.

**Scientific boundary:** Stability is agreement among seed-specific training
partitions, not agreement with a biological label. No biological or technical
label was opened, and no training-alignment or held-out-generalization endpoint,
method ranking, robustness claim, gene topology, fusion, new-data result,
optimization, or manuscript claim was produced in MV5-Q.

## 13. MV5-R prediction-locked clustering-outcome prefreeze evidence

**Observed:** The prefreeze binds 18 sources at accepted MV5-Q commit
`f16321c`, including `joined_metadata_cellcounts.csv` by SHA-256
`e2b6f1e7d7dc05ba8e80627aad4ea1677c3b718a844397acd06596d7e22a18d0`
without copying or tracking the external file. Its structural audit confirms
124 source rows/18 studies, 90 candidate samples/15 studies, five tissues, two
approaches, and an exact match to the MV5-Q sample axis.

**Observed:** The registry freezes PAM as primary, average linkage as
sensitivity, six descriptive training ARI/NMI endpoints, and two training-only
plurality held-out balanced-accuracy endpoints. The immutable execution queue
has 2,400 unique units: 150 groups by two algorithms by eight endpoints. Twelve
validations, ten abort rules, aggregation, uncertainty, missingness, public-
safety, and resource rules are prospective and explicit.

**Observed:** Three studies contain both sequencing approaches. The contract
therefore preserves tissue-stratified study-block resampling for tissue but
uses global study-block resampling for approach, avoiding false splitting of
mixed studies. All 18 source hashes validate; 8/8 regenerated public artifacts
are byte-identical; focused tests pass 18/18; the full suite passes 803 checks
with two expected skips; and the staged-tree standard package check is
`Status: OK`.

**Scientific boundary:** No real ARI, NMI, balanced accuracy, p-value,
confidence interval, method ranking, tissue ranking, or biological/technical
interpretation was computed or opened. Held-out study prediction is forbidden;
fold-specific cluster IDs cannot be pooled; labels may be used only under the
frozen training-only prediction rule in a separately authorized execution.

## 14. MV5-S prediction-locked clustering-outcome evidence

**Observed:** Prospective engine commit `c3f8da0` and execution-freeze commit
`4f7f73c` bind the exact max-NMI, ARI, training-only plurality prediction,
seed/fold aggregation, tissue-stratified/global-study bootstrap, atomic resume,
and public-safety implementations before outcomes. All 2,400 units complete,
producing 9,000 training seed metrics, 18,000 private held-out predictions,
3,000 held-out seed summaries, and 40 macro contexts with zero p-values and
zero method selection.

**Observed:** Held-out tissue balanced accuracy spans 0.0103 to 0.2578 across
the 20 fixed representation-distance-algorithm contexts; held-out approach
balanced accuracy spans 0.4925 to 0.5000. Training fold-mean ranges are tissue
ARI -0.0066 to 0.4796/NMI 0.0295 to 0.5395, study ARI -0.0176 to 0.3815/NMI
0.0240 to 0.4676, and approach ARI -0.0339 to 0.1976/NMI 0.0052 to 0.1078.
Training folds overlap and remain descriptive.

**Observed:** Independent validation passes 12 categories, including every
training metric, held-out prediction, bootstrap, and aggregation. Clean repeat
reproduces 2,400/2,400 private outcomes and 8/8 outcome tables byte-for-byte.
Resume reuses 2,400/2,400 units while 4,800 artifact/status files remain
unchanged. First-pass resources are 51.888 unit-seconds, 0.570 seconds maximum
unit time, and 740,687,872 bytes peak RSS. Final verification parses 5/5 MV5-S
scripts, passes 831 repository checks with two expected skips, and yields a
clean staged-tree package `Status: OK`.

**Scientific boundary:** These are complete secondary clustering outcomes, not
a method-ranking or manuscript-claim exercise. No p-value, winner, refit,
post-outcome tuning, spectral promotion, robustness result, gene topology,
fusion, new-data result, or optimization claim was produced.

## 15. MV5-T selection-resistant robustness/gap-gate evidence

**Observed:** The prospective gate at `7f6784d` freezes 14 committed sources
plus all 150 accepted private SCT/integrated coordinate-file hashes. Three
minimum/median/maximum fold-seed pairs expose 270 sample views per
representation; all paired 384-by-30 shapes, sample/cell identities, nested
192-in-256 selections, first-20 coordinates, finite values, and nonzero norms
pass.

**Observed:** Ten candidate families are scored without consulting MV5-S values.
Nested 192/256-cell counts, 20-PC truncation, and cosine-chord geometry are
admitted as four one-factor-at-a-time configurations. The unopened resource
queue has exactly 24 groups. Ten gate artifacts reproduce byte-for-byte;
focused tests pass 16/16 and the full suite passes 847 checks with two expected
skips.

**Observed:** Prior measurements project 15.542 worker-hours and 10.18 GB for
the four-setting full grid. Full execution is denied. Only a 24-group, one-
worker, two-hour, 2-GiB resource admission is authorized next.

**Scientific boundary:** MV5-T calculates no robustness diagram, distance,
retrieval endpoint, cluster, label outcome, ranking, or claim. All admitted work
is post-outcome secondary sensitivity and cannot rescue or replace completed
results.

## 16. MV5-U bounded robustness resource-admission evidence

**Observed:** All 24 frozen label-closed units complete across three folds, two
representations, and four one-factor-at-a-time cell-count/PC/geometry
configurations. The run builds 2,160 views, 1,069,189 finite H0/H1 intervals,
4,320 landscape summaries, 1,536 exact landscape pair rows, and 768 matched
energy rows.

**Observed:** The revised dissertation-aligned landscape definition is enforced:
H0/H1 separate, essential H0 excluded, all active consecutive levels, no level
cap, no uniform grid, and exact critical-pair L2 integration. Independent
validation passes 15/15 categories, including all-view MST oracles, analytic
square H1, analytic exact landscape, sampled energy, source identities, and
configuration isolation.

**Observed:** Maximum unit time is 55.508 seconds, total measured unit time is
895.449 seconds, peak process-tree RSS is 622,227,456 bytes, and private
production storage is 288,635,915 bytes. A clean repeat matches 168/168
scientific artifacts byte-for-byte. Resume reuses 24/24 units with all 240
private files unchanged.

**Verification:** The focused MV5-U correction tests pass 26 expectations, the
full repository suite passes with two established CRAN skips, and a build/check
of the exact staged index finishes with package `Status: OK`.

**Correction:** Initial independent validation exposed a Python `True`/`False`
versus R logical parsing defect in the validator. Commit `8bc2718` adds strict
cross-language boolean parsing and tests. The frozen production implementation
digest is unchanged, no scientific artifacts were modified, and corrected
validation passes.

**Scientific boundary:** This is resource and semantic admission evidence, not
a robustness outcome. Labels and outcomes remain closed and full robustness is
explicitly unauthorized. The next eligible action is only a prospective
streamed full-robustness execution gate.

## 17. MV5-V streamed full-robustness prefreeze evidence

**Observed:** The exact label-closed scope contains 600 groups, 54,000 views,
282,800 heldout-training biological pairs, 565,600 H0/H1 exact-landscape rows,
2,880 deterministic subchunks, 282,800 energy rows, and 1,131,200 assembled
four-method rows. All 150 private coordinate identities join one-to-one to the
accepted pair scopes; no label or outcome field enters the queue.

**Observed:** Accepted historical PH, landscape, and retrieval-input stages
project 69,381.174 seconds (19.273 worker-hours) across four configurations.
The prospective program cap is one worker, 30 worker-hours, 4 GiB RSS/group,
and 16 GiB private storage, with configuration-stratified stop decisions.

**Scientific boundary:** The calculation primitives are bound, but the full-
group orchestration engine is deliberately unbound. All 600 queue rows retain
`execution_authorized=FALSE`; MV5-W runner implementation/binding and one real
label-closed smoke are required before any configuration launch. MV5-V computes
no robustness outcome, comparison, ranking, or claim.

## 18. MV5-W full-robustness launch-readiness evidence

**Observed:** The prospectively selected integrated PC20 smoke completes 90
views, 425 directed pairs, 850 exact all-active H0/H1 landscape rows in four
subchunks, 425 energy rows, and 1,700 label-closed method rows in 70.085 seconds
at 541,970,432 bytes peak RSS and 15,559,362 private bytes.

**Observed:** Thirteen independent validation categories pass. All eight
deterministic artifacts reproduce byte-for-byte, and validation-only resume
reuses the group with all 11 files unchanged.

**Scientific boundary:** Only the first PC20 configuration's 150 groups are
eligible for a separate prospectively bound execution. No label, retrieval
outcome, robustness comparison, rank, winner, or claim was calculated; the
other three configurations remain unauthorized.

## 19. MV5-X PC20 configuration execution evidence

**Observed:** All 150 prospectively authorized PC20 groups complete on the
first pass with 13,500 views, 70,700 heldout-training biological pairs, 27,000
landscape summaries, 141,400 exact all-active H0/H1 landscape rows, 70,700
energy rows, and 282,800 four-method rows. The other 450 queue rows remain
unauthorized and unexecuted.

**Observed:** Total measured group time is 11,163.624 seconds (3.101 worker-
hours), maximum group time is 185.348 seconds, peak process-tree RSS is
638,365,696 bytes, and final private storage is 2,487,457,825 bytes. These are
well within the frozen 30-hour program, 600-second/4-GiB group, and 16-GiB
storage caps.

**Validation:** Fifteen independent categories pass across all 150 groups and
1,350 manifest-declared artifacts. All 13,500 H0 MST checks pass; all 282,800
method rows reconstruct from H0/H1 landscape and energy inputs; 16/16
prospectively marked repeat artifacts match byte-for-byte; and all 1,650 files
remain unchanged through a full 150-group validation-only resume.

**Scientific boundary:** This accepts the PC20 calculation as a robustness
input, not a robustness outcome. Labels, ranks, winners, comparisons, and
biological outcomes remain closed. A separate prospective outcome-evaluation
contract is required before label access; the other three configurations
remain unauthorized.

## 20. MV5-Y PC20 robustness-outcome prefreeze evidence

**Observed:** All 150 accepted PC20 groups and all eight representation/family
axes pair exactly to the accepted 30-coordinate retrieval inputs on fold, seed,
held-out query, training reference, and mapped method. The paired scope is
70,700 biological pairs and 282,800 method rows, with zero missing or excess
rows. The structural 90-sample/15-study external key join passes while label
values and accepted outcome tables remain unopened.

**Frozen:** The complete retrieval-only panel has 24 estimands: 16 direct
PC20-minus-PC30 changes and eight H0/H1 topology-increment DID changes over MRR
and fixed-1-NN. Every estimand receives paired tissue-stratified study-block
bootstrap uncertainty; exactly four MRR DID tests form the sole Holm-adjusted
family. No equivalence/noninferiority margin or claim is authorized.

**Validation:** A clean second assembly reproduces all 13 generated contract
ledgers byte for byte.

**Boundary:** MV5-X contains zero within-training PC20 pairs, so clustering
sensitivity is not identifiable without 525,350 new label-closed
within-training biological pairs across both representations and new frozen
cluster artifacts. All 150 MV5-Y execution flags remain false; ranks, labels,
outcomes, comparisons, and selection remain zero. MV5-Z may later execute only
the frozen retrieval contract.

## 21. MV5-Z PC20 retrieval-robustness execution evidence

**Prediction firewall:** A committed prediction lock (`c16f2b2`) binds 282,800
canonical PC20 ranking rows after all 178 sources and all 150 groups validate.
Independent reconstruction of every rank and tie passes before tissue access.

**Observed:** All 24 fixed PC20-minus-PC30 estimands complete. Primary MRR
topology-increment changes are -0.02036 (integrated H0), -0.06444 (integrated
H1), -0.02567 (SCT H0), and 0.00537 (SCT H1). Integrated H1 and SCT H0 paired
bootstrap intervals are below zero, but all four Holm-adjusted sign-flip
p-values are at least 0.352. The result is heterogeneous sensitivity evidence,
not equivalence, robustness success, superiority, or a default-setting result.

**Validation:** Fifteen independent categories reconstruct 282,800 ranks,
3,600 query-method outcomes, 10,800 query estimands, all aggregations, exact
bootstrap/sign matrices, 24 intervals, four raw/Holm p-values, and 150 private
hashes. All 150 private artifacts and 16 deterministic public files repeat
byte-for-byte; 300 private result/status files and all 17 runner public files
remain unchanged through resume.

**Boundary:** PC20 clustering and the other three robustness configurations are
not authorized. A selection-resistant continuation gate must decide the next
calculation using the pre-existing configuration order and complete evidence.

## 22. MV5-AA selection-resistant continuation evidence

**Order and decision:** The canonical MV5-V configuration order is PC20,
cosine-chord, nested 192 cells, then nested 256 cells. A committed rule binds
the complete PC20 analysis without estimate- or subgroup-driven selection and
authorizes only the unique next cosine configuration.

**Scientific basis:** PC20 changes coordinate count under Euclidean geometry;
cosine-chord tests radial-scale dependence while holding 384 cells and 30
coordinates fixed. The first sensitivity cannot resolve the second named
alternative, so continuation is scientifically identifiable and one-factor-at-
a-time regardless of favorable or unfavorable PC20 subgroups.

**Scope and validation:** The later label-closed scope is 150 groups, 13,500
views, 70,700 pairs, 141,400 landscape requests, 70,700 energy rows, and
282,800 method rows. Twelve standalone validation categories pass; two complete
11-ledger assemblies are byte-identical.

**Boundary:** No cosine calculation, rank, label, outcome, clustering, or
nested-cell work occurred. A bound cosine-only execution sprint is next.

## 23. MV5-AB label-closed cosine calculation evidence

**Observed calculation:** Exactly 150 cosine-chord groups complete: 13,500
row-normalized views, 141,400 exact H0/H1 landscape rows, 70,700 energy rows and
282,800 method rows. Runtime is 2.608 worker-hours; all resource caps pass.

**Independent geometry evidence:** A helper-independent validator renormalizes
all raw coordinates and matches all 5,170,500 H0 MST deaths (maximum error
`5.11e-15`) plus 30 stratified direct energy values (maximum error `4.88e-15`).
Fifteen artifact categories and two byte-identical validation assemblies pass.

**Reproducibility:** Sixteen deterministic repeat artifacts match exactly. All
1,650 private paths, sizes, hashes and timestamps remain unchanged through a
150-group validation-only resume.

**Boundary:** This accepts a label-closed calculation, not a cosine biological
result. Rankings, labels, outcomes, clustering and nested-cell settings remain
closed; a separate prediction-locked outcome prefreeze is next.

## 24. MV5-AC cosine retrieval-outcome prefreeze evidence

**Exact pairing:** All eight representation/family axes pair exactly between
the accepted 384-cell/30-coordinate Euclidean and cosine-chord calculations.
All 282,800 keys match on representation, fold, seed, query, reference, and
mapped method, with zero missing or excess rows.

**Frozen analysis:** The complete panel has two endpoints, 16 direct geometry
changes, eight topology-increment DIDs, paired tissue-stratified study-block
bootstrap intervals, and four Holm-adjusted MRR H0/H1 sign-flip tests. The
prediction lock orders by distance then canonical sample ID and must precede
tissue-only access. No equivalence claim is authorized.

**Validation and boundary:** 187 sources and all 150 private manifests bind;
11 independent categories pass; 14/14 ledgers repeat byte-identically. Labels,
cosine ranks, endpoints, outcomes, and selection remain zero. Clustering is not
identifiable because all 525,350 required within-training biological pairs are
absent; nested-cell settings remain unauthorized.

## 25. MV5-AD cosine retrieval-robustness execution evidence

**Prediction firewall:** Engine `c22d667` generated 282,800 canonical cosine
ranks. Independent reconstruction passed before the complete prediction lock
was committed at `2e2f9f3`; tissue values and accepted endpoint rows remained
unopened until after that commit.

**Observed:** All four topology-increment cosine-minus-Euclidean MRR changes
were negative. Integrated H0/H1 were -0.07216/-0.10106 (Holm 0.2734/0.1052).
SCT H0/H1 were -0.11833/-0.10088 (Holm 0.0224 each). This is a bounded geometry
sensitivity result, not equivalence, universal Euclidean superiority, or a
default-setting claim.

**Validation:** Fifteen independent categories reconstructed 282,800 ranks,
3,600 cosine outcomes, 3,600 Euclidean pairings, 10,800 query estimands, all
aggregations, exact bootstrap/sign matrices, 24 intervals, four raw/Holm
p-values, and 150 private hashes. Twenty-two deterministic repeat comparisons
pass; all 600 prediction/outcome unit files survive resume unchanged.

**Boundary:** Cosine clustering and both nested-cell settings remain closed. A
selection-resistant continuation gate must decide whether nested 192 cells is
the next justified calculation without subgroup or estimate selection.

## 26. MV5-AE selection-resistant nested-cell continuation evidence

**Decision:** The canonical order fixes nested 192 after completed PC20 and
cosine. Both prior 24-estimand/four-test panels are bound without slicing; the
decision interface consumes no representation, dimension, tissue, endpoint,
seed, estimate, interval, p-value, or method ranking.

**Distinctness and feasibility:** Nested 192 changes only cell depth while
retaining 30 PCs, Euclidean geometry, and deterministic inclusion within the
frozen 384 cells. Six real admissions and both complete prior configurations
support a one-worker, six-hour, 4-GiB later envelope.

**Scope and validation:** Exactly 150 later label-closed groups are authorized:
13,500 views, 70,700 pairs, 141,400 landscape rows, 70,700 energy rows, and
282,800 method rows. Ten independent categories and 11/11 byte repeats pass;
calculation/ranking/labels/outcomes/clustering/nested-256 remain zero.

## 27. MV5-AF label-closed nested-192 calculation evidence

**Observed calculation:** Exactly 150 nested-192 Euclidean groups complete:
13,500 views, 141,400 exact H0/H1 landscape rows, 70,700 energy rows, and
282,800 method rows. Runtime is 1.233 worker-hours; maximum unit time is 53.302
seconds, peak RSS is 454,893,568 bytes, and private storage is 1,254,700,479
bytes. Every resource cap passes.

**Exact nesting and independent evidence:** The frozen selection is the first
192 cells in a sample/seed/cell-ID SHA-256 order, not source rows 1–192. A
helper-independent validator reconstructs every selected identity, proves
192-within-256 inclusion, matches all 2,578,500 H0 MST deaths, 60 direct H1
diagrams, 60 exact landscape distances, and 30 energy distances. All 1,350
manifest hashes and 282,800 method rows validate.

**Reproducibility and boundary:** Sixteen group-repeat artifacts and 11/11
clean validator ledgers match byte-for-byte. All 1,650 private paths, hashes,
sizes, and timestamps remain unchanged through a 150-group resume. This accepts
only the label-closed calculation; ranks, labels, outcomes, clustering, nested
256, gene/fusion/new-data/Rust/default/claim work remain closed.

## 28. MV5-AG nested-192 retrieval-outcome prefreeze evidence

**Exact comparison:** The fixed contrast changes only cell depth from the
accepted 384 cells to its deterministic nested-192 subset. All 150 SCT and
integrated groups cross-link to the exact accepted 384-cell coordinate source;
30 coordinates, Euclidean geometry, folds, seeds, pair axes, landscape
semantics, and method families remain fixed.

**Pairing and analysis freeze:** All eight method axes pair exactly over
282,800 baseline and 282,800 nested-192 rows, with zero missing, excess, or
duplicate keys. The future distance/sample-ID tie order, two endpoints, 16
direct changes, eight topology-increment DIDs, blocked bootstrap, four MRR
sign-flip tests, and Holm family are frozen while tissue and endpoint values
remain unread.

**Validation and boundary:** The prefreeze binds 188 sources, passes eight
explicit acceptance criteria and 12 independent categories, and reproduces 15
production plus one validator ledger byte-for-byte. Retrieval is identifiable;
clustering is rejected because 525,350 required within-training pairs are
absent. Ranks, labels, outcomes, selection, clustering, and nested 256 remain
zero.

## 29. MV5-AH prediction-locked nested-192 outcome evidence

**Prediction lock:** All 282,800 nested-192 ranking rows were constructed and
independently reconstructed across 150 groups with all 188 sources valid, then
committed at `1a197a8` before tissue access. Clean prediction repeats reproduce
the complete gzip and all 150 private ranking payloads; all 300 private files
remain unchanged through resume.

**Complete results:** Tissue-only execution produced 7,200 endpoint rows and
all 24 fixed estimands. The four primary topology-increment MRR changes are
integrated H0/H1 -0.01080/-0.01115 and SCT H0/H1 +0.00131/-0.01366. All four
intervals cross zero and all Holm p-values equal 1.0. This is no detected
change, not equivalence or default evidence.

**Validation and correction:** Fifteen independent outcome categories
reconstruct endpoints, baseline pairing, all aggregation, bootstrap/sign
matrices, nulls, p-values, Holm adjustment, and private hashes. A resume-only
inference-file rewrite was corrected prospectively at `41bc7c7`; fresh accepted
runs preserve the same scientific values, reproduce 17 deterministic public
ledgers plus 150 private outcomes and inference matrices, and leave all 301
private paths/hashes/sizes/timestamps unchanged through resume. Clustering and
nested 256 remain closed.

## 30. Change history

| Version | Date | Summary |
|---|---|---|
| 0.3.32 | 2026-08-11 | Completed the prediction-locked nested-192 retrieval analysis: independently reconstructed and committed 282,800 ranks before tissue access; executed 7,200 endpoint rows and all 24 fixed estimands; found small primary topology-increment MRR changes with all intervals crossing zero and Holm p=1 without equivalence claims; passed 10 prediction and 15 outcome validation categories; audited a resume-only inference-file correction; reproduced rankings, 17 outcome ledgers, 300 private payloads and inference matrices; and preserved all 300 prediction plus 301 outcome files through immutable resumes while keeping clustering and nested 256 closed |
| 0.3.31 | 2026-08-11 | Completed the nested-192 retrieval-outcome prefreeze: bound 188 sources, proved exact 384-source identity for all 150 groups and exact pairing of all eight 282,800-row method axes, froze 24 complete-reporting estimands and four primary tests before label access, rejected clustering for 525,350 missing within-training pairs, passed eight criteria plus 12 independent categories, and reproduced 16/16 ledgers byte-for-byte while keeping ranks, labels, outcomes and nested 256 closed |
| 0.3.30 | 2026-08-11 | Completed the exact 150-group label-closed nested-192 calculation in 1.233 worker-hours; clarified and reconstructed the SHA-256 nested-cell order, proved all 13,500 192-within-256 selections, matched 2,578,500 H0 MST deaths plus 60 H1, 60 landscape and 30 energy oracles, validated all 282,800 method rows, reproduced 16 group and 11 validator artifacts byte-for-byte, and preserved all 1,650 paths/hashes/sizes/timestamps while keeping ranks, labels, outcomes, clustering and nested 256 closed |
| 0.3.29 | 2026-08-11 | Completed the post-cosine selection-resistant continuation gate: bound both complete result panels without slicing, preserved nested 192 as canonical position three, excluded nine subgroup/result selection inputs, authorized exactly 150 later label-closed groups under six-admission/two-full-precedent resource evidence, passed 10 independent categories and 11 byte repeats, and kept calculation, rankings, labels, outcomes, clustering, and nested 256 closed |
| 0.3.28 | 2026-08-11 | Completed the committed prediction-locked cosine retrieval analysis: independently reconstructed and committed 282,800 ranks before tissue access; executed 7,200 endpoint rows, 24 intervals and four Holm tests; retained negative topology-increment changes including SCT H0/H1 Holm 0.0224; passed 15 independent categories, 22 repeat comparisons, and unchanged 600-file prediction/outcome resumes without equivalence, clustering, nested-cell, or default claims |
| 0.3.27 | 2026-08-11 | Completed the cosine retrieval-outcome prefreeze: bound 187 sources, proved all eight 282,800-row Euclidean/cosine axes pair exactly, froze a durable distance/sample-ID prediction lock plus 24 complete-reporting estimands and four primary tests, validated 11 categories, repeated 14 ledgers byte-identically, and excluded clustering because 525,350 within-training pairs are absent while keeping rankings, labels and outcomes closed |
| 0.3.26 | 2026-08-11 | Completed the exact 150-group label-closed cosine-chord calculation in 2.608 worker-hours; independently renormalized all 13,500 views, reconstructed 5,170,500 H0 MST deaths and 30 stratified energy values, validated all 282,800 method rows and 1,350 manifest hashes, repeated 16 artifacts byte-identically, and preserved all 1,650 paths/hashes/sizes/timestamps through resume while keeping rankings, labels and outcomes closed |
| 0.3.25 | 2026-08-11 | Froze and applied the selection-resistant post-PC20 continuation rule, bound all 24 estimands/four primary tests and the canonical configuration order without subgroup selection, authorized exactly 150 later label-closed cosine groups, independently validated 12 categories, and reproduced two 11-ledger assemblies byte-for-byte while executing no new calculation or outcome |
| 0.3.24 | 2026-08-11 | Completed the committed prediction-locked PC20 retrieval robustness analysis: 282,800 canonical ranks, 7,200 endpoint rows, 24 paired estimands/intervals and four Holm-adjusted tests; independently reconstructed all values, repeated 150 private and 16 public artifacts byte-identically, and preserved 300 private files through resume while retaining heterogeneous negative sensitivity without equivalence or selection claims |
| 0.3.23 | 2026-08-11 | Froze the retrieval-only PC20 robustness evaluation after proving all eight 30-PC/PC20 prediction axes pair exactly; bound 178 sources, 24 complete-reporting estimands, paired blocked uncertainty and four-test multiplicity, reproduced 13/13 ledgers byte-identically, and formally rejected clustering from directed-only MV5-X artifacts while keeping ranks, labels, outcomes, and execution at zero |
| 0.3.22 | 2026-08-11 | Completed and independently validated the exact 150-group PC20 robustness calculation: 13,500 views, 141,400 H0/H1 landscape rows, 282,800 reconstructed method rows, 16 byte-identical repeat artifacts, and an unchanged 1,650-file full resume in 3.101 worker-hours, while keeping all outcomes and other configurations unauthorized |
| 0.3.21 | 2026-08-10 | Completed MV5-W launch readiness: implemented/bound the full-pair atomic runner, executed one prospectively selected real 90-view/425-pair PC20 smoke, passed independent numerical validation, eight-artifact byte repeat and unchanged 11-file resume within caps, and authorized only a later 150-group PC20 execution with outcomes closed |
| 0.3.20 | 2026-08-10 | Completed MV5-V streamed full-robustness prefreeze: bound 176 sources/calculation primitives, derived the exact 600-group/54,000-view and 565,600-landscape-row scope, froze deterministic subchunks/atomic resume/repeat/validation/abort contracts and a 30-hour/16-GiB envelope, while leaving orchestration and all execution/outcomes unauthorized |
| 0.3.19 | 2026-08-10 | Completed MV5-U bounded label-closed robustness resource admission: 24 units/2,160 views, exact dissertation-aligned H0/H1 landscapes, independent numerical validation, 168 byte-identical clean-repeat artifacts, zero-rebuild 240-file resume, passing resource caps, and an audited validator-only cross-language boolean correction; full robustness remains unauthorized |
| 0.3.18 | 2026-08-10 | Completed the MV5-T selection-resistant robustness/gap gate: froze 164 sources, validated paired coordinate readiness, admitted four one-factor-at-a-time cell-count/PC/geometry configurations, emitted a 24-group no-outcome admission queue, rejected immediate full execution on 15.54-hour/10.18-GB projections, and retained all outcomes closed |
| 0.3.17 | 2026-08-10 | Executed and independently validated the complete MV5-S prediction-locked clustering outcome contract: 2,400 units, 9,000 training metrics, 18,000 held-out predictions, blocked bootstrap intervals, byte-identical clean repeat, zero-rebuild resume, and complete secondary reporting without p-values or method selection |
| 0.3.16 | 2026-08-10 | Completed the MV5-R prediction-locked clustering-outcome prefreeze: bound 18 source identities and the external label hash, verified the exact 90-sample join without copying labels, froze two algorithm roles/eight endpoints/aggregation/uncertainty/aborts/resources, and emitted a 2,400-unit no-outcome execution queue |
| 0.3.15 | 2026-08-10 | Completed and independently validated MV5-Q label-closed clustering artifacts for all 150 frozen analysis groups, including all candidate PAM grids, one-SE selections, canonical PAM/average partitions, held-out assignments, ten-alias audit, byte-identical maximum-fold repeat, and zero-rebuild full resume while keeping all outcomes closed |
| 0.3.14 | 2026-08-10 | Completed and independently validated all frozen MV5-P label-closed training-distance inputs, matrices, exact oracles, maximum-group repeats, all-unit resume, and observed resource/storage evidence while keeping clustering and outcomes closed and explicitly correcting the scope of the prefreeze storage estimate |
| 0.3.13 | 2026-08-10 | Completed the MV5-O label-closed production prefreeze: froze exact sources, group/chunk/baseline queues, executable landscape and baseline implementations, resource/storage/abort/resume rules, and independent/byte-repeat/real-runner evidence while keeping full production, clustering, labels, and outcomes at zero |
| 0.3.12 | 2026-08-10 | Completed the MV5-N label-closed clustering and resource gate: froze training-only PAM and deterministic held-out assignment, retained the dissertation-aligned exact all-active-level H0/H1 landscape definition, reproduced the complete training-pair scope, passed bounded independent exact admissions and matched baselines, and kept full production and all outcomes closed |
| 0.3.11 | 2026-08-09 | Completed 6,750 corrected integrated cells-as-observations H0/H1 records with independent identity/file/diagram validation, stored and fresh MST evidence, exact complete-group repeat, immutable resume, and measured authorization for dissertation-aligned exact all-active-level integrated landscapes |
| 0.3.10 | 2026-08-09 | Completed all 75 fixed-D1-panel label-closed integrated coordinate groups: 6,750 views and 450 mappings passed 675 independent checks, exact complete-group repeat, zero-rebuild resume, and cap-passing measured reprojection; integrated cell PH is separately authorized but not executed |
| 0.3.9 | 2026-08-06 | Completed MV5-A/MV5-B immutable outcome-label-closed LOSO manifests, exact matched-baseline implementations and fixtures, and deterministic synthetic Seurat reference mapping; advanced only to MV5-C feasibility |
| 0.3.8 | 2026-08-06 | Froze the prospective MV-05 sample-level LOSO benchmark, matched baselines, label firewall, clustering/uncertainty/multiplicity rules, existing-data feasibility, and transductive-integration exclusion without computing outcomes |
| 0.3.7 | 2026-08-05 | Completed MV-04 immutable all-level exact H0/H1 landscape distances, eligible R-oracle and deterministic validation, fit-scope-bound normalization, contribution/resource profiling, bounded diagram-distance sensitivities, and the negative Rust gate; advanced only to MV-05 statistical-plan design |
| 0.3.6 | 2026-08-05 | Completed MV-03 representation-native residual recovery, paired extraction, frozen feature panels, matched subsampling, fit-scope scaling/PCA, 132 monitored corrected H0/H1 jobs, and five-seed feasibility summaries; advanced MV-04 only |
| 0.3.5 | 2026-08-05 | Implemented and analytically validated orientation-safe typed cell/gene topology objects, corrected-only PH dispatch, explicit legacy stamping, immutable cache/provenance identities, and full-shape versus fixture eligibility separation |
| 0.3.4 | 2026-08-05 | Froze the MV-01 dual-view scientific contract, deterministic existing-data pilot manifest, distance/clustering hierarchy, leakage and provenance rules, resource envelope, and stop/go criteria |
| 0.3.3 | 2026-08-05 | Recorded independent landscape-oracle results, the Persim sign-crossing failure, and hard invalidation of historical diagrams/results caused by feature-as-point PH orientation |
| 0.3.2 | 2026-08-05 | Added exact/adaptive reference-engine validation, analytical and representative benchmark evidence, established-library assessment, and the negative current Rust decision |
| 0.3.1 | 2026-08-05 | Recorded project-owner approval of the dissertation-aligned all-level, exact/error-controlled landscape target while preserving implementation, eligibility, and broader author-team gates |
| 0.3.0 | 2026-08-05 | Added exact level-depth, cap, grid-energy, clustering-stability, H0/H1 contribution, runtime, and memory evidence; rejected a universal fixed grid/count and proposed error-controlled all-level L2 |
| 0.2.0 | 2026-08-05 | Reconciled persistence-landscape definitions against code and legacy artifacts; recorded grid/orientation/level contradictions, H0 dominance, proposed contract, and invalidation boundary |
| 0.1.1 | 2026-08-03 | Confirmed `scPHcompare` as the canonical evolved repository, recorded contributor-credit requirements, and noted confidential review-material handling |
| 0.1.0 | 2026-08-03 | Initial document-, code-, and current-literature-grounded evidence baseline |
