# scPHcompare Project Evidence and Scientific Audit Baseline

## Document control

| Field | Value |
|---|---|
| Status | Draft evidence baseline; findings require verification before manuscript use |
| Version | 0.3.9 |
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

- GSE120221, 25 bone-marrow samples.
- Used as a technical/homogeneity check.
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

## 11. Change history

| Version | Date | Summary |
|---|---|---|
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
