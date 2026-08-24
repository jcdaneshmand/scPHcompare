# scPHcompare Modernization and Resubmission Plan

## Document control

| Field | Value |
|---|---|
| Status | Draft for owner and author-team review |
| Version | 0.4.19 |
| Created | 2026-08-03 |
| Canonical repository | `jcdaneshmand/scPHcompare` |
| Historical source repository | `jcdaneshmand/PH_ClusteringApp` |
| Plan owner | Jonah Daneshmand |
| Scientific contributors to preserve | Jonah Daneshmand, Dr. Eric Rouchka, Dr. Akshitkumar Mistry, Dr. Julia Chariker, and other contributors identified through repository history |
| Evidence baseline | `docs/PROJECT_EVIDENCE.md` |
| Confidential reviewer matrix | `docs/private/REVIEWER_RESPONSE_MATRIX.md` (Git-ignored; do not publish) |

This is the canonical auditable roadmap for the public `scPHcompare` project. The earlier plan in `PH_ClusteringApp` is retained as a historical preliminary draft because it was written before the repository lineage, dissertation, preprint, and reviewer reports were available.

## 1. Project outcome

The project will evolve from a dissertation-era analysis into a reproducible and defensible sample-level framework for studying how single-cell integration changes biological and technical structure.

The preferred provisional framing is:

> `scPHcompare` is a sample-level topological auditing framework for measuring how integration changes biological conservation and technical mixing across collections of scRNA-seq point clouds.

This framing remains provisional until the redesigned analysis establishes which claims are supported.

The final project should:

- preserve the historical analysis and contributor credit;
- provide one tested, documented, installable R package and analysis workflow;
- compare methods at a common observation level;
- distinguish biological conservation from batch removal;
- quantify sensitivity to outliers, cell counts, sampling, metrics, filtrations, and integration choices;
- establish a practical use case and external validation;
- use current literature to state novelty precisely;
- regenerate all quantitative figures and tables from versioned code;
- optimize only measured bottlenecks;
- produce a citable release aligned with the revised manuscript.

## 2. Non-negotiable controls

1. **Preserve before changing.** Historical code, manuscript outputs, and expensive artifacts must be inventoried before cleanup.
2. **Correctness before speed.** No Rust or deep optimization work begins until the primary statistical comparison is valid and reproducible.
3. **Same unit for primary comparisons.** Sample-level PH methods must be compared with sample-level non-topological baselines.
4. **Labels may evaluate, not silently design, the method.** Oracle `k` results must be labeled as such and supplemented with non-oracle analyses.
5. **Blocked validation.** Study, tissue, technology, and sample nesting must be respected in resampling and validation.
6. **Claims have provenance.** Every substantive claim must be tied to project evidence, an external source, or an explicit hypothesis.
7. **Figures are reproducible.** Quantitative figures must be generated from source data and scripts; manual/AI-assisted elements must have provenance.
8. **Confidential review stays private.** Reviewer reports and detailed response text remain under Git-ignored `docs/private/` unless the editors authorize publication.
9. **Credit is explicit.** Dr. Eric Rouchka and Dr. Akshitkumar Mistry will receive project credit. Exact author order and CRediT roles are decided by the author team.
10. **Negative results remain visible.** The project will report conditions under which PH adds no value or is less effective.

## 3. Current evidence summary

### Repository state

- `scPHcompare` is the owner-confirmed public evolution of `PH_ClusteringApp`.
- It is already structured as an R package with `DESCRIPTION`, `NAMESPACE`, `renv.lock`, documentation, tests, license, and citation metadata.
- At inspection, local `main` was two cached commits behind `origin/main` and had untracked `docs/` and `example_run.r` content.
- No GitHub Actions workflow was present.
- The full postprocessing path calls `assignRandomGroup()`, but no tracked implementation was located; existing tests mock it.

### Scientific state

- Primary analysis: 124 samples from eight tissues across Raw, SCT Individual, SCT Whole, and Seurat-integrated representations.
- Secondary homogeneous collection: 25 GSE120221 bone-marrow samples.
- PH: Euclidean Vietoris-Rips H0/H1, with Betti, Euler, landscape, bottleneck, spectral, and landscape-distance representations.
- Reported result: persistence-landscape clustering improves relative to conventional clustering after integration.
- Highest-priority risks:
  - cell-level conventional methods and replicated sample-level PH assignments are compared with the same cell-level metric machinery;
  - the `k=25` bone-marrow result is a technical distinctness check;
  - tissue, study, sequencing technology, retained cell count, and sample imbalance are confounded;
  - persistence-landscape definitions differ between preprint, dissertation, and code;
  - null models, curve tests, p-values, and filtration truncation require formal review;
  - there is no independent heterogeneous biological validation task.

These are audit findings to verify, not final judgments.

## 4. Success criteria

### Repository success

- Clean installation from a fresh checkout.
- `R CMD check` and automated CI pass.
- A deterministic toy workflow runs without mocked production functions.
- Dependencies, inputs, configuration, outputs, and compute requirements are documented.
- No confidential or oversized artifacts are accidentally tracked.

### Scientific success

- Primary comparisons operate on the same sample units.
- At least one external multi-study validation task is held out from method design.
- Biological conservation and technical mixing are assessed separately.
- Results are robust or explicitly sensitive to outliers, cell count, sampling seed, distance, filtration, and integration method.
- Direct uncertainty is reported for differences between methods.
- Main conclusions survive study-blocked evaluation or are narrowed honestly.

### Manuscript success

- Novelty is defined relative to verified prior work.
- The practical task is clear in the title, abstract, introduction, and methods.
- Every figure/table is reproducibly generated and traceable.
- References are verified for bibliographic accuracy, support, corrections, and retractions.
- Author roles and acknowledgements are approved by the contributor team.

## 5. Work phases

Detailed task-level execution controls, stable task IDs, evidence requirements, gate checklists, and the append-only work log are maintained in [`docs/plans/README.md`](docs/plans/README.md). The phase summaries below remain the strategic source of truth; the subplans govern execution and audit evidence.

## Phase 0 — Preserve, synchronize, and establish provenance

### Objective

Create a safe, canonical starting point without losing historical work or confidential material.

### Tasks

- Record current local and remote branch state before synchronizing.
- Review the two commits currently ahead on cached `origin/main` and preserve local untracked files.
- Inventory tracked, modified, deleted, ignored, and untracked content.
- Classify large artifacts as raw input, intermediate, reference output, publication result, cache, or disposable.
- Create checksums and a manifest for irreplaceable/expensive artifacts.
- Link `PH_ClusteringApp` history to `scPHcompare` in public documentation without merging repositories blindly.
- Preserve confidential editorial correspondence under `docs/private/` and verify it remains ignored.
- Identify the exact configuration and outputs used for BioRxiv v2.
- Record provenance for dissertation and manuscript figures, including manual and AI-assisted construction.

### Deliverables

- Repository and artifact inventory.
- Historical lineage note.
- Protected-artifact manifest.
- BioRxiv v2 run/output manifest.
- Figure provenance ledger.

### Gate G0

- No important local work is unclassified.
- Confidential materials are protected from publication.
- The author team agrees on the canonical historical baseline.

## Phase 1 — Audit scientific validity and implementation fidelity

### Objective

Determine whether the current code implements the documented method and whether the comparison supports the manuscript claims.

### Tasks

#### Observation-unit audit

- Trace the row represented by every cluster assignment and metric: cell, sample, study, or group average.
- Quantify the impact of copying sample-level PH clusters to all cells.
- Recalculate PH metrics once per sample and compare with current cell-weighted results.
- Determine exactly how PCA and integration spaces are shared across datasets.

#### Representation audit

- Resolve the disagreement among first-level landscapes, all-level L2 aggregation, and code-level `rowMeans()`.
- Verify matrix orientation and default `TDA::landscape()` level behavior for the locked package version.
- Verify BDM, SDM, and LDM formulas with analytical toy examples.
- Verify common grids, scaling, infinite deaths, threshold censoring, and H0/H1 combination.

#### Statistical audit

- Reconstruct every reported metric and p-value from source outputs.
- Audit random-group creation, cluster-size preservation, study blocking, seeds, and Monte Carlo resolution.
- Replace exact-zero Monte Carlo p-values with finite-sample estimates.
- Define multiple-testing families.
- Verify Wasserstein and KS assumptions for Betti, Euler, and landscape curves.
- Determine whether direct method-difference inference supports words such as “significantly outperforms.”

#### Data audit

- Regenerate cell counts and metadata tables.
- Resolve cases where reported post-filter counts exceed pre-filter counts.
- Verify all accessions, tissue labels, study nesting, and sequencing approach.
- Quantify group imbalance and cell-count dependence.

#### Execution audit

- Resolve missing/implicitly sourced functions such as `assignRandomGroup()`.
- Run unit tests and package checks in the locked environment.
- Add failures discovered during the audit as regression tests before fixing them.

### Deliverables

- Methods-to-code traceability matrix.
- Statistical audit report.
- Dataset integrity report.
- Minimal reproducible failures and regression fixtures.
- Claim disposition table: supported, narrowed, untested, or contradicted.

### Gate G1

- The author team knows which existing results remain valid.
- The primary analysis unit and landscape definition are approved.
- No performance optimization begins before G1.

## Phase 2 — Establish a reproducible baseline and public repository health

### Objective

Create a trustworthy executable baseline for the corrected method.

### Tasks

- Define typed inputs and outputs for every pipeline stage.
- Create a tiny analytical fixture with known PH behavior.
- Create a realistic benchmark subset.
- Lock package and system dependencies.
- Add deterministic seed policy and provenance metadata.
- Add a dependency-aware workflow, preferably evaluating `targets`.
- Add clean-environment CI for parsing, tests, documentation, package check, and toy smoke run.
- Improve `.gitignore` and size checks without exposing private materials.
- Make logs structured, bounded, and machine-readable.
- Add checkpoint validation based on content/configuration hashes rather than file existence alone.
- Document hardware, runtime, memory, and restart behavior.

### Deliverables

- Reproducible toy pipeline.
- Corrected reference outputs and tolerances.
- CI workflow.
- Dependency-aware pipeline prototype.
- Public repository audit and remediation list.

### Gate G2

- A clean clone can run the toy analysis.
- Corrected reference outputs are stable under repeated runs.
- The public repository contains no confidential correspondence.

## Phase 3 — Literature, novelty, reference, and figure audit

### Objective

Rebuild the scholarly foundation and visual evidence without relying on unverified legacy or LLM-generated language.

### Tasks

#### Literature and novelty

- Verify the prior-work identifiers raised during external review.
- Review scTDA, scGeom, topological/persistent PCA, persistent Laplacian, Hodge/differential-topology approaches, persistent community methods, and recent single-cell TDA reviews.
- Review modern integration benchmarking, especially the separation of biological conservation and batch removal.
- Review evidence that integration can remove biological signals or introduce distortion.
- Maintain a novelty matrix comparing problem, analysis unit, topology, integration, validation, biological interpretation, and software contribution.
- Reframe novelty as empirical/framework/software novelty unless a genuinely new algorithm is developed and validated.

#### Reference verification

- Match every title, author list, venue, date, page/article number, DOI, PMID/PMCID, and publication status to an authoritative record.
- Check corrections, expressions of concern, and retractions.
- Verify sentence-level support, not only citation metadata.
- Replace preprints with peer-reviewed versions where available.
- Correct known preliminary issues, including mismatched early citations and the malformed DOI for the Wang et al. single-cell topology article.

#### Figure audit

- Classify every figure as code-generated, manually assembled, AI-assisted, or unknown.
- Rebuild all quantitative figures from tidy versioned outputs.
- Regenerate the synthetic PH schematic from deterministic code.
- Remove or clearly mark degenerate comparisons.
- Simplify dense multi-panel plots and use accessible palettes, labels, and captions.
- Confirm every visual statement against underlying values.

### Deliverables

- Literature evidence table.
- Novelty matrix.
- Verified reference ledger.
- Figure provenance and regeneration ledger.
- Revised contribution statement.

### Gate G3

- Every citation and major claim has verified support.
- Every quantitative figure is reproducible.
- Novelty claims are precise and defensible.

## Phase 4 — Redesign the primary scientific comparison

### Objective

Build a fair sample-level benchmark that separates biology from technical variation.

### Primary estimand

Determine whether PH-derived sample representations provide stable, generalizable information about biological identity and integration-induced distortion beyond matched non-topological sample representations.

### Tasks

#### Matched sample-level baselines

- Pseudobulk expression distances.
- Sample centroids in a shared PCA or integrated latent space.
- Cell-type composition distances where annotations are available.
- Distributional or optimal-transport distances between sample embeddings.
- Graph-derived sample summaries.
- The same clustering algorithms applied across PH and non-PH distance matrices.

#### Validation design

- Study-blocked cross-validation or leave-one-study-out testing.
- Balanced tissue sensitivity analyses.
- Cell-count-matched repeated subsampling.
- Use the existing project datasets for the audit and initial redesigned benchmark.
- Evaluate a new external dataset only if the existing data leave a decision-relevant evidence gap; document the trigger, expected value, comparability, access, and cost before acquisition or processing.
- Controlled simulations with known biological and batch structure.
- Separate biological-conservation and batch-removal endpoints.
- Direct uncertainty for differences between methods.
- Oracle-`k` results labeled separately from selected/stability-based `k`.

#### Robustness controls

- Outlier injection/removal.
- Cell-identity shuffle and geometry-preserving controls.
- Composition-preserving and composition-changing controls.
- Repeated cell/landmark subsampling.
- Sensitivity to feature number, reduced dimension, and random seed.
- Common versus normalized filtration scales.

### Deliverables

- Prespecified statistical analysis plan.
- Matched sample-level benchmark.
- External-validation design.
- Robustness/control suite.

### Gate G4

- The primary benchmark is fair, blocked, and reproducible.
- The practical question is defined before broad method expansion.

## Phase 5 — Expand integration, topology, and clustering methods

### Objective

Test generalizability without creating an unstructured method zoo.

### Integration methods

Select a justified subset representing different paradigms, potentially including:

- Seurat RPCA/CCA;
- Harmony;
- LIGER;
- Scanorama;
- FastMNN;
- scVI/scANVI where environment and data permit.

Use native, documented output spaces. Do not force methods into an inappropriate common representation merely for convenience.

### Topological inputs and representations

- Euclidean, correlation, and cosine distances.
- Full expression versus a shared reduced space.
- H0/H1 primary analysis with a small H2 feasibility study.
- Correctly defined persistence landscapes.
- Betti/Euler summaries.
- Persistence images or other vectorizations if justified by the benchmark.
- Landmark, witness, sparse-Rips, alpha, or newer scalable complexes where their assumptions match the data.

### Clustering and related tasks

- Hierarchical linkage variants.
- Spectral clustering.
- PAM/k-medoids.
- Graph community detection such as Louvain/Leiden on sample graphs.
- Density-based clustering only when the sample count and geometry support it.
- Continuous retrieval/classification and distance-based tests that do not require a known `k`.

### Guardrails

- Prespecify the primary comparison set.
- Separate exploratory methods from confirmatory results.
- Control multiplicity and report all prespecified results.
- Prefer interpretable ablations over exhaustive parameter fishing.

### Deliverables

- Integration-by-topology benchmark matrix.
- Clustering/task comparison report.
- Sensitivity and ablation analyses.
- Generalizability conclusion with explicit boundaries.

### Gate G5

- Claims identify which integration/topology combinations work, fail, or remain uncertain.
- The Seurat-specific result is either generalized or explicitly labeled specific.

## Phase 6 — Establish biological interpretability and practical relevance

### Objective

Show what a topological signal means and how a researcher would use it.

### Tasks

- Select a biological case with known trajectories, states, perturbations, or cell-cycle structure.
- Avoid assigning biological meaning to H1/H2 solely because a loop/void exists.
- Link topology to annotations through controlled association, perturbation, and sensitivity analyses.
- Test whether topology improves integration selection beyond established metrics.
- Prototype a diagnostic report that identifies:
  - residual study signal;
  - loss of biological conservation;
  - integration-induced local/global distortion;
  - unstable results under subsampling or outliers.
- Test the proposed diagnostic on held-out data.
- Decide whether topology-aware parameter selection is supported or should remain future work.

### Deliverables

- Biological case study.
- Practical integration-audit workflow.
- Interpretation guide with explicit limitations.
- Decision on topology-aware optimization claims.

### Gate G6

- A user can state what decision `scPHcompare` helps make.
- Biological interpretations are supported rather than metaphorical.

## Phase 7 — Profile and optimize

### Objective

Reduce runtime and memory after the corrected analysis is stable.

### Measurements

- Wall-clock time and CPU utilization per stage.
- Peak memory, dense conversions, and object copies.
- Serialization/checkpoint I/O.
- Worker startup and data-transfer overhead.
- Scaling with cells, genes, dimensions, threshold, samples, and workers.
- Accuracy/stability versus subsampling or approximation.

### Optimization order

1. Remove invalid or unnecessary computation.
2. Preserve sparse data and avoid dense copies.
3. Compute once and reuse through dependency-aware caching.
4. Bound and configure parallelism.
5. Improve algorithms and data structures.
6. Use mature optimized libraries.
7. Evaluate GPU or alternative complex implementations.
8. Introduce Rust only for a remaining measured hotspot.

### Rust gate

A Rust component is approved only if:

- profiling identifies a stable dominant hotspot;
- an existing Ripser/GPU/sparse-complex solution is insufficient;
- a narrow input/output contract exists;
- equivalence and error tolerances are defined;
- benchmarks show a material end-to-end benefit;
- installation and cross-platform maintenance remain acceptable.

Likely candidates are bounded transformation or distance kernels, not an assumed wholesale rewrite of PH already delegated to compiled Ripser code.

### Deliverables

- Benchmark suite and bottleneck report.
- Before/after runtime and memory results.
- Approximation/stability trade-off report.
- Architecture decision record for any Rust component.

### Gate G7

- Optimization has measured benefit and preserves approved behavior.

## Phase 8 — Rewrite the manuscript and response package

### Objective

Build the paper around validated evidence rather than patching the previous narrative.

### Tasks

- Use the confidential response matrix internally to verify every concern is addressed.
- Rewrite the title, abstract, introduction, and contribution statement after analyses stabilize.
- Define dataset-level clustering and observation units immediately.
- Define integration versus harmonization consistently.
- Report H0/H1 scope and any H2 feasibility result.
- Report outlier and subsampling sensitivity.
- Separate results, interpretation, and hypotheses.
- Generate all manuscript values, figures, and tables from the workflow.
- Create a claim-to-code/data/configuration map.
- Draft responses supported by exact analyses and figure/table references.
- Conduct an author-team CRediT review, preserving credit for Dr. Rouchka and Dr. Mistry and recording all contributors accurately.

### Deliverables

- Revised manuscript.
- Supplement and reproducibility appendix.
- Confidential response-to-reviewers draft.
- Claim-to-evidence matrix.
- Approved authorship/CRediT statement.

### Gate G8

- Every concern has evidence, a correction, or a justified limitation.
- Manuscript and released workflow agree exactly.

## Phase 9 — Release and archive

### Objective

Publish a durable, citable research artifact.

### Tasks

- Clean-machine reproduction.
- GitHub release and semantic version.
- DOI archive through an approved repository.
- Data/artifact manifests and checksums.
- Runtime/hardware documentation.
- Reproducible example and full-analysis instructions.
- Final privacy, confidentiality, license, and citation audit.
- Confirm that private reviewer materials are absent from the release.

### Deliverables

- Versioned public release.
- Citable archive.
- Reproduction report.
- Final software/data availability language.

### Gate G9

- The release is installable, reproducible, citable, and aligned with the manuscript.

## 6. Immediate work packages

## WP-001 — Canonical-state and artifact audit

- **Phase:** 0
- **Status:** Proposed
- **Scope:** Synchronization plan, local-change preservation, artifact classification, BioRxiv run identification, confidential-file checks.
- **Scientific behavior change:** None.
- **Acceptance:** No unclassified irreplaceable files; exact paper outputs/configuration located or documented as missing.

## WP-002 — Observation-unit and PCA comparability audit

- **Phase:** 1
- **Status:** Proposed; highest priority
- **Scope:** Trace cell/sample units, PH assignment replication, PCA spaces, metric weighting, oracle `k`, and bone-marrow validation interpretation.
- **Scientific behavior change:** None during diagnosis.
- **Acceptance:** Reproduce current metrics and quantify how sample-level recalculation changes them.

## WP-003 — Persistence-landscape specification audit

- **Phase:** 1
- **Status:** In progress; definition, reference oracle, independent cross-check, and historical-diagram eligibility complete
- **Scope:** Reconcile paper, dissertation, TDA defaults, matrices, aggregation, distance calculation, and grids.
- **Scientific behavior change:** None until specification is approved.
- **Acceptance:** One mathematical specification, analytical fixtures, and a list of outputs requiring regeneration.

## WP-004 — Reference and figure provenance ledger

- **Phase:** 3
- **Status:** Proposed
- **Scope:** Verify citations and identifiers, check sentence-level support, classify all figures, locate source data/scripts.
- **Acceptance:** Every citation and figure has a provenance/status record.

## WP-005 — Statistical analysis plan

- **Phase:** 1 and 4
- **Status:** Proposed
- **Scope:** Estimands, blocked resampling, null models, Monte Carlo p-values, multiplicity, functional tests, uncertainty, external validation.
- **Acceptance:** Approved before expensive reruns.

## WP-006 — Dual-view and multiview topology framework

- **Phase:** 1, 4, 5, and 7
- **Status:** MV-06 complete with negative fusion disposition; MV6-H exact-commit outcome execution passes 15/15 independent and 13/13 byte-repeat gates; report cell/gene views separately and advance to MV-07 synthesis
- **Scope:** Corrected cell topology, deliberately specified gene topology, H0/H1 landscapes and diagram distances, matched clustering/baselines, and staged cell/gene fusion.
- **Scientific behavior change:** None until the MV-01 contract and subsequent implementation gates are approved.
- **Acceptance:** Each view passes independent definition, correctness, eligibility, and feasibility gates before fusion; a frozen confirmatory configuration precedes any full biological rerun.

## 7. Risk register

| ID | Risk | Likelihood | Impact | Mitigation | Status |
|---|---|---:|---:|---|---|
| R-01 | Current headline result changes under same-unit evaluation | High | Critical | Run WP-002 before manuscript or optimization work | Open |
| R-02 | Landscape outputs were generated from an inconsistent definition | High | Critical | Run WP-003 and regenerate dependent outputs | Open |
| R-03 | Tissue recovery reflects study/technology/cell count | High | Critical | Blocked validation, matched subsampling, external datasets | Open |
| R-04 | Confidential reviewer reports are accidentally published | Medium | Critical | Git-ignore `docs/private/`, CI secret/private-path check, release audit | Mitigated locally; verify continuously |
| R-05 | Old figures or text contain unverifiable AI/manual content | High | High | Provenance ledger and full regeneration/claim audit | Open |
| R-06 | Integration expansion becomes computationally unmanageable | Medium | High | Prespecify a representative method subset and use cached fixtures | Open |
| R-07 | Rust consumes effort without solving the bottleneck | Medium | High | Enforce G7 Rust gate | Open |
| R-08 | Novelty remains incremental after comparison with prior work | Medium | High | Reframe as validated framework/software contribution or narrow target venue | Open |
| R-09 | External validation data are unavailable or too confounded | Medium | High | Identify candidate datasets early and use simulation only as complement | Open |
| R-10 | Contributor roles are misrepresented during rewrite | Low | High | Author-team CRediT review and written approval | Open |

## 8. Decision log

| ID | Date | Decision | Basis | Status |
|---|---|---|---|---|
| D-001 | 2026-08-03 | `scPHcompare` is the canonical public project; `PH_ClusteringApp` is historical lineage | Owner confirmation and manuscript/dissertation repository links | Confirmed |
| D-002 | 2026-08-03 | Scientific-validity audit precedes optimization | Document/code mismatches and external review priorities | Confirmed |
| D-003 | 2026-08-03 | Rust is an evidence-gated optimization, not a predetermined rewrite | PH already uses compiled Ripser; bottleneck is unmeasured | Confirmed |
| D-004 | 2026-08-03 | Confidential reviews remain under Git-ignored `docs/private/` | Explicit editorial restriction | Confirmed |
| D-005 | 2026-08-03 | Both Dr. Rouchka and Dr. Mistry will receive project credit | Owner instruction | Confirmed |
| D-006 | 2026-08-03 | Use sample-level matched baselines for the primary comparison | Dissertation/code show a cell-versus-sample comparability risk | Proposed pending audit |
| D-007 | 2026-08-03 | Provisional paper framing is integration auditing rather than universal PH superiority | Better fit to evidence, literature, and external feedback | Proposed pending results |
| D-008 | 2026-08-03 | Use existing project data first; add new data only after a documented evidence-gap decision | Owner instruction; controls scope and avoids premature data expansion | Confirmed |
| D-009 | 2026-08-03 | Keep the dissertation and preprint PDFs local and out of Git; preserve `example_run.r` untracked until its usefulness is assessed | Owner instruction and uncertainty about the example's historical role | Confirmed |
| D-010 | 2026-08-03 | Fast-forward local `main` to verified `origin/main` at `3910b15` and preserve the audit foundation on a `codex/` branch | Fresh divergence/collision verification and owner authorization | Completed |
| D-011 | 2026-08-03 | Treat `example_run.r` as a likely scratch launcher on hold, not a canonical example or paper-run command | No tracked history, redundant README coverage, and unresolved/placeholder inputs | Confirmed pending any contrary provenance |
| D-012 | 2026-08-04 | Preserve `MIN_CELLS = 250` for historical reproduction, identify its exclusions as pre-PH QC dispositions, and require structured provenance before optimization | Code/log/artifact reconciliation plus need for a behavioral contract for future Rust or other optimization | Implemented; threshold sensitivity remains future work |
| D-013 | 2026-08-05 | Use all active persistence-landscape levels with exact or error-controlled integration as the corrected primary target; retain `k=1` only as paper sensitivity; do not activate yet | Eight-stratum audit rejected a universal cap/fixed uniform point count; owner confirmed that the dissertation-aligned definition should govern corrected work | Project-owner approved; eligible pilot production validation complete; no package-default or biological activation |
| D-014 | 2026-08-05 | Adopt `landscape_reference_v1` as a non-default correctness oracle; do not activate it or begin Rust work | Analytical and representative benchmarks show exact/adaptive agreement and compact output, while large-diagram certification and independent exact-distance validation remain open | Implemented locally; production engine decision pending |
| D-015 | 2026-08-05 | Classify all nine audited historical persistence-diagram artifacts and every derived distance/result as scientifically ineligible; allow only labeled historical reproduction or performance stress use | Current and historical code pass feature-by-cell matrices directly to ripserr's row-as-point interface, contradicting the dissertation's cell-by-cell intent | Confirmed by code, documentation, dissertation, and artifact signatures; corrected rerun required |
| D-016 | 2026-08-05 | Retain the R reference as oracle; reject Persim 0.3.8's built-in norm; advance corrected Persim critical-pair integration only as a production prototype; keep fixed-grid GUDHI and Rust out of the primary path | R exact, corrected Persim, and SciPy agree to floating precision; built-in Persim fails sign-crossing cases; controlled scaling identifies a promising route | Accepted for immutable eligible MV-04 pilot bundles only; no public-default activation |
| D-017 | 2026-08-05 | Expand the prospective method into distinct cell-topology and gene-topology views; keep cell topology confirmatory and gene topology secondary until independently validated | Project-owner direction following discovery of the legacy feature-as-point orientation; the two orientations represent different scientifically meaningful objects | Confirmed; MV-02 typed implementation complete, scientific eligibility pending MV-03 |
| D-018 | 2026-08-05 | Evaluate multiple methods in a staged three-layer hierarchy: within-sample geometry, between-sample topological distance, then sample clustering | Prevents conflation of distances and an unstructured method/parameter search | Confirmed planning rule |
| D-019 | 2026-08-05 | Admit cell/gene fusion only after both views pass independent gates; begin with normalized equal-weight distance fusion and require comparison against both components before advanced multiview methods | Preserves interpretability, controls scope, and prevents outcome-informed fusion tuning | Confirmed planning rule; fusion remains exploratory |
| D-020 | 2026-08-05 | Freeze `cell_topology_v1` for the pilot as 384 matched cells represented in 30 shared PCs from a common 500-gene panel, with Euclidean cell geometry | Matches dissertation intent, controls cell count, avoids raw high-dimensional Euclidean geometry, and keeps both views on one source matrix | MV-02 implementation/tests complete; MV-03 eligible-data feasibility required |
| D-021 | 2026-08-05 | Freeze `gene_topology_v1` for the pilot as the same 500 genes under Euclidean correlation-chord geometry across the same 384 cells | Creates an intentional metric gene-coexpression view while isolating orientation from feature/cell changes | MV-02 implementation/tests complete; secondary view only; MV-03 feasibility required |
| D-022 | 2026-08-05 | Use complete Vietoris-Rips H0/H1 for the primary pilot; prohibit sample-specific corrected thresholds; retain all-level H0/H1 landscape L2 and PAM with label-free stability-selected `k` as the primary downstream hierarchy | Resolves filtration comparability, preserves the approved landscape definition, and prevents arbitrary-distance k-means/oracle-label tuning | G-MV1 technical freeze; resource-gated in MV-03 through MV-05 |
| D-023 | 2026-08-05 | Freeze a deterministic 14-sample existing-data feasibility manifest, five 384-cell seeds, and staged time/memory/storage stop rules | Covers eight tissues, historical cell-count extremes, and bone-marrow assay approaches without presenting the pilot as a biological validation set | Manifest generated and hashed; fresh-QC reconciliation required before PH |
| D-024 | 2026-08-05 | Admit only typed `cell_topology_v1`/`gene_topology_v1` objects to the corrected PH entry point; keep reduced fixtures scientifically ineligible and stamp the current pipeline as `legacy_gene_view_v0` | MV-02 analytical, orientation, metric, provenance, cache, permutation, and PH fixtures pass; typed and legacy result-key families are disjoint | G-MV2 technical implementation complete; no biological activation |
| D-025 | 2026-08-05 | Use the corrected persistent Persim critical-pair batch engine for immutable MV-04 all-level exact H0/H1 pilot distances; keep components separate, raw combination secondary, normalization fit-scope-bound, and Rust rejected | 408 primary distances pass matrix, eligible R-oracle, analytical, deterministic, and capped resource checks; sensitivity methods are complete only where feasible and otherwise technically excluded | G-MV4 technical gate complete; advance only to MV-05 statistical-plan design |
| D-026 | 2026-08-06 | Freeze the prospective MV-05 benchmark at the sample unit with LOSO study blocks, separate biological/technical axes, cell H0/H1 versus matched cell-energy, secondary gene topology versus gene-correlation, label-free PAM stability, and blocked uncertainty/multiplicity | Existing metadata support 90 candidate samples across five multi-study tissues; the pilot supports no cross-study tissue endpoint; legacy evaluation is label-informed and cell-replicated; current integration is transductive | Plan-freeze subgate accepted; biological execution and G-MV5 remain open until inductive integration and staged benchmark gates pass |
| D-027 | 2026-08-06 | Implement outcome-label-closed immutable LOSO manifests, exact matched baseline formulas, and a separate reference-only Seurat transfer-projection route; never fall back to transductive integration | All 18 folds validate; analytical and scientific-shape baseline fixtures pass; two synthetic inductive mappings each find 111 anchors and produce identical embedding hashes while preserving reference identity | MV5-A/MV5-B technical gate complete; advance only to one-tissue/all-seed MV5-C feasibility |
| D-028 | 2026-08-06 | Preserve training-only feature selection while separating held-out cell-projection eligibility from gene-topology eligibility; retain exact all-level H0/H1 and prohibit naïve MV5-D execution until caching, chunking, and pair scope pass resource gates | Ten real fold/seed jobs complete cell topology in both representations; SCT gene is complete in 5/10 jobs and structurally unavailable in the one-reference fold; 750 exact landscape distances and 85 matrices are label-closed; naïve full projection exceeds both caps | MV5-C technical feasibility complete; G-MV5 open; advance only to MV5-C2/MV5-D engineering specification |
| D-029 | 2026-08-06 | Cache SCT once per sample–seed, retain exact all-active-level H0/H1, and narrow primary execution to deterministic held-out-query/training-reference pair chunks; defer tasks requiring complete matrices | Six real caches and both cached LOSO directions reproduce MV5-C exactly; 250/250 chunked distances match exactly and 50/50 chunks resume immutably; full all-view lower bound is 25.82 h while SCT cell-primary is 18.68 h | MV5-C2 complete; conditionally advance only to label-closed MV5-D0 SCT cell normalization, then stop and reproject; G-MV5 open |
| D-030 | 2026-08-07 | Require per-sample source shards, runtime-complete v2 normalization identities, and matrix-only SCT payloads before MV5-D0 production; retain legacy v1 caches only as historical evidence | Monolithic source deserialization reaches ~45.74 GiB RSS; six bounded SCT entries pass at 3.361 GiB maximum and project to 3.009 GiB storage; current-runtime v2 rebuilds are byte-identical while legacy exact equivalence fails despite identical inputs | MV5-D0 production blocked pending compliant 90-sample source migration; landscape definition unchanged; G-MV5 open |
| D-031 | 2026-08-07 | Resolve MV5-D0 source migration with the existing per-sample Seurat files, accept 90 canonical sparse raw shards and 450 runtime-complete matrix-only SCT caches, then stop before folds | 53/53 recovered monolithic comparisons are exact; 37 missing samples recovered; 450/450 caches independently validate at 2.562 worker-hours, 1.81 GiB peak RSS, and 2.992 GB storage; SCT cell-primary lower bound falls to 10.525 h | MV5-D0 stage 1 complete; owner review required before MV5-D1 label-closed fold execution; outcomes and G-MV5 remain closed |
| D-032 | 2026-08-07 | Execute MV5-D1 with a strictly training-derived 500-gene panel, training-only scaling/PCA, and training-mean mapping for held-out absent features; stop before PH | Four corrected pilot folds pass, one repeats byte-identically, and 75/75 production caches with 6,750 cell views independently validate at 2.376 worker-hours, 1.95 GiB peak RSS, and 0.895 GB; zero downstream jobs | MV5-D1 complete; 8.510 h is an incomplete known-components lower bound because production cell PH is unmeasured; MV5-D2 bounded PH profiling next; outcomes and G-MV5 remain closed |
| D-033 | 2026-08-07 | Profile corrected SCT cell PH on a deterministic 30-view label-closed pilot spanning all folds, seeds, held-out/training roles, and mapping strata; require full-view MST, exact-repeat, and reduced cross-engine checks; stop before full PH and landscapes | 30/30 diagrams, 30/30 H0 MST checks, 5/5 byte-identical repeats, and 10/10 Ripserr/GUDHI H0/H1 checks pass; 110.554 worker-seconds, 228.3 MiB peak RSS, and 0.907 MB storage; projected full PH is 6.752/7.135/10.187 h under median/P90/maximum assumptions | MV5-D2 complete; combined SCT cell-primary projection is 15.262/15.645/18.697 h, below the 21.6-h planning cap; owner review required before a separate full 6,750-view PH cache sprint; landscapes, outcomes, and G-MV5 remain closed |
| D-034 | 2026-08-07 | Execute the complete 6,750-view SCT cell-PH cache in immutable 90-view fold-seed groups; require all-record validation, stored and independently recomputed H0 MST evidence, exact full-group repeat, safe resume, and zero downstream counters | 75/75 groups and 6,750/6,750 H0/H1 records pass; 6,750 stored and 75 fresh MST checks pass with zero error; 90/90 repeated records are object- and byte-identical; measured PH is 1.047 worker-hours, 273.1 MiB peak RSS, and 196.3 MB | MV5-D3 complete; measured-plus-projected SCT cell-primary total is 9.556 h with 12.044 h cap margin; separately specify the exact all-level landscape-distance stage; outcomes and G-MV5 remain closed |
| D-035 | 2026-08-07 | Execute only frozen held-out-query to training-reference SCT cell landscape pairs with separate H0/H1, all active levels, exact critical-pair integration, explicit essential-H0 exclusion, immutable chunks, R-oracle admission, and timing-separated scientific output | 35,350 biological pairs / 70,700 rows in 360 chunks pass; four eligible R exact checks agree within 1.42e-14; 4/4 repeated v3 distance files are byte-identical; 1.165 worker-hours, 349.9 MiB peak RSS, and 65.1 MB private output; superseded timing-embedded pass is scientifically identical for all rows | MV5-D4 complete; all cell-primary precomputation is measured at 7.150 h with 14.450 h cap margin; label-closed training-fitted bundle/baseline assembly next; outcomes and G-MV5 remain closed |
| D-036 | 2026-08-08 | Assemble immutable SCT cell retrieval inputs with separate raw H0/H1, descriptive raw composite, matched same-cell energy, and same-panel training-standardized pseudobulk; do not fit topology component scales from held-out-query rows | 75/75 groups, 35,350 biological pairs, 176,750 unique ranked rows, and 375/375 method groups pass; every topology row and 450 baseline oracles validate within 1.14e-13; admission, 150-file resume, and six-file public assembly are byte-stable; 267.9 MiB peak RSS | MV5-D5 complete; prediction-locked MV5-E retrieval endpoint evaluation next; clustering, integration, gene topology, fusion, and G-MV5 remain closed |
| D-037 | 2026-08-08 | Accept the prediction-locked SCT cell retrieval evaluation and retain its null/negative primary result without tuning, replacement, or tissue selection | The frozen D5 hash, 75 bundles, and 375/375 completions passed before label opening; 2,250/2,250 endpoints were estimable; H0−energy MRR was −0.00396 (95% blocked CI −0.0748 to 0.0627) and H1−energy was −0.00819 (−0.1085 to 0.0569), with Holm p=1 for both; independent reconstruction and a byte-identical repeat passed | MV5-E complete; corrected SCT topology does not outperform matched energy on the frozen primary endpoint; tissue heterogeneity is diagnostic only; next stage must be separately specified without outcome-driven method changes |
| D-038 | 2026-08-08 | Reconstruct inductive integration from accepted raw shards using exact D1 panels and D0-selected cells; use query-specific active subsets without replacement; require reference immutability and a measured resource gate before full execution | Four outcome-blind structural folds completed 360 coordinate views and 38 mappings at 2.67–3.18 GiB; all independent checks and a byte-identical repeat passed; reserved complete projection is 14.963 worker-hours and 1.537 GB | MV5-F complete; separately authorize only the 75-group integrated coordinate cache next; full execution, PH, outcomes, clustering, gene topology, fusion, and new data remain closed |
| D-039 | 2026-08-09 | Execute all 75 fixed-D1-panel integrated coordinate groups under the label-closed MV5-F contract; require independent validation, exact complete-group repeat, immutable resume, and measured reprojection before PH | 75/75 groups, 6,750 coordinate views, and 450 mappings pass 675/675 independent checks; repeat RDS SHA is `eacde0ed…7d7`; resume rebuilt zero groups; 9.107 measured coordinate hours and 12.379 reserved total hours pass all caps | MV5-G complete; separately authorize integrated cell H0/H1 PH only; landscapes, distances, retrieval, outcomes, clustering, gene topology, fusion, and new data remain closed |
| D-040 | 2026-08-09 | Give MV5-G coordinates their own typed integrated cell-view identity and execute complete Euclidean VR H0/H1 before any landscape construction | 75/75 groups and 6,750/6,750 records pass all independent identity/file/diagram/stored-MST checks plus 75 fresh MSTs; a 90-view group is byte-identical; resume rebuilt zero groups; measured PH is 1.098 h, 274.9 MiB peak, and 179.9 MB | MV5-H complete; separately authorize exact all-active-level integrated H0/H1 landscape distances; outcomes/clustering/gene/fusion/new data remain closed |
| D-041 | 2026-08-09 | Apply the dissertation-aligned exact all-active-level persistence-landscape definition to the corrected MV5-H integrated cell diagrams, with separate H0/H1, explicit essential-H0 exclusion, no universal cap or grid, and directed label-closed query-to-training pairs only | 70,700/70,700 rows in 360 chunks pass all-record identity/numerical/level accounting; four exact R oracles agree within 7.11e-14; 14/14 maximum-group files repeat byte-identically; all 720 production files survive zero-rebuild resume; measured landscape work is 0.576 h, 243.9 MiB peak, and 1.255 GB including staging | MV5-I complete; separately authorize label-closed integrated retrieval-input assembly; outcomes/clustering/gene/fusion/new data remain closed |
| D-042 | 2026-08-10 | Assemble a separate immutable integrated cell retrieval-input family with raw H0/H1, descriptive raw composite, energy recomputed on the exact MV5-G coordinates, and the shared training-standardized pseudobulk context; prohibit held-out-derived component scaling and all evaluation | 75/75 groups, 35,350 biological pairs, 176,750 unique ranking IDs, and 375/375 method groups pass; all topology rows and 450 direct baseline oracles validate within 9.81e-13; representative/maximum groups, all 150 resume files, and all six public artifacts are byte-stable; 287.8 MiB peak RSS | MV5-J complete; authorize only a separately specified prediction-locked integrated retrieval evaluation; comparison with SCT, clustering, gene topology, fusion, new data, and manuscript claims remain closed |
| D-043 | 2026-08-10 | Accept the prospectively locked integrated-cell retrieval evaluation without reading the SCT outcome, and retain its null/negative H0/H1 results without tuning, replacement, or tissue selection | All 2,250 endpoints were estimable; H0-energy macro MRR was -0.005159 (95% blocked CI -0.140977 to 0.164341) and H1-energy was -0.040203 (-0.165125 to 0.062618), with Holm p=1 for both; 11/11 independent checks and all 15 byte-identical repeat comparisons passed after a documented conservative numerical-boundary correction | MV5-K complete; only a separately prespecified comparison of the already locked SCT and integrated results may proceed next; clustering, gene topology, fusion, new data, and claim promotion remain closed |
| D-044 | 2026-08-10 | Compare accepted SCT and integrated retrieval through paired topology-minus-energy difference-in-differences, with an explicit known-marginals/pre-joint timing disclosure and exact shared-pseudobulk identity control | The contract was committed at `b3f7e28` before endpoint joining; 2,250/2,250 pairs and 450/450 pseudobulk controls pass; H0 DID MRR was -0.001195 (95% CI -0.136724 to 0.167949) and H1 DID was -0.032016 (-0.114822 to 0.052640), Holm p=1 each; 11/11 independent categories and 13/13 byte repeats pass | MV5-L complete; no favorable representation DID supported; audit the highest-value remaining benchmark axis before authorizing clustering, technical mixing, robustness, gene topology, fusion, or new data |
| D-045 | 2026-08-10 | Select the next benchmark axis with a frozen no-outcome gate over scientific value, reviewer relevance, validity, artifact readiness, resources, and outcome-selection safety | Label-free clustering contract/resource gating scored 45, robustness 43, external validation 40; technical mixing failed the current validity gate because study is nested within tissue and LOSO queries lack same-study comparators; 262,675 training pairs and 525,350 H0/H1 rows per representation are missing; all seven independent checks and seven byte repeats pass | MV5-M complete; authorize only MV5-N training-only PAM/held-out-medoid contract, exact pair inventory, and bounded label-closed resource admission; no full matrices or outcomes |
| D-046 | 2026-08-10 | Freeze label-closed clustering as training-only PAM with five-seed one-SE stability selection, deterministic canonical labels, immutable held-out medoid assignment, and average linkage as the sole sensitivity; retain separate exact all-active-level H0/H1 landscapes and matched energy/pseudobulk baselines | Exact inventories reproduce 262,675 training pairs and 525,350 H0/H1 rows per representation; 384/384 admitted landscape rows, 12/12 independent R exact oracles, 12/12 byte repeats, immutable resume, and 384/384 admitted baseline rows pass; conservative total production projection is 16.117 worker-hours including reserve | MV5-N complete; full matrices remain closed pending a separate prospective MV5-O execution specification; all biological and technical labels, cluster outcomes, method selection, and claims remain closed |
| D-047 | 2026-08-10 | Prefreeze complete training-matrix production as 150 deterministic fold-seed-representation groups, 4,340 group-cached landscape chunks, 150 energy units, and 75 shared-pseudobulk units, all bound to immutable source/implementation roots, atomic status/output schemas, caps, validation, repeat/resume, and non-retrying abort rules | 19/19 independent checks and 7/7 byte-repeat artifacts pass; the exact landscape runner reproduces 32/32 accepted real-diagram distances, one independent R oracle, clean byte repeat, and immutable resume; 16.117 h and conservative 1.278 GB projections pass 21.6-h/10-GiB caps; production/clustering/outcomes remain zero | MV5-O prefreeze complete; authorize a separate goal to execute only the frozen label-closed distance queues, stopping before production clustering or label opening |
| D-048 | 2026-08-10 | Accept complete label-closed training-distance production while preserving separate exact all-active-level H0/H1 landscapes, representation-specific energy, shared pseudobulk, immutable source/unit identities, and zero downstream clustering/outcomes | All 150 groups, 4,565 units, 1,838,725 values, and 525 matrix components pass; 12/12 exact R oracles, 66/66 frozen repeat outputs plus one supplemental pseudobulk output, and all 4,565 immutable resume units pass; observed work is `12.044379` h with `492163072` B maximum RSS and `4570070656` B private root. The 1.278-GB output-focused prefreeze estimate omitted complete group-local interval staging, but a disclosed `6.011`-GB reserve forecast and final observed storage remain below 10 GiB | MV5-P complete; authorize only a separately scoped label-closed clustering-artifact sprint under the frozen MV5-N contract, stopping again before labels, ARI/NMI, or any biological/technical outcome |
| D-049 | 2026-08-14 | Accept published `main` at `d0192d35` as the canonical audited continuation baseline while preserving all nine publication branches and keeping releases, binaries, DOI actions, defaults, new calculations, and manuscript claims closed | PRs #111-#119 merged in the authorized order; exact final tree `188e24a7`; main R run `31768108363` passed in 22m19s; private/excluded-artifact guard passed | Publication merge technically complete; phase gates still require their own scientific/owner dispositions |
| D-050 | 2026-08-14 | Begin MV-06 only as a bounded label-closed four-stratum fusion feasibility sprint; require a separate matched gene-view scale-up/resource gate before blocked evaluation | MV-04 supplies matched cell/gene H0/H1 matrices for two 10-sample and two 4-sample strata, while the 90-sample MV-05 benchmark lacks full gene-view matrices | MV6-A authorized; blocked outcomes, advanced fusion, new data, defaults, claims, and release actions remain closed |
| D-051 | 2026-08-14 | Accept MV6-A as technical evidence that corrected cell/gene pilot geometries are nondegenerate and nonredundant enough to justify a prospective scale-up/resource gate, without claiming utility | All eight source hashes, 16 scales, 102 pairs, 510 weight rows, 12 correlations, 84 neighbor rows, 44 matrix hashes, 12 independent categories, and 11 byte-identical repeat files pass; 10-sample composite correlations are near zero with limited neighbor overlap | MV6-B resource/eligibility inventory authorized; full gene calculation, blocked outcomes, advanced fusion, defaults, claims, and release remain closed |
| D-052 | 2026-08-14 | Stop MV6-B before gene PH or blocked fusion because neither accepted representation supplies a complete strict 6,750-instance matched gene axis; require the owner to select a revised gene estimand or explicitly narrow MV-06 to pilot-only gene/fusion evidence | All 150 accepted private records verify; SCT has 71 incomplete held-out views across 31/75 groups plus 379 unresolved held-out variance checks; MV5-G stores cell coordinates rather than an integrated gene-expression payload; 12/12 independent categories pass | G-MV6 remains open; global-core rerun or cell-focused disposition are the two recommended directions; MV6-05/MV6-06 and MV-07 freeze await owner input |
| D-053 | 2026-08-14 | Pursue the owner-selected global-core direction with one label-closed, explicitly transductive technical 500-gene panel across all existing samples/seeds; retain training-only fold transforms and stop before PH | All 450 accepted caches verify; 2,536 unique genes pass presence/finite/nonconstant/technical rules, margin 2,036; seed top-500 overlap is 491–493 with rank correlations >0.999; independent validation passes 12/12 | MV6-C passes; authorize only a separately prefrozen bounded matched-SCT source/PH/landscape resource profile; full 6,750-view rerun, blocked fusion, outcomes, and G-MV6 remain closed |
| D-054 | 2026-08-14 | Retain Option A and both topology orientations, but require production-style batched/streamed landscape acceleration before a full matched rerun; do not optimize PH | Five fold sources, 20 cell/gene PH records, ten all-level landscape pairs, 7/7 byte repeats, and 12/12 independent categories pass; maximum PH projections total 5.86 h while independent per-pair landscapes inflate the complete maximum to 837.83 h; prior grouped cell execution shows the projection is not intrinsic | MV6-D completes as `revise_for_targeted_acceleration`; authorize only a separately prefrozen MV6-E batched landscape admission/resource sprint; full rerun, fusion, outcomes, and G-MV6 remain closed |
| D-055 | 2026-08-14 | Admit grouped Persim and the accepted Rust exact landscape kernel; prefer Rust only for explicit hash-verified private WSL production and retain grouped Persim as the canonical portable fallback | Both engines pass 20 R references, 180 cross-engine rows, 20 reverses, 40 self zeros, byte repeats, resume, and resource guards; maximum landscape projections are 42.814 h versus 0.873 h (49.015x) | MV6-E completes as `admit_both_with_rust_preferred_private`; authorize only a separately prefrozen MV6-F full matched source/PH/landscape production stage; public defaults, binary distribution, blocked fusion, outcomes, and G-MV6 remain closed |
| D-056 | 2026-08-14 | Freeze the full matched dual-view production as 75 atomic fold-seed groups with streaming views, 13,500 typed PH jobs, 141,400 exact four-component Rust rows, no silent fallback, and staged admission | All input/Rust hashes pass; the queue independently reconstructs; 14 sources bind to one implementation root; 10/10 validations, 19/19 focused expectations, and 6/6 byte repeats pass with zero production | MV6-F prefreeze completes as `prefreeze_pass_stage1_only`; authorize only an externally monitored maximum-group run, repeat, oracles, and resume before the other 74 groups, blocked fusion, or outcomes |
| D-057 | 2026-08-14 | Accept MV6-F maximum-group stage 1 and authorize only the remaining 74 frozen label-closed production groups under the unchanged atomic/resource/abort contract | Primary and repeat finish in 344.799/356.902 s at 3.045/3.107 GB peak process-tree RSS; 3/3 scientific artifacts repeat byte-identically; 5/5 resume identities remain unchanged; R and grouped Persim each pass 12/12 balanced cell/gene H0/H1 depth-stratified oracles with no level cap | MV6-F stage 1 completes as `pass_stage2_label_closed_only`; fusion, clustering, outcomes, public Rust adoption, release, and claims remain closed until stage 2 is complete and independently validated |
| D-058 | 2026-08-14 | Rebind MV6-F to one remediated implementation root before stage two, because the accepted P1 integrity fixes changed two frozen dependencies and the runner correctly rejects mixed roots | The original queue/Rust identities remain fixed; a corrected 23-file root `5a1258e8…8d292`, serial live-cap monitor, independent validator, 23/23 focused expectations, deterministic artifacts, and complete package-aware suite pass; an earlier root was quarantined after its resume checker caught unquoted WSL paths before execution | Authorize only a remediated maximum-group reexecution, repeat, scientific equivalence, R/Persim oracle, and resume gate; the other 74 groups remain closed until that admission CSV passes |
| D-059 | 2026-08-14 | Admit the exact remaining 74 MV6-F groups under the corrected remediated root and serial fail-closed monitor | Primary/repeat finish in 337.429/349.147 s at 3.100/3.030 GB peak tree RSS; 180 diagrams and 6,500 landscape rows are scientifically identical to the parent root; 3/3 byte repeats, R 12/12, Persim 12/12, resume 5/5, and 8/8 rebind categories pass | Stage two is authorized label-closed only; fusion, clustering, outcomes, public Rust adoption, release, and claims remain closed until complete-production validation and resume pass |
| D-060 | 2026-08-14 | Stop stage two after the first pending group breaches the frozen 8-GiB RSS cap; preserve the failure and prefreeze one unchanged-runner diagnosis at the already reserved 12-GiB aggregate ceiling | Group order 2 reached 8,747,204,608 B after 228.451 s and was killed with no published artifact; no later group launched; the exact monitor/runner/queue/root/Rust identities are bound and automatic retry remains false | Authorize only the one-group 12-GiB diagnostic after quarantining partial state; if it fails, optimize/profile PH before further production; H1 and all scientific scope remain unchanged |
| D-061 | 2026-08-14 | Accept the exact-group resource diagnosis and replace the 8-GiB per-group execution ceiling with a separately bound serial 12-GiB policy for the unchanged 74-row stage-two queue | The diagnosed group completed atomically in 423.050 s at 9,575,215,104 B peak process-tree RSS; its full directory validates under the unchanged queue `f5471633…10bb5`, implementation root `5a1258e8…8d292`, scientific runner, and Rust library `51d3fca4…160d`; the prospective policy reproduces all 74 rows in exact order, binds its driver/monitor, retains one worker/no retry/1,800-s/10-GiB storage caps, and keeps labels and all downstream jobs closed | Authorize serial stage-two execution only; stop at the first group above 12 GiB or any contract/artifact failure; complete-production validation, immutable resume, fusion, clustering, outcomes, and claims remain closed |
| D-062 | 2026-08-14 | Accept complete MV6-F label-closed production under the serial 12-GiB policy, while requiring a separately committed whole-corpus immutable-resume gate before fusion | All 74 stage-two rows completed or were validated reuse; actual monitored work was 21,538.531 s (5.983 h), peak RSS was 9,575,215,104 B, 11 groups exceeded the former 8-GiB ceiling, and retained private state was 624,237,551 B; the original monitor validated/reused 74/74 rows; an independent validator passed 6/6 categories across 75 groups, 35,350 biological pairs, 141,400 component rows, and balanced cell/gene H0/H1 totals | Authorize only the prefrozen 375-artifact plus canonical-metrics immutable-resume check; fusion, clustering, labels, outcomes, defaults, release, and claims remain closed until it passes |
| D-063 | 2026-08-14 | Accept MV6-F whole-corpus immutable resume and close matched dual-view production before any fusion design | The committed checker reran the original canonical monitor and preserved SHA-256, byte size, and modification time for all 375 scientific group artifacts plus the canonical metrics (376/376 rows); no group rebuilt, no partial state appeared, and labels/downstream jobs remained closed | MV6-F is complete; authorize only a separately committed complete-corpus label-closed MV6-G fusion prefreeze that fixes components, fold-local scaling, weight panel, blocked estimands, nulls, multiplicity, and label firewall before opening outcomes |
| D-064 | 2026-08-14 | Freeze complete-corpus fusion around training-only median component scales, nine fixed rankings, equal cell/gene fusion as the sole primary, two comparator contrasts, blocked inference, and a prediction-before-label firewall | The accepted 75-group corpus reconstructs 262,675 training pairs, 1,050,700 component rows, 300 scales, and 318,150 query rankings; implementation root `ab7039f1…31db` binds five sources; 12/12 independent categories, 16/16 focused expectations, and 10/10 byte repeats pass under queue root `f5471633…10bb5`; labels and downstream jobs remain zero | MV6-G prefreeze completes as `prefreeze_pass_stage1_training_scale_sentinel_only`; authorize only the maximum-group label-closed scale/ranking sentinel under 1,800-s, 12-GiB, 5-GiB, one-worker, no-retry caps before full production |
| D-065 | 2026-08-14 | Bind the maximum-group MV6-G sentinel to one atomic scientific runner, external resource monitor, independent formula/scale reconstruction, R/Persim oracles, clean repeat, and immutable resume before execution | Ten sources form stage-one root `883bbe32…16e2`; exact scope is 2,080 training pairs, 8,320 component rows, four scales, 1,625 query pairs, and 14,625 rankings; 11/11 launch categories, 17/17 focused expectations, Python parsing, and 2/2 launch byte repeats pass; production remains zero | Authorize monitored primary plus one clean repeat only; full label-closed production remains closed until every numerical, deterministic, resource/projection, and resume gate passes |
| D-066 | 2026-08-15 | Invalidate the first MV6-G launch root after its resume checker mishandles a spaced WSL path; preserve unchanged attempt artifacts and bind an argument-safe `processx` correction before reexecution | Attempt-one primary/repeat and all scientific/oracle gates pass, but the resume child exits on argument count before reuse validation; all five hashes/bytes/mtimes remain unchanged; corrected root `6a76a11d…ce82` changes only the checker plus its regression and passes 11/11 launch categories and 2/2 byte repeat | Quarantine attempt one; authorize a clean primary/repeat rerun under the corrected root and require the complete original admission ladder again; full production and labels remain closed |
| D-067 | 2026-08-15 | Accept the clean corrected-root MV6-G maximum-group sentinel and authorize only a separately prefrozen serial label-closed completion policy | Primary/repeat complete in 221.397/227.835 s at 166.314/166.302 MB peak RSS; 3/3 scientific artifacts repeat byte-identically; all four scales and 14,625 formulas/ranks reconstruct; R and Persim pass 12/12 depth-stratified cell/gene H0/H1 oracles spanning 19–1,856 active levels; 10/10 categories and 5/5 immutable resume pass; conservative total projection is 4.747 h and ~0.605 GB | MV6-G stage one completes as `pass_stage2_label_closed_only`; prefreeze the unchanged remaining 74 groups before execution; metadata, outcomes, advanced fusion, clustering, and claims remain closed |
| D-068 | 2026-08-15 | Freeze a dynamic MV6-G group runner for the 65–89-training-sample range, but authorize only maximum-group byte equivalence before serial production | The exact remaining workload is 74 groups, 260,595 training pairs, 1,042,380 component rows, 33,725 query pairs, and 303,525 rankings; root `9bf8614d…2a71c` binds eight sources; 9/9 categories and 3/3 prefreeze repeats pass under unchanged caps and label firewall | Execute only the accepted maximum group through the general runner and require 3/3 scientific byte identity plus resource pass; remaining production and labels stay closed |
| D-069 | 2026-08-15 | Admit the dynamic MV6-G runner after exact maximum-group rebind, while keeping the remaining 74 groups closed until a prospective serial execution/resume policy is committed | The monitored general runner exited zero in 226.340 seconds at 169,119,744 B peak RSS and 8,064,620 B private output; `training-distances.csv`, `scales.csv`, and `rankings.csv` reproduce the accepted sentinel byte-for-byte with labels closed and zero outcome/fusion jobs | Prefreeze the serial driver, aggregate caps, complete validator, and immutable-resume checker next; do not launch remaining production until that gate passes |
| D-070 | 2026-08-15 | Authorize the remaining 74 MV6-G groups only through a separately hash-bound serial driver with complete-corpus and immutable-resume gates | Execution root `38440b86…0091` preserves scientific root `9bf8614d…2a71c`; exact workload and 12-GiB/1,800-s/5-GiB/12-worker-hour caps pass 10/10 independent categories, 3/3 byte repeats, 7/7 focused expectations, and the 1,592-test complete suite; production state is empty | Execute in exact order with one worker and no retry; stop on the first failed cap/identity/artifact; then require 75-group validation and a separately committed 445-file immutable-resume pass before opening labels |
| D-071 | 2026-08-15 | Preserve the first MV6-G production failure and rebind only the execution monitor to retain child stdout/stderr before one diagnostic rerun | Execution order 2 exited 1 after 211.699 s at 169,197,568 B with no scientific output or cap breach; all source records contain finite H0 and H1; old pipes discarded the child error; corrected execution root `0c75e854…ee2b` preserves scientific root `9bf8614d…2a71c`, binds both log hashes, and again passes 10/10 validation, 3/3 repeat, 7/7 focused, and 1,592-test complete gates | Authorize exactly one unchanged first-group diagnostic rerun; no later group may launch until the captured cause is assessed and a prospective continuation decision is committed |
| D-072 | 2026-08-15 | Classify the reproduced first-group failure as an execution-queue schema mismatch and add an exact parent-runner pair-count alias without changing scientific code | Captured stderr identifies an `NA` expected-pair guard; the MV6-G workload supplies `query_biological_pairs` while the admitted runner consumes `biological_pairs`; root `97ffd576…3b38` makes them exact integers and enforces equality in builder, policy, independent validator, monitor, and regression tests; 10/10 categories, 3/3 repeats, 10/10 focused, and the 1,595-test suite pass | Authorize one clean execution-order-2 rerun; require complete scientific validation before authorizing execution orders 3–75; labels and outcomes remain closed |
| D-073 | 2026-08-15 | Accept the corrected execution-order-2 MV6-G group and authorize the remaining exact serial queue | The group completed atomically in 218.503 s at 168,214,528 B peak RSS; 2,080 training pairs, 8,320 component rows, four scales, 1,625 query pairs, and 14,625 rankings validate; independent scale/formula errors are `3.553e-15`/`9.326e-15` below `1.745e-13`; 7/7 categories and the label firewall pass | Run orders 3–75 serially with validated reuse of order 2 and no retry; then stop for independent 75-group validation and a separately committed immutable-resume gate |
| D-074 | 2026-08-15 | Stop after nine complete MV6-G groups when execution order 11 exposes a stale pre-sort block-index invariant; preserve all state and require corrected-root equivalence before restarting | The stopped group exits at the final ranking-axis check after 311.165 s and 184,246,272 B RSS; read-only reconstruction shows 81/81 blocks pass before sorting and falsely fail after stale indices are applied; the correction only rebuilds blocks after canonical sorting, with 72/81-sample regressions; root `8b0a1e42…cec0c` passes 9/9 prefreeze, 3/3 repeat, 15 focused, and 1,601 complete tests | Existing nine groups are quarantined and not reused; authorize only a clean maximum-group three-artifact byte-equivalence rebind before any production restart; labels and outcomes remain closed |
| D-075 | 2026-08-15 | Admit the corrected rank-block scientific root after maximum-group byte equivalence, while keeping production closed pending a rebuilt completion policy | Root `8b0a1e42…cec0c` completes the accepted sentinel in 207.573 s at 162,803,712 B RSS and reproduces training distances, scales, and rankings byte-for-byte (3/3); labels and all downstream counters remain closed/zero | Rebuild the serial completion policy against the corrected root, require its independent/repeat/focused/full-suite gates, and do not reuse the quarantined nine-group attempt |
| D-076 | 2026-08-15 | Authorize a clean MV6-G 74-group restart under corrected scientific root `8b0a1e42…cec0c` and execution root `deb03fbc…f745` | The rebuilt policy retains exact workload and 12-GiB/1,800-s/5-GiB/12-worker-hour caps; 10/10 independent categories, 3/3 byte repeats, 10/10 completion-focused expectations, 15/15 production-focused expectations, and the 1,601-test complete suite pass; active production state is empty | Execute all 74 groups from clean state with one worker/no retry, stop on first failure, then require independent 75-group validation and a separately committed 445-file immutable-resume pass |
| D-077 | 2026-08-15 | Accept corrected-root MV6-G complete production and independent 75-group validation, while keeping prediction lock and labels closed pending immutable resume | All 74 restart groups complete in 5.823640 worker-hours at 186,503,168-B peak RSS and 531,251,838-B final private state; the accepted sentinel plus restart corpus reconstructs 262,675 training pairs, 1,050,700 component rows, 35,350 query pairs, and 318,150 rankings; 7/7 complete checks pass; bounded cross-device publication and validation-queue schema remediations changed no scientific artifact | Authorize only the prefrozen 445-file immutable-resume checker; require hash/byte/mtime preservation before prediction-manifest prefreeze or any label access |
| D-078 | 2026-08-15 | Close MV6-G computation after the prefrozen complete-corpus immutable-resume pass and advance only to prospective prediction-manifest prefreeze | The unchanged completion driver reused 74/74 groups; all 445 bound artifacts retained identical SHA-256, byte size, and modification time; the checker exited zero in 387.3 seconds with labels closed and no biological outcomes | Authorize a separately committed prediction/outcome manifest and independent validation before any metadata read; fusion evaluation, clustering interpretation, and claims remain closed |
| D-079 | 2026-08-15 | Accept the MV6-H complete prediction lock and authorize only its fixed blocked outcome evaluation after the lock becomes durable | Root `c752408f…f4fd` binds 75 groups, 375 group files, 318,150 rankings, nine methods, two endpoints, two primary MRR contrasts, and all inference/implementation identities; independent validation passes 13/13, repeat 9/9, focused tests 21/21, and full suite 1,580 pass/4 established skips; metadata was hash-bound but never opened | Commit the lock, then execute at that exact commit with receipt-before-read and post-label rehash; advanced fusion, clustering, defaults, release, and claims remain closed |
| D-080 | 2026-08-15 | Close G-MV6 with equal-weight fusion rejected and carry cell/gene topology forward as separate views | Exact-commit blocked evaluation gives F0.5−cell MRR −0.01019 [−0.06484, 0.03601], raw p 0.835/Holm 1, and F0.5−gene +0.02686 [−0.08472, 0.08101], raw p 0.561/Holm 1; the required both-positive rule fails; independent validation passes 15/15 and all 13 outcome artifacts repeat byte-identically | Do not admit advanced/learned/tissue-specific fusion or select descriptive F0.25 despite its observed lead; advance to MV-07 robustness/confounding synthesis with complete negative reporting |
| D-081 | 2026-08-15 | Begin MV-07 with a selection-resistant evidence map and run existing-artifact confounding diagnostics before any new PH or data | Twenty-three immutable scientific/implementation sources rehash; 14 robustness and 10 confounding axes are classified; the corrected all-active-level H0/H1 landscape contract is unchanged; independent validation passes 10/10 and all 10 synthesis artifacts repeat byte-identically | Authorize only a separately committed MV7-B no-new-PH influence, retained-cell-count, and approach-stratification prefreeze; gene-panel/metric reruns, external data, method selection, defaults, and claim promotion remain closed |
| D-082 | 2026-08-15 | Narrow comparative cell-versus-gene claims after fixed no-new-PH influence/confounding diagnostics | The composite MRR contrast changes sign after deleting study SRA716608 or colon; study and tissue influence exceed 0.05, while retained-cell and mixed-study approach flags do not trigger; independent validation passes 7/7 and all 12 production artifacts repeat byte-identically | Advance to MV7-C existing-data synthesis; do not select a view, infer causal technology effects, authorize gene reruns/new data, or promote manuscript claims before the author gate |
| D-083 | 2026-08-15 | Treat existing data as sufficient for a methods-focused paper but insufficient for external-generalization or causal technology claims | MV7-C reconciles corrected landscapes, complete cell robustness, conditional gene evidence, negative fusion, and material study/tissue influence into ten explicit claim boundaries and four expansion options | Set G-MV7 to owner/author-team decision required; recommend manuscript/figure claim mapping before more PH, or a prospective external-data audit if generalization is the chosen ambition |
| D-084 | 2026-08-15 | Separate “full data” by estimand: keep the corrected 90-sample cross-study benchmark primary, admit the 34 retained single-study samples only to a prospectively specified 124-sample descriptive topology analysis, and retain the three below-250 candidates as a separate threshold sensitivity | The complete flow is 127 candidates, three explicit pre-PH exclusions, 124 retained, and 90 primary samples; all 34 descriptive-only samples support 384 cells and six depth-extreme source/SCT sentinels pass independently, while 16 approach labels disagree across metadata sources | Advance to MV7-E metadata-provenance and descriptive-fit prefreeze; do not recalculate the primary 90, run expanded PH/outcomes, change the landscape contract, or admit the three low-cell samples without later gates |
| D-085 | 2026-08-15 | Resolve approach provenance from official accession methods and activate the prespecified 124-sample global-core panel fallback before expanded PH | `Approach.y` matches public metadata 124/124 and all 16 disputed GEO records; the corrected primary 90 are all scRNA, making the old mixed-approach diagnostic not estimable; 33/34 added samples contain all 500 accepted features and one lacks only KLF2; five transductive per-seed fits, separate cell/gene views, complete VR H0/H1, 152,520 landscape components, caps, and label firewall pass 17/17 independent and 17/17 byte-repeat gates | Authorize only MV7-F's 34 raw shards and 170 SCT caches, then freeze the exact all-124 global-core panel before PCA or PH; primary 90, landscapes, labels, outcomes, and external data remain otherwise unchanged/closed |
| D-086 | 2026-08-17 | Recommend a prospectively frozen exact-500 HCA raw-read path before any external topology claim; retain common-475 only as a secondary harmonization sensitivity | The HCA Cell Ranger 3.0.0 Ensembl93 axes are exact but exclude 25 target genes by historical biotype policy; corrected common-475 PH and all-active-level H0/H1 landscapes pass complete validation; the paired comparison independently reproduces `material_panel_sensitivity`, localized to cell H1 (median Spearman 0.916 and top-10 overlap 0.702 below 0.95/0.80) | Authorize MV8-H planning only. Require owner approval for the 79.194-GiB/48-FASTQ resource commitment and a prospective custom-reference/software/QC/resource/firewall/validation contract before download or processing; labels and outcomes remain closed |
| D-087 | 2026-08-17 | Authorize the exact MV8-H 48-file HCA acquisition under a separately frozen, resumable, checksum-verified contract while keeping all raw processing and outcomes closed | Owner authorized including the download; the official Azul manifest resolves 48 unique files, eight exact units, and 85,034,239,918 bytes; the minimal target-complete Ensembl93 annotation contains 33,563 genes and all exact 500 targets; 13/13 independent prefreeze checks pass; 2.132 TB free passes the 1.5-TiB post-download reserve | Commit the prefreeze, run one small file sentinel, then resume through all 48 only after exact sentinel size/SHA-256; require complete download closure and separate exact Cell Ranger 3.0.0/reference validation before `mkref`, `count`, labels, or outcomes |
| D-088 | 2026-08-17 | Accept the exact MV8-H file-order-1 sentinel and authorize the committed downloader to resume through all 48 FASTQs | The 394,373,114-byte HCA_BM_001 L002 I1 file atomically publishes with exact official SHA-256 `4464d4ea…1dd302`; independent validation confirms size, hash, gzip identity, receipt, zero partials, cache cap, and the 1.5-TiB post-remaining reserve | Run the remaining 47 files serially/resumably; require an independent 48-file checksum closure before any reference build or raw processing; labels and outcomes remain closed |
| D-089 | 2026-08-17 | Close the exact MV8-H 48-FASTQ acquisition and authorize only frozen Ensembl-93 primary-assembly input acquisition and identity validation | The private cache contains 48/48 final files and exactly 85,034,239,918 bytes with zero partials; the downloader completion marker and empty error log pass; an independent validator reopens and recomputes all 48 SHA-256 digests, byte sizes, gzip identities, and receipt joins with no mismatch while preserving the 1.5-TiB reserve | Acquire and hash-bind the official Ensembl release-93 primary-assembly FASTA; separately validate exact Cell Ranger 3.0.0 and the custom reference before `mkref` or `count`; labels, outcomes, downstream topology, and deletion remain closed |
| D-090 | 2026-08-17 | Close the exact Ensembl release-93 primary-assembly input identity and stop before any Cell Ranger acquisition, reference build, or raw processing | The 881,214,682-byte archive has local SHA-256 `2a27436d…2406aa0`; the live official checksum file and independent local BSD `sum` both return `07119 860562`; a complete gzip stream read passes and the public bundle contains identities only | Require the owner to provide or explicitly install/accept access for exact Cell Ranger 3.0.0, then separately validate the runtime and custom reference before `mkref`; raw processing, topology, labels, outcomes, and deletion remain closed |
| D-091 | 2026-08-17 | Prospectively amend only the MV8-H counter runtime to the installed Cell Ranger 8.0.1 while retaining the admitted 3.0.0 matrices as historical comparators | The exact installed distribution (`aafd39e2…fdf3`), launcher, STAR, samtools, Python, and required CLI controls are content-bound; Ensembl-93 FASTA/source GTF, the independently repeated 33,563-gene target-complete GTF, exact-500 panel, and common-475 panel remain unchanged; six builder artifacts repeat byte-identically, independent validation passes 12/12, focused tests pass 178/178, and the complete suite passes 2,629 with four established skips | Authorize only construction of the frozen custom reference; require complete-tree/reference-feature validation before prefreezing one complete-unit count sentinel; all count, QC, topology, landscape, label, outcome, and deletion operations remain closed |
| D-092 | 2026-08-18 | Accept the conservatively built Cell Ranger 8.0.1 custom reference and authorize only a prospective one-complete-unit count-sentinel prefreeze | The owner requested lower resources where needed; the four-core/32-GiB command was committed before execution. `mkref` completed in 13,603 s with 450 monitor samples, 30.494-GiB peak RSS, 27,819,306,242-B peak run tree, and no disk/free-space breach. Independent validation passes 15/15: the retained 19-file/20,765,871,518-B reference has tree SHA `5e2aff9e…99f1c`, exact embedded Ensembl-93 FASTA/GTF, 33,563 genes, 2,565,751 feature records, all exact 500 and common 475 targets, and zero partials. An RSS composite-Boolean omission is disclosed; all samples independently pass the intended 32-GiB gate and the launcher is corrected forward-only. The locked library verifies 265/265 exact records; focused tests pass 254/254 and the complete suite passes 2,705 with two established skips | Prefreeze exactly one complete-unit count sentinel under explicit SC3Pv2/exon-only/BAM-disabled/no-secondary/resource/firewall gates; do not execute `count`, open matrices, run QC/PCA/PH/landscapes, access labels/outcomes, process remaining units, or delete artifacts yet |
| D-093 | 2026-08-18 | Close the MV8-H count-sentinel prefreeze evidence and retain all raw processing closed until owner-tokened execution | Independent validation now passes 14/14 across selection/fastrq/reference/command/resource/future contract/firewall/repeat/report/runner and non-execution gates; primary and repeat artifacts are byte-identical for all gated public files and decision/firewall text is aligned with explicit forward-only, 4 cores/32 GiB, and no universal level cap language | Continue to keep `count` execution, matrix access, remaining units, QC/PCA/PH, labels, outcomes, and deletion forbidden until explicit execution authorization and a separate owner-approved run record is added |
| D-094 | 2026-08-18 | Close the owner-authorized one-unit MV8-H Cell Ranger 8.0.1 execution as a metadata-only run record while preserving the prefreeze as a historical pre-execution record | The frozen HCA_BM_002 pipestance reports explicit successful completion at 17:47:13, whitespace-only stderr, no `_failed`/`_error` markers, exact Cell Ranger 8.0.1/SC3Pv2/exon-only/no-BAM/no-secondary/4-core/32-GiB controls, peak RSS 20,141,236,224 B, 14,348 s elapsed, 161,962,028 B final workspace, and 1,903,052,161,024 B free space; execution closure passes 12/12 and the seven public closure artifacts repeat byte-identically; matrix/QC contents were not opened | Keep matrices, expression/QC values, barcodes, labels, outcomes, landscapes, remaining units, and deletion closed; require a separate owner-approved structural/QC review and remaining-unit decision |
| D-095 | 2026-08-18 | Reconcile the strict original no-count objective with the later owner-authorized execution without rewriting history | The prefreeze remains independently closed at 14/14 with execution false; the later execution is separately closed at 12/12 with a 7/7 deterministic repeat; the reconciliation audit explicitly marks the literal no-count requirement as not achieved after the later authorization | Do not mark the original no-count goal complete; preserve the historical prefreeze and execution records separately, and treat package-build hygiene plus any structural/QC review as later gates |
| D-096 | 2026-08-18 | Bind clean-export package verification without treating the environment dependency gate as a biological or sentinel failure | Git archive at `afcebc9` contains zero `tmp/` entries and zero PDFs; canonical `R CMD build` produces `scPHcompare_0.1.0.tar.gz` (436,297 B); `R CMD check` and force-suggests-disabled check both stop at missing declared R imports/suggests before code checks; no packages were installed | Retain the clean export and dependency-gate record; resolve dependency availability in a dedicated release environment before claiming full package-check completion; keep the strict no-count goal unresolved |
| D-097 | 2026-08-18 | Revalidate the exact existing one-unit execution after owner resumed the goal and authorized count if needed, without launching a second count | Existing closure builder passes 12/12 again; all eight public execution-closure artifacts match the published hashes; resume record explicitly says second count started `no` and downstream firewall remains closed | Treat the existing one-unit record as the only execution; do not run another count unless a future plan explicitly changes the exactly-one-sentinel contract |
| D-098 | 2026-08-18 | Record the owner clarification that the no-count clause was unintended and is superseded | Owner explicitly authorized Cell Ranger count if needed and clarified that the no-count clause was not intended; the exact one-unit execution already exists and is independently revalidated 12/12 with 8/8 artifact equality | Accept exactly one completed sentinel execution as the current contract; preserve the historical prefreeze and keep matrix/QC/label/outcome/landscape/remaining-unit/deletion gates closed |
| D-099 | 2026-08-18 | Admit the completed HCA_BM_002 sentinel to metadata-only structural/QC review | New admission builder passes 8/8: seven expected outputs with no BAM/secondary, raw/filtered HDF5 shapes/dtypes, filtered 5,037 columns matching aggregate metrics, parseable Cell Ranger QC metadata, execution resource provenance, and downstream firewall; primary/repeat public artifacts match 7/7 | Next gate is an owner decision on matrix-content/QC review; labels, outcomes, landscapes, remaining units, and deletion remain closed |
| D-100 | 2026-08-18 | Prepare the allowlisted matrix-content/QC review prefreeze without opening matrix content | Prefreeze binds filtered/raw H5 plus aggregate metrics as the only potential post-approval inputs; molecule/web payloads, labels, outcomes, landscapes, remaining units, and deletion are forbidden; eight validation checks pass and eight public artifacts repeat byte-identically | Require explicit owner approval before opening any matrix/QC content; stop on drift, resource breach, unauthorized access, or nonreproducible repeat |
| D-101 | 2026-08-18 | Execute the approved one-unit matrix/QC review while preserving topology and biological firewalls | Allowlisted filtered/raw H5 values and aggregate metrics opened for HCA_BM_002 only; 8/8 checks pass; 5,037 filtered cells, 4,614 pass the frozen 384-cell QC-depth rule, median 801 detected features/cell and 2,545 counts/cell; resource cap passes; deterministic public artifact repeat is 7/7 | Next gate is an owner decision on topology review or stop; no barcodes, labels, outcomes, PCA/PH/landscapes, other units, or deletion are authorized |
| D-102 | 2026-08-18 | Bind the corrected HCA dual-view topology contract and test all five prespecified 384-cell seeds before PH execution; do not substitute or zero-pad the frozen 500-gene panel | The topology prefreeze binds 384 cells, the immutable 30-PC projection, separate Euclidean cell and correlation-chord gene views, separate H0/H1, and the dissertation-aligned all-active-level landscape definition. Label-closed SCT panel screening passes the five-seed/QC axis but retains only 497/500, 499/500, 497/500, 499/500, and 499/500 genes for seeds 20260805–20260809; every seed is therefore blocked from the exact-500 topology contract without a new panel decision | Preserve the exact-500 block and all five public seed hashes; do not run PH/landscapes or open labels/outcomes/other units. Owner decision is required on whether HCA validation should proceed under the already-planned ordered common-475 secondary panel or remain blocked pending another exact-500 strategy |
| D-103 | 2026-08-18 | Execute the approved common-475 secondary HCA topology and all-active-level landscape feasibility gate for HCA_BM_002 only | The frozen 384-cell/30-PC cell view and 475-gene correlation-chord gene view complete separate Ripserr H0/H1 PH; private diagrams are retained only in ignored `tmp`. Exact critical-breakpoint integration publishes 1,225 active levels across four view/dimension profiles (cell H0 383, cell H1 81, gene H0 474, gene H1 287), with no fixed grid or level cap; primary/repeat aggregate CSVs are byte-identical and the focused firewall test passes | Treat this as label-closed methods feasibility, not an external validation or biological result. Preserve the exact-500 block; keep fusion, labels, outcomes, remaining HCA units, and deletion closed. Next gate is a separate owner decision on external common-475 validation or an exact-500 strategy |
| D-104 | 2026-08-19 | Close the approved eight-unit common-475 HCA external technical-validation gate with separate cell/gene topology and exact all-active-level landscapes | All eight manifest-bound raw H5 inputs rehash to the frozen MV8-D identities; each unit deterministically supplies 384 cells and the ordered 475-gene panel. The run publishes 32 separate Ripserr profiles (cell/gene x H0/H1) and 32 exact critical-breakpoint landscape profiles with every active consecutive level retained. Primary/repeat scientific hashes agree after excluding only nondeterministic elapsed time; landscape artifacts are byte-identical; independent validation passes 11/11; no labels, outcomes, fusion, exact-500 recovery, or deletion was opened | Treat as label-closed external technical replication evidence, not a biological or clustering claim. Preserve exact-500 as blocked and keep fusion, labels, outcomes, and deletion closed. Next gate is an external common-475 result interpretation/claim map or an explicitly approved exact-500 recovery strategy |
| D-105 | 2026-08-19 | Bind the MV8-J provisional claim-to-evidence map before authorizing exact-500 recovery or new external analyses | The map reconciles the prior MV7-J claim inventory with MV8-I's eight-unit technical ranges: cell H0 383, cell H1 282–383, gene H0 474, gene H1 409–1,055 positive intervals; cell H1 37–66 and gene H1 85–288 active landscape levels. It classifies 14 claims as supported method/technical, conditional descriptive, not supported/evaluated, decision-required, or blocked; the independent claim-map test passes | Draft a methods-focused manuscript around corrected landscapes, separate views, reproducibility, and explicit panel/generalization limits. Run exact-500 HCA recovery only if a specific manuscript claim or figure is shown to require it; keep labels, outcomes, fusion, clustering, and deletion closed |
| D-106 | 2026-08-19 | Close the bounded HCA_BM_002 exact-500 Cell Ranger feasibility gate without opening downstream biology | Direct read-only decompression of the existing Cell Ranger 8 filtered feature/barcode axes finds 33,563 unique features, all 500 exact-panel IDs, all 475 common-panel IDs, and 5,037 filtered cells; the approved matrix/QC summary independently reports 4,614 cells passing the frozen 384-cell minimum. Existing execution metadata records 14,348 s, 20,141,236,224-B peak RSS, four cores, and 32 GiB; independent validation passes 10/10 and the public artifact test passes | Exact500 is technically feasible for HCA_BM_002. Do not process the remaining seven units or run exact500 PH/landscapes unless a later claim-specific gate authorizes it; labels, outcomes, fusion, manuscript work, and deletion remain closed |
| D-107 | 2026-08-20 | Close the owner-authorized seven-unit MV8-H exact-500 raw-read recovery as a label-closed input/runtime validation, without opening topology or biology | All seven remaining HCA units complete Cell Ranger 8.0.1 under the frozen SC3Pv2/exon-only/BAM-disabled/no-secondary/4-core/32-GiB contract. Corrected independent validation passes 68/68: every unit has a 33,563-feature axis, all 500 exact-panel IDs, all 475 common-panel IDs, at least 384 QC-eligible cells, required outputs, success markers, clean stderr, terminal performance evidence, and deterministic input bindings; the public closure records terminal `_perf`/`_timestamp` evidence where live samples were sparse and publishes no matrices, barcodes, labels, outcomes, PH, landscapes, or fusion | Treat this as exact-500 recovery feasibility and execution evidence only. Keep all matrix-derived biology, topology, landscapes, labels, outcomes, fusion, manuscript claims, and deletion closed until a separate claim-specific gate |
| D-108 | 2026-08-21 | Close the HCA_BM_002 exact-500 transform-contract audit before any exact-500 topology scale-up | Standard SCT retains 497/500 exact-panel features. Low-count SCT (`min_cells=1`) retains 500/500 and has finite standardized values plus a finite immutable 30-PC projection, but three standardized genes have zero variance across the frozen 384-cell selection, so the required 500-gene correlation-chord view is undefined. Two fresh low-count runs agree on the source-matrix SHA-256 and all gate results; no PH, landscapes, labels, outcomes, fusion, or clustering was computed | Do not run exact-500 topology under this contract. Do not silently drop, pad, or alter the three features. A future exact-500 route requires a separately approved scientific contract for cell selection, transform, or topology geometry; common475 remains the admitted external technical-validation view |
| D-109 | 2026-08-21 | Close MV8-L's prospective all-QC-cell gene-side exact-500 feasibility test before any mixed-scope topology contract | All 4,614 cells passing frozen QC were used in deterministic barcode order without panel-aware selection. Both standard SCT and `min_cells=1` retain all 500 panel features and are finite, but each leaves one exactly zero-variance standardized gene, so the unchanged correlation-chord gene distance remains undefined. Independent repeats agree; no other HCA matrix is opened | Do not adopt a view-specific exact-500 topology contract or inspect other units. Common475 remains the admitted external dual-view analysis. Any future exact500 route requires a separately approved representation or association-geometry change |
| D-110 | 2026-08-21 | Close MV8-M's prospective exact-500 gene-representation audit and identify representation-contract drift before any topology rerun | On the identical ordered 500 genes and all 4,614 frozen-QC HCA_BM_002 cells, later-production SCT `data` reproduces one constant gene and invalid exact-500 correlation-chord geometry, while SCT Pearson residuals recovered from the fitted model retain finite nonzero variance for all 500 and form a valid exact-500 distance; log-normalized RNA also passes as a diagnostic only. Corrected primary/repeat scientific summaries and common475 stability metrics are byte-identical; 12/12 checks pass in 559.74/567.16 s at 3,879,828/4,289,784 KiB. Common475 SCT-data versus residual distances are globally similar but locally non-interchangeable (Spearman 0.9268; mean top-10 overlap 0.6813) | Recommend restoring the original corrected MV01/MV03 Pearson-residual definition as the candidate gene representation, but do not adopt or silently replace MV05-D0/MV07/HCA evidence. Require owner approval and a prospective internal-124 plus external-eight migration/sensitivity prefreeze before any PH, landscapes, clustering, fusion, outcomes, or claims |
| D-111 | 2026-08-21 | Prospectively freeze the owner-approved MV8-N internal-124/external-eight Pearson-residual migration sensitivity before any rerun | Independently rehashed all 124 internal raw-count shards and eight current exact-reference HCA matrices; froze 132 all-QC standard-SCTransform fits evaluated on unchanged 384-cell axes, 1,280 proposed gene views, 1,272 planned PH records, 28 exact all-active-level landscape groups, and 40 descriptive paired comparison strata. The design separates fit-scope, representation-layer, and net migration effects while leaving cell topology immutable; 18/18 metadata-only checks pass | After commit, authorize only MV8-O's internal minimum/maximum plus HCA_BM_002 source/view sentinel under serial 12-GiB/1,800-second caps, with repeats for the maximum internal and external bridge. Keep the remaining 129 fits, PH, landscapes, comparisons, labels, outcomes, clustering, fusion, claims, and default adoption closed |
| D-112 | 2026-08-21 | Close the bounded MV8-O Pearson-residual source/view sentinel before any full migration production | Five serial all-QC standard-SCTransform source runs pass 11/11 independently rehashed closure checks: internal minimum primary; internal maximum primary/repeat; HCA_BM_002 primary/repeat. Required internal exact500 and external common475/exact500-residual correlation-chord geometries are finite and nondegenerate; HCA SCT-data exact500 remains the expected diagnostic-only invalid view. Eligible matrix/distance hashes match both repeats. Every run stays below 1,800 s and 12 GiB; the maximum internal repeats peak at 12.708 and 12.703 GB, leaving a narrow but positive margin. The documented glmGamPoi-native fallback is classified as a performance diagnostic, not an error or changed transform | Require a separate full-source-production prefreeze that explicitly decides how to manage the narrow maximum-source memory margin before the remaining 129 fits. PH, landscapes, comparisons, clustering, fusion, labels, outcomes, claims, and default adoption remain closed |
| D-113 | 2026-08-21 | Prospectively freeze serial full source production after MV8-O, before opening topology | MV8-P passes 8/8 metadata-only checks. It binds the 129 remaining sources (122 internal, seven external) to all-QC standard SCTransform plus Pearson residuals on unchanged selected-384 axes, ordered by increasing frozen-QC cell count. The largest remaining source has 9,071 cells, below the 11,475-cell completed maximum sentinel. The conservative one-worker/zero-automatic-retry/1,800-s/12-GiB policy remains; the measured maximum sentinel leaves a 176.7-MB margin | After commit, launch only the 129 serial source fits and private cache construction. Rehash every cache into a source-production closure before considering PH, landscapes, comparisons, clustering, fusion, labels, outcomes, claims, or default adoption |
| D-114 | 2026-08-22 | Implement and verify the resumable MV8-P serial source-production runner and prospective MV8-Q closure before launch | The runner consumes exactly the frozen 129-row queue in ascending cell-count order, resolves and SHA-256 checks all 122 internal and seven external sources, uses one monitored child with the existing 1,800-s/12-GiB caps, writes an atomic public resource ledger after every job, and permits resume only from a rehashed completed strict prefix. The worker prospectively freezes the seven external selected-384 axes. The MV8-Q builder independently rehashes all private caches/audits and requires 132-source coverage while preserving all topology and outcome firewalls. Focused checks pass 22/22 and the complete package suite passes 3,305 with two established optional skips | After this implementation commit, launch MV8-P into fresh private/public roots. Monitor it without competing with unrelated WSL workloads. Stop on the first source, geometry, stderr, elapsed, or RSS failure; otherwise build MV8-Q only after 129/129 completion |
| D-115 | 2026-08-23 | Preserve MV8-P's evidence-backed stop and recover jobs 124-129 through a fresh bounded overlay | MV8-P accepted 123/129 fits before job 124 stopped after 23.13 s at 2.12 GB because a 654.06-MiB `future` export exceeded the runtime's 500-MiB default; neither production cap was approached. The owner approved a method-neutral 2-GiB `future.globals.maxSize` allowance, exactly one explicit second attempt for job 124, first attempts for jobs 125-129, immutable original evidence, and merged independent closure. One worker plus the 1,800-s/12-GiB caps remain | Commit the recovery prefreeze/runner/closure support before launch. Write only fresh overlay roots, stop on the first new failure, and keep PH, landscapes, comparisons, clustering, fusion, labels, outcomes, claims, cleanup, and deletion closed |
| D-116 | 2026-08-23 | Preserve MV8-PR's accepted job 124 and second stopped receipt, then recover jobs 125-129 under a prospective 14-GiB cap | MV8-PR job 124 passed in 424.9 s at 11.52 GB with valid required geometry. Job 125 then stopped at 13.025 GB, only 140.3 MB above the frozen 12-GiB cap, without an accepted cache or audit. The owner authorized a fresh five-job overlay, one explicit retry only for job 125, first attempts for jobs 126-129, and a 14-GiB monitored cap. The one-worker, 1,800-s, 2-GiB `future` allowance, transform, axes, and geometry policies remain | Commit the second amendment before launch; preserve both prior roots, stop on the first new failure, and merge only through an independently rehashed 132-source closure. All topology and outcome stages remain closed |
| D-117 | 2026-08-23 | Close full Pearson-residual source production across the preserved MV8-P, MV8-PR, and MV8-PS roots | MV8-PS completed jobs 125-129 under the 14-GiB amendment; maximum peak was 13.306 GB at job 125. MV8-Q passes 17/17 checks after independently rehashing all 129 private cache/audit pairs, preserving both stopped receipts, confirming required internal exact500 and external common475/exact500-residual geometries, and binding the three MV8-O primary fits for 132/132 source coverage. No PH, landscape, clustering, fusion, label, outcome, or claim computation occurred | Source representation is closed. The next sprint must prospectively freeze topology production, including the PH queue/backend fallback, exact all-active-level landscape contract, resource staging, and comparison firewall before computing any diagram |
| D-118 | 2026-08-23 | Prospectively freeze corrected full topology production across all 132 source fits before computing any new diagram | MV8-R passes 20/20 metadata-only checks and binds 1,272 source-produced views, excluding eight invalid external SCT-data exact500 diagnostics. It prospectively corrects 16 common475 panel-metadata hashes without changing matrix/distance evidence and identifies that all eight historical MV8-I selected-cell axes differ from the current exact-reference axes. The corrected queue contains 1,272 gene plus eight same-axis external cell PH records; 1,240 internal cell/gene comparators remain immutable. Ripserr-primary/GUDHI-resource-fallback policy and the exact all-active-level, H0/H1-separate, no-grid/no-cap streamed landscape contract are frozen. No PH or landscape was executed | After commit, authorize only eight same-axis external baseline units and the bounded 23-record PH sentinel. Keep the remaining 1,257 PH records, all 28 landscape groups, 40 comparisons, clustering, fusion, labels, outcomes, adoption, and manuscript claims closed. Landscapes additionally require rebuilding and hash-rebinding the accepted Rust kernel against canonical R oracles |
| D-119 | 2026-08-23 | Prospectively freeze MV8-S same-axis external baselines and the bounded 23-record PH sentinel before execution | MV8-S independently rehashes all 16 current HCA filtered/raw H5 inputs, both actual source caches, the ordered common-475 panel, and its 30-PC reference. The 23-record queue contains exactly 8 cell and 15 gene views, including seven source-produced gene sentinels. Eight baseline primaries/repeats, 23 PH primaries/repeats, full-view H0 MST oracles, three subset and one full Ripserr/GUDHI checks, one-worker monitoring, exact fallback rules, exact execution-head publication, a 126,600-s aggregate allowance, and an 8-GiB private cap are bound; 16/16 prefreeze checks and 29/29 focused assertions pass with zero execution | Commit and push the exact execution machinery, then run only the 66 monitored units into fresh ignored/private and public-ledger roots. Require MV8-T independent rehash/reconstruction closure before any full-PH prefreeze. The remaining 1,257 PH records, all landscapes/comparisons, clustering, fusion, labels, outcomes, adoption, and claims remain closed |
| D-120 | 2026-08-23 | Preserve MV8-S's first-child implementation stop and prospectively authorize one helper-only fresh-root replacement | Commit `218af286…` stopped on HCA_BM_001 before cell selection because the child omitted `R/toy_baseline.R`, which defines `.with_preserved_seed`. The immutable receipt records exit 1 after 32.74 s at 931,418,112 B, below both caps, with no baseline output and zero PH/landscape work. MV8-SA binds the exact failed ledger/log hashes and the helper-only remediation; 12/12 recovery checks and 42/42 focused assertions pass, including independent implementation/manifest rehash | Commit and push MV8-SA, then launch exactly one complete replacement into fresh v2 roots with the amendment path and exact new execution head bound. Preserve the failed roots. Keep the original science, resources, zero within-run retries, landscape definition, and every downstream firewall unchanged |
| D-121 | 2026-08-23 | Preserve the invalidly expanded MV8-SA launch head and require mechanically derived exact-head provenance for one fresh v3 replacement | The v2 launch supplied `ecc01ad1…` instead of actual Windows Git HEAD `ecc01ada…`; it was immediately interrupted and its exact orphan child terminated. The v2 roots contain only two zero-byte logs, with no public ledger, baseline, PH, or landscape artifact. MV8-SB passes 12/12 checks and the combined focused suite passes 54/54 with historical MV8-SA helper hashing and current MV8-SB runner/closure hashing separated correctly | Commit and push MV8-SB. In one PowerShell launch command, capture the full head directly from `git rev-parse HEAD` and pass that same value to the environment and runner; use fresh v3 roots and bind both recovery audits. Preserve v1/v2 roots and keep science, resources, landscapes, and downstream firewalls unchanged |
| D-122 | 2026-08-23 | Close the mechanically Git-bound MV8-S sentinel through independent MV8-T reconstruction | Execution head `38e2bdeb…` completed 66/66 monitored units in 1,277.54 aggregate child-seconds: eight current-input baselines plus byte-identical repeats, 23 Ripserr PH records plus byte-identical repeats, and four Ripserr/GUDHI checks yielding eight exact H0/H1 comparisons. No fallback or unexpected stderr occurred; peak RSS was 5,465,161,728 B and private storage 20,533,642 B. MV8-T independently rehashes 20 inputs and all 31 primary/repeat artifact pairs, reconstructs every typed view, reruns all 23 full-view H0 MST oracles, binds MV8-SA/MV8-SB, and passes 20/20 checks | The bounded sentinel is closed. A separate full-PH production prefreeze may now be planned, but no remaining PH, landscape, comparison, clustering, fusion, label, outcome, adoption, or manuscript-claim execution is yet authorized. Landscape production still requires Rust rebuild/hash rebind and canonical R oracle admission |
| D-123 | 2026-08-23 | Prospectively freeze the exact 1,257-record full-PH production remainder after MV8-T | MV8-U passes 23/23 public-receipt checks and subtracts the 23 independently closed MV8-T identities from the immutable 1,280-row MV8-R queue. All remaining records are source-produced gene views: 1,236 internal and 21 external; 625 SCT-data and 632 Pearson-residual; 131 accepted source caches. A label-blind serial order runs 14 common475 records, 625 exact500 residual records, then 618 exact500 SCT-data records. All four strata have sentinel measurements. The measured-max projection is 15.25 h and 164 MiB; a two-times 30.49-h planning projection fits within a 72-h aggregate stop and 1-GiB private cap. Current WSL memory/disk exceed the prospective launch gate. Ripserr remains primary at 1,800 s/8 GiB, with one exact GUDHI attempt only after an RSS stop at 1,800 s/12 GiB. Entry, runner, strict-prefix resume, and independent 1,280/1,280 closure machinery are hash-bound; no PH or landscape was executed | After commit and a fresh launch recheck, MV8-V may run exactly the 1,257 jobs with one worker and zero automatic retries. Stop and preserve on timeout, ordinary failure, ambiguous partial evidence, or aggregate cap. MV8-W must reconstruct/revalidate all 1,257 records and bind the 23 MV8-T records before any landscape sprint. Landscapes retain the approved all-active-level/no-grid/no-cap/exact-or-error-controlled streamed definition but remain closed pending full-PH closure and Rust rebuild/hash-rebind/R-oracle admission; comparisons, clustering, fusion, labels, outcomes, adoption, and claims remain closed |
| D-124 | 2026-08-23 | Preserve MV8-V's valid first PH child and prospectively recover the parent helper-load omission without retrying science | At execution head `79c2147…`, production job 1 completed Ripserr PH in 1.86 s at 96.1 MB with empty child stderr, an atomic 40,518-byte record, and a completed resource-ledger row. The parent then stopped before selected-metric publication because it had not loaded the already-bound `mv08s_validate_ph_record_v1` helper. MV8-VA independently reconstructs the frozen typed view, validates the record and full-view H0 MST oracle, rehashes six stopped-run files, and proves the two prior wrapper-launch attempts produced no R process, output root, or scientific record. All 20/20 recovery checks pass; the original runner/root evidence remains immutable | After the MV8-VA commit, byte-copy—not recompute—the admitted job-1 record into fresh private/public roots, publish an exact one-record completed prefix, and use the separately hash-bound helper-loaded runner with `--resume` at job 2. Zero retries or recomputations are authorized. Preserve all original/empty-launch evidence. PH science, one-worker caps/fallback/aggregate policies, landscapes, comparisons, clustering, fusion, labels, outcomes, adoption, and claims remain unchanged and closed as applicable |
| D-125 | 2026-08-23 | Preserve MV8-VA's zero-output bootstrap type check and bind a value-equality-only remediation | The committed bootstrap reopened and validated evidence but stopped before creating fresh roots because `file.info()` byte counts were doubles while CSV byte counts were integers and the script used storage-type-sensitive `identical()`. A first logged reproduction also carried a mistyped manually expanded head; a second used the exact mechanically read Windows Git head. Both produced only a 66-byte runtime notice and 54-byte validation error, with zero roots, copies, recomputations, retries, or science. MV8-VB passes 12/12 checks and confirms all MV8-VA audit artifacts plus every non-bootstrap implementation remain immutable | Commit the one-line numeric value-equality amendment before any new bootstrap. Then run exactly one mechanically Git-bound corrected bootstrap into the still-fresh v2 roots, require a byte-identical one-record prefix with zero retries/recomputation, and start only the MV8-VA recovery runner with `--resume` at job 2. All science, resources, landscapes, comparisons, clustering, fusion, labels, outcomes, adoption, and claims remain unchanged |
| D-126 | 2026-08-23 | Preserve MV8-VB's zero-output bootstrap self-hash stop and explicitly admit its amendment chain | The mechanically Git-bound corrected invocation reached the original MV8-VA implementation guard and stopped because the bootstrap did not yet traverse the MV8-VB binding for its own changed hash. The interactive invocation created no raw log file and MV8-VC records that limitation explicitly; both fresh v2 roots remain absent, with zero copies, recomputations, retries, or science. MV8-VC passes 14/14 checks and binds the immutable MV8-VA audit, immutable MV8-VB audit, and current bootstrap hash | After the MV8-VC commit, run one corrected bootstrap with the explicit manifest-verified amendment-prefreeze path. Require a byte-identical one-record prefix with zero retries/recomputation, then resume only at production order 2. Science/resources and every downstream firewall remain unchanged |
| D-127 | 2026-08-23 | Bind the successful one-record recovery prefix and the runner's amendment-chain admission before resume | The MV8-VC-bound bootstrap accepted exactly the independently validated job-1 PH record by byte copy: 40,518 bytes, SHA-256 `50398bd6…`, H0 MST passed, zero recomputation, zero retries, and resume order 2. Prelaunch inspection found that the helper-loaded recovery runner still compared the amended bootstrap against its historical MV8-VA hash. MV8-VD passes 16/16 checks and binds both the current bootstrap and the runner's manifest-verified amendment traversal before execution | After the MV8-VD commit, launch only the committed recovery runner with explicit MV8-VA and MV8-VD prefreeze paths, one worker, and `--resume`. It must retain the strict one-record prefix and begin new PH at order 2. MV8-U science/resource/fallback/aggregate policy is unchanged; landscapes and all downstream work remain closed |
| D-128 | 2026-08-24 | Close the exact MV8-V full-PH queue and bind its recovery provenance to independent MV8-W reconstruction | MV8-V completed 1,257/1,257 records with one worker, zero retries, zero fallback, empty launch stderr, zero partials, 4,421.03 aggregate attempt-seconds, 4,367,720,448-B peak child RSS, and 75,605,218 B of private evidence. MV8-W independently rehashes all 131 source caches, reconstructs all 1,257 typed views, reopens every PH artifact, reruns every full-view H0 MST oracle, and combines the immutable 23 MV8-T records for exact 1,280/1,280 coverage with 20/20 checks. MV8-WX binds the complete MV8-VA→MV8-VD recovery chain, the byte-identical job-1 bootstrap, the exact committed recovery head, zero recomputation, and zero retry with 26/26 aggregate-only checks | Full PH is closed. No landscape, comparison, clustering, fusion, label, outcome, adoption, or manuscript-claim job is authorized by this closure. The next prospective stage may only rebuild and hash-bind the accepted Rust landscape kernel, admit it against canonical R oracles, and freeze the exact all-active-level streamed landscape queue before any landscape calculation |
| D-129 | 2026-08-24 | Prospectively freeze the post-MV8-W Rust landscape-kernel rebuild and canonical-R admission gate | MV8-X passes 20/20 metadata-only checks and selects one deterministic immutable stress pair for each of 28 future landscape groups: 14 H0 and 14 H1 across internal124/external8, common475/exact500, all three representations, and both cell/gene topology views. The approved dissertation-aligned definition remains all finite positive-persistence intervals, essential-H0 exclusion, all consecutive active levels, exact or error-controlled squared L2, separate H0/H1 primary outputs, no universal grid or level cap, and streamed/chunked evaluation. The dependency-free source, canonical R oracle/wrapper, PH inputs, Rust 1.97.1 installer identity, two clean builds, two independent oracle runs, one-worker/zero-retry resource limits, privacy rules, and 17 acceptance gates are hash-bound. No landscape distance was calculated | Commit and push MV8-X before acquiring the isolated project-local toolchain. Then rebuild twice, require byte-identical binaries, native C-ABI and strict unit/format/lint checks, run the 28 real-PH oracle pairs plus nine analytic/fallback fixtures twice, and publish MV8-Y independent admission closure. Production landscapes and every comparison, clustering, fusion, label, outcome, adoption, and manuscript-claim stage remain closed |
| D-130 | 2026-08-24 | Preserve the first opaque MV8-X oracle stop and prospectively correct only its active-level diagnostic | Both clean Rust 1.97.1 builds are byte-identical to each other and to the previously accepted 442,272-byte candidate (`51d3fca…`); format, four unit tests, strict Clippy, and the native C ABI pass. Oracle run A then completed all 28 rows and stopped at the combined gate without publishing its root. An aggregate-only Rust diagnostic passes engine, reverse, and both self-zero gates 28/28 but reproduces the old active-level assertion only 14/28. Actual-code inspection shows the assertion incorrectly equated finite-interval count with the number of nonzero landscape levels; the kernel correctly reports critical-landscape depth, and existing HCA review code defines this as maximum overlap depth. MV8-XA binds the failed attempt, exact old/new runner and closure hashes, a canonical endpoint-sweep active-depth oracle, per-pair checkpoints, unchanged selection/binary/formula/R tolerances/resources, and 16/16 recovery checks | Commit and push MV8-XA. Run exactly one fresh observable replacement A at the new committed head; run B only if A passes. Preserve any checkpoint/failure evidence. No automatic retry, production landscape, comparison, clustering, fusion, label, outcome, adoption, or manuscript-claim execution is authorized. MV8-Y must traverse both MV8-X and MV8-XA before admission |
| D-131 | 2026-08-24 | Admit the hash-bound private Rust kernel through independent MV8-Y closure without opening production landscapes | The observable replacement A passes 28/28 frozen real-PH pairs and 9/9 analytic/fallback fixtures in 1,855.096 s at 262,983,680-B peak RSS; conditional run B independently reproduces the scientific and fixture tables byte-for-byte in 1,911.537 s at 264,826,880 B. Both remain below the 3,600-s/12-GiB one-worker caps and bind amendment head `2e8ae125…`. MV8-Y traverses immutable MV8-X plus MV8-XA, rehashes the official isolated toolchain, dependency-free source, duplicate 5.80/5.96-s builds, exact prior-candidate identity, native ABI/dependencies, 28 group oracles, maximum-overlap active depths, reverse/self invariants, nine fixtures including missing/corrupt-library fallback, A/B determinism, private-storage/privacy rules, and downstream firewalls with 37/37 checks | The explicit SHA-256-verified private WSL candidate is admitted only for a separate future landscape-execution prefreeze. R remains the canonical oracle; grouped Persim fallback remains available. The binary is not published or defaulted. No production landscape, comparison, clustering, fusion, label, outcome, adoption, or manuscript-claim job is yet authorized |
| D-132 | 2026-08-24 | Prospectively freeze exact streamed production-landscape execution and a maximum-burden sentinel without opening the full queue | MV8-Z passes 31/31 metadata-only checks and traverses MV8-R, MV8-W, MV8-X/MV8-XA, and MV8-Y. The immutable workload contains 28 H0/H1-separate groups and 152,744 unordered dimension-specific pairs, partitioned exactly into 628 hash-bound chunks of at most 250 pairs. Private group axes are label-blind and bind all 2,544 group-unit memberships; public evidence exposes no unit, job, accession, or path. The admitted 442,272-byte Rust candidate remains explicit-hash/private/not-default. The label-blind maximum-burden group is internal exact500 residual H1; its top pair has 3,214 and 2,702 finite intervals, so its canonical R check is prospectively routed to the admitted adaptive error-certified oracle rather than the 500-interval exact guard. Strict atomic resume, zero retry, one worker, no mixed-engine fallback, resources, schema, privacy, and downstream firewalls are frozen; no landscape was executed | After this commit, authorize only one fresh 250-pair Rust chunk, one separate fresh repeat, and one canonical-R error-certified check of its maximum-burden pair. MV8-ZA must independently rehash, compare scientific bytes, validate active depth/oracle/resources/privacy, and close the sentinel. The remaining 152,494 unique production pairs and all comparisons, clustering, fusion, labels, outcomes, adoption, and manuscript claims remain closed pending a new prospective production decision |
| D-133 | 2026-08-24 | Bind live elapsed/RSS enforcement around the unchanged MV8-Z sentinel workers before execution | MV8-ZA passes 12/12 metadata-only checks and binds the unchanged committed chunk/oracle workers plus a parent process-tree monitor. The monitor requires a mechanically supplied exact Git head and fresh roots, measures recursive RSS every 0.25 seconds, terminates only its child tree on the already frozen 3,600-second or 4-GiB cap, classifies expected completion stderr, writes aggregate public receipts, and permits exactly two 250-pair Rust children plus one R oracle with one worker and zero retry. No estimand, input, engine, queue, fallback, or downstream contract changes | Commit/push MV8-ZA, then launch exactly the three authorized children. Preserve any stop. If all pass, MV8-ZB independently rehashes private outputs, compares normalized scientific bytes, validates adaptive R agreement, resources, privacy, and firewalls before any full-production decision |
| D-134 | 2026-08-24 | Preserve the first monitored sentinel helper-load stop and bind a helper-only fresh-root recovery | At committed head `0afc8a7…`, the primary child stopped before any pair calculation because both workers omitted `R/landscape_contract.R`, which defines `finite_landscape_diagram()`. The monitor preserved one failed ledger row and exact logs after 1.59 s at 75.1 MB; zero pair outputs were written and the repeat/oracle never started. MV8-ZB passes 14/14 checks, binds exact old/new worker and monitor hashes, adds only the already accepted helper source plus amendment traversal, and retains unchanged science, caps, one worker, zero automatic retry, and all downstream firewalls | Commit/push MV8-ZB, then execute exactly one complete three-child replacement into fresh v2 roots. Preserve v1. Stop on any new failure; if clean, MV8-ZC independently closes determinism, adaptive oracle agreement, resources, privacy, and firewalls |

## 9. Status dashboard

| Phase | Status | Gate |
|---|---|---|
| 0. Preserve and establish provenance | Audited nine-slice publication stack merged; private PDFs/reviews, generated evidence bundle, binaries, and `example_run.r` remain excluded; baseline owner/author review remains | G0 not evaluated |
| 1. Scientific/implementation audit | Corrected landscape contract and orientation-safe cell/gene definitions validated; historical diagrams invalidated; consolidated claim/author review remains | G1 not evaluated |
| 2. Reproducible baseline/repository health | Published main CI passes privacy, exact restore, package check, and installed realistic H0/H1 fixtures; clean-clone user route and complete policy review remain | G2 not evaluated |
| 3. Literature/reference/figure audit | Initial literature sample only | G3 not evaluated |
| 4. Redesign primary comparison | Prediction-locked retrieval, clustering, complete cell robustness panels, corrected matrices, and null/negative results have auditable evidence; MV7-A maps remaining gene-view and confounding gaps without outcome selection | G4 not evaluated |
| 5. Expand methods | MV-06 complete: equal-weight fusion fails the both-component rule; outcome validation 15/15 and repeat 13/13 pass; cell and gene views continue separately into MV-07 | G-MV6 closed negative; G5 not evaluated |
| 6. Biological/practical validation | HCA exact-annotation/common-475 sensitivity complete; MV8-H exact 48-file acquisition, Ensembl-93 inputs, Cell Ranger 8.0.1 runtime/reference, one-unit count execution, matrix/QC review, eight-unit common-475 dual-view topology/all-active-level validation, and MV8-J claim map are auditable. MV8-M isolates the raw-read exact-500 failure to later-production SCT `data`; MV8-N freezes the Pearson-residual migration; MV8-O validates source construction; MV8-Q closes all 132 source fits; MV8-R freezes the corrected 1,280-record topology queue; MV8-T closes the Git-bound 66-unit sentinel; MV8-V completes the exact remainder; and MV8-W/MV8-WX independently close 1,280/1,280 PH records plus the complete recovery chain. MV8-X/MV8-XA/MV8-Y admit the private Rust engine, and MV8-Z now freezes the exact 28-group/152,744-pair streamed queue without calculating production landscapes | G6 not evaluated; only the bounded maximum-burden MV8-Z sentinel is authorized after commit. Full landscapes/comparisons, clustering, fusion, labels, outcomes, claims, and adoption remain closed |
| 7. Profile and optimize | The dependency-free Rust landscape kernel rebuilds byte-identically under isolated Rust 1.97.1, exactly matches the prior accepted binary, passes format/unit/lint/native-ABI/dependency gates, agrees with canonical R on 28/28 real-PH stress pairs plus 9/9 fixtures twice, and is now bound to 628 resumable chunks under the all-active-level/no-grid/no-cap definition. MV8-Z identifies a substantially heavier H1 maximum-burden sentinel and correctly routes its R oracle to error-certified adaptive integration. R remains canonical; grouped Persim remains a separately gated recovery route | G7 not evaluated; run and independently close only the fresh/repeat 250-pair sentinel plus one R oracle before any full production or binary distribution/default adoption decision |
| 8. Rewrite manuscript | Not started | G8 not evaluated |
| 9. Release/archive | Not started | G9 not evaluated |

## 10. Revision history

| Version | Date | Summary |
|---|---|---|
| 0.4.19 | 2026-08-24 | Preserved the first MV8-ZA child's 1.59-s/75.1-MB missing-helper stop with zero pair outputs and zero later children; MV8-ZB passes 14/14 and prospectively binds only `landscape_contract.R` loading plus amendment traversal for one fresh-root three-child replacement under unchanged science/resources/firewalls. |
| 0.4.18 | 2026-08-24 | Added the 12/12 MV8-ZA monitor-only prefreeze: recursive process-tree elapsed/RSS enforcement and exact-head/fresh-root guards are bound around the unchanged MV8-Z workers before the three-child maximum-burden sentinel; science and all full-production/downstream authorizations remain unchanged and closed. |
| 0.4.17 | 2026-08-24 | Prefroze MV8-Z with 31/31 metadata-only checks: the exact 28-group/152,744-pair H0/H1-separate all-active-level workload is bound into 628 atomic resumable chunks, the private admitted Rust candidate and all PH/group/pair axes are hash-bound, and only a fresh/repeat 250-pair maximum-burden H1 sentinel plus one adaptive error-certified canonical-R oracle are authorized; full production and all downstream work remain closed. |
| 0.4.16 | 2026-08-24 | Closed MV8-Y Rust admission with 37/37 checks: duplicate clean builds exactly reproduce the prior accepted candidate; independent A/B runs pass all 28 frozen real-PH groups and nine analytic/fallback fixtures with byte-identical scientific tables, correct maximum-overlap active depth, reverse/self invariants, and bounded resources/privacy; the private hash-explicit engine is admitted only for a separate landscape-execution prefreeze, with zero production or downstream authorization. |
| 0.4.15 | 2026-08-24 | Preserved the first opaque MV8-X oracle-gate stop; proved the rebuilt Rust candidate exactly matches the prior accepted binary and passes build/native checks; isolated the stop to an incorrect raw-interval-count assertion while engine/reverse/self gates pass 28/28; and prefroze one observable active-depth-corrected replacement with 16/16 MV8-XA checks under unchanged pairs, formula, R tolerances, resources, and downstream firewalls. |
| 0.4.14 | 2026-08-24 | Prefroze the post-MV8-W Rust landscape admission with 20/20 checks, one immutable stress pair for each of 28 future groups, exact dissertation-aligned all-active-level H0/H1-separate semantics, a dependency-free Rust 1.97.1 duplicate-build gate, two canonical-R oracle runs, nine fixtures, bounded resources/privacy, and zero production landscapes or downstream authorization. |
| 0.4.13 | 2026-08-24 | Closed MV8-V at 1,257/1,257 with one worker, zero retries/fallback/partials, empty stderr, and all caps/firewalls retained; MV8-W independently reconstructs and rehashes the complete 1,280/1,280 PH inventory with 20/20 checks, while aggregate-only MV8-WX binds the full recovery provenance with 26/26 checks and leaves landscapes plus all downstream science closed. |
| 0.4.12 | 2026-08-23 | Bound the successful byte-identical one-record MV8-VA prefix and prelaunch runner amendment traversal; MV8-VD passes 16/16, authorizes only the 1,256-record resume at order 2 with one worker/zero retries, and leaves the scientific/resource contract and all downstream firewalls unchanged. |
| 0.4.11 | 2026-08-23 | Preserved MV8-VB's zero-output self-hash-guard stop with its interactive capture limitation explicit; MV8-VC passes 14/14 and binds the MV8-VA→MV8-VB→current bootstrap implementation chain before one fresh-root, zero-retry, zero-recomputation bootstrap. |
| 0.4.10 | 2026-08-23 | Preserved MV8-VA's storage-type-sensitive byte-count stop and an invalid manually expanded reproduction head; a mechanically head-bound repeat confirms zero roots/copies/recomputations/retries, and MV8-VB passes 12/12 while binding the sole numeric value-equality remediation under unchanged recovery science and firewalls. |
| 0.4.09 | 2026-08-23 | Preserved MV8-V's valid 1.86-s/96.1-MB job-1 Ripserr record and parent-only missing-helper stop; MV8-VA passes 20/20 checks after independent typed-view/PH/MST validation, binds both empty wrapper-launch histories, authorizes a byte-identical fresh-root one-record prefix with zero retries/recomputations, and resumes only at job 2 under unchanged science/resources/firewalls. |
| 0.4.08 | 2026-08-23 | Prefroze the exact 1,257-record full-PH remainder with 23/23 checks, immutable source/view hashes, four sentinel-measured strata, a label-blind serial risk order, 15.25-h measured-max and 30.49-h two-times planning projections inside a 72-h stop, 164-MiB projected output inside a 1-GiB cap, strict-prefix resume, exact RSS-only GUDHI fallback, and prospective independent 1,280/1,280 closure while landscapes and all downstream science remain closed. |
| 0.4.07 | 2026-08-23 | Closed the mechanically Git-bound MV8-S v3 sentinel at 66/66 with eight exact baseline repeats, 23 Ripserr PH records and exact repeats, all H0 MST oracles, eight exact cross-engine H0/H1 checks, no fallback/unexpected stderr/partials, and low measured resource use; MV8-T independently rehashes/reconstructs the complete evidence and passes 20/20 while all later PH and landscapes remain separately gated. |
| 0.4.06 | 2026-08-23 | Preserved the first MV8-SA replacement's invalid manually expanded head, terminated its exact orphan child, proved v2 contains only two empty logs and zero scientific/public outputs, and prefroze one fresh v3 launch whose full head must be derived directly from Windows Git; MV8-SB passes 12/12 and the combined focused suite 54/54 without changing science, resources, landscapes, or downstream firewalls. |
| 0.4.05 | 2026-08-23 | Preserved MV8-S commit `218af286` stopping before cell selection because its baseline child omitted the preserved-seed helper import; bound the 32.74-s/0.93-GB/no-output receipt, added only `R/toy_baseline.R`, required immutable failed roots plus one fresh replacement, and passed 12/12 recovery checks and 42/42 focused assertions without changing science, resources, landscapes, or downstream firewalls. |
| 0.4.04 | 2026-08-23 | Prefroze MV8-S with fresh hashes for 16 current HCA H5 inputs, two source caches, the common-475 panel/reference, exactly 8 cell plus 15 gene PH records, deterministic baseline/PH repeats, full-view MST checks, reduced/full Ripserr-GUDHI checks, exact execution-head publication, serial caps, and prospective MV8-T independent closure; 16/16 audit checks and 29/29 focused assertions pass before execution. |
| 0.4.03 | 2026-08-23 | Prefroze MV8-R full topology production with 20/20 metadata checks: reconciled 16 common475 panel receipts, required a current-input same-selected-cell external baseline after detecting eight historical axis mismatches, bound the corrected 1,280-record PH queue and immutable 1,240-record internal comparator, retained separate H0/H1 and exact all-active-level streamed landscapes without a grid or level cap, and authorized only the 23-record PH sentinel after commit. |
| 0.4.02 | 2026-08-23 | Closed MV8-Q full Pearson-residual source production: jobs 125-129 completed under the authorized 14-GiB overlay, all 129 production caches/audits were independently rehashed across three immutable roots, both stopped receipts were preserved, 17/17 closure checks passed, and the three MV8-O primaries complete 132/132 source coverage while topology and outcomes remain closed. |
| 0.4.01 | 2026-08-23 | Preserved accepted MV8-PR job 124 and the immutable job-125 12-GiB stop; prefroze a fresh jobs-125-129 overlay with one explicit job-125 retry, first attempts for jobs 126-129, a prospective 14-GiB monitored cap, and unchanged one-worker/1,800-s/2-GiB-future/scientific policies. |
| 0.4.00 | 2026-08-23 | Preserved MV8-P's 123-fit accepted prefix and job-124 stopped receipt; prefroze a fresh jobs-124-129 overlay with a bounded 2-GiB `future` export allowance, one explicit retry only for job 124, unchanged one-worker/1,800-s/12-GiB caps, immutable original evidence, merged 132-source closure, and all downstream firewalls retained. |
| 0.3.99 | 2026-08-22 | Implemented and verified the resumable MV8-P serial source-production runner plus prospective MV8-Q closure: exact queue/source hashes, one-worker caps, strict-prefix resume, atomic ledgers, external selected-axis preflight, independent private-cache rehash, and topology/outcome firewalls are enforced; focused checks pass 22/22 and the full suite passes 3,305 with two established optional skips. |
| 0.3.98 | 2026-08-21 | Froze MV8-P full Pearson-residual source production after the successful MV8-O sentinel: 129 remaining sources (122 internal, seven external) are bounded to a single worker, zero automatic retries, 1,800 s, and 12 GiB, in ascending frozen-QC cell count. The largest remaining source is below the completed maximum sentinel; source fitting is authorized only after commit, while all topology and biological analysis remain closed. |
| 0.3.97 | 2026-08-21 | Closed MV8-O's bounded all-QC Pearson-residual source/view sentinel: internal minimum primary, internal maximum primary/repeat, and HCA_BM_002 primary/repeat all passed rehashed resource, geometry, deterministic-repeat, and scope-firewall checks. The maximum source remained under 12 GiB but with a narrow margin, so full production requires a separate resource prefreeze; no PH, landscapes, clustering, fusion, labels, outcomes, claims, or default adoption was opened. |
| 0.3.96 | 2026-08-21 | Completed MV8-N's metadata-only Pearson-residual migration prefreeze: independently hash-bound 124 internal raw shards and eight current exact-reference HCA matrices; froze 132 all-QC model fits, 1,280 proposed gene views, 1,272 planned PH records, 28 exact all-active-level landscape groups, and 40 paired descriptive strata; passed 18/18 checks and authorized only a three-unit source/view sentinel after commit. |
| 0.3.95 | 2026-08-21 | Closed MV8-M's label-closed exact-500 representation audit: reconciled the original Pearson-residual contract against later SCT-data cache drift; on identical 500 genes and 4,614 HCA_BM_002 cells, SCT data remains invalid with one constant gene while Pearson residuals and diagnostic log-normalized RNA form valid exact-500 correlation-chord geometries. Corrected primary/repeat summaries and common475 stability metrics are byte-identical, resources pass, and the residual candidate remains owner-gated pending a prospective internal/external migration sensitivity prefreeze. |
| 0.3.94 | 2026-08-21 | Closed MV8-L's all-QC-cell exact-500 feasibility test: expanding gene-side observations to all 4,614 frozen-QC cells restores full retention but not nonzero transformed variance; both documented SCT configurations leave one constant standardized gene and correctly block correlation-chord topology. Repeats and resource gates pass; cross-unit computation, PH, landscapes, labels, outcomes, fusion, clustering, manuscript work, and deletion remain closed. |
| 0.3.93 | 2026-08-21 | Closed MV8-K's exact-500 HCA_BM_002 transform-contract audit: standard SCT retains only 497/500 panel genes; low-count SCT retains 500/500 and projects to the frozen 30-PC model but leaves three exactly zero-variance standardized genes, correctly blocking the separate gene correlation-chord view. The repeat source matrix is deterministic, resources are recorded, and no topology, landscape, labels, outcomes, fusion, clustering, manuscript work, or deletion was opened. |
| 0.3.92 | 2026-08-20 | Closed the owner-authorized seven-unit MV8-H exact-500 raw-read Cell Ranger 8.0.1 recovery: all units completed under the frozen 4-core/32-GiB contract; corrected independent validation passed 68/68 with 33,563-feature axes, exact500=500, common475=475, 384-cell eligibility, required outputs, clean stderr, and terminal performance evidence; public audit/test artifacts remain label-, topology-, landscape-, fusion-, manuscript-, and deletion-firewalled |
| 0.3.91 | 2026-08-19 | Closed the bounded HCA_BM_002 exact-500 Cell Ranger feasibility gate: existing filtered features contain 33,563 unique IDs with all 500 exact-panel genes and the common475 subset, 5,037 filtered cells with 4,614 passing the frozen 384-cell QC minimum, and the existing run's 14,348-second/20.141-GB/4-core/32-GiB resource record passes; independent validation is 10/10, with remaining units and all downstream biology still closed |
| 0.3.90 | 2026-08-19 | Bound MV8-J's provisional claim-to-evidence map to the completed eight-unit common-475 HCA validation: reconciled 14 method, technical, descriptive, blocked, and decision-required claims; recorded cross-unit cell/gene H0/H1 ranges and active landscape-level ranges; recommended exact-500 recovery only for a claim-specific manuscript need, while keeping labels, outcomes, fusion, clustering, and deletion closed |
| 0.3.89 | 2026-08-19 | Closed MV8-I eight-unit common-475 HCA technical validation: all eight raw H5 identities rehashed, 32 separate cell/gene H0/H1 topology profiles and 32 exact critical-breakpoint all-active-level landscape profiles completed, primary/repeat scientific hashes agreed after excluding elapsed time, landscape artifacts were byte-identical, independent validation passed 11/11, and labels/outcomes/fusion/exact-500 recovery/deletion remained closed |
| 0.3.88 | 2026-08-18 | Executed the owner-authorized common-475 secondary HCA_BM_002 topology gate: separate 384-cell/30-PC and 475-gene correlation-chord views, separate H0/H1 Ripserr PH, and exact critical-breakpoint all-active-level landscapes with no fixed grid or level cap; 1,225 active levels across four profiles, byte-identical primary/repeat aggregates, focused firewall test, and exact-500 block preserved |
| 0.3.87 | 2026-08-18 | Bound the MV8-H corrected dual-view topology prefreeze and executed a five-seed label-closed SCT panel-retention gate: all five seeds retain 384 QC-eligible cells but the frozen exact-500 panel retains only 497/500, 499/500, 497/500, 499/500, and 499/500 genes, so exact-500 cell/gene PH is blocked without substitution; all H0/H1 landscape clauses and labels/outcomes/other-unit/deletion firewalls remain unchanged |
| 0.3.86 | 2026-08-18 | Executed the approved one-unit MV8-H matrix/QC review: allowlisted filtered/raw H5 and aggregate metrics only, 8/8 checks, 5,037 filtered cells with 4,614 passing frozen QC depth, median 801 features and 2,545 counts, resource cap and 7/7 repeat passing, topology/labels/outcomes/remaining/deletion firewall preserved |
| 0.3.85 | 2026-08-18 | Prepared the MV8-H allowlisted matrix-content/QC review prefreeze without opening matrix content: filtered/raw H5 plus aggregate metrics only, molecule/web/labels/outcomes/landscapes/remaining/deletion firewall, eight validation checks, deterministic 8/8 public artifact repeat, and explicit owner-approval gate |
| 0.3.84 | 2026-08-18 | Admitted the completed HCA_BM_002 sentinel to metadata-only structural/QC review: 8/8 validation, exact seven-output/no-BAM structure, raw/filtered HDF5 metadata, 5,037 filtered cells matching aggregate metrics, resource provenance, deterministic 7/7 public artifact repeat, and downstream firewall preserved |
| 0.3.83 | 2026-08-18 | Recorded the owner clarification superseding the unintended no-count clause; accepted exactly one completed MV8-H sentinel execution, preserved historical prefreeze evidence, revalidated closure 12/12 with 8/8 artifact equality, and kept all downstream gates closed |
| 0.3.82 | 2026-08-18 | Revalidated the exact existing MV8-H one-unit execution after owner resumed the goal: closure 12/12, eight public artifacts 8/8 byte-identical, no second count, downstream firewall preserved, and strict no-count history retained |
| 0.3.81 | 2026-08-18 | Bound clean-export package verification for MV8-H: zero `tmp/`/PDF entries, successful 436,297-byte source build, dependency-gate-limited `R CMD check` without installing packages, and explicit preservation of the strict no-count objective as unresolved |
| 0.3.80 | 2026-08-18 | Reconciled the strict original no-count objective with the later owner-authorized MV8-H execution; preserved the 14/14 historical prefreeze and 12/12 separate execution closure, explicitly recorded that the literal no-count condition is not achieved, and isolated the remaining package-build hygiene issue caused by the large ignored `tmp/` tree |
| 0.3.79 | 2026-08-18 | Closed the owner-authorized MV8-H HCA_BM_002 Cell Ranger 8.0.1 sentinel as a metadata-only execution record: explicit pipestance success, no diagnostic stderr/failure markers, exact frozen command and input bindings, 20.141-GB peak RSS, 14,348-second elapsed runtime, 12/12 closure checks, 7/7 deterministic repeat, focused MV8-H tests passing, and matrix/QC/label/outcome/landscape/remaining-unit/deletion firewall preserved |
| 0.3.78 | 2026-08-18 | Completed MV8-H Cell Ranger 8.0.1 count-sentinel prefreeze evidence as byte-identical primary/repeat across 14 independent checks: 12/12 contract/firewall closure, 2/2 deterministic artifact equality, and complete runner/report closure; launcher and builder now harden forward-only controls, retain no-kill/no-delete semantics, and preserve explicit 4 cores/32 GiB and no universal level cap language for public auditability |
| 0.3.77 | 2026-08-18 | Prospectively lowered MV8-H `mkref` to four cores/32 GiB, completed the exact Cell Ranger 8.0.1 custom reference in 13,603 s without resource breach, independently closed its 19-file/20.766-GB tree, exact Ensembl-93 FASTA/GTF, 33,563-gene and exact-500/common-475 axes, disclosed and corrected the RSS composite-Boolean omission, passed 265/265 dependency, 254/254 focused, and 2,705-pass complete verification, retained the cell/gene separate-H0/H1 all-active-level landscape contract, and opened only a count-sentinel prefreeze |
| 0.3.76 | 2026-08-17 | Prospectively amended MV8-H to the exact installed Cell Ranger 8.0.1 distribution, froze explicit SC3Pv2/exon-only/BAM-disabled/no-secondary commands, rebound unchanged Ensembl-93/exact-500/common-475 identities, passed deterministic/independent/focused/complete verification, preserved separate cell/gene H0/H1 all-active-level landscapes, and opened only custom-reference construction before a separately validated count sentinel |
| 0.3.75 | 2026-08-17 | Closed Ensembl release-93 primary-assembly input identity with exact 881,214,682-byte size, local SHA-256, matching official/local BSD `07119 860562`, and full gzip integrity; stopped before Cell Ranger acquisition, EULA, `mkref`, or raw processing pending owner-provided exact 3.0.0 runtime validation |
| 0.3.74 | 2026-08-17 | Closed the MV8-H exact 48-FASTQ acquisition after independent size/SHA-256/gzip/receipt validation of all 85,034,239,918 bytes with zero partials and preserved storage reserve; authorized only Ensembl-93 primary-assembly input acquisition and identity validation before a separate Cell Ranger/runtime/reference gate |
| 0.3.73 | 2026-08-17 | Accepted the independently validated 394,373,114-byte MV8-H FASTQ sentinel with exact official SHA-256, gzip/receipt/zero-partial/storage gates, and authorized only resumable completion of the remaining 47 files before a separate checksum closure |
| 0.3.72 | 2026-08-17 | Prefroze the owner-authorized MV8-H exact 48-file/79.194-GiB HCA acquisition, minimal 33,563-gene exact-500 annotation, Cell Ranger 3.0.0/SC3Pv2 target, unchanged QC/384-cell/dual-view/all-active-level landscape contract, resource/firewall/stop rules, and 13/13 independent checks; only a one-file sentinel then checksum-verified download are open |
| 0.3.71 | 2026-08-17 | Completed exact HCA annotation reconciliation and the full 124-sample paired 500-versus-475 source/PH/all-active-level landscape/comparison chain; all independent gates pass and material cell-H1 sensitivity supports an exact-500 raw-read prefreeze, while the 79.194-GiB download and processing remain owner-gated |
| 0.3.70 | 2026-08-15 | Completed MV7-E: resolved all 16 approach conflicts against official GEO methods, superseded the non-estimable mixed-approach diagnostic without changing topology, activated the fixed all-124 global-core fallback after one KLF2 absence, froze five transductive per-seed cell/gene estimands and the unchanged all-active-level H0/H1 landscape contract, and authorized MV7-F upstream caches only |
| 0.3.69 | 2026-08-15 | Completed MV7-D full-corpus reconciliation: resolved 127 candidates into 124 corrected-descriptive and 90 primary samples, demonstrated 6/6 omitted-source/SCT feasibility, exposed 16 approach-label provenance conflicts, preserved the landscape contract, and added a prospective MV7-E–MV7-J expansion plan before any new PH |
| 0.3.68 | 2026-08-15 | Completed MV7-C claim/data-sufficiency synthesis and author dossier; recommend a methods-focused manuscript before more PH, while reserving external validation for a consciously chosen generalization claim and requiring author-team approval |
| 0.3.67 | 2026-08-15 | Completed MV7-B no-new-PH diagnostics: study/tissue influence and cell-versus-gene sign instability require narrower comparative claims, while retained-cell and approach flags do not trigger; advance to MV7-C synthesis without authorizing reruns or new data |
| 0.3.66 | 2026-08-15 | Completed MV7-A selection-resistant evidence mapping: preserved the corrected landscape contract, classified 14 robustness and 10 confounding axes, reconciled measured compute, and authorized only a no-new-PH MV7-B prefreeze before gene reruns or external data |
| 0.3.65 | 2026-08-15 | Executed and independently validated exact-commit MV6-H blocked outcomes; equal fusion failed the required both-positive MRR rule, so G-MV6 closes negative, advanced fusion is rejected, and cell/gene views advance separately to MV-07 |
| 0.3.64 | 2026-08-15 | Froze and independently validated the complete MV6-H prediction lock over 75 groups and 318,150 rankings; all nine lock artifacts repeat byte-identically, metadata remains unopened, and only exact-commit blocked outcome execution is authorized |
| 0.3.63 | 2026-08-15 | Closed MV6-G computation after all 445 accepted corpus artifacts passed SHA-256, byte-size, and mtime preservation under the unchanged resume driver; authorized only prediction-manifest prefreeze before label access |
| 0.3.62 | 2026-08-15 | Completed clean corrected-root MV6-G production in 5.823640 worker-hours and passed independent 75-group validation 7/7; retained closed labels pending the 445-file immutable-resume gate |
| 0.3.61 | 2026-08-15 | Rebuilt and admitted corrected-root MV6-G serial execution at root `deb03fbc…f745` after independent, byte-repeat, focused, and 1,601-test gates; authorized a clean 74-group restart |
| 0.3.60 | 2026-08-15 | Admitted corrected MV6-G root `8b0a1e42…cec0c` after the maximum-group sentinel reproduced all three accepted scientific artifacts byte-for-byte; authorized only serial-completion re-prefreeze |
| 0.3.59 | 2026-08-15 | Preserved MV6-G's nine complete groups and rank-axis stop, diagnosed stale pre-sort indices at 81 training samples, corrected post-sort block validation at root `8b0a1e42…cec0c`, and authorized only maximum-group byte equivalence before restart |
| 0.3.58 | 2026-08-15 | Accepted the schema-corrected MV6-G first group after atomic completion, cap pass, and independent reconstruction of four scales, all nine formulas, and deterministic ranks; authorized the remaining 73 serial groups |
| 0.3.57 | 2026-08-15 | Captured MV6-G's first-group error as a `query_biological_pairs`/`biological_pairs` queue mismatch, added and independently enforced the exact runner alias at root `97ffd576…3b38`, and authorized one clean group rerun |
| 0.3.56 | 2026-08-15 | Preserved MV6-G's cap-passing first-group child failure, verified nonempty H0/H1 source diagrams, rebound private child-log capture at execution root `0c75e854…ee2b`, and authorized exactly one diagnostic rerun before later groups |
| 0.3.55 | 2026-08-15 | Prefroze the separately hash-bound MV6-G serial completion, complete validator, and 445-file resume checker under exact workload and aggregate caps; authorized only the remaining 74 label-closed groups |
| 0.3.54 | 2026-08-15 | Admitted the dynamic MV6-G runner after its monitored maximum-group rebind reproduced all three accepted scientific artifacts byte-for-byte; authorized only a separately prefrozen serial completion and immutable-resume policy |
| 0.3.53 | 2026-08-15 | Prefroze the dynamic MV6-G complete-production runner and exact 74-group workload, while authorizing only a monitored maximum-group three-artifact equivalence gate before any remaining group |
| 0.3.52 | 2026-08-15 | Accepted the corrected-root MV6-G maximum-group sentinel after byte repeat, independent scale/rank reconstruction, 12/12 R and 12/12 Persim oracles, resource/projection gates, and 5/5 immutable resume; only a separately prefrozen 74-group label-closed completion policy may proceed |
| 0.3.51 | 2026-08-15 | Preserved the first MV6-G sentinel's successful scientific gates but invalidated its launch root after the resume checker mishandled a spaced WSL path; rebound an argument-safe checker and regression at root `6a76a11d…ce82` before clean reexecution |
| 0.3.50 | 2026-08-14 | Prefroze MV6-G maximum-group execution under a ten-source implementation root with atomic output, 12-GiB/1,800-s/5-GiB caps, independent scale/ranking reconstruction, R/Persim oracles, clean repeat, and immutable resume; only primary and repeat sentinel runs are authorized |
| 0.3.49 | 2026-08-14 | Completed MV6-G complete-corpus fusion prefreeze: accepted queue root preserved; 262,675 training pairs, 1,050,700 component rows, 300 training-only scales, nine fixed rankings, blocked inference, and a strict prediction-before-label firewall pass independent validation and byte repeat; only the maximum-group sentinel is authorized next |
| 0.3.48 | 2026-08-14 | Closed MV6-F after the committed whole-corpus resume gate preserved SHA-256, bytes, and mtimes for 375 scientific artifacts plus canonical metrics; authorized only prospective label-closed MV6-G fusion prefreeze next |
| 0.3.47 | 2026-08-14 | Completed all MV6-F stage-two production under the serial 12-GiB policy; 74/74 canonical reuse metrics and 6/6 independent complete-production categories pass across 75 groups, and only a committed 375-artifact immutable-resume check is authorized next |
| 0.3.46 | 2026-08-14 | Accepted the successful 9.575-GB exact-group diagnosis and prefroze a separately hash-bound serial 12-GiB policy for all 74 unchanged MV6-F stage-two rows, with no scientific-root change, no retry, and downstream work still closed |
| 0.3.45 | 2026-08-14 | Preserved the first stage-two group's fail-closed 8.747-GB RSS stop and prefroze exactly one unchanged-science 12-GiB diagnostic under the existing serial aggregate allowance; no later group or downstream job launched |
| 0.3.44 | 2026-08-14 | Admitted MV6-F stage two after the corrected root reproduced all 180 diagrams and 6,500 landscape rows scientifically, repeated 3/3 artifacts byte-identically, passed 12/12 R and 12/12 Persim oracles plus 5/5 immutable resume and 8/8 rebind checks |
| 0.3.43 | 2026-08-14 | Completed MV6-F stage-two rebind prefreeze: preserved the queue and Rust identities, bound corrected root `5a1258e8…8d292`, added a serial live-cap/checkpoint monitor and complete validator, quarantined a pre-execution WSL resume-path defect, passed deterministic/independent/focused/full-suite gates, and authorized only remediated stage-one equivalence |
| 0.3.42 | 2026-08-14 | Completed MV6-F maximum-group stage 1: 180 dual-view PH records and 6,500 exact all-active-level landscape rows pass monitored caps, byte repeat, immutable resume, 12/12 R oracles, and 12/12 grouped-Persim oracles; authorized only the other 74 label-closed groups |
| 0.3.41 | 2026-08-14 | Completed MV6-F production prefreeze: exact 75-group/13,500-PH/141,400-component queue, bound streaming atomic runner, Rust/input roots, ten fail-closed rules, 10/10 independent validations, 19/19 focused expectations, six byte-repeat files, and zero production before a stage-1-only gate |
| 0.3.40 | 2026-08-14 | Completed MV6-E exact landscape acceleration admission: both grouped Persim and Rust pass R/cross-engine/invariant/determinism/resume/resource gates; Rust is preferred only for hash-verified private WSL production, reducing the maximum landscape projection to 0.873 worker-hours before a separate MV6-F production prefreeze |
| 0.3.39 | 2026-08-14 | Completed MV6-D bounded matched-SCT profiling: five fold sources, 20 cell/gene PH records, ten all-level landscape pairs, exact stage-1 repeat, and independent validation pass; PH is inexpensive, while per-pair landscape execution triggers a targeted batched/streamed acceleration gate before production or fusion |
| 0.3.38 | 2026-08-14 | Completed MV6-C global-core feasibility: all 450 accepted SCT caches verify, 2,536 genes satisfy the strict all-cache rule, the frozen 500-gene panel is stable across seeds and independently reconstructed, and only a bounded matched-SCT source/PH resource profile is authorized next |
| 0.3.37 | 2026-08-14 | Completed MV6-B prospective scale-up gate at a structural stop: 150 accepted records verify, 71/450 held-out views across 31/75 groups lack the full frozen panel, integrated artifacts contain no corrected gene-expression payload, and blocked fusion/MV-07 freeze now require an explicit owner estimand decision |
| 0.3.36 | 2026-08-14 | Completed MV6-A four-stratum label-closed fusion feasibility: all frozen sources/components/weights/neighbor diagnostics independently reconstruct and repeat byte-identically; low pilot redundancy justifies a separate scale-up/resource gate but no biological or blocked claim |
| 0.3.35 | 2026-08-14 | Reconciled the audited P01-P09 publication merge and exact main CI, corrected phase-dashboard drift, preserved the dissertation-aligned landscape contract, and authorized only a four-stratum label-closed MV6 fusion feasibility sprint before any larger matched gene-view or blocked evaluation |
| 0.3.34 | 2026-08-10 | Completed MV5-P label-closed distance production: 150 groups, 4,565 units, 1,838,725 values, 525 complete matrices, 12 exact R oracles, two 33-output maximum-group repeats, supplemental pseudobulk repeat, all-unit zero-rebuild resume, observed cap reconciliation, and explicit correction of the output-focused storage underestimate without fitting clusters or opening labels |
| 0.3.33 | 2026-08-10 | Completed MV5-O production prefreeze: bound 18 sources and four production implementations; froze 150 groups and 4,565 resumable units; added corrected 1.278-GB storage accounting, atomic/resume/abort rules, 12-oracle and maximum-group repeat plans; passed independent/repeat and real landscape-runner fixtures without executing production or outcomes |
| 0.3.32 | 2026-08-10 | Completed MV5-N label-closed clustering/resource gate: froze training-only PAM and deterministic held-out assignment, retained average linkage as sensitivity, instantiated exact complete-pair identity inventories, passed bounded exact landscape and matched-baseline admissions with independent oracles/repeats/resume, and projected full production at 16.117 worker-hours without authorizing it or opening outcomes |
| 0.3.31 | 2026-08-10 | Completed MV5-M no-outcome benchmark-gap gate: audited nine axes, independently reconstructed eight weighted candidate scores and exact fold pair scope, selected MV5-N clustering contract/resource gating, retained robustness second, blocked current technical mixing on identifiability, and authorized no full matrices or outcomes |
| 0.3.30 | 2026-08-10 | Completed MV5-L locked SCT-versus-integrated retrieval comparison: pre-join freeze, explicit known-marginals limitation, 2,250 exact pairs, 450 exact pseudobulk identity controls, paired H0/H1 topology-minus-energy DIDs, blocked uncertainty/Holm inference, 11-category independent reconstruction, and 13 byte-identical production artifacts |
| 0.3.29 | 2026-08-10 | Completed MV5-K prediction-locked integrated retrieval evaluation: 2,250 estimable endpoints, prespecified blocked uncertainty and multiplicity, null/negative H0/H1-versus-energy results, documented numerical-boundary correction, 11/11 independent validation categories, 15 byte-identical repeat artifacts, and zero SCT outcome access |
| 0.3.28 | 2026-08-10 | Completed MV5-J immutable label-closed integrated cell retrieval inputs: separate H0/H1, descriptive raw composite, exact-coordinate energy, shared pseudobulk context, 176,750 unique rankings, all-row/formula validation, byte-stable admission/resume/public assembly, and authorization for a separate prediction-locked integrated evaluation |
| 0.3.27 | 2026-08-09 | Completed MV5-I exact all-active-level integrated cell landscape distances: 70,700 independently validated H0/H1 rows, complete level accounting, four exact R oracles, byte-identical maximum-group repeat, 720-file immutable resume, complete storage projection, and authorization for separate label-closed retrieval-input assembly |
| 0.3.26 | 2026-08-09 | Completed MV5-H integrated cell PH: 6,750 typed cells-as-observations H0/H1 records, full independent identity/diagram/MST validation, byte-identical 90-view repeat, zero-rebuild resume, and cap-passing authorization for exact all-active-level landscape distances |
| 0.3.25 | 2026-08-09 | Completed MV5-G full fixed-D1-panel label-closed integrated-coordinate production: 75 groups, 6,750 views, independent validation, byte-identical complete-group repeat, immutable resume, and a cap-passing measured authorization for separate integrated cell PH |
| 0.3.24 | 2026-08-08 | Completed MV5-F fixed-D1-panel label-closed integration-induction resource gate: four structural real-data folds, immutable references, independent mapping validation, byte-identical repeat, and cap-passing 75-group mapping-plus-downstream projection without full execution or outcomes |
| 0.3.23 | 2026-08-08 | Completed MV5-E prediction-locked SCT cell retrieval evaluation: 2,250 complete query endpoints, five-seed/tissue macro summaries, paired study-block uncertainty and Holm inference, independently reproduced null H0/H1-versus-energy primary results, full tissue heterogeneity, and byte-identical production repeat |
| 0.3.22 | 2026-08-08 | Completed MV5-D5 immutable label-closed SCT cell retrieval inputs: separate H0/H1, descriptive raw composite, matched energy and pseudobulk context baselines, canonical rankings, independent formula validation, byte-stable resume/public assembly, and explicit no-leakage component-scale disposition |
| 0.3.21 | 2026-08-07 | Completed MV5-D4 exact all-active-level SCT cell landscape distances: 70,700 independently validated H0/H1 rows, explicit essential-H0 exclusion, eligible R-oracle agreement, byte-identical v3 group repeat, timing-separation correction, component assembly, and fully measured cap-passing cell-primary precomputation |
| 0.3.20 | 2026-08-07 | Completed MV5-D3 immutable full SCT cell-PH production: 6,750 independently validated H0/H1 records, complete stored MST evidence, 75 fresh MST checks, a byte-identical 90-view repeat, measured resources, and a cap-passing pre-landscape decision |
| 0.3.19 | 2026-08-07 | Completed MV5-D2 bounded label-closed cell-PH profiling: 30 representative H0/H1 diagrams, full-view MST correctness, byte-identical repeats, reduced Ripserr/GUDHI agreement, and a complete cap-passing SCT cell-primary projection without launching full PH or landscapes |
| 0.3.18 | 2026-08-07 | Completed MV5-D1 training-only SCT cell-fold coordinates: exact corrected pilot, 75 independently validated caches, 6,750 typed views, explicit held-out missing-feature mapping, and honest incomplete PH feasibility gate |
| 0.3.17 | 2026-08-07 | Completed MV5-D0 stage 1 from existing per-sample sources: 90 raw shards, 450 frozen selections, 450 independently validated v2 SCT caches, and mandatory post-cache reprojection |
| 0.3.16 | 2026-08-07 | Stopped MV5-D0 at monolithic-source and legacy-identity gates; added monitored matrix-only SCT caching, runtime-complete v2 identities, six-entry resource evidence, and source-migration requirement |
| 0.3.15 | 2026-08-06 | Completed MV5-C2 immutable sample–seed SCT caches, cached-fold equivalence, exact query-to-training chunks and resume, 90-sample execution plan, and conditional SCT cell-primary/no-go full-panel resource decision |
| 0.3.14 | 2026-08-06 | Completed MV5-C label-closed one-tissue/all-seed execution, exact H0/H1 landscapes, matched baselines, immutable neighbor/clustering artifacts, conditional gene-view disposition, and negative naïve full-run resource gate |
| 0.3.13 | 2026-08-06 | Completed MV5-A/MV5-B immutable label-closed LOSO manifests, analytical and scientific-shape matched baselines, and deterministic synthetic reference-to-query Seurat mapping without computing outcomes |
| 0.3.12 | 2026-08-06 | Froze the MV-05 sample-level LOSO statistical benchmark contract, executable registries/validators, existing-data feasibility, legacy-method exclusions, and inductive-integration gate without computing biological outcomes |
| 0.3.11 | 2026-08-05 | Completed MV-04 immutable exact all-level H0/H1 distance bundles, eligible-oracle/determinism checks, normalization and contribution contracts, bounded sensitivities, production profiling, and negative Rust decision |
| 0.3.10 | 2026-08-05 | Completed MV-03 representation-native corrected dual-view pilot diagrams, seed sensitivity, and resource feasibility |
| 0.3.9 | 2026-08-05 | Completed MV-02: added typed orientation-safe source/PCA/cell/gene/PH contracts, explicit legacy stamping, analytical and cache/metric fixtures, and a clean package gate without running biological data |
| 0.3.8 | 2026-08-05 | Completed MV-01: froze paired cell/gene mathematical contracts, method hierarchy, deterministic pilot manifest, provenance/leakage rules, resource envelope, and failure hypotheses without changing production PH |
| 0.3.7 | 2026-08-05 | Added the seven-sprint dual-view/multiview roadmap, preserving cell topology as confirmatory, gene topology as secondary, and fusion as gated exploratory work |
| 0.3.6 | 2026-08-05 | Independently cross-checked landscape distances, invalidated historical diagrams due feature-as-point orientation, profiled production candidates, and selected a non-activated corrected critical-pair prototype |
| 0.3.5 | 2026-08-05 | Added and benchmarked the non-default all-level exact/adaptive landscape reference oracle, recorded ADR-001, and kept scientific activation and Rust gated |
| 0.3.4 | 2026-08-03 | Fast-forwarded and verified main, documented test-environment limits, assessed the scratch example, and opened the BioRxiv v2 reconstruction trace |
| 0.3.3 | 2026-08-03 | Recorded local-only/Git-ignore policy for source PDFs and placed the uncertain example script on a preserved untracked hold |
| 0.3.2 | 2026-08-03 | Refreshed origin, confirmed fast-forward topology and no untracked/confidential collisions, and recorded a conditional integration decision |
| 0.3.1 | 2026-08-03 | Recorded Phase 0 canonical-state sprint: cached divergence assessment, artifact inventory, protected hashes, and safe synchronization procedure |
| 0.3.0 | 2026-08-03 | Added auditable phase subplans and work log; recorded the existing-data-first policy and conditional new-data decision gate |
| 0.2.0 | 2026-08-03 | Canonical evidence-driven plan incorporating repository lineage, dissertation/preprint findings, current literature, confidential external-review workstreams, and contributor-credit requirements |
