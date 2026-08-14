# scPHcompare Modernization and Resubmission Plan

## Document control

| Field | Value |
|---|---|
| Status | Draft for owner and author-team review |
| Version | 0.3.40 |
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
- **Status:** MV6-E admits both exact landscape engines after canonical-R, cross-engine, determinism, resource, and resume gates; Rust is preferred only for explicit hash-verified private WSL production, and a separate MV6-F full matched-production prefreeze is required before blocked fusion
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

## 9. Status dashboard

| Phase | Status | Gate |
|---|---|---|
| 0. Preserve and establish provenance | Audited nine-slice publication stack merged; private PDFs/reviews, generated evidence bundle, binaries, and `example_run.r` remain excluded; baseline owner/author review remains | G0 not evaluated |
| 1. Scientific/implementation audit | Corrected landscape contract and orientation-safe cell/gene definitions validated; historical diagrams invalidated; consolidated claim/author review remains | G1 not evaluated |
| 2. Reproducible baseline/repository health | Published main CI passes privacy, exact restore, package check, and installed realistic H0/H1 fixtures; clean-clone user route and complete policy review remain | G2 not evaluated |
| 3. Literature/reference/figure audit | Initial literature sample only | G3 not evaluated |
| 4. Redesign primary comparison | Prediction-locked retrieval, clustering, complete robustness panels, corrected matrices, and null/negative results have auditable evidence; author-level synthesis remains | G4 not evaluated |
| 5. Expand methods | MV6-E admits exact batched landscape production: grouped Persim is the portable fallback and hash-verified Rust is preferred privately; the full matched workload now projects below cap, so MV6-F production prefreeze is next before blocked fusion | G5 not evaluated |
| 6. Biological/practical validation | Not started | G6 not evaluated |
| 7. Profile and optimize | Rust landscape kernel passed complete numerical equivalence and nonpublishing four-platform candidate CI; R remains canonical/default and distribution/adoption remains open | G7 not evaluated |
| 8. Rewrite manuscript | Not started | G8 not evaluated |
| 9. Release/archive | Not started | G9 not evaluated |

## 10. Revision history

| Version | Date | Summary |
|---|---|---|
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
