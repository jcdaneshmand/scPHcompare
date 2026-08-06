# scPHcompare Modernization and Resubmission Plan

## Document control

| Field | Value |
|---|---|
| Status | Draft for owner and author-team review |
| Version | 0.3.8 |
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
- **Status:** MV-01 complete; pilot contracts/methods/manifest frozen; MV-02 implementation next
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
| D-013 | 2026-08-05 | Use all active persistence-landscape levels with exact or error-controlled integration as the corrected primary target; retain `k=1` only as paper sensitivity; do not activate yet | Eight-stratum audit rejected a universal cap/fixed uniform point count; owner confirmed that the dissertation-aligned definition should govern corrected work | Project-owner approved; reference validation complete; diagram eligibility, production error policy, and broader author-team confirmation pending |
| D-014 | 2026-08-05 | Adopt `landscape_reference_v1` as a non-default correctness oracle; do not activate it or begin Rust work | Analytical and representative benchmarks show exact/adaptive agreement and compact output, while large-diagram certification and independent exact-distance validation remain open | Implemented locally; production engine decision pending |
| D-015 | 2026-08-05 | Classify all nine audited historical persistence-diagram artifacts and every derived distance/result as scientifically ineligible; allow only labeled historical reproduction or performance stress use | Current and historical code pass feature-by-cell matrices directly to ripserr's row-as-point interface, contradicting the dissertation's cell-by-cell intent | Confirmed by code, documentation, dissertation, and artifact signatures; corrected rerun required |
| D-016 | 2026-08-05 | Retain the R reference as oracle; reject Persim 0.3.8's built-in norm; advance corrected Persim critical-pair integration only as a production prototype; keep fixed-grid GUDHI and Rust out of the primary path | R exact, corrected Persim, and SciPy agree to floating precision; built-in Persim fails sign-crossing cases; controlled scaling identifies a promising route but no eligible diagrams yet exist | Prototype decision accepted; no default activation |
| D-017 | 2026-08-05 | Expand the prospective method into distinct cell-topology and gene-topology views; keep cell topology confirmatory and gene topology secondary until independently validated | Project-owner direction following discovery of the legacy feature-as-point orientation; the two orientations represent different scientifically meaningful objects | Confirmed as roadmap scope; implementation pending MV-01/MV-02 |
| D-018 | 2026-08-05 | Evaluate multiple methods in a staged three-layer hierarchy: within-sample geometry, between-sample topological distance, then sample clustering | Prevents conflation of distances and an unstructured method/parameter search | Confirmed planning rule |
| D-019 | 2026-08-05 | Admit cell/gene fusion only after both views pass independent gates; begin with normalized equal-weight distance fusion and require comparison against both components before advanced multiview methods | Preserves interpretability, controls scope, and prevents outcome-informed fusion tuning | Confirmed planning rule; fusion remains exploratory |
| D-020 | 2026-08-05 | Freeze `cell_topology_v1` for the pilot as 384 matched cells represented in 30 shared PCs from a common 500-gene panel, with Euclidean cell geometry | Matches dissertation intent, controls cell count, avoids raw high-dimensional Euclidean geometry, and keeps both views on one source matrix | G-MV1 technical freeze; MV-02 tests required before data execution |
| D-021 | 2026-08-05 | Freeze `gene_topology_v1` for the pilot as the same 500 genes under Euclidean correlation-chord geometry across the same 384 cells | Creates an intentional metric gene-coexpression view while isolating orientation from feature/cell changes | G-MV1 technical freeze; secondary view only |
| D-022 | 2026-08-05 | Use complete Vietoris-Rips H0/H1 for the primary pilot; prohibit sample-specific corrected thresholds; retain all-level H0/H1 landscape L2 and PAM with label-free stability-selected `k` as the primary downstream hierarchy | Resolves filtration comparability, preserves the approved landscape definition, and prevents arbitrary-distance k-means/oracle-label tuning | G-MV1 technical freeze; resource-gated in MV-03 through MV-05 |
| D-023 | 2026-08-05 | Freeze a deterministic 14-sample existing-data feasibility manifest, five 384-cell seeds, and staged time/memory/storage stop rules | Covers eight tissues, historical cell-count extremes, and bone-marrow assay approaches without presenting the pilot as a biological validation set | Manifest generated and hashed; fresh-QC reconciliation required before PH |

## 9. Status dashboard

| Phase | Status | Gate |
|---|---|---|
| 0. Preserve and establish provenance | Main synchronized; P0-02 complete; P0-03/P0-04/P0-05 in progress | G0 not evaluated |
| 1. Scientific/implementation audit | Landscape target/oracle independently validated; historical diagrams invalidated; MV-01 cell/gene contracts technically frozen; implementation, other distances, statistics, and claim disposition remain | G1 not evaluated |
| 2. Reproducible baseline/repository health | Full testthat suite green; analytical fixture and structured provenance implemented; clean-clone/package-check validation remains | G2 not evaluated |
| 3. Literature/reference/figure audit | Initial literature sample only | G3 not evaluated |
| 4. Redesign primary comparison | MV-01 dual-view contract/method/pilot freeze complete; MV-02 implementation next | G4 not evaluated |
| 5. Expand methods | Staged distance/clustering/fusion hierarchy planned; execution waits for eligible dual-view diagrams | G5 not evaluated |
| 6. Biological/practical validation | Not started | G6 not evaluated |
| 7. Profile and optimize | Corrected critical-pair prototype is promising on controlled cell-scale workloads; eligible-diagram batch profiling remains; Rust gated off | G7 not evaluated |
| 8. Rewrite manuscript | Not started | G8 not evaluated |
| 9. Release/archive | Not started | G9 not evaluated |

## 10. Revision history

| Version | Date | Summary |
|---|---|---|
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
