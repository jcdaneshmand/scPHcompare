# scPHcompare merge-readiness decision dossier

Date: 2026-08-13
Scope: draft PRs #111 through #118
Status: owner decision package; no merge, retarget, release, or push authorized

## Executive recommendation

The eight-PR publication stack is ready for an owner-authorized, sequential
merge procedure. The live heads match the audited ledger, all retained R checks
pass, the P08 four-host Rust candidate certification passes, every PR remains a
draft, `main` remains unchanged, and no private document, PDF, result binary,
credential, or local example entered the stack.

Recommended decisions:

1. Accept the dual-view scientific architecture and corrected persistence-
   landscape contract summarized below.
2. Accept P08's Windows PH-child launcher repair as a cross-platform transport
   fix that preserves the PH call and scientific parameters.
3. Merge the Rust source and certification work only as an optional,
   explicitly selected candidate. R remains canonical and default; artifact
   release, installation, and production adoption remain closed.
4. Keep the final remote-publication audit, run ledger, this dossier, and the
   eventual merge ledger out of P08. Publish them in a separate documentation
   PR after the eight-slice merge sequence.
5. Preserve the current public-preprint author order and existing software
   credit. Begin a low-friction credit discussion now, but make final CRediT
   roles, ORCIDs, affiliations, creator/contributor classification, and author
   approval mandatory before a DOI deposit, tagged release, or revised-
   manuscript submission—not before repository-cleanup merges.

Approval of this dossier is not itself authorization to merge. The owner must
separately say to begin the merge sequence.

## Decisions requested from the project owner

| ID | Decision | Recommendation | Consequence of approval |
|---|---|---|---|
| D1 | Is the scientific architecture faithful to the intended dissertation modernization? | Approve | Cell and gene topology remain distinct views; corrected landscapes become the shared comparison object. |
| D2 | Does the eight-PR organization provide an acceptable public history? | Approve | Preserve the dependency-ordered stack and merge it with merge commits, one PR at a time. |
| D3 | Should P08's Windows PH launcher fix be accepted? | Approve | Keep the same `ripserr::vietoris_rips()` call and parameters, but pass paths and numbers as child-process arguments. |
| D4 | Should the Rust candidate enter `main`? | Approve as candidate only | Source, wrapper, tests, and certification enter the repository; R remains canonical/default and no binary is released. |
| D5 | Where should the final audit and ledger be published? | Separate P09 documentation PR | Avoid changing P08's accepted head and rerunning its expensive matrix; add final merge outcomes after the stack lands. |
| D6 | When should coauthor roles be discussed? | Start now; finalize before release/DOI/manuscript | Repository cleanup may proceed while formal CRediT and archival metadata remain an explicit later gate. |
| D7 | May the sequential merge procedure begin? | Decide only after reviewing D1-D6 | A separate explicit instruction triggers the merge sequence. |

Suggested owner response:

> I approve D1-D6 and authorize the sequential merge procedure in D7.

The owner may approve D1-D6 while withholding D7 if more time is desired.

## Scientific architecture requiring human judgment

### Observation views

- `cell_topology_v1` is the corrected confirmatory view: cells are points in a
  shared cell-PC coordinate system, so topology describes cell-population
  organization.
- `gene_topology_v1` is a prespecified secondary view: the common named genes
  are points under a controlled correlation-chord geometry, so topology
  describes within-sample gene-expression relationships.
- `legacy_gene_view_v0` remains historical reproduction/stress material. It is
  not treated as a valid implementation of the new gene view.
- Cell and gene results must stand independently. Fusion is secondary and
  exploratory until its weights and validation are prediction-locked.

This resolves the historical observation-orientation defect without declaring
the dissertation valueless: it preserves gene-side topology as a meaningful,
newly controlled view and adds the intended cell-side view as the confirmatory
object.

### Corrected persistence-landscape contract

The accepted target is `full_l2_error_controlled_v1`:

- use every finite, positive-persistence interval;
- count essential/invalid intervals in provenance, but exclude the essential H0
  interval from landscape evaluation;
- retain every consecutive active landscape level with zero padding when
  depths differ;
- calculate H0 and H1 separately and report them as primary components;
- integrate squared level differences exactly or with independently validated,
  error-controlled integration over filtration scale;
- impose no universal uniform grid or landscape-level cap;
- use streamed/chunked computation instead of unnecessary dense
  grid-by-level arrays;
- permit `sqrt(d_H0^2 + d_H1^2)` only as a secondary descriptive combination,
  accompanied by both components and the H1 squared-distance contribution;
- treat `k=1`, fixed-grid, scale-normalized, or capped landscapes only as
  explicitly named sensitivities or historical compatibility calculations.

This is aligned with the dissertation's latest all-level L2 intent while making
the grid, dimension, provenance, and error policy explicit. It is the contract
implemented by the corrected R engine and evaluated by the optional Rust
kernel.

### What CI cannot approve

CI proves syntax, package health, frozen numerical checks, deterministic build
properties, fixtures, and candidate-accelerator absence/presence behavior. It
cannot decide that the biological interpretation, view hierarchy, H0/H1
reporting hierarchy, authorship, or future manuscript claims are appropriate.
Those human decisions are represented by D1, D4, and D6.

## Live GitHub snapshot

Snapshot taken 2026-08-13 using authenticated account `jcdaneshmand`.

- Repository: public, active, default branch `main`.
- Live `main`: `3910b15ce6d3f88197847a8aba94b184e7d9c034`.
- GitHub merge commits are enabled. Squash and rebase merges are also enabled,
  but are prohibited for this audited stack.
- `main` has no GitHub branch-protection rule. GitHub does not enforce an
  external approval; human review is a project-governance and scientific gate.
- All eight PRs are open drafts, report `MERGEABLE`, and have no recorded
  review decision.

| Slice | PR | Base -> head | Current head | Retained checks |
|---|---:|---|---|---|
| P01 | #111 | `main` -> P01 | `826dde01e6baaca653115ea9f8e0683e176e86f7` | R success, run 31747769123 |
| P02 | #112 | P01 -> P02 | `aeee92d7fa13f1c167a5f838199b109ef9951a7b` | R success, run 31747768694 |
| P03 | #113 | P02 -> P03 | `8feb1659f88aa04223d400f741ac5d2ef23d5f33` | R success, run 31747769836 |
| P04 | #114 | P03 -> P04 | `1fac88fe19a05796a498d4271bf6462d016e94ea` | R success, run 31747768121 |
| P05 | #115 | P04 -> P05 | `36164a9544777a41028c61124a72e357b03ebb00` | R success, run 31747768255 |
| P06 | #116 | P05 -> P06 | `625910ee07ce50234ce9b93fca6231cc9e3a2005` | R success, run 31747769410 |
| P07 | #117 | P06 -> P07 | `17030cc2fe060a9abe661e5ea2408f56b9fd3771` | R success, run 31747768559 |
| P08 | #118 | P07 -> P08 | `98b3cc1f6b835d60bf92fdcd7cb273841a48855e` | R success, run 31750612007; five Rust jobs success, run 31750611952 |

P02-P07 also display cancelled checks from superseded stacked-parent updates.
Each has a later retained successful run. The cancellations are concurrency
diagnostics, not failures of the current candidate heads.

## One-page review summaries

### P01 / PR #111 — foundation

Purpose: establish package health, dependency locking, reproducible toy and
realistic fixtures, PH provenance, H0/H1 execution observability, the corrected
landscape specification/reference engine, diagram eligibility findings, and
the dual-view topology roadmap/contract.

Primary risk: this is the broad conceptual foundation (138 files). A mistaken
view or landscape definition would propagate through the entire stack.

Recommendation: approve. The landscape target matches the latest dissertation
intent and explicitly corrects the prior grid, level, and observation ambiguity.

### P02 / PR #112 — typed topology and reproducible execution contracts

Purpose: turn cell/gene view definitions into typed constructors and validators;
complete the corrected pilot, topological-distance contract, statistical
benchmark plan, resource-safe scaffold, and normalization/cache identity gate.

Primary risk: axis identities, fitted transformations, or cache identities
could silently mix views or leak information.

Recommendation: keep cell topology confirmatory and gene topology secondary;
they answer different questions and should remain independently interpretable.

### P03 / PR #113 — label-closed cell pipeline

Purpose: implement label-closed cell folds, bounded profiling, full cell PH,
cell-landscape distances, retrieval inputs/evaluation, integration resource
gating, and integrated coordinate production.

Primary risk: outcome leakage and sample/cell identity drift across folds,
representations, or projections.

Recommendation: approve label-closed evaluation as the standard for revised
work; it is stronger than choosing settings with known biological labels.

### P04 / PR #114 — integrated topology and benchmark specification

Purpose: produce integrated cell PH and landscape distances, retrieval inputs
and evaluation, compare representations under locked rules, and define the
post-retrieval benchmark and clustering resource gates.

Primary risk: treating integration-induced distance change as biological
improvement without separating technical and biological effects.

Recommendation: keep integrated and unintegrated representations as separate
strata and compare them through locked evaluations.

### P05 / PR #115 — full distance, clustering, and PC20 execution

Purpose: produce label-closed full distance matrices and clustering artifacts,
freeze and execute prediction-locked outcomes, establish robustness admission,
and complete the PC20 workload.

Primary risk: selecting a favorable distance, cluster count, or configuration
after viewing outcomes.

Recommendation: approve the execution history as methods/audit evidence, not as
authorization of a final biological claim.

### P06 / PR #116 — prediction-locked robustness extensions

Purpose: evaluate PC20, cosine, and nested-192 configurations through frozen
continuation gates and prediction-locked outcome engines.

Primary risk: robustness analyses becoming post-hoc model search.

Recommendation: retain null, negative, and stopped extensions; the sequential
prefreeze/continuation trail strengthens the resubmission.

### P07 / PR #117 — corrected landscape matrices and consumers

Purpose: finish nested-192/256 work, reconcile the public landscape contract,
implement corrected public APIs, repair adaptive numerical integration,
validate realistic compatibility, add corrected artifacts and consumers, and
accept the complete corrected landscape matrix corpus.

Primary risk: this is the key mathematical slice. A silent fallback to legacy
`k=1`, fixed-grid, mixed H0/H1, or uncontrolled integration would invalidate
downstream distances and clustering.

Recommendation: approve the all-active, separate-H0/H1, exact/error-controlled
contract and keep every legacy behavior explicitly labeled.

### P08 / PR #118 — optional acceleration and publication evidence

Purpose: freeze stability resampling, benchmark a corrected-Persim oracle, add
and fully evaluate a bounded Rust landscape kernel, specify cross-platform
distribution, certify the candidate on four hosts, and define the lean-code /
separate-evidence publication architecture.

Primary risk: confusing numerical acceptance with production adoption.
Secondary risk: the late Windows launcher repair touches core PH orchestration.

Recommendation: merge as transparent optional research infrastructure. R stays
default, no binary enters the repository, and no release is authorized.

## Focused P08 review

### Rust candidate

Accepted evidence:

- one pair and one homology dimension are processed per call;
- all finite intervals and all consecutive active levels are used;
- H0 and H1 remain separate calls/results;
- exact piecewise-linear squared-L2 integration has no grid or level cap;
- there are no external Rust crates; Rust 1.97.1 is pinned;
- all 408 existing-data references passed: 318/318 exact and 90/90
  adaptive-certified;
- 408 reversed pairs were bit-identical and 112 self-comparisons exactly zero;
- four native host rows passed format, unit, strict-Clippy, two-clean-build
  identity, C ABI, R fixtures, dependencies, and R checks;
- R package checks passed with the candidate absent and present;
- no compiled binary is tracked or published.

Still closed: Rust-by-default, binary/release publication, a network installer,
the Linux glibc-2.17 release baseline, signed/attested release manifests,
cache/API adoption, and removal of canonical R.

Recommendation: merge P08 as research infrastructure and an optional candidate,
not as production adoption.

### Windows PH-child launcher repair

The pre-existing function wrote Windows paths directly inside an `Rscript -e`
source string. A path beginning with `C:\Users` can be parsed by R as a `\U`
escape before `ripserr` runs. Commit
`98b3cc1f6b835d60bf92fdcd7cb273841a48855e` changes only
`R/PH_Functions.R` (14 insertions, 7 deletions):

- dataset path, output path, `DIM`, and threshold become process arguments;
- the child reads them with `commandArgs(trailingOnly = TRUE)`;
- the same `ripserr::vietoris_rips(dataset, max_dim, threshold,
  return_format = "mat")` call remains;
- timeout, threshold selection, persistence-diagram serialization, monitoring,
  and provenance behavior remain unchanged.

The WSL real-child regression passed 27/27, and the final Windows host job
passed accelerator-absent and accelerator-present package checks.

Recommendation: accept. This is a portability correction at the process
boundary, not a change to the topological calculation.

## Private and generated artifact exclusion

The exact diff from live `main` to P08 adds 699 files. Checks passed:

- zero added path matches `docs/private`, `tmp`, `results`,
  `example_run.r`, PDF, Word, RDS/RData, archive, executable, or
  shared-library extensions;
- `git diff --numstat` reports no binary-diff row;
- the largest added blob is the 448,408-byte derived-evidence member manifest,
  not the 161,963,473-byte evidence ZIP;
- no GitHub token, AWS access-key pattern, or private-key header was found;
- `docs/Dissertation_SubmissionReady_October.pdf` and
  `docs/Jonah-BioRxiv_v2.pdf` exist locally and are Git-ignored;
- reviewer correspondence is not tracked;
- the local worktree contains only the intentionally untracked
  `example_run.r`;
- no generated Rust DLL, dylib, or shared object is tracked or released.

The 1,136 generated evidence members remain in a separate deterministic
161,963,473-byte ZIP with SHA-256
`0cd4fa24c407cfe3611c1e424889f38c15bb66d1784039c3d1291325155afa62`,
pending a separate DOI/deposit decision.

## Audit publication decision

Recommendation: use a separate P09 documentation PR after P08 lands.

1. The remote ledger and CI-remediation audit describe events after P08's
   scientific candidate was frozen.
2. Adding them to P08 would change its accepted head and trigger another
   expensive R/Rust matrix without changing science.
3. P09 can include final merge commit IDs, retarget/check outcomes, this owner
   decision record, and the closed-action statement.
4. P08 stays focused on accelerator/evidence architecture.

Until P09 is created, the audit, run ledger, and this dossier remain local on
`codex/phase-0-audit-foundation`; they are not public GitHub artifacts.

## Sequential merge and rollback procedure

### Non-negotiable rules

- Use GitHub merge commits only. Do not squash, rebase-merge, force-push, or
  delete a candidate branch during the sequence.
- Process exactly one PR at a time.
- Do not merge a child while its base names an unmerged parent branch.
- After each base change, wait for current checks and revalidate the slice.
- Stop on any unexpected file, head change, conflict, failed required check,
  altered landscape contract, or unexpected diff.

### Procedure

1. Re-read live `main`, PR bases/heads, draft states, bodies, and checks.
2. Mark only #111 ready and merge it into `main` with a merge commit. Record
   the pre-merge main SHA and returned merge commit. Do not delete P01.
3. Confirm the new `main` tree equals P01's accepted tree
   `bc581aa0ee1c4a6ddff50a6bcd726f322dbb729f`.
4. Retarget #112 to `main`. Confirm its head and accepted slice are unchanged,
   wait for successful checks, then merge with a merge commit.
5. Repeat retarget, exact-diff validation, check, and merge for #113 through
   #118 in order.
6. Confirm final `main` tree equals P08's accepted tree
   `1fbf3350a4f15fe0209cc0a9eb1c70e8ef3fd7f6` and rerun the
   sensitive-path/closed-action audit.
7. Create a separate draft P09 documentation PR with final audit, publication
   ledger, merge ledger, and owner decision.
8. Keep tags, releases, binaries, DOI/Zenodo actions, defaults, new
   calculations, and manuscript claims closed.

### Rollback points

- Before merge: leave the PR open/draft and stop.
- After retarget but before merge: retarget the child back to its recorded
  parent if needed, record why, and stop.
- After an erroneous merge: do not rewrite `main`. Prepare a forward
  `git revert -m 1 <merge-commit>` for explicit owner approval; never execute
  it automatically.
- Preserve candidate branches until the complete sequence and P09 are accepted.

## Credit and authorship checklist

Current public metadata preserves the bioRxiv order:

1. Jonah Daneshmand
2. Julia H. Chariker
3. Akshitkumar Mistry
4. Eric C. Rouchka

`CITATION.cff` lists all four as software authors and cites bioRxiv DOI
`10.1101/2025.07.24.666637`. `DESCRIPTION` identifies Jonah as R package
author/maintainer and Julia, Dr. Mistry, and Dr. Rouchka as contributors. The
Zenodo evidence metadata deliberately marks roles/order as provisional.

Before a tag, DOI, release, or manuscript submission, obtain from every author:

- preferred published name, ORCID, and affiliations;
- software-author/contributor classification;
- evidence-dataset creator/contributor classification;
- revised-manuscript order;
- CRediT roles (Conceptualization, Methodology, Software, Validation, Formal
  analysis, Investigation, Data curation, Visualization, Supervision, Project
  administration, Resources, Funding acquisition, Writing – original draft,
  and Writing – review & editing);
- approval of licenses, availability text, and claims tied to their work.

Do not infer Dr. Mistry's or Dr. Rouchka's final roles from job titles. Both
retain visible project credit now; exact roles require author-team confirmation
before archival release or resubmission.

## Exact scientific review inventory

The companion `merge-readiness-file-inventory-v1.csv` is the canonical exact
per-slice inventory. It records every current-base-to-current-head path, change
status, review class, and whether the path is landscape-key.

| Slice | Core human-review entry points |
|---|---|
| P01 | `R/landscape_contract.R`; `R/landscape_convergence.R`; `R/landscape_reference.R`; `R/realistic_fixture.R`; `R/PH_Calculation.R`; `R/PH_Functions.R`; `R/unified_pipeline.R` |
| P02 | `R/dual_view_topology.R`; `R/topological_distance_contract.R`; `R/mv05_benchmark_contract.R`; `R/mv05_benchmark_execution.R`; `R/mv05_inductive_mapping.R`; `R/mv05_resource_safe_execution.R` |
| P03 | `R/mv05d1_post_fold_projection.R`; `R/mv05d2_ph_profiling.R`; `R/mv05d3_ph_production.R`; `R/mv05d4_landscape_production.R`; `R/mv05d5_retrieval_inputs.R`; `R/mv05e_retrieval_evaluation.R`; `R/mv05f_integration_gate.R`; `R/mv05g_coordinate_production.R` |
| P04 | `R/mv05h_integrated_ph_production.R`; `R/mv05i_integrated_landscape_production.R`; `R/mv05j_integrated_retrieval_inputs.R`; `R/mv05k_retrieval_evaluation.R`; `R/mv05l_representation_comparison.R`; `R/mv05m_benchmark_gap_gate.R`; `R/mv05n_clustering_gate.R` |
| P05 | `R/mv05o_production_prefreeze.R`; `R/mv05p_distance_production.R`; `R/mv05q_clustering_artifacts.R`; `R/mv05r_outcome_prefreeze.R`; `R/mv05s_outcome_execution.R`; `R/mv05t_robustness_gate.R`; `R/mv05u_robustness_admission.R`; `R/mv05v_robustness_prefreeze.R`; `R/mv05w_launch_readiness.R`; `R/mv05x_configuration_execution.R` |
| P06 | `R/mv05y_robustness_outcome_prefreeze.R`; `R/mv05z_outcome_execution.R`; `R/mv05aa_robustness_continuation_gate.R`; `R/mv05ab_cosine_execution.R`; `R/mv05ac_cosine_outcome_prefreeze.R`; `R/mv05ad_outcome_execution.R`; `R/mv05ae_robustness_continuation_gate.R`; `R/mv05af_nested_execution.R`; `R/mv05ag_nested192_outcome_prefreeze.R`; `R/mv05ah_outcome_execution.R` |
| P07 | `R/landscape_reference.R`; `R/landscape_public_api.R`; `R/corrected_landscape_workflow.R`; `R/corrected_landscape_consumers.R`; `R/mv05an_landscape_public_contract_prefreeze.R`; `R/mv05ap_realistic_compatibility_gate.R`; `R/mv05apr1_realistic_rerun.R`; `R/mv05ar_opt_in_integration_prefreeze.R` |
| P08 | `R/landscape_rust_prototype.R`; `R/PH_Functions.R`; `rust/scph_landscape_kernel/src/lib.rs`; `rust/scph_landscape_kernel/include/scph_landscape_kernel.h`; `.github/workflows/rust-accelerator-certification.yml` |

The companion inventory also includes orchestration scripts, tests, scientific
specifications/audits, publication manifests, package metadata, and workflows
as separate review classes so no decision relies on an informal shortlist.

## Final readiness conclusion

No technical blocker is visible. The remaining gate is informed owner judgment,
followed by an explicitly authorized, observable one-PR-at-a-time merge
sequence. The recommended answer is to approve D1-D6. D7 should be authorized
only when the owner is ready for public `main` to advance.
