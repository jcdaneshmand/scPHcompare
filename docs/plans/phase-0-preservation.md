# Phase 0 Subplan — Preservation and Provenance

## Objective and boundary

Establish a loss-resistant canonical baseline. This phase changes no scientific behavior and performs no blind pull, merge, deletion, or artifact relocation.

## Task ledger

| ID | Task | Required evidence/output | Acceptance test | Status |
|---|---|---|---|---|
| P0-01 | Record repository state | Branch, HEAD, remotes, status, tracked/untracked/ignored inventory, submodules/LFS if any | Snapshot is dated, repeatable, and contains no confidential file contents | `needs_review` |
| P0-02 | Assess local/remote divergence | Commit graph and per-commit summaries for local and current remote tips; collision list for untracked files | A written synchronization procedure identifies every overwrite/conflict risk before mutation | `complete` — owner-authorized fast-forward verified at `3910b15` |
| P0-03 | Classify artifacts | Ledger labeling each large/important file as source, raw input, intermediate, reference output, publication result, cache, or disposable | Every irreplaceable or expensive artifact has an owner, retention class, and location | `in_progress` — visible groups classified; historical run/figure artifacts remain open |
| P0-04 | Protect artifacts | SHA-256 manifest, size, modification date, provenance, backup status | Checksums verify; confidential paths are ignored; no raw source was modified | `in_progress` — initial public/private manifests created; backup status unresolved |
| P0-05 | Reconstruct BioRxiv v2 run | Dataset/configuration/package/output/figure map with known gaps | Each reported result maps to an output and configuration, or is explicitly marked unrecoverable | `in_progress` — initial trace identifies count/path/configuration gaps |
| P0-06 | Document lineage | Public note explaining `PH_ClusteringApp` → `scPHcompare`, without copying obsolete artifacts blindly | Owner confirms canonical repository and historical relationship | `not_started` |
| P0-07 | Inventory figure provenance | Figure/panel ledger: code-generated, manual, AI-assisted, mixed, or unknown; source path and regeneration status | Every paper/dissertation figure has a status and next action | `not_started` |
| P0-08 | Evaluate G0 | Gate record with evidence links and author-team decision | All G0 checklist items have explicit pass/fail/defer disposition | `not_started` |

## Execution details

- Capture read-only state before synchronization.
- Do not treat file existence as provenance; use checksums and configuration identifiers.
- Keep the original PDFs immutable.
- Store only generalized reviewer-derived tasks publicly; keep original correspondence in `docs/private/`.
- Record missing provenance as a finding rather than reconstructing it from memory.

## Gate G0 checklist

- [ ] No important local work is unclassified.
- [ ] Irreplaceable/expensive artifacts have verified checksums and retention decisions.
- [ ] Confidential materials are excluded from Git and release inputs.
- [ ] BioRxiv v2 outputs/configurations are located or gaps are documented.
- [ ] Figure provenance is recorded.
- [ ] Author team confirms the canonical historical baseline.

## Required gate record

Record date, participants, evidence reviewed, exceptions, decision (`approved`, `conditional`, or `not approved`), and next tasks in `WORK_LOG.md` and the main decision log.
