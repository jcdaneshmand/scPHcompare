# P0-03 Initial Artifact Inventory — 2026-08-03

## Classification policy

| Class | Meaning | Default treatment |
|---|---|---|
| Source/code | Versioned implementation or documentation source | Track and review |
| Raw/input | Input not derived by this repository | Preserve; document access and rights |
| Reference output | Derived output used to validate or support prior results | Preserve until reproducibility is established |
| Publication source | Manuscript/dissertation source artifact | Preserve; decide distribution rights |
| Planning/evidence | Auditable modernization documentation | Track after owner review |
| Confidential | Material restricted from public distribution | Keep ignored and outside releases |
| Cache/staging | Reproducible transient material | Deletion only after explicit review/authorization |
| Unknown | Provenance or use not yet established | Preserve until classified |

## Inventory

| ID | Path/group | Class | Size/count | Integrity/provenance | Retention/action | Status |
|---|---|---|---:|---|---|---|
| A-001 | `R/`, `man/`, `tests/`, package root metadata | Source/code | 77 tracked files total across repository | Git history at local HEAD | Preserve; reconcile cached upstream before further edits | Classified |
| A-002 | `inst/extdata/inputs/metadata_MultiTissueAnalysis.csv` | Raw/input metadata | 19,820 bytes | SHA-256 in protected manifest; tracked | Preserve; map to source datasets and data dictionary during P1-08 | Classified, provenance incomplete |
| A-003 | `inst/extdata/inputs/custom_iteration_inputs.R` and templates/config | Source/configuration | 3 configuration files plus toy placeholder | Tracked | Preserve; test against reconstructed paper run | Classified |
| A-004 | `Result Tables/*.csv` | Reference output/publication result | 4 files; 1,188,636 bytes total | SHA-256 in protected manifest; tracked | Preserve until BioRxiv v2 reconstruction and replacement outputs are approved | Classified, run provenance incomplete |
| A-005 | `docs/Dissertation_SubmissionReady_October.pdf` | Publication source | 6,137,261 bytes | SHA-256 recorded | Preserve immutable and local; explicitly Git-ignored by owner decision | Classified/protected |
| A-006 | `docs/Jonah-BioRxiv_v2.pdf` | Publication source | 5,137,625 bytes | SHA-256 recorded | Preserve immutable and local; explicitly Git-ignored by owner decision | Classified/protected |
| A-007 | `PROJECT_PLAN.md`, `docs/PROJECT_EVIDENCE.md`, `docs/plans/`, `docs/audits/` | Planning/evidence | Untracked public documentation | Created from source review and owner decisions | Review, then commit intentionally as one documented scope | Classified |
| A-008 | `docs/private/` | Confidential | 2 files; 15,884 bytes | Ignored; private hashes recorded separately | Preserve locally; verify exclusion continuously; never release without authority | Classified/protected |
| A-009 | `renv.lock`, `renv/activate.R`, `.Rprofile` | Reproducibility control | Tracked; lockfile 523,197 bytes | Git history | Preserve; validate on clean environment in Phase 2 | Classified |
| A-010 | `renv/staging/` | Cache/staging | 18 ignored entries; approximately 1.2 MiB | Incomplete install staging appearance | Retain for now; candidate for later explicit cleanup | Classified |
| A-011 | `example_run.r` | Source/example candidate | 333 bytes; untracked | References nonexistent `metadata.csv`; historical purpose uncertain | Preserve untracked on hold; assess provenance and usefulness before correcting, publishing, or deleting | Classified/hold |
| A-012 | Full single-cell expression objects and paper-run intermediate PH objects | Unknown/not present in repository | Not found in current repository inventory | Location and availability unknown | Locate during P0-05; do not infer they are lost | Open |
| A-013 | Original figure source data/scripts and manually assembled assets | Unknown/partial | Not established | Dissertation/preprint figures inventoried conceptually only | Locate and classify during P0-07 | Open |

## Initial conclusions

1. The repository is small and source-oriented; the two untracked PDFs account for most local non-Git storage.
2. The four tracked result tables are important legacy reference outputs and should not be replaced until their run provenance is reconstructed.
3. Raw single-cell matrices and the complete historical intermediate-output set are not visible in this repository snapshot. Their external locations remain a Phase 0 question.
4. The ignored `renv/staging` tree looks disposable, but this audit does not authorize deletion.
5. The current example script may be disposable, but its purpose has not been established. It remains untracked and preserved pending a provenance/usefulness check; no deletion is authorized.

## P0-03 disposition

`in_progress`: all currently visible top-level artifact groups are classified, but P0-05 and P0-07 must locate or explicitly mark missing the paper-run inputs/intermediates and figure sources before P0-03 can be considered exhaustive.
