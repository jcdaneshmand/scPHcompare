# P0-01 Canonical Repository State — 2026-08-03

## Audit identity

| Field | Recorded value |
|---|---|
| Repository | `E:\Repositories\Jonah\PH Pipeline Repo\scPHcompare` |
| Canonical public repository | `https://github.com/jcdaneshmand/scPHcompare.git` |
| Audit date | 2026-08-03, America/New_York |
| Audit mode | Read-only Git/filesystem inspection |
| Branch | `main` |
| Local HEAD | `64ec4fd5e871e7a0d0c30db8146b766c5579726f` |
| Description | `v0.9.9-60-g64ec4fd` before adding the current untracked planning documents |
| Cached upstream | `origin/main` at `3910b15ce6d3f88197847a8aba94b184e7d9c034` |
| Cached divergence | Local-only commits: 0; cached-upstream-only commits: 2 |
| Cached upstream last updated | 2026-01-05 16:27:09 -0500 |

## Working-tree state

The tracked working tree had no modified or deleted tracked files. Git reported:

| Class | Count | Notes |
|---|---:|---|
| Tracked files | 77 | Includes R sources, package metadata, tests, documentation, result tables, and `renv` controls |
| Public untracked files | 18 | Three root files and fifteen files below `docs/` |
| Ignored files | 20 | Two confidential documentation files and eighteen `renv/staging` entries |
| Declared submodules | 0 | No `.gitmodules` file |
| Git LFS inspection | Unavailable | `git-lfs` was not installed; no LFS declaration was observed during this audit |

### Public untracked groups

| Group | Files | Classification/action |
|---|---:|---|
| Planning and evidence | 14 | New auditable project documentation; review before intentional commit |
| Source PDFs | 2 | Publication/dissertation source documents; protect and decide public distribution policy before commit |
| `.gitignore` | 1 | Confidentiality control; review and commit intentionally |
| `example_run.r` | 1 | Candidate example, currently invalid because it references missing `inst/extdata/inputs/metadata.csv` |

The actual metadata file present is `inst/extdata/inputs/metadata_MultiTissueAnalysis.csv`. The example must be corrected and exercised before publication.

### Ignored groups

| Group | Files | Classification/action |
|---|---:|---|
| `docs/private/` | 2 | Confidential editorial/reviewer material; retain locally, keep ignored, exclude from release and public logs |
| `renv/staging/` | 18 | Incomplete/package-install staging cache; reproducible/disposable candidate, but no deletion is authorized by this audit |

## Tracked-file distribution

| Top-level group | Tracked files |
|---|---:|
| Root | 9 |
| `.github` | 1 |
| `R` | 15 |
| `man` | 34 |
| `tests` | 5 |
| `inst` | 6 |
| `renv` | 3 |
| `Result Tables` | 4 |

## Repository storage

Git reported 970 loose objects using approximately 3.31 MiB and no packs or garbage. The largest non-Git files are the dissertation PDF (5.853 MiB), preprint PDF (4.900 MiB), tracked result tables (0.639 and 0.463 MiB), `renv.lock` (0.499 MiB), and ignored `renv/staging` package files.

## Controls and limitations

- No fetch, pull, merge, reset, checkout, deletion, stage, commit, or push was performed.
- Git required a per-command `safe.directory` override because filesystem ownership differs from the execution identity. Global Git configuration was not changed.
- The upstream comparison is explicitly **cached**, not current. The cached ref is approximately seven months old and must be refreshed before synchronization.
- File counts are point-in-time evidence. Re-run the commands in the synchronization plan immediately before any mutation.

After writing the five public Phase 0 audit files, the public-untracked count increased from 18 to 23; adding the ignored confidential hash manifest increased the ignored count from 20 to 21. This expected post-audit state does not alter the captured pre-documentation snapshot.

## P0-01 disposition

`needs_review`: the state capture is complete and reproducible, but owner review is required before the subplan status becomes `complete`.
