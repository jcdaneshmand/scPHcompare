# `example_run.r` Provenance and Usefulness Assessment — 2026-08-03

## Identity

| Field | Value |
|---|---|
| Path | Repository root `example_run.r` |
| Size | 333 bytes |
| Created/modified | 2026-01-05 10:22:22 local filesystem time |
| SHA-256 | `D205E5EAAFEA27D4AED51146DB1C20CDEC62E56506CDEA239BF4727EDB736185` |
| Git state | Untracked and not ignored |

## Provenance search

- No commit on any local or fetched remote ref contains `example_run.r` at that path.
- Direct path history is empty.
- The file was created locally on the same day as several January 2026 package-maintenance commits, but no commit establishes its author or purpose.
- Its `run_unified_pipeline()` call resembles the README quick start rather than a unique analysis configuration.

## Validity/usefulness

- The script references `inst/extdata/inputs/metadata.csv`, which does not exist.
- The repository instead contains `metadata_MultiTissueAnalysis.csv`, but substituting it would not make the run reproducible because none of its 127 `.RData` paths resolve on this machine.
- The custom-iteration file referenced by the script exists but contains placeholder `path/to/...` values.
- The README already contains maintained quick-start and toy-data examples.
- The script requests Seurat and Harmony together plus all post-processing modules; this may be useful as a memory aid, but it is not currently a validated example or the recovered paper command.

## Decision

Classify as **likely scratch launcher; hold**.

- Preserve untracked for now; do not stage or publish.
- Do not delete without explicit owner approval.
- Do not repair it speculatively by changing the metadata filename.
- If Phase 2 produces a clean toy smoke command, prefer a tested maintained example derived from that workflow rather than reviving this file unless new provenance is found.
