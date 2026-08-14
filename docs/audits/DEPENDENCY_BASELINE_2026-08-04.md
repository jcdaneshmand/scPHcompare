# Dependency baseline checkpoint

Date: 2026-08-04

Scope: authoritative `renv` lockfile selection, isolated Ubuntu restoration,
strict package validation, and installed-package persistent-homology smoke test

Branch: `codex/phase-0-audit-foundation`

Preceding package-health checkpoint: `4a1028f`

Dependency baseline commit: `81ca5a99cd442af241476cff7533a973c3e39922`

No push or remote mutation was performed.

## Outcome

The tracked lockfile, rather than the drifted local project library, is the
authoritative dependency baseline. Its existing scientific package versions
were preserved while its package set was reconciled to the complete declared
dependency closure in `DESCRIPTION`, including `Suggests` needed for strict
checks and optional supported workflows.

The resulting baseline has 265 package records and SHA-256:

```text
73EA4688D292661BE7847419919E12057B9F7188AF739D6BB938B2E35C6A0150
```

An exact isolated restore of this lockfile completed under Ubuntu WSL and R
4.4.1. A strict source build and `R CMD check` then reported:

```text
Status: OK
0 errors | 0 warnings | 0 notes
```

The installed-package PH/provenance fixture also completed with Harmony 1.2.4,
`ripserr` 0.3.0, three H0 features, one H1 feature, and two provenance CSVs.

## Authority decision

The preceding package-health audit found six version/source mismatches between
the tracked lockfile and the working project library and 29 stale packages in
the library. Updating the lock from that library would therefore have silently
changed the reproducibility target. Instead, this sprint:

1. restored the preceding lock exactly into an isolated library;
2. verified that scientific stack with a strict package check and real PH
   fixture;
3. derived the full dependency closure from `DESCRIPTION` against that verified
   environment; and
4. wrote a new lockfile without changing any overlapping package version,
   source, repository, or GitHub SHA.

The old lock contained 232 records. The new lock contains 265 records:

- 54 records were added to make testing and declared optional workflows
  reproducible;
- 21 unused records were removed; and
- 0 of the records present in both locks changed version or source identity.

The principal additions are the `testthat`/`mockr` toolchain and the exact
closures for `batchelor` 1.20.0, `BiocSingular` 1.20.0, `rliger` 2.1.0, and
`scCustomize` 3.0.1. The key scientific versions remain Harmony 1.2.4 from
GitHub commit `200114cfad09bee95543b478b3642c133d59159c`, `ripserr` 0.3.0,
Seurat 5.3.0, SeuratObject 5.0.2, TDA 1.9.4, and TDAstats 0.4.1.

The removed records were:

```text
backports, blob, broom, cellranger, conflicted, DBI, dbplyr, dtplyr,
gargle, googledrive, googlesheets4, haven, ids, modelr, readxl, rematch,
reprex, rvest, selectr, tidyverse, uuid
```

These were remnants of the removed umbrella `tidyverse` dependency and are not
in the current complete declared dependency closure.

## Restore and validation evidence

Environment:

- Ubuntu 22.04.4 LTS under WSL;
- R 4.4.1;
- Bioconductor 3.19;
- `renv` 1.1.5; and
- isolated libraries beneath the Git-ignored `tmp/` validation workspace.

Preceding-lock verification:

1. The 232-record lockfile with SHA-256
   `A2E6C4EFEF59A1E449F09AA43C6AE843596A611F1B9CBD814FACCD0ACDA48D2D`
   was restored nontransactionally into a clean isolated library in 253.3
   seconds.
2. Comparison with the restored and R base/recommended libraries found no
   missing package and no version/source mismatch.
3. The complete declared `Suggests` closure was hydrated at its exact cached
   versions without upgrading the runtime stack.
4. A vanilla strict check with `_R_CHECK_FORCE_SUGGESTS_=TRUE` completed in
   1,238 seconds with `Status: OK`.
5. A separately installed package completed the real four-point PH/provenance
   fixture with `h0=3 h1=1 provenance_files=2`.

Final-lock verification:

1. The 265-record final lock was restored into a second clean isolated library
   in 299.4 seconds.
2. Across that target plus R base/recommended libraries, all 264 expected
   non-`renv` records were present at exact locked versions. There were no
   unexpected target packages and none of the 21 removed legacy packages was
   present. The running exact `renv` 1.1.5 supplied the excluded `renv` record.
3. A vanilla source build and strict `R CMD check --no-manual` with
   `_R_CHECK_FORCE_SUGGESTS_=TRUE` completed in 1,562.8 seconds with `Status:
   OK`; all examples and tests passed.
4. A fresh source installation followed immediately by the installed-package
   fixture reported:

   ```text
   scPHcompare=0.1.0
   harmony=1.2.4 ripserr=0.3.0
   ph_status=completed h0=3 h1=1 provenance_files=2
   ```

## Resolver findings and operational cautions

Setting `renv/settings.json` to traverse `Suggests` globally was tested and
rejected. With `renv` 1.1.5 it caused resolver failures involving
`BiocVersion` and then `BiocManager`, even though the records and repositories
were present. The settings file was restored exactly to its previous state.
Optional/test packages are therefore explicit lockfile records, but future
lock refreshes must derive and pass the declared package set deliberately
rather than enabling global `Suggests` traversal until that resolver behavior
is understood or a verified `renv` upgrade is adopted.

The first final restore required seeding the exact cached `BiocVersion` record
before `renv` completed resolution. The resulting library was then compared
record by record and passed the strict check. This cache-seed requirement is
environment bootstrap debt and should be tested in a future cold-cache/CI
restore before claiming one-command reconstruction on a new machine.

`R CMD build` copies the source tree before applying `.Rbuildignore`. Because
the disposable libraries were placed beneath the repository, the final build
spent substantial time copying ignored `tmp/` content. Future isolated restore
libraries should be placed beside or outside the repository; `tmp/` remains
Git-ignored only as a safety net for local validation artifacts.

## Preserved and excluded artifacts

- `docs/Dissertation_SubmissionReady_October.pdf` remains Git-ignored and
  untracked.
- `docs/Jonah-BioRxiv_v2.pdf` remains Git-ignored and untracked.
- `docs/private/` and ad hoc `provenance/` output remain Git-ignored.
- `example_run.r` remains deliberately untracked and unchanged, with SHA-256
  `D205E5EAAFEA27D4AED51146DB1C20CDEC62E56506CDEA239BF4727EDB736185`.
- Disposable validation libraries and scripts under `tmp/` are not part of the
  baseline and are removed after their evidence is recorded.

## Next action

Add a cold-cache dependency restore to continuous integration or a documented
bootstrap procedure, then begin the Phase 2 reproducible-baseline work on
determinism, a maintained toy dataset, and benchmarkable stage-level runtime
measurements. Rust feasibility work should begin only after those measurements
identify a stable computational bottleneck and numerical-equivalence contract.
