# Package health checkpoint

Date: 2026-08-04

Scope: local repository checkpoint, Ubuntu R package validation, namespace and
documentation repair, clean installation, and installed-package PH/provenance
smoke test

Branch: `codex/phase-0-audit-foundation`

Code commit: `d431295ed2260e2d5f2e1115c8f3daa2bebe80e9`

Preceding provenance checkpoint: `1113fea`

No push or remote mutation was performed.

## Outcome

The package builds, installs, loads, runs its examples, and passes its complete
test suite in Ubuntu WSL with R 4.4.1. The final non-CRAN package check reports:

```text
Status: OK
0 errors | 0 warnings | 0 notes
```

A separate clean-library smoke test installed `scPHcompare` into an empty
temporary library, loaded that installed copy, ran a real four-point
`ripserr::vietoris_rips()` child process through `process_and_monitor()`, and
verified:

- attempt status `completed` with exit status 0;
- three finite H0 features and one H1 feature;
- a nonempty persistence diagram;
- one sample-flow CSV and one PH-attempt CSV written by the installed package.

The temporary library and runtime artifacts were deleted automatically when
the validation process exited.

## Baseline and resolution trace

The first full check reported 1 error, 6 warnings, and 2 notes. The main causes
and resolutions were:

| Finding | Resolution |
|---|---|
| Broad imports produced 26 namespace masking messages | Removed the umbrella `tidyverse` import, added missing base imports, and explicitly excluded conflicting exports from remaining package imports |
| `PercentageFeatureSet()` was passed unsupported `perl = TRUE` | Selected hemoglobin features with `grep(..., perl = TRUE)` and passed the exact feature vector; QC thresholds and filtering defaults were unchanged |
| `LICENSE` was not in the DCF form required by `MIT + file LICENSE` | Converted `LICENSE` to DCF and retained the full GitHub-facing text as `LICENSE.md` |
| Non-ASCII source characters triggered portability findings | Replaced punctuation-only characters with ASCII equivalents |
| Seurat v5 layer cleanup called nonexistent `RemoveLayer()` | Used the supported `Layers()` and `LayerData() <- NULL` interfaces |
| Static analysis found undeclared base functions and data-mask symbols | Added explicit base imports and a narrow `globalVariables()` declaration |
| Patchwork helpers were unqualified | Added the dependency and qualified the calls |
| Cross-iteration comparison had lost its initialization and loop scaffold | Restored input validation, output creation, group discovery, logging, reference-cell handling, and iteration traversal from the historical implementation without reintroducing script-level `library()` calls |
| Landscape helper defaults referenced a free `grid_points` symbol | Replaced the invalid defaults with the documented 500-point grid |
| Spectral clustering referenced an undefined `similarity_matrix` | Constructed a Gaussian similarity kernel from the distance matrix using its median positive distance, added input validation, and added a regression test |
| Rd usage, parameter, and example declarations disagreed with code | Reconciled source and Rd documentation; internal helper examples now identify the internal namespace without expanding the public API |
| One Euler example attempted eight worker processes during check | Marked this genuinely heavy example `\dontrun`; computational behavior remains exercised by the test suite |

No scientific default for QC, integration, persistence dimension, filtration,
permutation count, or bootstrap count was changed in this checkpoint.

## Verification evidence

Environment:

- Ubuntu 22.04.4 LTS under WSL;
- R 4.4.1;
- repository `renv` library available for dependencies.

Checks:

1. `devtools::check(document = FALSE, manual = FALSE, cran = FALSE)` completed
   in 3 minutes 22 seconds with status OK and 0/0/0.
2. The complete `testthat` suite passed within the package check.
3. The source-loaded post-processing test file passed 19 expectations after
   adding the spectral regression.
4. The clean-library installed-package PH/provenance smoke test passed with
   `ph_status=completed h0=3 h1=1 provenance_files=2`.
5. `git diff --check` passed before the code commit.

## Preserved and excluded artifacts

- `docs/Dissertation_SubmissionReady_October.pdf` remains Git-ignored and
  untracked.
- `docs/Jonah-BioRxiv_v2.pdf` remains Git-ignored and untracked.
- `docs/private/` remains Git-ignored.
- ad hoc `provenance/` runtime output remains Git-ignored.
- `example_run.r` remains deliberately untracked and unchanged, with SHA-256
  `D205E5EAAFEA27D4AED51146DB1C20CDEC62E56506CDEA239BF4727EDB736185`.

## Declared residual environment debt

`renv::status()` is not synchronized even though package checking and clean
installation succeed. The lockfile retains 29 installed/recorded packages that
are no longer detected as used, primarily the removed `tidyverse` dependency
closure. Six installed versions differ from the lockfile:

- `curl` 7.0.0 recorded versus 6.2.2 installed;
- `fs` 1.6.6 versus 1.6.5;
- GitHub `immunogenomics/harmony` versus CRAN 1.2.3;
- `miniUI` 0.1.2 versus 0.1.1.1;
- `Rcpp` 1.0.14 versus 1.1.1;
- `rprojroot` 2.1.1 versus 2.0.4.

The lockfile was intentionally left unchanged because choosing whether the
lockfile or the current library is authoritative is a separate dependency
policy decision. A blind snapshot here would mix a large environment rewrite
into a code-health repair and could silently change the reproducibility target.

## Recommended next action

Run a dedicated dependency-baseline sprint: choose the authoritative R and
Bioconductor versions, restore the current lockfile in a clean library, execute
the package check and installed PH fixture there, then update the lockfile only
from that verified environment. This is the narrowest next step before runtime
profiling or a Rust feasibility prototype.
