# Realistic fixture, CI, and observability audit — 2026-08-05

## Decision

This sprint establishes a redistributable production-route smoke fixture and a
local CI definition without changing the normal scientific defaults, pushing a
branch, publishing an artifact, or adding new biological data. The fixture is
appropriate for software and reproducibility validation. It is not evidence
for a biological claim and does not replace a representative benchmark cohort.

## Scope and protected state

- Repository: `scPHcompare`, branch `codex/phase-0-audit-foundation`.
- Execution environment: Ubuntu 22.04 under WSL, R 4.4.1, the locked `renv`
  library, Seurat 5.3.0, Harmony 1.2.4, and ripserr 0.3.0.
- Existing project data remained the only scientific data in scope. No external
  dataset was downloaded or substituted.
- `docs/private/`, all PDFs, and `example_run.r` remain excluded from public
  package/CI inputs. `example_run.r` remains untracked and was not modified.
- The GitHub remote was not pushed or otherwise mutated. The workflow is local
  configuration until the owner chooses to publish the branch.

## Actual production route inspected

The implementation was based on direct inspection and execution of:

- `R/PH_Calculation.R`: sparse loading, Seurat construction, fixed QC,
  SCTransform, merge, representation extraction, Harmony/Seurat integration,
  PH dispatch, reconciliation, and output assembly;
- `R/PH_Functions.R`: subprocess execution, polling, intermediate persistence,
  progress logging, batching, and shared-file behavior;
- `R/ph_utils.R`: stable sample naming and PH batch reconciliation;
- `R/unified_pipeline.R`: wrapper boundaries and post-processing dispatch;
- `R/Integration_flexible.R`: the full Seurat-integration defaults and its
  200-cell/3,000-feature production assumptions.

The older `generate_toy_data()` route was not reused: its 100 features cannot
pass the production minimum of 500 detected genes per cell.

## Fixture contract

`run_realistic_fixture()` generates two deterministic sparse matrices. Each
sample contains 520 genes and 20 cells, with 60 ribosomal genes, five
mitochondrial genes, five hemoglobin genes, and 450 ordinary synthetic genes.
Every gene is detected in more than three cells and every cell has at least 500
detected genes. Thus the fixture exercises all current fixed gates:

- `nFeature_RNA >= 500` and `<= 9000`;
- mitochondrial percentage `<= 20`;
- ribosomal percentage `> 5`;
- gene prevalence in more than three cells;
- post-QC sample count `>= MIN_CELLS`, with the fixture explicitly setting the
  smoke threshold to 20.

The run then executes individual and merged SCTransform, Harmony integration,
and H0 persistent homology for the Harmony representation. The new
`ph_representations` selector makes this focused smoke route possible. Its
default is every historical representation, so ordinary scientific runs retain
their prior behavior.

The packaged reference is
`inst/extdata/realistic_fixture_reference.csv`. For seed `20260805` it requires:

- two loaded and two PH-eligible samples;
- loaded features `520;520` and post-QC cells `20;20`;
- one `Harmony Integration` iteration with matrices `520x20;520x20`;
- finite fraction at least `0.999999999999`, with numeric tolerance `1e-8`;
- two completed PH attempts;
- input SHA-256
  `51b2da0ab384c8a6d4e4a0905a6158753db0f399c8124d9a5cf48a9c373cff70`;
- rounded-Harmony SHA-256
  `79f61fe938182b061ab63900eecf377ec87f2d80d8a17689bf10b31f8457f6eb`.

## Production defects exposed and repaired

The realistic run exposed four defects that mocked wrapper tests did not:

1. `orig.ident` was a factor but a merge-prefix `vapply()` required character.
   The value is now explicitly coerced.
2. Seurat reduced underscore-containing cell/sample names to a shared prefix,
   collapsing two fixture batches into one. Seurat construction now explicitly
   preserves the full stable sample ID in `project.name` and `orig.ident`.
3. Harmony 1.2.4 returned cells by features while the pipeline assigned
   feature-by-cell dimnames. The code now detects the two supported
   orientations, transposes when required, and rejects any unexpected shape.
4. Parallel PH workers wrote one shared RDS and progress CSV concurrently. One
   successful result could be lost nondeterministically. Workers now return
   their results without shared writes; the parent commits results and progress
   serially in stable batch order.

## Runtime and memory evidence

`run_unified_pipeline()` now records `metadata_load`, `ph_processing`, optional
`postprocessing`, and `pipeline_total` stages. Each row includes timestamps,
elapsed seconds, status, error text, process RSS before/after, and RSS delta.

Two installed-package fixture repetitions produced identical stable manifests
and hashes:

| Run | PH processing (s) | Total (s) | PH RSS delta (bytes) | Total RSS delta (bytes) |
|---|---:|---:|---:|---:|
| 1 | 65.3698 | 65.8798 | 96,153,600 | 97,079,296 |
| 2 | 65.3648 | 65.8016 | 96,329,728 | 97,300,480 |
| Final packaged-reference validation | 65.4000 | 65.9134 | 95,674,368 | 96,604,160 |

The legacy subprocess poll interval is 60 seconds, so it imposes a visible
runtime floor even though both PH children finish earlier. RSS is sampled for
the parent R process at stage boundaries; it is neither peak memory nor an
aggregate of child-process memory. A sampling profiler and a shorter adaptive
poll interval remain separate performance work.

## CI definition

`.github/workflows/r-package-ci.yml` uses:

- read-only repository permissions and concurrency cancellation;
- Ubuntu 22.04 and R 4.4.1;
- the verified `scripts/restore_locked_dependencies.R` bootstrap;
- a lockfile/activation-script keyed `actions/cache@v4` cache;
- package build and `R CMD check --no-manual`;
- a clean installed-package realistic fixture run;
- upload of only the dependency manifest and synthetic fixture evidence;
- an early failure if Git tracks `docs/private/`, a PDF, or `example_run.r`.

Action choices were checked against the maintained r-lib/actions repository
(https://github.com/r-lib/actions) and GitHub's dependency-cache documentation
(https://docs.github.com/en/actions/concepts/workflows-and-actions/dependency-caching).
The YAML parsed successfully with the locked Ubuntu `yaml` package. The
workflow has not run on GitHub because this sprint intentionally did not push.

## Validation

- Realistic installed-package run 1: passed the complete specified route.
- Realistic installed-package run 2: passed with identical stable contract and
  identical input/Harmony hashes after the shared-write repair.
- Final exact-source installed-package run: passed the packaged reference
  validator with the same contract and hashes.
- Final source test suite: 8 contexts, 134 expectations, 0 failures, 0 warnings,
  0 skips.
- Final `R CMD build .`: passed.
- Final `R CMD check --no-manual scPHcompare_0.1.0.tar.gz`: `Status: OK`,
  including installation, syntax, namespace, documentation, examples, and
  tests.
- An earlier `--as-cran` diagnostic completed with policy warnings/notes,
  including the non-mainstream `SeuratDisk` dependency and the large import
  surface. The actionable documentation and package-boundary findings were
  repaired; dependency-surface reduction remains planned optimization work.

## Limitations and next evidence

- Only the supported Harmony integration path is covered end to end. The full
  Seurat integration route retains production-scale assumptions and needs a
  separately designed representative fixture/benchmark.
- The fixture checks H0 only to keep CI bounded. The analytical square/line
  baseline remains the H1 numerical correctness contract.
- Synthetic counts deliberately satisfy fixed QC gates and are unsuitable for
  biological interpretation or manuscript performance claims.
- The rounded Harmony hash is a locked Ubuntu/R/dependency reference. Cross-OS
  or dependency-upgrade acceptance should use the invariant/tolerance fields
  and deliberately approve any new reference hash.
- CI is statically validated locally but remains unobserved on a hosted runner
  until a branch is pushed.
- Rust optimization remains unjustified by the measured PH kernel; the 60-second
  orchestration floor and realistic data scaling should be addressed first.
