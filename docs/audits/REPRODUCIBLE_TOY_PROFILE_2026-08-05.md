# Reproducible toy baseline and preliminary profile

Date: 2026-08-05

Scope: cold-cache dependency bootstrap, deterministic analytical fixture,
stage-level timing, strict package validation, and preliminary Rust gate

Branch: `codex/phase-0-audit-foundation`

Dependency baseline commit: `81ca5a9`

Implementation commit: `d9af01deee5d81c5e6292d5fb5551adef06c58c1`

No push or remote mutation was performed.

## Outcome

The committed 265-record dependency baseline can be reconstructed from an
empty cache under Ubuntu WSL and R 4.4.1. The final verifier found 265 exact
version/source-commit matches and zero mismatches. A maintained public
`run_toy_baseline()` route now runs real persistent homology on two analytical
point clouds, preserves the caller's RNG state, enforces known birth/death
values with a `1e-10` tolerance, writes stable scientific digests, and emits
structured stage timings.

The complete source build and strict package check passed:

```text
Status: OK
0 errors | 0 warnings | 0 notes
```

Three repeated installed-package toy runs produced identical scientific
manifests. The tiny-workload profile does **not** support beginning a Rust
rewrite: the installed ripserr kernel is already native code, each PH stage
took less than 0.45 seconds, and cold namespace/orchestration load dominated
the first invocation.

## Cold-cache bootstrap

Configuration:

- Ubuntu 22.04.4 LTS under WSL;
- R 4.4.1 and Bioconductor 3.19;
- `renv` 1.1.5 bootstrapped by the committed `renv/activate.R`;
- lockfile SHA-256
  `73EA4688D292661BE7847419919E12057B9F7188AF739D6BB938B2E35C6A0150`;
- initially absent cache, library, and renv-state directories beneath
  `/var/tmp/scphcompare-cold-bootstrap-20260805`; and
- no access to the project library or ordinary renv cache.

The standard Bioconductor 3.19 archive no longer serves
`BiocVersion_3.19.1.tar.gz`. The bootstrap script therefore verifies the lock
record, clones the official Bioconductor `RELEASE_3_19` branch, checks out exact
commit `99e637d62c373025b9a757047a144a03c32905cd`, installs version 3.19.1, and
writes a commit marker because `R CMD INSTALL` does not retain that package's
`git_last_commit` field. The verifier requires that marker and checks retained
Git commit metadata for other Bioconductor/GitHub records, including Harmony.

The definitive source-only restore completed in 4,043.6 seconds (67 minutes
23.6 seconds). A follow-up synchronized scan took 5.0 seconds and reported:

```text
lock_records=265
exact_records=265
mismatched_records=0
```

The complete record-level evidence is
[`dependency-bootstrap-manifest-2026-08-05.csv`](dependency-bootstrap-manifest-2026-08-05.csv).
Its SHA-256 is
`C33B19A554A9924F1B22533F779F5B790C0A2AA0BDABC22E3FE8AE23E60C1628`.
This proves reconstruction on the tested Ubuntu/R platform, not cross-platform
portability. Source-only CI should cache the verified renv library; rebuilding
this full optional/test closure on every job would be unnecessarily expensive.

## Analytical correctness contract

`run_toy_baseline(output_dir, seed = 20260805L)` rotates two fixed point clouds
using a seed-selected orthogonal transform. Rotation exercises seed handling
without changing pairwise distances. The function saves and restores the
caller's prior `.Random.seed`.

The enforced reference values are:

| Fixture | H0 | H1 | Expected birth/death values |
|---|---:|---:|---|
| Rotated unit square | 3 | 1 | H0 deaths `(1, 1, 1)`; H1 `(1, sqrt(2))` |
| Rotated four-point line | 3 | 0 | H0 deaths `(1, 1, 1)` |

The accepted absolute numerical tolerance is `1e-10`. The verified run had
maximum absolute reference error `2.22044604925031e-16` and produced:

```text
input_sha256=e29bf7edb5f60ae9cffdca0819668465587aa281fc99aad007f08d1f889ce46a
persistence_sha256=7578999ae0b670219e64d6db152873f6848c25c087255c31eba0494c7bed5a6d
```

Stable outputs:

- `toy_inputs.rds`;
- `persistence_diagrams.rds`;
- `baseline_manifest.csv`;
- `ph_attempt_log.csv`; and
- `stage_timings.csv`.

Wall-clock timestamps and elapsed times are intentionally excluded from the
scientific manifest. Output directories must be empty, preventing stale files
from being mistaken for a new baseline.

The older `generate_toy_data()` helper remains available for sparse-matrix and
metadata-loading demonstrations. Its 100-feature matrices cannot satisfy the
full pipeline's default 500-feature-per-cell QC threshold. The README now says
so explicitly and no longer presents that helper as an end-to-end pipeline
example.

## Verification evidence

1. The focused installed-package toy test passed 28 expectations with no
   failures, warnings, or skips under `NOT_CRAN=true`. It exercised both real
   ripserr child processes, RNG isolation, topology/reference values, output
   round trips, timing schema, and stale-output rejection.
2. A source build and strict `R CMD check --no-manual` with
   `_R_CHECK_FORCE_SUGGESTS_=TRUE` completed in 310.1 seconds wall time with
   `Status: OK`; all examples and tests passed.
3. Three installed-package profile repetitions used seed `20260805`, R 4.4.1,
   ripserr 0.3.0, scPHcompare 0.1.0, and the exact cold-restored library. All
   three scientific manifests were identical.
4. `git diff --check` passed before the implementation commit.

## Preliminary stage profile

The compact tracked summary is
[`toy-baseline-profile-2026-08-05.csv`](toy-baseline-profile-2026-08-05.csv).
Its SHA-256 is
`34F09166D495A5A6523907CF220E46B2AD37FF3586F2C51E2A23890566F55989`.
Times are seconds across three repetitions on the same 40-logical-core WSL
host:

| Stage | Sample | Median | Min | Max |
|---|---|---:|---:|---:|
| End to end | all | 0.748 | 0.696 | 40.093 |
| Namespace/orchestration overhead | all | 0.011 | 0.010 | 39.247 |
| Persistent homology | square | 0.342 | 0.342 | 0.443 |
| Persistent homology | line | 0.393 | 0.343 | 0.394 |
| Generate inputs | all | 0.0002 | 0.0002 | 0.0073 |
| Write artifacts | all | 0.0013 | 0.0013 | 0.0017 |

The 40.093-second first invocation contains a 39.247-second interval before the
instrumented package stages begin. The same interval is about 0.01 seconds on
warm repetitions. By execution boundary, this is attributable to first
namespace/dependency loading plus wrapper orchestration; it is not a sampled
call-stack profile and should not be overinterpreted as a single function.

The actual installed matrix method
`ripserr:::vietoris_rips.matrix()` calls `ripserr:::ripser_cpp()`, which invokes
the compiled `_ripserr_ripser_cpp` symbol through `.Call` in `ripserr.so`.
Therefore the current PH kernel is already native code. On this analytical
fixture, rewriting it in Rust cannot provide material end-to-end benefit.

## Performance and Rust decision

Rust gate status: **not met; no Rust implementation authorized by this
sprint**.

Evidence-backed next performance questions are:

1. reduce or defer broad namespace imports and remeasure cold package load;
2. profile representative sample sizes and the complete preprocessing,
   integration, PH, distance, and clustering stages;
3. measure peak memory, dense conversions, subprocess startup, and scaling
   curves; and
4. evaluate batching/caching and mature native PH alternatives before any new
   language boundary.

The toy profile is a correctness and startup baseline, not a representative
scientific workload benchmark. It does not establish how PH scales with cell
or feature count and does not support claims about large-sample bottlenecks.

## Preserved and excluded artifacts

- `docs/Dissertation_SubmissionReady_October.pdf` and
  `docs/Jonah-BioRxiv_v2.pdf` remain Git-ignored and untracked.
- `docs/private/` and ad hoc `provenance/` output remain Git-ignored.
- `example_run.r` remains deliberately untracked and unchanged, with SHA-256
  `D205E5EAAFEA27D4AED51146DB1C20CDEC62E56506CDEA239BF4727EDB736185`.
- No real dissertation dataset was modified or used.
- Disposable cold libraries, check directories, and profile outputs are
  removed after the compact manifests and audit are retained.

## Next action

The next reproducibility sprint should create a small redistributable fixture
that passes the actual Seurat/QC/integration route, then add CI using the
verified bootstrap script and a persistent renv cache. Only after that fixture
exists should stage timing and peak-memory instrumentation expand to the full
pipeline.
