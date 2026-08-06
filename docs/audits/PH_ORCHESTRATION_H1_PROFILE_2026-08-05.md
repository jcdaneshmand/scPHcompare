# PH Orchestration and H1 Profile Audit — 2026-08-05

## Decision

H1 persistent homology is retained as a first-class, explicitly requested
pipeline dimension and is now required in the synthetic installed-package CI
route. It is not yet promoted as a biological result: the fixture proves
software behavior, bounded cost, invariants, and repeatability only.

## Goal and protected scope

This sprint inspected the actual `processx`/`Rscript`/`ripserr` path, removed
the known orchestration delay, added sampled process-tree memory evidence, and
measured H0 versus H1 on the existing realistic synthetic fixture. PDFs,
`docs/private/`, and the untracked `example_run.r` were not added to Git. No
remote operation was performed.

## Root cause and implementation

`process_and_monitor()` previously called `Sys.sleep(60)` after starting each
Ripser child. A child that completed immediately was therefore not observed
until the full interval expired. The monitor now uses the process handle's
bounded wait, which returns immediately on child exit, with a default 0.25 s
sampling interval and a separate 60 s progress-message interval. The public
controls are `ph_poll_interval`, `ph_progress_log_interval`, and
`ph_max_time_per_sample`.

Every PH attempt now records the polling configuration and poll-boundary
samples for monitor RSS, child RSS, descendant RSS, aggregate process-tree RSS,
sample count, and maximum observed tree process count. “Peak” in these field
names means the maximum sampled value, not a continuous kernel-level maximum.
Appending to a legacy attempt CSV fills new additive fields with missing values
while retaining required identifier/status schema checks.

No PH input matrix, distance, threshold, Ripser algorithm, or numerical
tolerance changed.

## Benchmark protocol

- Ubuntu WSL, R 4.4.1, locked project environment.
- Seed `20260805`; two synthetic samples, each 520 genes by 20 cells.
- Real sparse-RData loader, fixed QC, Seurat, SCTransform, Harmony, and Ripser.
- Harmony representation only; infinite Ripser threshold (`-1`).
- Two repetitions each for maximum homology dimension 0 and 1.
- The profiler aborts on incomplete samples, input/Harmony hash drift,
  within-dimension PH hash drift, H0 drift under H1, or absent H1 intervals.
- Reproduction command:
  `Rscript scripts/profile_realistic_ph_dimensions.R <empty-dir> 2`.

The exact summary is versioned in
`docs/audits/ph-dimension-profile-2026-08-05.csv`.

| Max dimension | Repeats | Median total (s) | Max total (s) | Median PH attempt (s) | Max sampled PH tree RSS | H0 counts | H1 counts |
|---:|---:|---:|---:|---:|---:|---|---|
| 0 | 2 | 10.4645 | 12.898 | 0.4772220 | 75,657,216 B | 519;519 | not requested |
| 1 | 2 | 6.1645 | 6.186 | 0.6635815 | 83,054,592 B | 519;519 | 922;856 |

The older installed H0 baseline was approximately 65.8 s. The new measured H0
maximum was 12.898 s, an approximately 80% reduction, while its rounded PH hash
remained `64b599...b80b7a8`. H1 added about 0.19 s to the median individual PH
attempt and about 7.4 MB to the maximum sampled child-tree RSS. End-to-end H1
was not slower in this small repeated sample because preprocessing/startup
variance dominated.

Both H1 repetitions produced the same rounded PH hash
`b0a747...44670fc`, the same H0 counts as H0-only runs, and H1 counts of 922 and
856. The locked realistic reference is now dimension-specific for H0 and H1.

## CI policy

The required local CI definition now installs the checked source and runs both
H0 and H1 realistic fixtures. The extra run is justified by sub-second PH work,
bounded sampled tree RSS near 83 MB, deterministic output, and direct coverage
of the dissertation/paper's H1 path. H2 remains outside required CI and needs a
separate scientific and scaling decision.

## Limitations and next action

Final validation used the clean built source in Ubuntu WSL:

- complete source suite: 151 expectations, 0 failures, 0 warnings, 0 skips;
- `R CMD check --no-manual`: `Status: OK`;
- installed-package H0/H1 profile: both packaged dimension references passed;
- GitHub Actions YAML: parsed successfully;
- tracked protected-artifact query: no PDFs, `docs/private/`, or
  `example_run.r` returned.

- This fixture treats expression-matrix rows as points, following current code;
  the scientific orientation and preprocessing definition still require the
  Phase 1 method audit.
- Synthetic H1 intervals have no biological interpretation.
- Poll-boundary RSS can miss a short-lived maximum.
- Only two repetitions and one small fixture were used; this is a regression
  and feasibility profile, not a publication performance claim.
- Hosted CI remains unobserved until the branch is deliberately pushed.

The next highest-value sprint is the scientific landscape-specification audit:
reconcile paper, dissertation, and code definitions for landscape levels,
grids, aggregation, H0/H1 combination, and distances before optimizing or
adding clustering methods.
