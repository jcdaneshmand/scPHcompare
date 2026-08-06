# Landscape oracle and diagram-eligibility audit — 2026-08-05

## Outcome

The dissertation-aligned landscape definition remains the corrected target, and `landscape_reference_v1` is independently validated as a correctness oracle. No historical persistence-diagram artifact is scientifically eligible for the corrected analysis, however, because the generating code supplied feature-by-cell expression matrices directly to a PH API that treats rows as points. The historical diagrams may be retained only as labeled performance-stress artifacts.

The most promising production prototype uses Persim's exact critical-pair construction with a corrected exact segment integral. It is not activated: Persim 0.3.8's built-in `p_norm()` is inaccurate for tested landscape differences that cross zero, the prototype has not been exercised on newly generated eligible diagrams, and efficient batched R/Python integration has not been designed. Rust remains gated.

## Reproducible scope

The audit used an isolated Ubuntu WSL Python environment, without modifying the locked R environment. The recorded versions are Persim 0.3.8, GUDHI 3.12.0, NumPy 2.2.6, SciPy 1.15.3, scikit-learn 1.7.2, and psutil 7.0.0. The isolated environment was removed after evidence generation.

The public scripts are:

- `scripts/profile_landscape_python_oracles.py`
- `scripts/landscape_reference_worker.R`
- `scripts/audit_landscape_oracle_and_diagrams.R`

Raw outputs are under `docs/audits/landscape-oracle-2026-08-05/`. Historical RDS inputs and source PDFs were read only and were not copied into the repository.

## Independent numerical cross-check

The frozen corpus contains single, translated, overlapping H0/H1, narrow, deep, and sign-changing landscape pairs. Three independent paths were compared: the R exact breakpoint-stream oracle; direct all-level evaluation with partitioned SciPy quadrature; and Persim exact critical pairs, evaluated both with its built-in norm and with the corrected exact linear-segment identity.

R exact and SciPy agreed to a maximum absolute difference of `8.88e-16`. The corrected Persim calculation agreed to floating-point precision throughout the corpus and scaling benchmarks. The new analytical regression uses `[0,2]` versus `[0.25,2.25]`; its squared H0 distance is exactly `7/64`.

Persim's built-in norm did not pass. In the mixed H0/H1 fixture it returned `0.322748612183951` instead of `0.330718913883074` for H0 and `0.152616076042685` instead of `0.140237893119751` for H1. Inspection localized the discrepancy to the sign-crossing branch used to integrate a piecewise-linear difference. This is a locally reproduced candidate defect, not an upstream issue report. Persim documentation describes `PersLandscapeExact` as exact and shows that subtraction can yield negative critical-pair values, which is precisely the case the regression exercises.

GUDHI's Python `Landscape` interface is a fixed-grid vectorizer rather than the documented C++ exact representation. At 1,000 points its relative error was about `3.7e-5` on the scaling corpus; at 250 points it was about `5.7e-4`. It remains useful for display or sensitivity analysis, not the primary exact/error-controlled distance.

## Observation-orientation audit

The dissertation states that each sample should yield a cell-by-cell Euclidean distance matrix (PDF page 54; printed page 40). Current code instead extracts Seurat assay matrices in feature-by-cell orientation in `R/PH_Calculation.R` and passes them without transposition through `process_and_monitor()` to `ripserr::vietoris_rips()` in `R/PH_Functions.R`. The historical scripts in `PH_ClusteringApp` use the same path. The ripserr contract states that a point-cloud row is a point and a column is a dimension.

This establishes a hard orientation conflict: the computed points are genes/features, while the dissertation intended cells. It is not repaired by the small-sample exclusion policy, by landscape recalculation, or by coincidental interval counts. Correcting it requires a new PH input contract and a clean diagram rerun.

All nine audited historical artifacts are therefore classified `ineligible_orientation_conflict` with allowed use `performance_stress_only`. Supporting symptoms include:

- integrated H0 artifacts with 2,999 intervals, consistent with 3,000 integration features acting as points;
- raw and per-sample SCT artifacts whose inferred H0 point lower bound exceeds the recorded cell count for every matched sample;
- the 122-sample retry artifact using positional rather than stable sample identifiers;
- threshold candidate files that exist but are not manifest-bound to specific diagram artifacts.

The detailed CSV records artifact SHA-256 digests, sample identifiers, cell-count matches, H0/H1 counts, inferred point bounds, threshold provenance, eligibility, and permitted use.

## Controlled production-scale profiling

Because no existing diagram passed eligibility, production candidates were benchmarked on deterministic nested-interval workloads sized to the existing post-filter cell distribution, not on mislabeled historical biology. Across 149 existing samples, filtered cell counts ranged from 396 to 11,475, with quartiles 1,177, 2,663, and 4,681. Boundary workloads used 500, 1,000, 2,500, and 5,000 intervals.

| Candidate | 500 intervals | 1,000 | 2,500 | 5,000 | Accuracy decision |
|---|---:|---:|---:|---:|---|
| Corrected Persim critical pairs | 0.026 s | 0.095 s | 0.447 s | 1.639 s | Agreed with SciPy to floating precision |
| SciPy direct quadrature | 0.622 s | 1.580 s | 6.058 s | 19.278 s | Independent correct comparator |
| R exact reference | 4.481 s | 17.038 s | timeout | timeout | Correct oracle, not production-scale |
| R adaptive reference | 8.004 s | 19.609 s | timeout | timeout | Correct in corpus; estimate is not certification |
| GUDHI grid, 1,000 points | 0.045 s | 0.074 s | 0.185 s | 0.428 s | Approximate; fixed-grid primary rejected |

Times are kernel medians for one boundary repetition; the timeout was 30 seconds. The scaling corpus used three repetitions and showed zero within-candidate distance drift. Full process-startup and peak-RSS measurements are in the raw files. The corrected Persim path used approximately 144–155 MB peak process-tree RSS through 5,000 intervals; GUDHI's dense 1,000-point path rose to roughly 335 MB at 5,000 intervals.

These workloads test numerical and computational feasibility only. They do not establish the geometry, depth, or filtration behavior of corrected H0/H1 diagrams, and they do not support biological claims.

## Decisions and next gate

- Keep `landscape_reference_v1` as the non-default R correctness oracle.
- Reject unmodified Persim 0.3.8 `p_norm()` as an oracle or production distance.
- Advance corrected exact critical-pair integration only as a prototype candidate.
- Keep fixed-grid GUDHI landscapes to display/sensitivity roles.
- Do not activate a scientific default, regenerate headline results, or start Rust work.

The next sprint must first define and test a cell-as-observation PH input contract, including the analysis coordinate space, dimensionality reduction, filtration policy, memory bounds, and provenance. It should then generate a small set of corrected eligible H0/H1 diagrams and batch-profile the candidate landscape engine on those diagrams. A persistent worker or equivalent batch interface is required; paying roughly 1.4 seconds of Python startup for every one of 7,626 sample pairs would dominate the 124-sample calculation.

## Validation

- Complete source suite: 220 expectations passed with zero failures, warnings, or skips.
- Built source package: `R CMD check --no-manual` completed with `Status: OK` on Ubuntu WSL R 4.4.1.
- Python profiler compiled successfully under the recorded isolated Python 3.10.12 environment.
- `git diff --check` passed before commit.

## Evidence files

- `persim-r-exact-crosscheck.csv`: R, SciPy, built-in Persim, corrected Persim, and GUDHI corpus results.
- `historical-diagram-eligibility-detail.csv` and `historical-diagram-eligibility-summary.csv`: artifact-level scientific-reuse decisions.
- `candidate-scaling-benchmark.csv` and summary: three-repetition scaling evidence.
- `candidate-production-boundary-benchmark.csv` and summary: 30-second boundary evidence.
- `corrected-cell-scale-summary.csv`: existing post-filter cell-count range used to size controlled workloads.
- `python-environment-manifest.csv` and `r-session-info.txt`: dependency and runtime provenance.

## Primary references

- ripserr point-cloud orientation: <https://cran.r-project.org/web/packages/ripserr/refman/ripserr.html>
- Persim exact-landscape behavior: <https://persim.scikit-tda.org/en/latest/notebooks/Persistence%20landscapes.html>
- Persim exact implementation: <https://persim.scikit-tda.org/en/latest/_modules/persim/landscapes/exact.html>
- GUDHI exact/grid landscape distinctions and scaling warning: <https://gudhi.inria.fr/doc/latest/group___persistence__representations.html>
