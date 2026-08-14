# Persistence-landscape reference-engine audit — 2026-08-05

## Outcome

The sprint produced and validated a non-default implementation of `full_l2_error_controlled_v1`. It confirms that the dissertation-aligned all-level H0/H1 definition is implementable as a correctness reference without dense landscape serialization. It does not activate the corrected scientific method.

## Reproducible command

From the canonical repository in Ubuntu WSL:

```text
Rscript --vanilla scripts/benchmark_landscape_reference.R \
  <historical-dir>/PD_list_after_retries_unintegrated_sctWhole.rds \
  docs/audits/landscape-reference-2026-08-05
```

The historical RDS was read only. It remains outside the public repository. The benchmark manifest records its filename and MD5 digest but does not copy or serialize any source diagram.

## Correctness corpus

The automated corpus covers empty diagrams, exclusion of infinite/invalid intervals, H0/H1 separation, the single-tent energy identity, translation and scaling laws, a width-0.002 narrow feature, overlapping and deep landscapes, symmetry, deterministic repetition, exact/adaptive agreement, the exact-size guard, non-activation provenance, and streamed all-level values against `TDA::landscape`.

Final validation: 26 focused reference-engine expectations and 219 complete-suite expectations passed with zero failures, warnings, or skips. The built source tarball passed `R CMD check --no-manual` with `Status: OK` under Ubuntu WSL R 4.4.1.

## Benchmark corpus and results

The benchmark used five analytical pairs, one complete pair chosen deterministically from the smallest eligible existing diagrams, and a 25-persistence-quantile-interval-per-dimension workload from each of the two largest diagrams in the selected artifact. The latter is a computational workload sample, not a biological subsample and not evidence for biological conclusions.

| Workload | Exact median s | Adaptive median s | Absolute combined-distance difference |
|---|---:|---:|---:|
| Single tent | 0.0030 | 0.0050 | 0 |
| Translated tents | 0.0030 | 0.0215 | 0 |
| Overlapping H0/H1 | 0.0055 | 0.0475 | 0 |
| Narrow feature | 0.0020 | 0.0055 | `1.12e-18` |
| Deep 25-level landscape | 0.0080 | 0.1180 | 0 |
| Complete existing small pair | 0.0030 | 0.0165 | `2.78e-17` |
| Existing large-diagram quantile workload | 0.0385 | 0.3975 | 0 |

Each case used one unrecorded warmup followed by two timed repetitions; every repetition delta was zero. Returned reference objects occupied approximately 9.0–9.3 KB. These object sizes show that output serialization is compact; they are not peak-process-memory measurements. Adaptive squared-distance error estimates ranged from zero to approximately `4.80e-13`, below the requested `1e-8` tolerance in this corpus.

Raw evidence:

- `landscape-reference-2026-08-05/landscape-reference-benchmark.csv`
- `landscape-reference-2026-08-05/landscape-reference-benchmark-manifest.csv`
- `landscape-reference-2026-08-05/landscape-reference-session-info.txt`

## Interpretation boundary

The benchmark validates numerical behavior on tractable and representative workloads only. It does not establish production runtime for H0 diagrams with thousands to tens of thousands of active levels, validate the biological eligibility of existing diagrams, compare cluster outcomes, or support manuscript claims.

## Open gates

1. Obtain approval before installing Persim/GUDHI in an isolated environment, then cross-check exact distances rather than only pointwise landscape values.
2. Determine whether an established exact critical-pair engine can be integrated without quadratic memory at project scale.
3. Establish a certified or conservatively bounded large-diagram error policy; the current adaptive path records estimates and refinement deltas.
4. Benchmark eligible full diagrams after the observation orientation and filtration policy audits.
5. Keep Rust gated until an algorithmically optimized established implementation is still an end-to-end bottleneck.
