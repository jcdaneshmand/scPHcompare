# Landscape Level/Grid Convergence and Feasibility Audit — 2026-08-05

## Outcome

Persistence landscapes remain central to scPHcompare, but neither a universal level cap nor a universal 250/500/1,000-point uniform grid is scientifically adequate. The proposed corrected primary contract is now **all active levels with exact or error-controlled integration**, with `k=1` retained only as a paper-compatibility sensitivity. The corrected contract remains unactivated pending implementation and the separate diagram-eligibility audits.

## Scope and safeguards

The audit used the eight existing dataset-by-representation diagram strata separately: large heterogeneous and bone-marrow collections, each under Integrated, Raw, SCT Individual, and SCT Whole preprocessing. Interval-depth inventory covered all 124 or 25 diagrams. Distance and clustering sensitivity used all 25 bone-marrow samples and a deterministic 16-sample large-data subset containing two metadata-matched samples from each of eight tissues.

These are methodological sensitivity results, not biological validation. Existing persistence diagrams remain conditional on the observation-unit and filtration audits. No public scientific default was changed, no historical artifact was overwritten, and no manuscript result was regenerated.

## Methods

- Levels evaluated initially: `1, 5, 10, 25, 100, 250, 500, 1,000`.
- Uniform dimension-specific grids: 250, 500, and 1,000 points over each stratum’s finite support.
- Exact-full extension: every active level on the 250-point grid, including H0 depths above 30,000.
- Distance agreement: Spearman correlation and relative Frobenius error against the reference distance matrix.
- Clustering agreement: ARI at `k=2,5,8` and cophenetic correlation for average and Ward.D2 linkage.
- H0/H1 contribution: pairwise H1 energy fraction `d1^2/(d0^2+d1^2)`.
- Grid fidelity: trapezoidal all-level energy compared with the exact identity `sum_i (death_i-birth_i)^3/12`.
- Convergence thresholds: Spearman at least 0.999, relative distance error at most 0.01, minimum clustering ARI at least 0.95, and omitted energy at most `1e-4` where the grid resolved that energy.

The accelerated diagnostic evaluator matched TDA exactly on analytical H0/H1 fixtures. It was used only for profiling and did not activate a new package path.

## Finding 1 — exact level depth is representation-dependent

| Stratum | Exact H0 depth | Exact H1 depth | Conditional cap on 250-point grid |
|---|---:|---:|---:|
| Large Integrated | 2,999 | 4 | 2,999 |
| Large Raw | 34,839 | 4,327 | 250 |
| Large SCT Individual | 32,935 | 725 | 32,935 |
| Large SCT Whole | 5,625 | 160 | 5,625 |
| Bone Integrated | 5,999 | 17 | 5,999 |
| Bone Raw | 16,789 | 1,428 | 100 |
| Bone SCT Individual | 16,789 | 1,428 | 100 |
| Bone SCT Whole | 11,944 | 105 | 11,944 |

The conditional caps in the last column cannot be promoted to primary settings because the 250-point grid under-resolves narrow features. Their variation nevertheless proves that one fixed `K` cannot represent every preprocessing condition honestly.

At `k=1`, correlation with the coarse-grid exact-full distance ranged from 0.546 to 0.989 and relative distance error from 0.081 to 0.974. At `k=25`, several integrated/SCT strata still had correlations of only 0.774–0.919 and relative errors of 0.705–0.874. First-level and five-level historical definitions therefore are not neutral approximations to full landscapes.

## Finding 2 — rank convergence concealed grid under-resolution

At the 1,000-level cap, 250-point grids produced combined-distance Spearman correlations of 0.999972–1.0 against 1,000-point grids, relative errors below `1.8e-4`, and identical tested cluster partitions. On rank and clustering criteria alone, 250 points appeared sufficient.

The analytic energy check contradicted that conclusion. At 1,000 points, only three of eight strata passed a one-percent worst-sample energy-error threshold in both dimensions. Examples of maximum relative error include:

- Large Raw: 0.215 in H1.
- Large SCT Individual: 0.023 in H0 and 1.0 in H1.
- Bone Raw: 0.062 in H1.
- Bone SCT Individual: 1.0 in H0 and 0.465 in H1.
- Large Integrated: 0.0053 in H1 despite excellent rank agreement.

Some narrow landscapes were therefore missed completely even on the dimension-specific 1,000-point grid. Agreement between two uniform grids is not sufficient evidence that either resolves the underlying functions.

## Finding 3 — H1 is retained but empirically small in the legacy diagrams

Under the coarse-grid exact-full calculation, median pairwise H1 energy fractions ranged from approximately `1.55e-8` to `1.33e-4`. H1 may still contain discriminating pairs or biological information, but the unweighted combined distance is numerically dominated by H0 in these diagrams. Corrected analyses must report H0 and H1 separately before any secondary combination.

## Feasibility and performance

- Original TDA benchmark on one large Raw H0 diagram: approximately 23–25 seconds at 500 points for `K=25`, `100`, or `250`; runtime was dominated by interval/grid evaluation rather than requested output depth.
- Accelerated exact-full H0 evaluation for all diagnostic samples in one stratum took approximately 0.14–1.77 seconds before distance accumulation.
- The complete first pass took about 297 seconds; the uniform exact-full pass took about 234–252 seconds.
- Largest recorded landscape profile object was approximately 1.11 GB; read-only process monitoring observed about 1.18 GiB RSS during an exact-full pass.

Dense landscape serialization is therefore wasteful, but streamed/chunked full-level distance accumulation is feasible. This is a credible isolated optimization target. It does not yet justify choosing Rust over an established exact landscape-distance implementation; implementation benchmarking should decide that separately.

## Decision recommendation

1. **Primary levels:** every active consecutive landscape level; no universal cap.
2. **Primary integration:** exact or adaptively refined shared evaluation with an explicit error target checked against analytical energy identities.
3. **Compatibility sensitivity:** `k=1` on a declared common domain, clearly labeled as the paper definition.
4. **Display grids:** regular grids are allowed for figures, but their point count does not define the primary distance.
5. **Comparison strata:** keep each dataset-by-representation analysis separate; do not pool filtration ranges implicitly.
6. **Dimensions:** report H0 and H1 distances separately; combined unweighted L2 is secondary with H1 energy fractions.
7. **Activation:** do not activate `full_l2_error_controlled_v1` until the exact/adaptive engine passes analytical and tractable independent-reference tests, records error/provenance, and diagram eligibility is resolved.

## Auditable outputs

All machine-readable outputs are under `docs/audits/landscape-convergence-2026-08-05/`. The consolidated decision file is `landscape-feasibility-decision-summary.csv`; input hashes, deterministic sample selection, interval depths, energy capture, grid fidelity, distance/clustering stability, H1 contribution, runtime, and memory are preserved in the accompanying CSVs. Reproduction scripts are under `scripts/`.

## Next technical gate

Implement and benchmark an exact or error-controlled streaming landscape-distance engine. Evaluate established implementations first; a Rust component becomes reasonable only if this isolated, now-measured kernel remains a bottleneck and can be validated bit-for-bit or within a declared tolerance against the reference implementation.
