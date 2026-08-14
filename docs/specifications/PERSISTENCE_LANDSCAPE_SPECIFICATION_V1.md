# Persistence-Landscape Specification v1

## Document control

| Field | Value |
|---|---|
| Status | Project-owner-approved target definition; not activated |
| Date | 2026-08-05 |
| Owner | Author team |
| Scope | Persistence landscapes for H0/H1, their summaries, distances, and compatibility artifacts |
| Approval gate | P1-05 / G1 |

This specification separates historical compatibility from the corrected method. On 2026-08-05, the project owner approved `full_l2_error_controlled_v1` as the target definition for future corrected work. It does not authorize reuse of legacy landscape-derived results and does not silently change the current public pipeline. Activation requires a validated reference engine, eligible input diagrams, broader author-team confirmation, and a clean rerun.

The target is aligned with the dissertation’s most recent conceptual definition: L2 aggregation across all available landscape levels within each homology dimension, followed by an optional L2 H0/H1 combination. This specification adds the numerical and provenance requirements that were previously missing.

The first non-default reference engine is implemented and independently validated on analytical and controlled scaling workloads (`landscape_reference_v1`). The diagram audit found that every historical diagram was generated with features rather than cells as points, so none is eligible for corrected analysis. Activation now requires freshly generated cell-oriented diagrams, production-engine validation on them, and broader author-team confirmation.

## 1. Mathematical object

For sample `s` and homology dimension `h`, let the finite positive-persistence intervals be

`D[s,h] = {(birth_i, death_i): birth_i < death_i < infinity}`.

The persistence landscape is the ordered family `lambda[s,h,k](t)`, where `k = 1, ..., K[h]`. Essential intervals and invalid or zero-persistence intervals are excluded from landscape evaluation and must be counted in provenance.

The canonical in-memory and serialized orientation is:

- rows: common filtration-grid values `t`;
- columns: consecutive landscape levels `k`;
- shape: `length(grid[h]) x K[h]`;
- metadata: specification ID, grid values, level values, dimension, software versions, input-diagram digest, and truncation provenance.

Serialized corrected artifacts must never rely on matrix-shape inference to recover orientation.

## 2. Filtration-domain evaluation

The primary evaluation domain is shared across every sample participating in one comparison stratum and is derived separately for H0 and H1 from the finite support of that stratum. “Stratum” means one explicitly identified dataset/representation/comparison analysis; pooling strata is not implicit.

Dimension-specific evaluation is required because H0 and H1 can occupy very different filtration ranges. A single grid dominated by H0 can undersample H1. Diagnostic profiling showed that 250/500/1,000 equally spaced points can preserve pairwise-distance ranks while still missing substantial or even complete landscape energy for narrow features. Therefore no fixed uniform point count is approved for primary distances.

Primary distances require either exact integration or a shared adaptive/error-controlled quadrature whose tolerance is tested against the exact all-level energy identity `sum_i (death_i - birth_i)^3 / 12`. The quadrature must record its domain, nodes, error target, achieved error bound, and refinement history. A regular grid may still be used for display after the distance calculation is complete.

Absolute filtration scale is primary. Any sample-normalized or dimension-normalized filtration analysis is a separately named sensitivity analysis, not a replacement applied without disclosure.

Unreconciled per-sample grids remain prohibited for pairwise landscape distances: equal array indices must represent equal filtration values. Adaptive evaluation may use pair-specific refinement only when both functions in that pair are evaluated on the same nodes and the integration error is controlled. A fixed `[0,1]` grid is permitted only when the underlying filtration has explicitly been transformed to `[0,1]` and that transformation is recorded.

## 3. Landscape levels

The proposed corrected primary representation uses every active consecutive level. Diagnostic exact-full depths ranged from 2,999 to 34,839 for H0 and from 4 to 4,327 for H1 in the profiled strata. No universal fixed cap preserved the corrected representation across preprocessing conditions. The implementation should therefore stream or chunk level contributions rather than serialize dense grid-by-level matrices.

If computation requires a cap, it must be named, for example `full_l2_error_controlled_K1000_sensitivity_v1`, and accompanied by an exact omitted-energy assessment plus distance-rank and clustering-stability sensitivity. A capped result must not be called “all levels.”

The first-level-only representation is retained as a prespecified paper-compatibility sensitivity, not the corrected primary analysis.

## 4. Summaries and distances

The standard dimension-specific L2 landscape distance is primary:

`d_h(i,j) = sqrt(sum_k integral (lambda[i,h,k](t) - lambda[j,h,k](t))^2 dt)`.

Integration must be exact or meet a declared, verified error tolerance. Trapezoidal integration is acceptable only on an adaptively refined shared node set whose error is verified. A raw Frobenius/discrete Euclidean norm is a legacy calculation because its magnitude depends on grid density and omits integration weights.

For a one-dimensional display curve within one homology dimension, the all-level summary is:

`A[s,h](t) = sqrt(sum_k lambda[s,h,k](t)^2)`.

The mean across levels is not an L2 summary and is not used in the corrected analysis. `lambda_1(t)` remains a separately labeled sensitivity curve.

H0 and H1 distances are reported separately as primary results. The unweighted combination

`d_combined(i,j) = sqrt(d_0(i,j)^2 + d_1(i,j)^2)`

is secondary and must be accompanied by the dimension components and, when summarized across pairs, the H1 energy fraction `d_1^2 / (d_0^2 + d_1^2)`. This avoids implying that an algebraic H0/H1 combination makes their empirical contributions balanced.

Pointwise H0/H1 curve aggregation is allowed only on an explicitly shared display grid. With dimension-specific primary grids, H0 and H1 curves remain separate.

## 5. Versioned contracts

| ID | Levels | Grid | Distance | Intended use |
|---|---|---|---|---|
| `legacy_k1_unit_grid_v0` | `k=1` | 100 points on `[0,1]` | Discrete Euclidean | Reproduce current public-package behavior only |
| `paper_k1_common_grid_v1` | `k=1` | Dimension-specific common absolute grid | Trapezoidal L2 | Paper-definition sensitivity |
| `full_l2_common_grid_v1` | Consecutive `1:K[h]` | Fixed-count dimension-specific grid | Trapezoidal L2 | Superseded proposal; convergence audit rejected a universal fixed point count |
| `full_l2_error_controlled_v1` | Every active consecutive level | Shared adaptive nodes or exact integration | Exact/error-controlled L2 over levels and filtration | Project-owner-approved corrected primary target; not activated |

Historical five-level artifacts are not assigned a corrected contract ID: they combine inconsistent matrix orientations, undocumented per-diagram grids, and index-wise Frobenius distances.

## 6. Invalidation and migration

The following historical outputs cannot be reused as corrected results:

- serialized persistence-landscape matrices without grid metadata;
- landscape distance matrices and clustering derived from them;
- aggregated landscape curves, tests, null distributions, and figures derived by `rowMeans()` or ambiguous orientation;
- manuscript values or claims whose evidence depends on those outputs.

The separate observation-unit audit has now invalidated all nine inspected historical persistence-diagram artifacts: feature-by-cell assay matrices were passed directly to a row-as-point PH API, contrary to the dissertation's cell-by-cell intent. Consequently bottleneck, spectral, and landscape results derived from those diagrams are also ineligible for corrected analysis. Historical artifacts remain available only for explicitly labeled reproduction or performance-stress work.

Migration requires regenerating landscapes from eligible diagrams, storing the full contract metadata, recomputing dimension-specific and combined distances, rerunning clustering/statistics, and regenerating every affected table and figure. Legacy and corrected artifacts must live under distinct specification IDs and must never share cache keys.

## 7. Required acceptance evidence

- Analytical fixtures distinguish first-level, mean-level, and all-level L2 behavior.
- Translation fixtures prove that per-sample grids are rejected.
- Out-of-range fixtures prove that fixed unit-grid truncation is detected.
- H0/H1 component distances agree with hand calculations before combination.
- Grid-resolution, analytic energy fidelity, and level-cap convergence are quantified on realistic data.
- Exact or adaptive integration agrees with analytical fixtures and an independent implementation on tractable diagrams. This is complete for the R exact oracle, corrected Persim critical-pair integral, and SciPy quadrature; unmodified Persim 0.3.8 fails the sign-changing-difference fixture.
- Every result records achieved integration error and whether all levels or a named cap were used.
- Software, grid, levels, dimensions, diagram digest, and specification ID are embedded in result provenance.
- The author team records whether `full_l2_error_controlled_v1` is approved, revised, or rejected.

## 8. Decisions still required

1. Obtain broader author-team confirmation of the project-owner-approved all-level L2 target and `k=1` sensitivity hierarchy.
2. Approve each dataset-by-representation combination as a separate comparison stratum unless a scientifically justified pooling rule is written.
3. Approve the production exact/error-controlled implementation after it is batch-benchmarked on freshly generated eligible diagrams.
4. Decide whether a scale-normalized sensitivity is scientifically useful in addition to absolute scale.
5. Approve how separate H0/H1 results and the secondary combined result enter the revised manuscript.

## 9. Reference implementation status

`landscape_reference_distance()` is an internal, explicitly non-default interface for `full_l2_error_controlled_v1`. It uses exact streamed breakpoint integration for tractable diagrams and partitioned adaptive quadrature with a tighter-repeat check above a provisional exact-complexity guard. The result records separate H0/H1 distances, the secondary combined distance, H1 squared-distance contribution, interval counts, numerical controls, and versioned provenance.

The provisional 200-interval exact guard is an implementation-safety threshold, not a level cap and not part of the scientific definition. Adaptive integration includes all active levels at every evaluation point. Its reported QUADPACK error and refinement delta are numerical error estimates, not a formal certified bound; primary-method activation therefore remains blocked on a large-diagram error-policy decision.

Architecture and benchmark evidence are recorded in ADR-001, ADR-002, `docs/audits/LANDSCAPE_REFERENCE_ENGINE_2026-08-05.md`, and `docs/audits/LANDSCAPE_ORACLE_AND_DIAGRAM_ELIGIBILITY_2026-08-05.md`.
