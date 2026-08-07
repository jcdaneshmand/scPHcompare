# MV-04 topological-distance and production-engine validation

| Field | Value |
|---|---|
| Date | 2026-08-05 |
| Scope | Immutable first-seed eligible MV-03 pilot diagrams |
| Primary contract | `full_l2_exact_critical_pairs_v1` |
| Production engine | `persim_0.3.8_corrected_critical_pairs_batch_v1` |
| Views | `cell_topology_v1`; `gene_topology_v1` |
| Homology | H0 and H1 retained separately |
| Biological interpretation | Prohibited |
| Gate | G-MV4 passes for technical advancement only |

## Outcome

MV-04 converted 56 unique corrected diagrams in eight frozen strata into eight
complete, immutable H0/H1 sample-distance bundles. The primary computation uses
all active consecutive landscape levels and exact integration over signed
piecewise-linear differences. No universal level cap or uniform grid was
introduced. Essential intervals were already excluded by the diagram-bundle
builder, and H0 and H1 remain separate in every persisted object.

The persistent Python worker built each landscape once per stratum and
dimension, then calculated all unordered sample pairs. It completed 408
dimension-specific distances in 2,833.29 internal seconds (2,865.9 seconds at
the capped wrapper), with peak process RSS 3,186,323,456 bytes. The finite
interval input was 8,898,315 bytes; the primary pair and build tables were
232,875 and 30,253 bytes. The run stayed below the frozen two-hour and 8-GiB
caps.

## Correctness and determinism

- Four analytical critical-pair tests pass, including the sign-crossing case
  that invalidated Persim 0.3.8's built-in norm.
- All eight H0 matrices and all eight H1 matrices are finite, nonnegative,
  complete, symmetric, exactly zero on the diagonal, and bound to the expected
  sorted sample and diagram identifiers.
- Four eligible R exact-breakpoint oracle cases span H0 and H1. Maximum absolute
  production-versus-oracle disagreement is `1.4210854715202e-14`.
- A second batched repetition of one four-sample stratum reproduced all 12 H0/H1
  distances and squared distances bit-for-bit.
- Every bundle contains its contract, method, input diagram IDs, matrix hashes,
  fitted-scale cache keys, and a bundle cache key. The final artifact manifest
  records SHA-256 hashes and byte sizes.

The R implementation remains the correctness oracle. The corrected Persim path
is the measured pilot production engine; Persim's built-in `p_norm()` remains
prohibited.

## Landscape scaling profile

Construction, rather than pair integration or interpreter startup, dominates
the eligible workload. The slowest single landscape build took 195.58 seconds;
the slowest pair difference took 23.48 seconds. The hardest diagram was a large
Seurat-integration gene-view H1 object with 2,802 finite intervals, 1,624 active
levels, and 2,239,722 critical points. Resource measurements by stratum and
dimension are retained in `mv04-landscape-resource-profile-2026-08-05.csv` and
`mv04-landscape-build-profile-2026-08-05.csv`.

This is a geometry-dependent critical-pair/object-growth bottleneck. Replacing
only the R/Python call boundary would not address it.

## H0/H1 combination and normalization

The raw dissertation-aligned Euclidean combination is retained only as a
secondary output. H0 dominates it in every pilot stratum: median H1 squared
contributions range from about 0.006% to 0.73%, and no observed stratum supports
calling the raw combination balanced.

For later cross-view work, each component also receives a separate median
off-diagonal scale fit. The fitted object records its fit-scope ID, sorted fit
sample IDs, matrix hash, scale, and cache key. Application to held-out samples
is separate from fitting, and zero/near-zero components fail rather than being
silently weighted. Current pilot scales use
`descriptive_all_pilot_samples`, contain no biological labels, and do not
authorize fusion. Future blocked evaluation must refit scales within each
training partition.

## Diagram-distance sensitivities

TDA 1.9.4 bottleneck and Wasserstein-p=1 were evaluated under explicit time and
memory caps. Complete matrices or complete matrix groups are retained only when
all required unordered pairs pass validation. Wasserstein-p=1 completed all 48
dimension-pair rows across the four four-sample strata, but its H0 pair work
alone was roughly an order of magnitude slower than bottleneck before reaching
the harder large gene-H1 diagrams. Full Wasserstein matrices are therefore
technically excluded, not treated as missing at random.

The full bottleneck attempt preserves every completed matrix group and records
any group that cannot finish under the 1,800-second/4-GiB cap as technically
excluded. It completed 183 of 408 dimension-pair rows in 11 complete matrix
groups; the first large SCT gene-H1 group did not finish before the cap. Exact
resource evidence, matrix checks, and
landscape-versus-sensitivity rank correlations are in the corresponding
`mv04-*-feasibility`, `mv04-sensitivity-matrix-validation`, and
`mv04-distance-rank-sensitivity` CSV files. Sensitivities do not replace the
primary landscape matrices.

## Rust gate

Rust remains rejected for this stage.

- A stable hotspot is identified, but only one full eligible primary repetition
  exists, so uncertainty and end-to-end improvement are not established.
- The hotspot is algorithmic construction and object growth, not simply a
  language boundary.
- Established or lower-object-overhead mathematically equivalent approaches
  have not been exhausted on the eligible bundle.
- Cross-platform R-package installation and maintenance ownership are not
  resolved.

The narrow typed contract and numerical oracle are now strong enough to test a
future alternative engine, but they are not sufficient to authorize a rewrite.

## Gate G-MV4

- All primary matrices pass mathematical and identifier checks: **pass**.
- Eligible production results agree with the independent R oracle: **pass**.
- H0 and H1 remain separately persisted and reported: **pass**.
- Sensitivities are complete where feasible and otherwise technically excluded
  with retained evidence: **pass**.
- Runtime, peak memory, I/O sizes, build time, and pair time are measured:
  **pass at frozen pilot scale**.
- Rust remains rejected because every gate criterion is not satisfied: **pass**.

| Question | Disposition |
|---|---|
| Scientific contract coherent? | approve for immutable corrected pilot distances |
| Correctness demonstrated? | pass |
| Computation feasible? | yes for primary pilot; bounded sensitivity exclusions apply |
| Biological interpretation permitted? | prohibited |
| Next action | advance to MV-05 statistical-plan and fair single-view benchmark design |

## Boundary

This sprint did not select clusters, use biological labels, combine cell and
gene views, learn fusion weights, add data, regenerate manuscript claims, change
the package's public scientific default, or begin Rust implementation.

## Package validation

The complete source test suite passes 399/399 expectations with no warnings or
skips. `R CMD build` succeeded, and `R CMD check --no-manual` completed with
`Status: OK` (0 errors, 0 warnings, 0 notes) under the locked Ubuntu R 4.4.1
environment. Python syntax compilation and the engine's four analytical
self-tests also pass.
