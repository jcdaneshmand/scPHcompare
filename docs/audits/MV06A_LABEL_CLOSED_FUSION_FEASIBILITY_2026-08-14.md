# MV6-A label-closed fusion feasibility

| Field | Value |
|---|---|
| Date | 2026-08-14 |
| Prefreeze commit | `b17b8f9` |
| Initial engine commit | `42a411c` |
| Identity correction | `cd0006e` |
| Serialization correction | `31f1d42` |
| Contract | `mv06a_label_closed_fusion_feasibility_v1` |
| Inputs | Eight frozen MV-04 distance bundles |
| Labels/outcomes | Zero |
| Gate | Technical feasibility passes; blocked evaluation remains closed |

## Outcome

All eight predeclared source hashes verified. Four matched cell/gene strata
were reconstructed with identical sample axes and separate H0/H1 components:

- two 10-sample large multi-tissue technical panels;
- two 4-sample bone-marrow technical panels;
- 16 nondegenerate median off-diagonal scale fits;
- 102 unordered sample pairs;
- 510 complete five-weight-grid rows;
- 12 cell-versus-gene correlation diagnostics;
- 84 deterministic per-sample neighbor diagnostics;
- 44 normalized/composite/fusion matrix hashes.

Every equal-weight result is exactly the frozen 0.25 convex contribution from
`cell-H0`, `cell-H1`, `gene-H0`, and `gene-H1`. H0 and H1 remain independently
identifiable. No raw landscape was averaged across views, no component was
dropped, and no weight was selected from results.

## Descriptive complementarity evidence

The two 10-sample panels provide the only informative nearest-neighbor
diagnostic in this sprint. Their cell-versus-gene composite geometry is weakly
correlated and has limited local-neighborhood overlap:

| Stratum | Composite Pearson | Composite Spearman | Mean shared neighbors of 3 | Mean Jaccard | Exact set fraction |
|---|---:|---:|---:|---:|---:|
| `large__sct_whole` | 0.00358 | 0.08867 | 1.5 | 0.38 | 0.10 |
| `large__seurat_integration` | -0.02615 | -0.04730 | 1.1 | 0.28 | 0.10 |

The H0/H1 component correlations in those panels also remain small (Pearson
from -0.1420 to 0.0958). This is evidence that the corrected cell and gene
views are not simple numerical duplicates on the pilot. It is not evidence
that their difference is biological, desirable, stable across resampling, or
predictively useful.

For each four-sample bone panel, frozen `k = min(3, n - 1)` equals three and
therefore includes every other sample. Its neighbor overlap is necessarily one
and carries no complementarity information. The six-pair correlations are
reported completely but are too small and unstable for scientific inference.

## Validation and determinism

The independent validator does not call the production fusion constructor. It
rehashes all eight sources and reconstructs:

1. sample axes and all four raw matrices;
2. all 16 median scale fits;
3. every normalized component and within-view composite;
4. all 510 weight-grid rows;
5. all correlations and deterministic neighbor orderings;
6. all 44 serialized matrix hashes;
7. the artifact manifest and prohibited-output counters.

All 12 validation categories pass within `1e-12`. A clean repeat reproduced all
11 output files with identical bytes and SHA-256 hashes. The temporary repeat
directory was removed after comparison; the accepted evidence is retained in
`docs/audits/mv06a-feasibility-evidence/`.

The first run stopped before output because legacy MV-04 stores the view in
`stratum_id` while `method_id` stores the distance engine. Commit `cd0006e`
corrected only that identity interpretation. The next run stopped before
output while flattening an R aggregate matrix; `31f1d42` corrected only stable
serialization and added a regression test. Neither correction changed the
prefrozen formula, weight grid, inputs, or scientific scope.

## Gate disposition

| Requirement | Result |
|---|---|
| Frozen source identities and axes | Pass |
| Four nondegenerate normalized components | Pass |
| Complete unselected weight grid | Pass |
| Deterministic component/contribution visibility | Pass |
| Independent numerical reconstruction | Pass |
| Clean byte-identical repeat | Pass |
| Biological or blocked utility | Not evaluated and not identifiable here |

MV6-A passes technical feasibility and authorizes only MV6-B prospective
matched gene-view scale-up/resource gating. The observed low redundancy makes
that question worth evaluating; it does not authorize the calculation by
itself. Blocked outcomes, learned fusion, clustering claims, new data, defaults,
manuscript claims, release, DOI, binaries, PDFs, reviewer material, and
`example_run.r` remain closed.
