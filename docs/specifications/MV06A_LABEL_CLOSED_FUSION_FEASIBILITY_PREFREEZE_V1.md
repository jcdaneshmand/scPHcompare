# MV6-A label-closed fusion feasibility prefreeze v1

| Field | Frozen value |
|---|---|
| Contract | `mv06a_label_closed_fusion_feasibility_v1` |
| Date | 2026-08-14 |
| Accepted base | `d0192d35a4ab52006aa83b0ad3b0ad6a19f066cb` |
| Scientific role | Technical/descriptive feasibility only |
| Labels/outcomes | Prohibited |
| Views | `cell_topology_v1`, `gene_topology_v1` |
| Dimensions | H0 and H1 retained separately |
| Primary fusion | Convex equal-weight fusion of four separately normalized components |
| Advanced fusion | Prohibited |

## Purpose

Test whether the existing corrected cell/gene pilot matrices are aligned,
nondegenerate, technically complementary, and suitable for a later
prospectively scaled fusion program. This sprint cannot establish biological
utility, blocked generalization, a preferred weight, a manuscript claim, or a
package default.

## Frozen inputs

Only these eight immutable MV-04 bundles are allowed:

| File | Bytes | SHA-256 |
|---|---:|---|
| `bone__integrated__cell_topology_v1.rds` | 1,374 | `b2988fbe2b69d350f1c219c5bf575b4883acaea563f097df34d30ac28bd31789` |
| `bone__integrated__gene_topology_v1.rds` | 1,381 | `cc79c500b1d2d2d5303913aae9ed68a2c3ccf1ac063d75ad2c4901e9ed2df7c8` |
| `bone__sct_whole__cell_topology_v1.rds` | 1,375 | `82991bd488d7dc33a383ede6de755e7d7a30e3b2f0b87ac977c2eb0961ca189a` |
| `bone__sct_whole__gene_topology_v1.rds` | 1,374 | `9616a928b9fbf6497c1d74a7e481ec9fc02064391f5e0e378134ca17af4588a2` |
| `large__sct_whole__cell_topology_v1.rds` | 3,798 | `0a4d57baef8c6fb7498d44ca00efa97ded9173522cffafc2aa238147f372db8e` |
| `large__sct_whole__gene_topology_v1.rds` | 3,794 | `fc3e4675e9b09287626e5174fd8f0c5f02eda078fe5167c364a0519bc2af0ac5` |
| `large__seurat_integration__cell_topology_v1.rds` | 3,795 | `7346d8b9571d5c997d5f673a3b38c47a0f9f0b64bf59a8493bca961246a37e30` |
| `large__seurat_integration__gene_topology_v1.rds` | 3,810 | `846a5efe764f4f444711af3e4bdaa7000a7530cfa2a82b0770c0bb2e5ae7bc1d` |

The four required paired strata are `bone__integrated`,
`bone__sct_whole`, `large__sct_whole`, and
`large__seurat_integration`. Cell and gene bundles must have identical sorted
sample axes within a stratum. H0 and H1 matrices must be complete, finite,
nonnegative, symmetric, named, and zero-diagonal.

## Frozen normalization and fusion

For component matrix `D_j` and declared fit sample set `T`, fit

`s_j = median{D_j[a,b] : a < b; a,b in T}`.

The scale must be finite and greater than `sqrt(.Machine$double.eps)`. Apply
`Z_j = D_j / s_j` without centering, clipping, transformation, or label access.
The pilot fit scope is all samples in the stratum and must be labeled
`mv06a_descriptive_all_pilot_samples`. This scope is descriptive only. A later
blocked evaluation must refit every scale on each training partition.

Within-view composites are:

- `C = 0.5 * Z_cell,H0 + 0.5 * Z_cell,H1`;
- `G = 0.5 * Z_gene,H0 + 0.5 * Z_gene,H1`.

For frozen gene weight `w`, fusion is:

`F_w = (1 - w) * C + w * G`.

The complete weight grid is `w in {0, 0.25, 0.5, 0.75, 1}`. The technical
primary is `w = 0.5`, which equals a 0.25 convex weight on each of the four
normalized H0/H1 components. No result may select, drop, or relabel a weight.
The arithmetic convex form is frozen because it preserves component
identifiability and does not reintroduce raw H0 scale dominance.

## Frozen diagnostics

For every stratum:

1. record all four fitted scales and source hashes;
2. retain every unordered pair's four normalized components, cell/gene
   composites, all five fusion distances, and equal-weight contribution;
3. calculate Pearson and Spearman correlations for cell versus gene H0, H1,
   and within-view composites over all unordered pairs;
4. calculate deterministic nearest-neighbor overlap separately for H0, H1,
   and within-view composites, using `k = min(3, n - 1)`, ascending distance,
   then ascending canonical sample ID for exact ties;
5. report mean overlap count, mean Jaccard index, exact-neighbor-set fraction,
   and per-sample values;
6. report complete weight-grid matrix hashes and distance ranges;
7. verify a clean repeat produces identical stable outputs.

Correlations and overlaps are descriptive geometry diagnostics. Tissue,
approach, study, manuscript endpoints, clustering outcomes, and labels may not
be read or joined.

## Independent validation

The validator must not call the production fusion constructor. It must:

- independently rehash every source;
- reconstruct all four median scales from raw H0/H1 matrices;
- independently reconstruct every component, composite, and weight-grid row;
- verify sample axes, pair completeness, correlations, neighbor ordering, and
  deterministic identities;
- reject missing, excess, duplicated, nonfinite, or label-bearing output;
- confirm all source files are unchanged.

## Gate and stop rules

MV6-A passes technical feasibility only if all source hashes and axes match,
all scales are nondegenerate, every frozen diagnostic is complete, independent
reconstruction agrees within `1e-12`, and stable outputs repeat exactly.

Passing authorizes only MV6-B prospective scale-up/resource inventory. It does
not authorize a full gene calculation, blocked outcomes, clustering claims,
Similarity Network Fusion, learned weights/kernels, new data, package defaults,
Rust adoption, manuscript claims, tag, release, DOI, PDFs, reviewer material,
or tracking `example_run.r`.

Abort without substitution if an input hash/axis differs, a component is
degenerate, any result is nonfinite, a label-bearing source is required, an
output omits a frozen weight/component, or independent reconstruction fails.
