# MV7-H full-landscape closure

Date: 2026-08-16

Status: complete; MV7-I descriptive prefreeze authorized

## Completed scientific artifact

The full 124-sample, five-seed, dual-view landscape corpus is complete:

- 20/20 seed × cell/gene × H0/H1 groups;
- 7,626 unordered sample pairs per group;
- 152,520 component-distance rows total; and
- 38,130 rows for each of cell-H0, cell-H1, gene-H0, and gene-H1.

Every row uses finite positive intervals, excludes essential H0, retains all
active consecutive landscape levels, and reports exact streamed squared-L2
distance without a fixed grid or level cap. H0 and H1 remain separate; no
dimension combination was computed in MV7-H.

## Resource outcome

The primary stress group used 663.772 seconds. The remaining 19 groups used
2,988.372 charged seconds, for 3,652.144 total landscape seconds against the
23,232.367-second cap. Peak process-tree RSS was 204,963,840 bytes against the
12-GiB cap. The 20 production group artifacts occupy 90,201,469 bytes.

Gene-H1 dominated at approximately 654–672 seconds per remaining seed. Gene-H0
used approximately 17 seconds per group; cell-H0 used approximately 15 seconds;
and cell-H1 used approximately 33–35 seconds. The accepted Rust kernel therefore
made the full all-level calculation comfortably feasible without changing the
scientific definition.

## Independent validation

The complete-corpus validator passes all 13/13 categories. It independently
verified every group hash, every canonical 124-sample pair axis, all 152,520
row contracts, component balance, Rust identity, resource limits, and the
closed label/downstream firewall.

All four view/dimension components have byte-identical repeats. Four frozen
maximum-depth R oracles pass:

| Component | Intervals | Absolute Rust/R error | Acceptance tolerance | Method |
|---|---:|---:|---:|---|
| cell-H0 | 383 + 383 | 1.273293e-11 | 8.338114e-07 | exact |
| cell-H1 | 533 + 492 | 1.053326e-10 | 1.062933e-09 | adaptive v4 |
| gene-H0 | 499 + 499 | 7.105427e-15 | 3.049324e-10 | exact |
| gene-H1 | 3,201 + 2,664 | 5.024166e-14 | 1.907064e-11 | adaptive v4 |

The cell-H1 reference used 6/27 coarse/fine roundoff-recovery splits while
preserving the original global error budget. The heavy fallback-seed gene-H1
oracle required no splits and agreed especially closely. The analytical H1
fixture passed. Whole-landscape immutable resume passed across 61 production
and checkpoint files with identical hashes, sizes, and mtimes.

## Next authorization

Only the MV7-I descriptive prefreeze is authorized. It may freeze H0/H1
component matrices, a clearly secondary unweighted H0/H1 composite, seed
aggregation, uncertainty summaries, label-free PAM and average-linkage
procedures, and the label-access firewall. It must not yet run clustering,
combine dimensions, access biological labels, compute outcomes, select a
favorable view/component/seed, or make result-dependent claims.
