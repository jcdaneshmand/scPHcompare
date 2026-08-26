# MV17-A source-inventory prefreeze attempt-2 failure

The exact-head `5d741ad289b521778cf48a07ea8ca1d7b0a07e79` MV17-A v2 builder
passed the recovered closure-manifest gate and exactly rehashed all 1,272 bound
MV8 PH artifacts. It then stopped before audit publication because the validator
incorrectly required every artifact in the accepted MV8 binding to have 475 or
500 gene points.

Aggregate-only diagnosis showed 1,264 `gene_topology_v1` artifacts and eight
adopted `cell_topology_v1` sentinel artifacts. The eight cell records have 384
points and are legitimate members of the historical MV8 closure, but MV17's cell
estimand is sourced from the independently closed MV13/MV14 corpus. File bytes,
file hashes, diagram hashes, schemas, label state, and outcome state all passed.

Recovery retains the complete mixed-closure rehash, then inventories only the
1,264 gene PH artifacts (2,528 H0/H1 dimension records) and the corresponding 26
gene-only landscape groups (626 chunks; 152,688 pairs). The eight historical
cell records and two cell-distance groups are excluded from the MV17 gene
estimand rather than misclassified or mutated. Attempt 2 remains preserved as an
empty local root. A distinct v3 root is required after this record and correction
are committed.
