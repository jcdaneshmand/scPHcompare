# MV5-AX larger corrected-matrix and speed gate

Date: 2026-08-12

## Outcome

All 56 accepted existing diagrams verify by result-file, diagram, stored
provenance, and eligibility identity. The frozen scope contains eight strata
and 204 within-stratum pairs: 114 exact-only pairs and 90 pairs with adaptive
H1. Bone strata have four samples and candidate `k=2:3`; large strata have ten
samples and candidate `k=2:9`.

The conservative planner projects 6.02 adaptive worker-hours, or 3.01 hours
wall when the two independent large gene strata run concurrently. Each process
remains internally single-worker and capped at 2 GiB. This scheduling-only
speedup preserves numerical order within each sidecar and all scientific
semantics.

Rust remains a justified future candidate for the adaptive repeated-sort
kernel, but it need not delay this bounded run. Any Rust implementation must be
validated against the current strict `1e-8` reference outputs before use.

MV5-AY is authorized to calculate all eight complete corrected matrices.
Partitions remain closed pending a separate stability design.
