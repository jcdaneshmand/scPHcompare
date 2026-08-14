# ADR-018: Bind corrected H0/H1 to integrated coordinates before landscapes

## Status

Accepted and executed on 2026-08-09.

## Context

MV5-G completed all fixed-panel inductive integrated coordinates and authorized
integrated cell PH, but no diagrams existed for that representation. Reusing the
SCT MV5-D3 identity would misstate the scientific input, while moving directly
to landscapes would blur diagram correctness with summary-distance behavior.

## Decision

Give each MV5-G coordinate matrix its own scientific-eligible typed integrated
`cell_topology_v1` identity. Preserve its 384 ordered cell rows and 30 integrated
coordinate columns, bind it to the accepted MV5-G file/group/payload/coordinate
hashes, and compute complete Euclidean Vietoris–Rips H0/H1 over field 2.

Retain one explicit essential H0 interval. Require 383 finite H0 deaths to match
the Euclidean MST, require all H1 intervals to be finite and positive, serialize
one immutable record per view, and group execution only for source-loading
efficiency. Admit two groups, then execute all 75 with two workers under the
900-second, 4-GiB, 12-hour, and 1-GB PH guards. Stop before landscapes.

## Evidence

All 75 groups and 6,750 views completed with zero failures. The cache contains
2,592,000 H0 intervals and 1,545,943 H1 intervals. Every coordinate identity,
record identity, file hash, diagram invariant, stored MST check, and scope check
passed independent reconstruction; 75 fresh full-view MST calculations passed.

One independently executed 90-view group is object- and byte-identical to its
production counterpart. A completed-queue rerun rebuilt zero groups and
preserved every view-audit and result-set hash.

Measured PH used 1.098 worker-hours, 274.9 MiB peak process-tree RSS, 54.64
seconds for the slowest group, and 179.9 MB of private storage. The measured
coordinate-plus-PH and reserved downstream total is 12.169 worker-hours and
1.269 GB, within all frozen caps.

## Consequences

A separate integrated landscape-distance sprint is authorized under the
project-owner-approved dissertation-aligned definition: H0/H1 separate, the
essential H0 class excluded at landscape input, all active consecutive levels,
and exact critical-pair L2 integration without a universal grid or level cap.
MV5-H itself ran zero landscape, distance, retrieval, clustering, gene-view,
fusion, new-data, or biological-outcome jobs.
