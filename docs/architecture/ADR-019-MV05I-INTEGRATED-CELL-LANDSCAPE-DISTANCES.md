# ADR-019: exact all-active-level landscapes for integrated cell topology

Date: 2026-08-09
Status: accepted for MV5-I execution

## Context

MV5-H completed the corrected integrated cells-as-observations H0/H1 diagrams.
The project now needs comparable sample-to-sample topology values without
reintroducing the dissertation-era ambiguity around fixed grids, arbitrary
landscape caps, or combined H0/H1 representations.

## Decision

MV5-I reuses the accepted MV5-D4 mathematical contract on the MV5-H records:
finite positive-persistence intervals only, essential H0 excluded, H0 and H1
computed separately, every active consecutive landscape level, and exact
piecewise-linear L2 integration over critical-pair segments. There is no
universal grid and no universal level cap.

Only held-out-query-to-training-reference pairs are created inside each frozen
LOSO fold-seed group. This gives 35,350 biological pairs and 70,700 dimension
rows. The scope is directed and rectangular, so it intentionally contains no
self pairs or reciprocal rows and does not claim square-matrix symmetry or
diagonal validation.

## Consequences

The resulting distances match the most mature dissertation definition while
making its operational details explicit and testable. H0 and H1 can be studied
separately; a secondary Euclidean component combination is allowed but cannot
replace them. Runtime and storage remain bounded through immutable chunks and
resume validation. Retrieval, biological evaluation, clustering, gene
topology, multiview fusion, and new data require later separately authorized
sprints.
