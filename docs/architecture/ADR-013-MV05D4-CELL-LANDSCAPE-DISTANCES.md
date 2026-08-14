# ADR-013: MV5-D4 exact cell landscape distances

## Status

Accepted 2026-08-07. The complete label-closed SCT cell landscape-distance
stage is finished; outcome evaluation remains a separate gate.

## Context

MV5-D3 produced all corrected cell H0/H1 diagrams, but the primary retrieval
analysis requires sample-to-sample distances under the dissertation-aligned
all-level landscape definition. Historical grid-based and first-level results
are not eligible, and full fold matrices would add unsupported work beyond the
frozen query-to-training endpoint.

## Decision

1. Calculate exactly 35,350 eligible query-to-training biological pairs as
   70,700 separate H0/H1 rows in 360 immutable chunks.
2. Exclude all 6,750 essential H0 intervals before landscape construction and
   admit only finite positive-persistence intervals.
3. Use every active level and exact critical-pair segment integration; retain
   the R exact breakpoint engine as an independent eligible-diagram oracle.
4. Keep H0 and H1 primary and separate; publish combined distance and H1
   squared contribution only as secondary components.
5. Bind every request to source record, diagram, serialized-file, manifest, and
   implementation identities; reject partial or stale resume states.
6. Admit two-worker production only after representative and maximum-size
   groups pass time, RSS, storage, exact-oracle, and scope gates.
7. Reject the first complete output set because it embedded nondeterministic
   timing in scientific CSVs, despite exact scientific equivalence; accept only
   timing-separated v3 outputs with a byte-identical full-group repeat.
8. Stop with outcomes closed and all clustering, integration, and gene-view
   counters at zero.

## Evidence

- 4,682,085 finite intervals staged and exactly 6,750 essential H0 intervals
  excluded;
- representative 850-row and maximum 3,250-row admission groups passed;
- four complete eligible H0/H1 R-oracle checks agree within `1.42e-14`;
- 75/75 groups, 360/360 chunks, and 70,700/70,700 dimension rows independently
  validate;
- 35,350 H0/H1 component pairs assembled with secondary combination and H1
  contribution;
- all four files in a fresh 850-row group repeat are byte-identical;
- all 70,700 v2/v3 scientific records are identical, while accepted v3 removes
  per-pair timing from scientific output;
- 1.165 worker-hours, 349.9 MiB peak process-tree RSS, and 65.1 MB private
  result/status storage; and
- total measured SCT cell-primary precomputation is 7.150 worker-hours, leaving
  14.450 hours below the 21.6-hour planning cap; and
- the complete current-source suite passed and the isolated R 4.4.1 package
  check reported `Status: OK`.

## Consequences

The exact cell-topology distance layer is complete and no longer projected.
The next eligible sprint is label-closed assembly of training-fitted component
scales, immutable query-to-training retrieval bundles, and matched baseline
distances. Labels and outcomes remain closed until those prediction artifacts
and their failure policy are immutable.
