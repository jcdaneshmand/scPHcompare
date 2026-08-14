# ADR-014: Freeze immutable label-closed SCT cell retrieval inputs

- Status: Accepted
- Date: 2026-08-08
- Decision scope: MV5-D5 only

## Context

MV5-D4 completed exact H0 and H1 landscape distances for every frozen
held-out-query-to-training-reference pair. The statistical benchmark still
required matched sample-level baselines and immutable neighbor rankings before
labels could be opened. The existing MV5-C prototype demonstrated the formulas
on one tissue but did not provide full-cohort production identities.

The earlier plan also allowed a training-fitted H0/H1 component scale. The
actual MV5-D4 scope contains no within-training topology distances, so such a
scale cannot be fitted without either leaking held-out queries or computing a
large new topology scope.

## Decision

1. Keep raw H0 and raw H1 as separate confirmatory retrieval methods.
2. Retain `sqrt(H0^2 + H1^2)` only as an explicitly descriptive secondary
   method.
3. Do not estimate a topology component scale from held-out query distances.
   Do not launch the 262,675-pair within-training topology expansion in MV5-D5.
4. Build the matched energy-distance baseline from the exact same 384 cells and
   training-fitted 30-PC coordinates.
5. Build the pseudobulk context baseline from the same training-derived,
   training-standardized 500-gene panel and the accepted held-out missing-
   feature mapping.
6. Rank by distance and then canonical sample ID, recording exact ties.
7. Bind every result to fold, mean-profile, MV5-D4, and implementation hashes;
   publish atomically and refuse stale cache reuse.
8. Keep labels and all outcome, clustering, integration, gene-view, and fusion
   work closed.

## Evidence

All 75 groups completed 35,350 biological pairs across five methods, yielding
176,750 unique ranked rows and 375 successful method-group records. Independent
validation reproduced every topological distance, checked 450 baseline pairs,
and found a maximum difference of `1.136868e-13`. Admission, public assembly,
and 150-file resume tests are byte-stable. Production used 1,461.45 aggregate
worker-seconds, 35.37 seconds maximum per group, 280,870,912 bytes peak
process-tree RSS, and 161,447,925 private result bytes.

## Consequences

The prediction-side artifacts required for the SCT cell retrieval endpoint are
immutable and can be evaluated without recomputation. The matched sample-level
baseline requirement is satisfied. H0/H1 interpretation remains transparent;
the project does not disguise H0 dominance through an outcome-influenced or
unsupported composite scale.

MV5-D5 does not support clustering, within-training distance contrasts,
integration comparisons, gene topology, or fusion. Those tasks remain gated.
The next stage may open frozen labels for prespecified retrieval evaluation only
after checking the MV5-D5 public hashes.
