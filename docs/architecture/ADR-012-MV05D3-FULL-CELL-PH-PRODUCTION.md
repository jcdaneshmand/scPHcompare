# ADR-012: MV5-D3 full cell-PH production

## Status

Accepted 2026-08-07. MV5-D3 is complete and ends at immutable SCT cell H0/H1
diagrams; landscapes and all biological analysis remain a separate gate.

## Context

MV5-D2 established that corrected 384-cell, 30-PC persistence diagrams were
numerically correct, deterministic, and feasible, but it produced only 30
pilot diagrams. The landscape analysis cannot be trusted until every intended
fold-seed view has an identity-bound, independently validated diagram. Running
one process per view would also repeatedly deserialize the same large fold
cache, wasting time without strengthening scientific isolation.

## Decision

1. Build exactly 6,750 diagrams from the 75 validated MV5-D1 caches and reject
   any source or view outside the frozen public manifest.
2. Group execution by fold-seed cache, but serialize and identify each of the
   90 diagrams independently.
3. Retain separate complete H0/H1 diagrams, including one explicit essential
   H0 class, under the validated MV5-D2 contract.
4. Admit full production only after two groups pass the 900-second, 4-GiB,
   storage, identity, and correctness guards.
5. Require atomic writes, immutable resume, per-group checkpoints, and refusal
   to overwrite stale or ambiguous results.
6. Independently validate every record and stored H0 MST result, recompute one
   full MST oracle in each group, and byte-compare a fresh 90-view group repeat.
7. Keep outcomes closed and stop with all landscape, distance, clustering,
   integration, gene-view, and outcome counters at zero.

## Evidence

- 75/75 groups and 6,750/6,750 views completed with zero failures;
- 2,592,000 H0 and 2,096,835 H1 intervals retained;
- all 6,750 files, identities, static diagrams, and stored MST checks passed;
- 75/75 independently recomputed full-view H0 MST checks passed with zero
  recorded maximum error;
- a fresh 90-view group repeat was object- and byte-identical in every record;
- immutable resume preserved all 90 first-group records and its checkpoint;
- 1.047 PH worker-hours, 273.1 MiB peak process-tree RSS, and 196.3 MB of
  private result storage; and
- measured-plus-projected SCT cell-primary work is 9.556 hours, leaving 12.044
  hours below the 21.6-hour planning cap; and
- the complete current-source suite passed and the isolated R 4.4.1 source-
  package check reported `Status: OK`.

## Consequences

The corrected SCT cell-primary analysis can proceed from a complete diagram
cache rather than from estimates or historical feature-oriented artifacts.
Grouped execution reduces repeated source-cache I/O while preserving
per-view provenance and restartability. The next eligible sprint is a
label-closed, bounded landscape-distance stage using the approved all-active-
level exact/error-controlled definition. It must continue to report H0 and H1
separately and may not open outcomes merely because the diagrams exist.
