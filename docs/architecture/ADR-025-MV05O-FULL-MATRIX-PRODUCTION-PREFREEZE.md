# ADR-025: prefreeze full training-matrix production before execution

- **Date:** 2026-08-10
- **Status:** accepted prospectively; production not executed
- **Scope:** corrected cell-topology training distances and matched baselines

## Context

MV5-N proved that complete training matrices are scientifically specified and
fit the resource envelope, but it intentionally executed only bounded
minimum/representative/maximum profiles. A full run contains 1,838,725 values
and must not begin from mutable scripts, implicit queues, or outcome-aware
choices.

The admission landscape runner is a correctness oracle but is hard-coded to
12 groups of at most 32 rows. Reusing it as a per-chunk production launcher
would also discard group-level landscape caching and invalidate the resource
projection.

## Decision

Freeze a distinct production family before execution:

- regenerate exact training pairs per group and require the accepted MV5-N
  group/chunk roots;
- run 150 deterministic fold-seed-representation groups with at most two
  workers;
- retain 4,340 resumable landscape chunks while constructing landscapes once
  per group;
- run 150 representation-specific energy units and 75 shared pseudobulk units;
- bind four implementation hashes and one 18-artifact source-freeze root;
- require atomic output/status pairs and hash-validated immutable resume;
- enforce 21.6 worker-hours, 900 seconds, 4 GiB, and 10 GiB guards;
- require 12 independent exact oracles, two maximum-group clean repeats, all
  formula checks, and an all-unit resume proof; and
- keep production, clustering, labels, outcomes, and downstream counters zero
  during this prefreeze.

## Consequences

The future matrix run is mechanically authorized only after this prefreeze is
committed. Any source or implementation change invalidates the authorization.
No complete matrix or cluster result is produced by this decision.

The projected private footprint is conservatively 1.278 GB including baseline
and 25% status/validation reserve, not merely the 619 MB landscape estimate.

Descriptive training-partition alignment and inductive held-out generalization
remain distinct later label-open estimands. Neither is part of matrix
production or may influence its configuration.
