# MV17-G full H0/H1 calibration production-design specification

Date: 2026-08-26
Status: prospective design only; no source selection, queue materialization, null generation, PH, or production execution authorized

## Scientific axis and precision

Use the immutable MV17-C exact500/seed-20260805 axis: 132 cell units and 132 gene units. H0 and H1 remain separate, essential H0 remains excluded, and the four closed summaries remain finite interval count, total persistence, maximum persistence, and the exact all-active-level landscape squared L2 norm.

Cell units retain the three admitted nulls: coordinate permutation, covariance Gaussian, and radial-density clouds. Gene units retain those three plus within-row axis shuffling. Use exactly 99 fixed null replicates per compatible unit/family and the plus-one greater-or-equal tail, yielding a minimum attainable probability of 0.01 and worst-case Monte Carlo standard error of 0.05. These are calibration summaries, not labels, biological tests, or multiplicity-adjusted manuscript claims.

## Deterministic queue

Order the 132 private units within each view by canonical identity token. The public seed formula depends only on view, private unit order, family order, and replicate. Seeds start at 174001 and are unique across all 91,476 primary null runs.

The primary workload is 91,476 null runs plus 264 observed runs, grouped into 924 unit/null-family children plus 264 observed children: 1,188 atomic children and 91,740 scientific runs. The repeat workload reuses the MV17-C private minimum/median/maximum burden units for each view with identical seeds: 21 null-family children plus six observed children, 27 atomic children and 2,085 scientific runs. Every repeated metric must match primary exactly.

Each null-family child receives one source matrix and its complete 99-seed contiguous block, returns 99 H0/H1 by four-summary records atomically, and retains the H0 MST oracle for every replicate. Observed children run once. No dimension, summary, family, unit, failure, or unstable result may be silently dropped.

## Measured resource envelope

MV17-C used 816.83 seconds for 195 primary runs and 296.40 seconds for 65 repeat runs. The conservative rate is therefore 4.56 seconds per scientific run. At that rate the fixed MV17-G workload projects to 118.845 serial hours; a 25% margin projects to 148.556 hours. A seven-day aggregate cap, one worker, zero automatic retries, 1,800-second and 8-GiB per-child caps, 12-GiB private cap, and 64-MiB public cap are the design boundaries.

The next implementation sprint must support atomic checkpoint adoption after interruption. Adoption may skip only a completed child whose source, seed block, implementation, result schema, scientific hash, empty streams, and resource receipt all validate exactly. Failed or partial children are preserved; they are not retried without a separately committed recovery gate.

## Firewalls and next gate

This design does not authorize a private 264-unit locator, real queue, null, PH calculation, full calibration, or real localization. A source-safe selector, grouped worker, serial checkpoint runner, independent closure, tests, and exact-head prefreeze must be implemented and committed before execution.

Labels, outcomes, tissues, clustering, view ranking, fusion, biological interpretation, manuscript claims, H2 substitution, cleanup, and deletion remain closed. Real localization is a separate MV17-G lane requiring its own source/output/resource prefreeze after full-calibration execution is safely gated.
