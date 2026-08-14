# MV5-AP realistic landscape compatibility and resource gate

Date: 2026-08-12  
Accepted base: `5f60a17`  
Gate implementation: `10b7eb2`  
Independent validator: `94ea8db`

## Outcome

MV5-AP stops safely with a major numerical-readiness issue. The new public
scientific landscape definition remains methodologically coherent, but its
current exact guard and adaptive implementation are not ready for realistic
workflow integration. No opt-in integration sprint, workflow default change,
or artifact migration is authorized.

The frozen full evaluation was intentionally not run after its first sentinel
pair triggered the predeclared abort rule. Continuing through the remaining
pair matrix would have generated knowingly incomplete or differently tuned
evidence after the gate had already failed.

## Frozen corpus and subset

The accepted MV-04 manifest contains 56 tracked, corrected persistence
diagrams in eight strata: both cell and gene views, bone and large cohorts, and
integrated, SCT, and Seurat-integration representations where available. H0 has
383–499 finite intervals per diagram; H1 spans 79–2,802.

The deterministic, outcome-independent rule freezes three diagrams per stratum:
minimum, middle ordered, and maximum H1 depth after ordering by H1 count,
sample ID, and diagram ID. This yields 24 diagrams. All 24 diagram hashes, file
hashes, stored provenance hashes, corrected result classes, and scientific-
eligibility flags verify exactly. The subset is fully frozen, but only the
first two-diagram sentinel was executed because it triggered the abort rule.

## Two reproduced blockers

The sentinel is the minimum/middle-depth pair from the bone integrated cell
stratum (383 H0 intervals and 130/137 H1 intervals).

1. The public exact default (`exact_max_intervals = 200`) rejects it immediately
   because realistic H0 contains 383 finite intervals. Raising the explicit
   computational guard to 500 succeeds in 5.318 and 5.516 seconds across the
   two clean builds.
2. Adaptive integration at the promised `1e-8` absolute/relative tolerance
   fails in both builds with `extremely bad integrand behaviour`, even with the
   subdivision limit raised from 200 to 1,000. It emits no partial result and no
   uncertified distance, which is the correct failure behavior.

The exact sentinel distances are H0 `54.1735947138941`, H1
`6.23387432164545`, and descriptive combined `54.5310879525003`.

## Diagnostic evidence, not a policy change

An explicit `1e-6` adaptive diagnostic succeeds in 22.019 and 22.451 seconds.
Its largest distance difference from exact is approximately
`1.005e-9`, well within the diagnostic comparison threshold. This supports the
landscape definition and points to numerical error policy/partitioning as the
problem. It does not authorize silently loosening tolerances or changing the
default.

Explicit legacy mode returns zero for H0, H1, and combined on this pair because
the historical fixed `[0,1]` grid does not cover the sentinel's informative
filtration range. This is descriptive compatibility evidence only—not a
method winner or biological result.

## Serialization and repeat

The raised-guard exact versioned result round-trips through RDS with identical
object, distances, and cache key in both builds. Seven scientific/provenance
ledgers and both 16-category independent validation ledgers are byte-identical.
All probe fields except wall time are identical, including statuses, error
messages, distances, error estimates, object sizes, source IDs, and cache keys.
Compressed RDS byte sizes differ slightly because runtime metadata is retained;
each file round-trips exactly within its build.

## Decision and repair boundary

MV5-AP authorizes only a numerical-engine remediation prefreeze (provisionally
MV5-AQ). That sprint should diagnose the H1 integration partition/crossing
behavior, establish a realistic guard policy that distinguishes computational
safety from landscape depth, and prove strict error control on the frozen
sentinel before rerunning MV5-AP. It must not silently loosen tolerance, fall
back to a grid, cap landscape levels, change workflow defaults, overwrite
legacy artifacts, or proceed to clustering/fusion/new data/Rust/claims.

Because this is a major issue, automatic sprint continuation stops here pending
project-owner confirmation.
