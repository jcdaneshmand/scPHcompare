# MV5-BA corrected-Persim equivalence and speed benchmark

Date: 2026-08-13
Status: accepted rejection; R retained; Rust prefreeze authorized

## Decision

Corrected Persim is mathematically supported but rejected as the production
replacement. The accepted R exact/error-controlled engine remains canonical.
The prospectively frozen Rust trigger is satisfied on throughput grounds, so
MV5-BB may define a Rust-kernel prefreeze. Rust implementation, additional seed
production, partitions, labels, and outcomes remain closed.

## Bound evidence

MV5-BA staged all 56 raw-hash-verified diagrams and all 408 accepted MV5-AY
dimension references. The adoption gate was evaluated first on a deterministic
worst-depth panel: the three largest interval-sum H1 pairs from each large gene
stratum. This six-pair panel contains 12 separate H0/H1 comparisons.

The project-controlled Persim calculation passed all three analytical fixtures,
including the sign-changing case that invalidates Persim's built-in norm. It
also passed all 12 worst-depth comparisons. Maximum absolute error in squared
distance was `3.92595352928515e-13`, within every accepted exact/adaptive
certificate. This supports corrected critical-pair integration as an
independent oracle; it does not support production adoption.

## Memory architecture trials

The first prototype retained landscapes for all 56 diagrams and exceeded the
frozen 2 GiB/process cap, reaching 3,072,057,344 bytes before termination. A
one-stratum-at-a-time design also breached the cap at 2,343,881,728 bytes in a
large gene stratum. The accepted timing design retained only one sample pair at
a time and explicitly released heap memory between pairs. It completed with a
peak of 1,960,210,432 bytes. Thus pair-bounded execution resolves the retention
failure, but with little memory margin.

## Matched timing

The same six pairs were recomputed from exact canonical RDS diagrams by the R
engine, with cache identities repeated exactly. R pair times were 259.224 to
314.626 seconds; corrected-Persim times were 300.169 to 502.519 seconds. Median
times were 295.515 seconds for R and 370.556 seconds for Persim. The candidate
speedup was `0.797491`: Persim was about 25% slower, and was slower on all six
pairs. Zero pairs met the frozen threefold replacement gate.

The full 408-dimension Persim corpus was not resumed after adoption became
impossible. Two retention designs had already breached memory, and the
pair-bounded design had failed the throughput gate on every prospectively
selected worst-depth pair. Continuing would add compute without changing the
frozen decision. The full accepted R corpus and its 204 matrices remain intact.

## Boundaries

- Persim 0.3.8 built-in `p_norm()` remains prohibited.
- Corrected Persim remains an independent tractable/worst-case oracle only.
- The accepted R engine remains the production reference.
- MV5-BB may freeze a narrow Rust kernel, FFI, equivalence corpus, resource
  limits, and fallback policy.
- MV5-BB may not implement Rust, calculate new seeds, cut partitions, or open
  biological labels or outcomes.

## Public evidence

- `mv05ba-speed-comparison-2026-08-13.csv`
- `mv05ba-equivalence-summary-2026-08-13.csv`
- `mv05ba-architecture-trials-2026-08-13.csv`
- `mv05ba-benchmark-summary-2026-08-13.csv`
- `mv05ba-continuation-decision-2026-08-13.csv`
- `mv05ba-independent-validation-2026-08-13.csv`

Private diagrams, references, critical-pair metrics, and timing traces remain
under ignored `tmp/` storage.
