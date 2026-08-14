# MV5-BB Rust landscape-kernel prefreeze

Date: 2026-08-13
Status: accepted; bounded prototype specified; toolchain authorization required

MV5-BB converts the MV5-BA throughput trigger into a narrow prospective Rust
contract without changing the accepted R implementation. R remains canonical
and owns input validation, all scientific provenance, H0/H1 separation,
certification, atomic artifacts, public APIs, and fallback. Corrected Persim
remains an independent oracle only.

The proposed kernel accepts two canonical finite birth/death arrays and one
homology dimension through a versioned C ABI. It returns exact piecewise-linear
squared-L2 distance plus level/event diagnostics and a stable status code. It
must use all consecutive active landscape levels, zero-pad unequal depth, avoid
uniform grids and universal caps, allocate one pair at a time, use one internal
thread, and remain deterministic. Rust cannot read project files/metadata,
generate PH, construct views, name caches, combine H0/H1 as primary, cluster,
select `k`, or access labels/outcomes.

The equivalence ladder contains 443 required results: three analytical cases,
20 tractable R exact cases, 12 MV5-BA worst-depth results, all 318 exact MV5-AY
dimensions, and all 90 adaptive-certified MV5-AY H1 results. Only the first
three tiers are needed to evaluate a prototype; all tiers and a later separate
decision are required for production adoption.

The prototype gate requires at least 3x median speedup on the same six-pair
panel, no pair slower than R, peak RSS at or below 1 GiB, two clean deterministic
build/runs, stable panic/error fallback, memory-safety tooling, and exact-staged
package checks both with and without Rust available. Any failure leaves the R
path unchanged.

No Rust toolchain is currently installed in Ubuntu. MV5-BC is therefore blocked
on owner authorization to obtain a pinned toolchain. No download, installation,
Rust source, added seed, partition, label, outcome, default, or legacy change
was made in MV5-BB.
