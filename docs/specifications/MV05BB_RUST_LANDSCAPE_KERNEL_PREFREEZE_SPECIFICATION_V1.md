# MV5-BB Rust landscape-kernel prefreeze v1

Date: 2026-08-13
Authorization: MV5-BA `4a0670b`

## Decision and boundary

MV5-BB specifies a narrow, optional Rust prototype for the corrected persistence-
landscape squared-L2 kernel. The accepted R engine remains canonical. Corrected
Persim remains an independent oracle. This sprint does not install a toolchain,
write Rust, alter an R default, calculate additional seeds, select `k`, cut a
partition, or access labels/outcomes.

MV5-BC may implement a bounded prototype only after explicit authorization to
install or otherwise provide a pinned Rust toolchain. Prototype success cannot
authorize production adoption; that requires a later full equivalence gate.

## Scientific ownership

R owns input validation, canonical row ordering, source identities, essential-
interval accounting, H0/H1 separation, result certification, cache identity,
atomic artifacts, fallback, and every public API. Rust receives only two
canonical finite birth/death arrays and one requested dimension. It returns a
squared distance plus diagnostic counts and a status code.

The Rust kernel must calculate the same dissertation-aligned object: all finite
intervals, every consecutive active landscape level, zero padding at unequal
depth, exact piecewise-linear squared-L2 integration, no uniform grid, and no
universal level cap. H0 and H1 are separate calls and remain separate results.

## Prototype architecture

- Stable versioned C ABI; borrowed read-only R buffers and R-owned scalar
  outputs.
- Pair-bounded allocation, one internal thread, and at most two independent
  processes scheduled by R.
- Deterministic event ordering and compensated or pairwise summation.
- Panic cannot unwind across FFI; stable nonzero status codes return control to
  R.
- R recomputes canonically on any error, nonfinite/invalid output, failed
  certificate, or unavailable shared library. No partial Rust artifact survives.
- Unsafe Rust is denied except a minimal reviewed FFI shim.

Rust is forbidden from reading files or metadata, generating PH, constructing
cell/gene views, selecting levels by approximation, naming caches, combining H0
and H1 as a primary result, clustering, selecting `k`, or accessing labels and
outcomes.

## Equivalence ladder

1. Three analytical/sign-crossing fixtures within `1e-12` squared-distance
   absolute error.
2. Twenty tractable exact R-oracle results within `1e-10` absolute plus
   `1e-10` relative error.
3. The 12 MV5-BA worst-depth H0/H1 certificates before any speed decision.
4. All 318 exact MV5-AY dimension results before production adoption.
5. All 90 adaptive-certified MV5-AY H1 results before production adoption.

Tier 3 is sufficient only to assess the prototype. Tiers 4 and 5 plus a
separate adoption decision are mandatory before the Rust result can be used in
production.

## Performance and resource gates

The prototype must achieve at least a threefold median speedup over R on the
same frozen six-pair panel, with no individual pair slower than R. Peak RSS is
limited to 1 GiB per pair-bounded process. Outputs and identities must repeat
exactly across two clean builds/runs. The R package must pass an exact-staged
check both when Rust is absent and when the prototype is available.

Failure of any equivalence, memory-safety, determinism, speed, resource,
fallback, or package gate rejects the prototype without changing the R path.

## Toolchain boundary

No Rust toolchain is currently present in the Ubuntu environment. MV5-BC needs
owner approval before downloading/installing a pinned toolchain. Its version,
installer/source hashes, crate lockfile, and build-command identity must be
frozen before compilation. No installation is authorized by this prefreeze.
