# MV5-AZ label-closed stability/resampling prefreeze v1

Date: 2026-08-12
Authorization: accepted MV5-AY production `a78f731` and continuing owner direction

## Purpose and boundary

MV5-AZ freezes the evidence needed to identify a non-oracle sample partition
from corrected persistence-landscape distances. It does not cut a tree, run
PAM, select `k`, access a biological label, or authorize the additional seed
calculation. The next executable scope is numerical-equivalence and speed
benchmarking only.

## Existing resampling evidence

Three existing asset families answer different questions and must not be
treated as interchangeable.

1. MV03 Stage C contains five seeds and both cell/gene views across eight
   strata, but only two samples per stratum. It measures technical diagram and
   distance reproducibility; it cannot identify a nontrivial `k` grid.
2. MV05C contains genuine six-sample/five-seed matrices within five fold-
   specific strata. It is useful as a method and acceleration-equivalence
   reference, but gene-view and representation coverage is incomplete and its
   fold axes are not the accepted ten-sample MV5-AY target.
3. MV5-AY supplies the corrected complete baseline matrices on the intended
   four- or ten-sample axes, but only for seed 20260805. It freezes geometry and
   axes; one seed cannot estimate partition stability.

No existing asset identifies the current primary `k`. Reusing MV03 or MV05C as
if it did would conflate technical reproducibility, a different fold target,
and current partition stability.

## Primary matched-axis target

The first adequately sized target is the existing large ten-sample panel,
analyzed as four independent representation/view strata:

- SCT whole, cell topology;
- SCT whole, gene topology;
- Seurat integration, cell topology; and
- Seurat integration, gene topology.

Each stratum must contain the same ten sample IDs in the same order at all five
frozen seeds. The accepted seed 20260805 contributes 40 existing diagrams and
180 existing unordered pairs across these four strata. A future launch would
add four seeds, 160 diagrams, and 720 unordered pairs. Under the accepted R
baseline, 360 added gene-view H1 pairs require the adaptive route. H0 and H1
remain separate. Bone four-sample panels remain descriptive/sensitivity scope;
they are not promoted as the first primary stability target.

## Frozen stability rule

For each representation, view, homology dimension, and candidate `k`:

1. candidate `k` is `2:min(10,n-1)`, hence 2 through 9 at `n=10`;
2. run PAM independently on each of the five seed-specific dissimilarity
   matrices, after a deterministic medoid and tie rule is frozen;
3. compute all ten pairwise adjusted Rand indices among the five partitions;
4. use their mean as the stability statistic;
5. estimate uncertainty by deleting one seed at a time and applying the
   delete-one-seed jackknife to the mean-pairwise-ARI statistic; and
6. select the smallest `k` whose mean stability is within one jackknife standard
   error of the maximum mean stability.

The Euclidean H0/H1 combination is descriptive only. Cross-view or
cross-representation fusion cannot select `k` in this stage. Biological labels,
outcomes, and oracle `k` remain closed until inputs, stability summaries,
selected `k`, and partitions are immutable.

## Landscape definition

Every future distance must preserve the accepted dissertation-aligned contract:
all finite intervals, all consecutive active levels, zero padding for missing
depth, no universal level cap, no uniform grid, separate H0/H1, and exact or
strict `1e-8` error-controlled L2 integration. Acceleration may change the
implementation but not this mathematical object.

## Acceleration before expansion

The first candidate is the already used corrected Persim 0.3.8 critical-pair
path with the project-controlled exact linear-segment integral. Persim's
built-in norm remains prohibited because it failed sign-changing fixtures.

Before adoption, the candidate must:

- agree with analytical/sign-changing fixtures within `1e-12` absolute error;
- agree on all 114 MV5-AY exact-only pairs within `1e-10` absolute plus `1e-10`
  relative squared-distance tolerance;
- agree on all 90 adaptive-H1 pairs within each accepted R achieved-error
  certificate plus a scale-aware floating allowance;
- repeat normalized results and identities exactly;
- stay below 2 GiB RSS per process with no more than two independent processes;
- be reproducibly dependency-locked; and
- achieve at least a threefold median wall-time speedup on a prospectively
  frozen high-depth gene-H1 benchmark panel to replace the R path.

If equivalence passes but speed does not, retain the accepted R implementation.
Rust begins only if the mature corrected engine misses throughput/packaging or
the projected full run breaches the resource cap. Rust must pass the same
equivalence corpus; language choice never relaxes scientific tolerance.

## Resource and abort rules

A future five-seed launch must project no more than 40 worker-hours, 12
wall-hours with at most two independent processes, 2 GiB RSS per process, and
20 GiB retained private artifacts. Stop on axis mismatch, source/hash drift,
landscape-contract drift, failed certificate/equivalence, nondeterminism,
resource breach, label/outcome access, or any attempt to consume a partition
before authorization.

## Next sprint

MV5-BA is authorized to benchmark corrected Persim equivalence and speed against
the accepted MV5-AY shards and selected public-safe aggregate fixtures. It may
recommend R, corrected Persim, or a Rust proof-of-concept gate. It may not
generate the four additional seeds or authorize partitions.
