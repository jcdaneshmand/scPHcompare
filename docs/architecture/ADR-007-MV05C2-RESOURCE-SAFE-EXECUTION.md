# ADR-007: MV5-C2 resource-safe execution and narrowed MV5-D boundary

## Status

Accepted on 2026-08-06 for label-closed SCT cell-primary precomputation. Full
MV5-D, outcome evaluation, gene-view execution, inductive integration,
full-matrix clustering, fusion, default changes, and manuscript claims are not
authorized.

## Context

MV5-C showed that the corrected definitions work on six samples but projected
the naïve 90-sample run at 151.46 representation worker-hours plus 91.05
landscape hours. Inspection of the actual worker found that every fold reran
`SCTransform` for every sample and seed even though the normalized sample
depends on the selected cells and seed, not the held-out study. The all-pairs
landscape stage also requested distances that are unnecessary for the frozen
primary held-out retrieval endpoint.

The statistical plan nevertheless contains secondary tasks that do require a
complete matrix. Reducing pair scope must therefore narrow the executable
endpoint set explicitly rather than claiming that query-to-training pairs can
support clustering or within-study contrasts.

## Decision

1. Cache SCT once for each sample–seed and reuse it across folds with immutable
   input/version identities and atomic publication.
2. Keep feature ranking, scaling, PCA, and any future integration reference fit
   inside each training partition. They are not moved into the sample cache.
3. Replace full within-fold landscape matrices with exact query-to-training
   requests for the primary continuous retrieval task only.
4. Keep H0/H1 separate and retain exact all-active-level integration. Pair
   reduction is not landscape approximation.
5. Run deterministic fold-local chunks with request and output hashes,
   structured resource guards, stale-cache refusal, and byte-immutable resume.
6. Defer full-matrix clustering and within-study pair contrasts. Do not fill
   missing matrix cells, symmetrize across incompatible fold fits, or infer
   them from another representation.
7. Reject full all-view MV5-D under the current cap. Its measured lower bound is
   25.82 hours before integrated reference mapping is added.
8. Conditionally permit a staged, label-closed SCT cell-primary execution whose
   current conservative projection is 18.68 worker-hours. Stop and re-project
   after normalization before launching all fold/seed diagrams.

## Evidence

- 450 cached normalizations replace 6,750 repeated operations (15-fold
  architectural reduction).
- Six real cache entries reproduce MV5-C SCT matrices exactly; build/resume are
  514.4/67.41 seconds at 4.72/0.49 GB peak RSS.
- Cached fold construction reproduces both MV5-C LOSO directions exactly in
  10.382–17.718 seconds per fold after cache loading.
- Query-to-training scope reduces landscape rows from 1,802,250 to 212,100.
- 250/250 real exact distances equal immutable MV5-C with zero difference.
- Resume reuses 50/50 chunks and preserves 100/100 artifact hashes.
- Full-view lower bound is 25.82 hours; SCT cell+gene is 22.25 hours without
  adequate reserve; SCT cell-primary is 18.68 hours.

## Consequences

The next safe action is an MV5-D0 SCT cell-primary label-closed precomputation,
starting with the 450 sample–seed normalization cache and a mandatory resource
reprojection. Gene topology remains scientifically valuable but conditional
and computationally deferred. Integrated cell topology needs its own cached
reference/mapping profile; the sample SCT cache alone does not prove it fits.

G-MV5 remains open because no biological or technical endpoint has been
computed. The resource architecture cannot be described as a completed
benchmark or evidence that topology works biologically.
