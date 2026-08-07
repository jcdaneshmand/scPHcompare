# ADR-006: MV5-C existing-data feasibility and cell-projection boundary

## Status

Accepted on 2026-08-06 for label-closed one-tissue feasibility. This decision
does not authorize outcome evaluation, method ranking, full-cohort execution,
fusion, a public-default change, or manuscript claims.

## Context

ADR-005 established synthetic inductive mapping but did not demonstrate that
real existing-data folds, complete H0/H1 diagrams, exact all-level landscapes,
and matched baselines were jointly feasible. The prospective bounded selection
rule chose the smallest eligible tissue with at least two studies and at least
384 filtered cells in every sample. It selected six samples from two liver
studies, with a real 5-versus-1 study imbalance.

The existing dual-view source contract required every panel gene to vary in
every sample. That condition is necessary for the Pearson-correlation gene
view but not for projecting held-out cells through a PCA fitted on training
samples. Applying it to both views incorrectly made a held-out constant gene
invalidate cell topology.

## Decision

1. Keep the prospectively selected six-sample/two-study pilot and five seeds;
   do not switch tissues after observing execution behavior.
2. Keep feature ranking, scaling, PCA, and integration reference fitting within
   the training study of each fold. Study remains the split variable; tissue
   and approach remain absent from execution artifacts.
3. Admit `scph_cell_projection_source_v1` only for training-fitted cell
   projection and pseudobulk means. It permits within-held-out-sample constant
   genes but remains ineligible for gene topology and gene-correlation
   baselines.
4. Preserve `scph_dual_view_source_v1` and its nonconstant-gene requirement for
   intentional gene topology. Do not select features using held-out variance
   merely to force a complete gene result.
5. Use `IntegrateEmbeddings()` only for inductive cell coordinates. Record
   integrated gene topology as structurally unavailable because this API does
   not produce a corrected held-out gene-expression matrix.
6. Retain H0 and H1 separately and use exact, all-active-level landscape L2.
   Do not introduce a landscape grid or level cap.
7. Build held-out neighbor rankings, PAM stability, one-standard-error `k`, and
   average-linkage partitions while labels remain closed. These are immutable
   predictions/partitions, not endpoint estimates.
8. Do not start MV5-D under the current implementation. A 90-sample projection
   fails both the 30-minute per-job and 24-worker-hour aggregate caps; cached
   normalization, chunked execution, and reduced pair scope are required first.

## Evidence

- Ten fold/seed jobs completed with zero canonical process failures in 4,760.71
  observed worker-seconds; median 484.685 seconds, maximum 556.26 seconds, and
  peak observed RSS 2,950,356,992 bytes.
- SCT and inductively integrated cell topology completed in all ten jobs.
- SCT gene topology completed in all five 5-reference/1-query jobs and was
  structurally unavailable in all five 1-reference/5-query jobs because the
  training-only panel contained constant genes in held-out samples.
- Integrated gene topology was structurally unavailable in all ten jobs by
  method definition.
- Thirty held-out mappings completed with 39–132 anchors and immutable query
  embedding/reference hashes.
- The exact landscape engine produced 750 finite distances over 150 diagrams
  and 25 completed strata in 136.41 internal seconds; the measured process run
  used 354,390,016 bytes peak RSS.
- Eighty-five immutable distance matrices, 425 held-out neighbor-ranking rows,
  and 17 complete five-seed clustering groups were built without outcomes.

## Consequences

The cell view is technically feasible in both LOSO directions. The gene view
is conditionally feasible and cannot yet be described as a complete
cross-study secondary view for this tissue. This is a useful negative result,
not permission to use held-out variance for feature repair. G-MV5 remains open
because no biological or technical endpoint has been evaluated and neither
view has a scientific works/fails/uncertain disposition.

The next authorized sprint is MV5-C2/MV5-D engineering specification: cache
normalization across reusable sample/seed scopes, define query-to-training pair
scope, stream/chunk distance construction, and re-project resources. Outcome
labels must remain closed until that gate passes.
