# MV17-C bounded H0/H1 calibration sentinel specification

Date: 2026-08-26  
Status: prospective implementation specification; no real source selection or null execution authorized

## Eligible primary axis

Use exactly one label-blind view per accepted unit and topology: exact500 panel,
seed 20260805, 132 units (124 internal plus eight external), cell coordinates as
384 selected cells by 30 shared PCs, and gene coordinates as 500 genes by 384
selected cells under the already closed Pearson-residual representation. H0 and
H1 remain separate. Existing diagrams and source matrices must be rehashed
against MV17-A/MV17-B closures before selection.

## Sentinel selection

For each of the cell and gene views, order the 132 eligible units by finite H1
interval count, then by a canonical private identity token. Select positions 1,
66, and 132 as minimum-, median-, and maximum-burden sentinels. Use those same
three view-specific units for H0 and H1 so dimension comparisons do not change
the source unit. Public evidence may publish roles and aggregate counts only;
unit IDs, paths, barcodes, and matrices remain in an ignored private locator.

## Null and summary contract

Run nine fixed replicates per compatible family. Cell uses coordinate
permutation, covariance Gaussian, and radial-density clouds; gene uses those
three plus within-row axis shuffling. This yields 189 complete H0/H1 PH jobs.
Seeds are fixed by `mv17c_null_seed_registry_v1()` before source selection.

For each dimension, retain finite positive intervals only and calculate:

1. finite interval count;
2. total persistence;
3. maximum persistence, or zero when empty; and
4. exact all-active-level landscape squared L2 norm, equal to the sum of
   persistence cubed divided by 12 across intervals.

Use the one-sided plus-one empirical tail `(1 + #{null >= observed}) / 10`.
With nine replicates, resolution is 0.10 and worst-case Monte Carlo standard
error is approximately 0.158. This sentinel is a feasibility/calibration-shape
gate, not inferential evidence; no multiplicity correction or significance
claim is permitted. A full calibration prefreeze must increase precision.

## Repeat, resources, and stops

Independently regenerate every maximum-burden job (63 jobs) and require exact
scientific agreement. Use one worker and zero automatic retries. Per child caps
are 1,800 seconds and 8 GiB RSS; aggregate caps are two hours, 2 GiB private,
and 32 MiB public. Stop without retry on source drift, singular covariance,
non-finite output, PH/MST failure, nondeterminism, stderr, or any cap breach.

Labels, outcomes, tissues, clustering, fusion, view ranking, biology, manuscript
claims, full-corpus calibration, real localization, H2, cleanup, and deletion
remain closed. Passing MV17-C authorizes only a separate MV17-G prefreeze after
the already closed MV17-D2 qualification is also bound.
