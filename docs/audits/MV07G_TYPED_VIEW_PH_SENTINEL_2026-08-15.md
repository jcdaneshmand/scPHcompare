# MV7-G typed-view and PH sentinel audit

Date: 2026-08-15

Status: complete; MV7-H full-PH/landscape prefreeze authorized

## Scope and boundary

MV7-G tested the exact full-124 descriptive transform and both corrected
topology orientations before full PH or landscape-distance calculation. The
sentinel used six prospectively frozen MV7-D sample identities across five
seeds. Tissue/depth information selected these identities earlier but was not
read during execution. Labels and outcomes remained closed.

For each seed, the runner reopened and identity-validated all 124 SCT caches,
used the exact MV7-FP 500-gene panel, standardized 47,616 equally sampled cells,
and fit one shared 30-component PCA. It then constructed six 384-cell Euclidean
PC views and six 500-gene Pearson-correlation chord-distance views. The exact
workload was five global-fit/source bundles and 60 complete Vietoris-Rips H0/H1
records. Landscape, distance, clustering, label, and outcome jobs were zero.

The implementation was frozen at `2ca771f`; explicit CSV-null normalization was
corrected at `0f3757a`; and the independently validated execution freeze was
committed at `b506046`.

## Production result

All 65 primary jobs completed serially with one worker and zero retries.

- five global fits averaged 192.147 seconds (185.963 to 194.822 seconds);
- maximum global-fit process-tree RSS was 2,198,278,144 bytes;
- 30 cell-PH jobs averaged 1.615 seconds and used at most 114,343,936 bytes;
- 30 gene-PH jobs averaged 2.283 seconds and used at most 157,941,760 bytes;
- aggregate primary plus representative-repeat worker time was 1,324.384
  seconds; and
- accepted private production state was 48,176,762 bytes.

Every cell diagram had 384 points and 383 finite H0 intervals. Every gene
diagram had 500 points and 499 finite H0 intervals. All 60 finite H0 death
multisets matched their independently calculated MST edge multisets with
maximum absolute error zero.

All H1 intervals were finite and had positive persistence. Cell views contained
319 to 442 H1 intervals; gene views contained 998 to 2,672 H1 intervals. This
supports retaining H1 as a first-class, separately reported topology component;
it does not yet establish a biological or clustering advantage.

## Independent and deterministic validation

The independent validator reopened all five source bundles and all 60 PH
records, revalidated typed orientation and identities, reran every full-view
MST comparison, checked finite-positive H1 structure, and passed 11/11
categories.

For the first seed, the first 32 ordered points of all 12 sentinel views were
recomputed with Ripserr and `TDA::ripsDiag`/GUDHI. All 24 H0/H1 comparisons
passed with maximum absolute error zero, including the explicit arbitrary
distance route for gene topology.

One complete global-fit bundle plus its 12 PH records repeated exactly:
13/13 serialized artifacts were byte-identical and SHA-identical.

## Full-corpus projection

The deliberately conservative projection scales the complete source bundle
linearly from six to 124 stored views even though the expensive global fit is
already performed once per seed. It projects:

- global fits and all typed views: 5.515 worker-hours and 759,215,627 bytes;
- 620 cell-PH records: 0.278 worker-hours and 18,688,557 bytes; and
- 620 gene-PH records: 0.393 worker-hours and 51,265,072 bytes.

The total is 6.187 worker-hours and 829,169,255 bytes, below the prospective
24-worker-hour and 8-GiB admission limits.

## Immutable-resume amendments

Before closure, whole-sentinel resume exposed three validation-only issues:

1. The original monitor refused to start when complete public evidence already
   existed, even though all private receipt-reuse logic was present. Commit
   `44860b8` added a no-overwrite public resume path and an 83-artifact
   hash/size/mtime validator.
2. Regenerated CSV bytes differed after type/precision round-tripping even when
   their tables were semantically equal. All 83 accepted artifacts were already
   proven unchanged 7/7. Commit `b84328b` changed only the regenerated-table
   comparison to strict column/row equality plus `1e-12` numeric tolerance;
   accepted public files remained read-only.
3. The regenerated decision initially counted the private resume snapshot in
   production-state bytes. Commit `d0ce1b5` restricted that accounting to the
   frozen production namespaces and exactly reproduced the accepted
   48,176,762-byte total.

The final whole-sentinel resume reused every receipt, semantically reproduced
all five production tables, and performed no overwrite. Across 78 private
scientific artifacts and five public production files, hashes, byte sizes, and
mtimes were unchanged with maximum mtime delta zero. The attempt and final
7/7 ledgers and the amendment register are retained publicly.

## Public and private separation

Public evidence contains aggregate source/PH metrics, hashes, interval counts,
MST errors, cross-engine results, repeat/resume ledgers, resource projection,
and decisions. It contains no expression matrices, typed-view payloads,
persistence diagrams, private paths, labels, outcomes, PDFs, or generated
binaries. Global transform/view bundles, PH records, logs, receipts, and the
resume snapshot remain under ignored private `tmp/` state.

## Regression gate

The complete source-loaded package suite passed after the production and
immutable-resume amendments: 362 reported test result groups, zero failures,
zero warnings, and the same four established environment/CRAN skips. No
accepted scientific artifact was rewritten by this verification.

## Decision

MV7-G passes. The exact cell-as-observation and gene-as-observation views are
both feasible and correctly separated. H0 and H1 are valid in both views, and
the measured resource projection authorizes only a separate MV7-H full-PH and
landscape prefreeze. It does not directly authorize full execution, label
access, clustering selection, or manuscript claims.

The revised landscape definition remains unchanged: finite positive intervals,
essential H0 excluded, all consecutive active levels, H0 and H1 separate,
exact or error-controlled squared-L2 distance, no universal fixed grid or
level cap, and streamed or chunked computation.
