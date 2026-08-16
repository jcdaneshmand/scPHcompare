# MV7-H exact PH resource-fallback prefreeze v1

Date: 2026-08-16

Status: prospectively specified after owner approval and before fallback
production

## Purpose

Complete the frozen 124-sample dual-view PH corpus without changing its
mathematical estimand when a gene-side Ripserr process exceeds the accepted
8 GiB resource cap.

## Invariants

The fallback does not change:

- the exact 124 samples, five seeds, 384 matched cells, or 500-gene panel;
- cell-side PCA coordinates or gene-side correlation-chord distances;
- complete Vietoris--Rips persistence through H1 over field 2;
- the threshold-minus-one full-filtration contract;
- separate H0 and H1 persistence diagrams;
- exclusion of essential H0 from landscapes;
- all active landscape levels, exact squared-L2 integration, no grid, no
  universal level cap, or streamed within-seed groups; or
- the closed state of labels, outcomes, clustering, dimension combination, and
  manuscript claims.

## Trigger and engine policy

Ripserr remains the primary PH engine for every view. GUDHI is eligible only
when all of the following are true:

1. the queued stage is `gene_ph`;
2. the typed view is `gene_topology_v1`;
3. the primary attempt ended specifically as `rss_cap_exceeded`;
4. no primary output was atomically published; and
5. no completed fallback receipt already exists.

Other failures remain fatal. Cell views cannot use the fallback. A qualifying
job receives one TDA/GUDHI attempt, one worker, a 1,800-second elapsed cap, and
a 12 GiB process-tree RSS cap. The original Ripserr receipt and the GUDHI
receipt are kept in separate immutable ledgers.

## Engine normalization

TDA/GUDHI is run on the exact 500-by-500 arbitrary distance matrix at its
maximum finite distance. GUDHI's capped essential H0 interval is removed and
the contract's single infinite essential H0 interval is restored. No other
interval is removed, truncated, rounded, rescaled, or imputed.

Each accepted fallback record must pass the existing typed-view provenance,
finite-H0 MST, finite-positive H1, payload-hash, cache-key, and atomic-output
validators. The selected engine and version are published per PH record.

## Repeat and independent validation

Every selected fallback artifact is repeated with the same engine even if it
is outside the original sentinel repeat set. The serialized PH record must be
byte-identical. Independent validation must prove that each fallback has:

- exactly one preceding primary RSS-cap failure;
- one completed GUDHI receipt below 12 GiB and 1,800 seconds;
- a valid typed PH record and MST oracle;
- a byte-identical repeat; and
- full H0/H1 interval equality wherever both engines complete.

The accepted feasibility evidence already includes zero-error full-scale H0
and H1 equality on a 500-gene, 2,221-H1-interval view and a successful GUDHI
calculation of the previously failing 3,222-H1-interval view.

## Resume and stop boundary

The five valid source bundles and fifteen valid PH records are reused only
after their receipt hashes and byte sizes are frozen and revalidated. The
failed primary output must remain absent. This amendment authorizes completion
and independent validation of the PH corpus only. Landscape production remains
closed until full-PH validation explicitly authorizes the one previously
selected stress group.
