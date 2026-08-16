# MV7-H gene-H1 resource gate

Date: 2026-08-16

Status: production paused; owner decision required for exact-engine fallback

## Gate result

The v4 full-PH run completed all five 124-sample source bundles, eight cell PH
records, and seven gene PH records. The next gene record, seed 20260805 sample
`SRA628554_SRS2664364`, exceeded the frozen 8 GiB process-tree RSS cap after
124.595 seconds. Its measured peak was 8,678,952,960 bytes. The process was
killed, no output was published, and every earlier atomic artifact remains
receipt-valid. Labels and outcomes stayed closed.

This is not the prior small-sample Ripserr discrepancy. The cell-side PH for the
same sample completed, while the failing object is the complete 500-gene
distance filtration through H1.

## Exact-engine feasibility

The already corroborating TDA/GUDHI engine was run on the exact frozen failing
view under a 12 GiB and 30-minute diagnostic cap. It completed in 404.091
seconds at 6,650,802,176 bytes peak RSS. After removal of GUDHI's capped
essential H0 interval, the diagram contained 499 finite H0 intervals and 3,222
finite positive H1 intervals. Its H0 deaths matched the independent MST oracle
with maximum absolute error 0.

A second full 500-gene GUDHI run targeted the heaviest gene view already
completed by Ripserr. GUDHI produced exactly the same 499 H0 and 2,221 H1
birth/death pairs as Ripserr, with maximum absolute error 0 in both dimensions.
This extends the existing 24/24 smaller Ripserr/GUDHI checks to a full heavy
gene diagram.

The operational tradeoff is material: on that completed view, Ripserr used
13.347 seconds and 577,155,072 bytes, while GUDHI used 399.543 seconds and
6,514,917,376 bytes. A universal GUDHI switch would therefore be needlessly
slow. GUDHI is useful specifically because its memory use stayed below the cap
on the view where Ripserr did not.

## Recommended decision

Approve a prospective exact resource-fallback amendment:

- keep Ripserr primary for every cell and gene view;
- after a recorded Ripserr RSS-cap failure only, run the same typed view,
  complete Vietoris--Rips filtration, H0/H1 dimensions, field, and maximum
  scale through TDA/GUDHI under a separate 12 GiB/30-minute cap;
- normalize only the engine-specific capped essential H0 representation;
- record the selected engine and both attempts in public aggregate evidence;
- require MST validation, finite-positive H1 validation, deterministic repeat
  for every fallback artifact, and full interval equality wherever both engines
  complete;
- keep landscapes unchanged: all active levels, exact squared-L2 integration,
  no grid, no level cap, separate H0/H1, and streamed groups; and
- retain all labels, outcomes, clustering, combination, and claims closed.

This changes implementation policy, not the mathematical estimand. It still
requires owner approval because it introduces mixed PH-engine provenance and a
higher fallback memory allowance. Raising the Ripserr cap alone is not
recommended because later memory demand is unknown. Using GUDHI for every gene
view is not recommended because of the measured runtime. Reducing, truncating,
or omitting gene H1 is rejected at this gate because it would change the central
scientific definition.
