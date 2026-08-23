# MV8-S same-axis external baseline and PH-sentinel prefreeze v1

**Status:** prospective execution contract

**Parent:** MV8-R full topology-production prefreeze, decision D-118

## Purpose

MV8-S is the bounded execution gate between completed Pearson-residual source
production and the 1,280-record corrected topology workload. It establishes the
missing current-input external SCT-data baseline on the exact selected-384 axes
and tests the frozen PH implementation before any full production commitment.

## Authorized work

- Reconstruct eight external common-475 baselines from each unit's current
  exact-reference filtered/raw H5 inputs, using the frozen seed-20260805
  selected-384 digest and selected-cell-only SCTransform fit.
- Build the frozen 30-PC cell view and correlation-chord gene view for each
  baseline, then repeat every baseline byte-for-byte.
- Compute exactly 23 H0/H1-separate PH records: all 16 external baseline views
  plus seven source-produced gene-view sentinels, followed by byte-identical
  repeats.
- Validate every full record against its typed-view H0 MST oracle.
- Compare Ripserr with TDA/GUDHI on three 32-point checks and one full external
  gene-view check, separately for H0 and H1 at tolerance `1e-6`.

## Execution and resource policy

Execution is serial: one runner, one child, one worker, and zero automatic
retries. Cell PH uses Ripserr under 600 seconds and 4 GiB. Gene PH uses Ripserr
under 1,800 seconds and 8 GiB. Only a monitored Ripserr RSS-cap stop on a gene
view permits the already approved exact TDA/GUDHI fallback under 1,800 seconds
and 12 GiB. That transition is a fallback, not a retry. The total worst-case
elapsed allowance is 126,600 seconds and private output is capped at 8 GiB.

All outputs are atomic and resumable only from rehashed completed records or a
recorded eligible Ripserr RSS-cap stop. Any other prior state fails closed.

## Landscape definition and firewall

MV8-S does not compute landscapes. The dissertation-aligned future definition
remains unchanged: H0 and H1 are primary separate outputs; include every finite
positive interval; exclude essential H0; retain every consecutive active level;
use exact or error-controlled squared-L2 integration; impose no universal grid
and no universal level cap; and stream or chunk rather than materializing dense
landscape matrices. Landscape execution additionally requires rebuilding and
hash-rebinding the Rust kernel and passing canonical R oracles.

The remaining 1,257 PH records, all landscapes, comparisons, clustering,
fusion, labels, outcomes, default adoption, cleanup, deletion, and manuscript
claims remain closed. A successful run still requires independent MV8-T closure
before a separate full-PH prefreeze can be considered.
