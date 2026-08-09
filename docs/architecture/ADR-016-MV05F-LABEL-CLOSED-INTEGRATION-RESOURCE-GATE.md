# ADR-016: Reconstruct integration from fixed D1 panels and gate full execution

## Status

Accepted and executed on 2026-08-08.

## Context

MV5-E found no primary SCT H0/H1 retrieval improvement over matched cell energy,
but this cannot tune or discard a prospectively planned integration comparison.
The accepted D0 caches do not retain Seurat SCT models, and the earlier MV5-C
pilot selected a new query-compatible panel rather than preserving each D1 panel.

## Decision

Rebuild each training reference from accepted raw shards using the exact D1
500-gene panel and D0-selected cells. Fit joint reference SCT and 30-PC PCA using
training samples only. Rebuild and map held-out queries independently, using
only frozen-panel features available in that query and never replacing missing
genes. Use transfer projection and integrated embeddings without label transfer.
Require the full reference digest to remain unchanged.

Select a four-fold real-data pilot solely from D1 structural workload evidence,
execute one seed with one worker under explicit caps, and project the complete
75-group mapping plus PH, landscape, baseline, and retrieval workload. Stop
before full production and all outcomes.

## Evidence

Four groups covering 1, 4, 8, and 25 held-out samples completed in 325.8 to
840.9 seconds. Peak process-tree RSS was 2.67 to 3.18 GiB. All 360 coordinate
views and 38 query mappings passed; references were unchanged. Queries retained
498 to 500 frozen-panel genes, with no replacements.

Independent validation passed seven categories for all bundles. Two clean
same-implementation admission builds produced byte-identical RDS SHA-256
`aa7d7c46ef0c4fa19855fb51c57529a4d33658990ddfa38b00d91ced2d98cd39`.

With a 25% reserve, the complete projection is 14.96 worker-hours, 1.54 GB
storage, 3.97 GiB peak RSS, and 1,051 seconds for the worst group. These pass
the 21.6-hour, 10-GB, 8-GiB, and 1,800-second caps.

## Consequences

A separately authorized full label-closed integrated cell-view execution is
technically feasible. It is not executed or biologically evaluated here. The
full run must preserve all frozen axes and must not use MV5-E outcomes for
tuning. Clustering, gene topology, fusion, new data, and claims remain gated.
