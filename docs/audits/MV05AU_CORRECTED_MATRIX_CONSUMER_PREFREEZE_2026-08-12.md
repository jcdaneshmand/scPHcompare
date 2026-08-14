# MV5-AU corrected landscape matrix consumer prefreeze

Date: 2026-08-12

Authorization: MV5-AT `9326726`; prospective contract `7db3e5c`.

## Outcome

MV5-AU passes and authorizes only MV5-AV: implement a verified read-only
corrected-sidecar loader and deterministic average-linkage dendrograms built
separately for H0 and H1.

The code audit distinguished true matrix consumers from sampled-function
consumers. Average-linkage is the unique safe first analytic consumer because
it is deterministic, label-free, needs no cluster count or kernel scale, and
retains the complete hierarchy. PAM and spectral clustering remain closed
pending a prospective `k` policy; Ward and other linkage methods remain closed
pending comparison. Betti/Euler, aggregated landscape, and cross-iteration
curve routines are schema-incompatible with a distance matrix.

## Frozen boundary

The future loader must verify the completion marker, pair shards, matrix class,
scientific contract, sidecar input and matrix keys, axes, and explicit
iteration/view identity. `view_id` must be supplied as `cell_topology_v1` or
`gene_topology_v1`; it cannot be inferred. The consumer returns separate H0
and H1 average trees. The combined matrix remains descriptive and cannot be
consumed. No partition, selected `k`, label, outcome, legacy alias, or source
mutation is permitted.

## Audit evidence

Ten candidate consumers were classified; only the verified loader and
average-linkage dendrogram were eligible. Two API contracts, six view/dimension
rows, six migration stages, 15 mandatory validations, and 14 abort rules were
frozen. Twelve source hashes reproduce, all 15 prohibited-change counters are
zero, nine public ledgers reproduce byte-for-byte across two clean builds, and
10 independent validation categories pass.

No corrected artifact was consumed and no source or workflow code changed.

## Next sprint

MV5-AV may implement and test only the loader and separate H0/H1 average trees,
including a bounded label-free smoke against existing MV5-AT artifacts.
Partitions, `k` selection, evaluation, visualization, defaults, legacy changes,
new data, claims, and optimization/Rust remain unauthorized.
