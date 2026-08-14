# MV5-AI selection-resistant nested-256 continuation gate

Date completed: 2026-08-11

Prospective engine HEAD: `87e59e095736e49debe8353f4d6c82b63ade12de`

Decision: authorize only a later label-closed calculation of
`nested_cells_256_pc30_euclidean_v1`.

## Evidence and selection firewall

The gate bound all 24 estimands, 24 intervals, and four primary tests from each
of the accepted PC20, cosine-chord, and nested-192 analyses. No result slice was
used. The decision helper received no representation, homology dimension,
tissue, endpoint, seed, estimate, interval, p-value, method rank, winner, or
subgroup value. Nested 256 remained the fourth and only unexecuted configuration
in the immutable MV5-V order.

The 43-source freeze was hash-exact. Independent validation passed all 11
categories, including reconstruction of 150 matched nested-192/nested-256 axes
with identical coordinate-source, pair-axis, and private-source identities.
The later calculation must verify actual per-view cell-ID inclusion before any
artifact is accepted.

## Authorized calculation scope

- 150 groups across 15 folds, five seeds, and two separate representations
- 30 coordinates and Euclidean point geometry
- 13,500 views and 70,700 directed biological pairs
- 141,400 separate exact all-active-level H0/H1 landscape rows
- 70,700 matched-energy rows and 282,800 four-method rows
- one worker; 600 seconds and 4 GiB RSS per group; six worker-hours and 4 GiB
  private storage for the configuration

The six real nested-256 admissions used at most 33.220 seconds and 493,096,960
bytes RSS. The full nested-192 execution used 1.233 worker-hours, at most 53.303
seconds per group, 454,893,568 bytes peak RSS, and 1,254,700,479 bytes private
storage. A conservative quadratic `(256 / 192)^2` projection gives 2.192
worker-hours, 94.760 seconds per group, 808,699,676 bytes RSS, and 2,230,578,629
bytes private storage, all within the frozen envelope.

## Reproducibility and boundary

Two clean builds produced byte-identical SHA-256 hashes for all 11 gate ledgers
and the independent-validation ledger. Calculation, rankings, labels, outcomes,
clustering, other configurations, and manuscript claims remained zero.

This authorization is a feasibility and design decision only. It is not a
robustness, equivalence, invariance, superiority, or default claim.
