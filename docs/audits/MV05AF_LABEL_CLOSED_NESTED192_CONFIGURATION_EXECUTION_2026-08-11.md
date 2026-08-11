# MV5-AF label-closed nested-192 configuration execution

Date: 2026-08-11

Authorization: MV5-AE commit `0dc460d`

Prospective engine: `08f8332`

Runtime/source binding: `e9ba7d2`

## Result

The exact authorized `nested_cells_192_pc30_euclidean_v1` calculation is
complete and accepted. All 150 groups published atomically under one worker.
The private result tree contains 13,500 192-by-30 views, 70,700 directed
heldout-training pairs, 27,000 H0/H1 landscape summaries, 141,400 exact
all-active-level landscape rows in 720 subchunks, 70,700 matched energy rows,
and 282,800 four-method rows.

The calculation used 4,439.486 monitored worker-seconds (1.233 hours). The
maximum group took 53.302 seconds, peak process-tree RSS was 454,893,568 bytes,
and final private storage was 1,254,700,479 bytes. These are below the frozen
600-second/4-GiB unit, six-hour configuration, and 4-GiB storage caps.

## Exact nesting definition

The code audit corrected an imprecise shorthand before launch. “First 192”
means the first 192 cells in the deterministic
`sha256_sample_seed_cell_nested_v1` ordering over sample ID, accepted seed,
and cell ID. It does not mean source matrix rows 1–192. Independent validation
reconstructed that ordering from every raw 384-cell view and proved that every
192-cell set is contained in the corresponding still-unexecuted 256-cell set.
All 30 frozen coordinates and Euclidean geometry were retained without refit or
normalization.

## Independent numerical validation

The validator avoided production scientific helpers and independently checked:

- all 13,500 selected point-ID and transformed-payload hashes;
- all 2,578,500 finite H0 MST deaths, with maximum absolute error
  `4.97e-14`;
- 60 stratified direct H1 diagrams over three fold strata, five seeds, and both
  representations;
- 30 stratified direct energy distances, with maximum absolute error
  `4.88e-15`;
- 60 stratified exact H0/H1 landscape distances using a separate Python
  implementation;
- all 1,350 artifact-manifest hashes and all 282,800 four-method rows;
- all aggregate cardinalities, exact all-active-level policies, resource
  limits, and label/ranking/outcome closure fields.

Fifteen aggregate validation categories passed. The two prospectively selected
minimum/maximum-fold group repeats reproduced all 16 deterministic artifacts
byte-for-byte. A clean second run of the geometry, landscape, and aggregate
validators reproduced 11/11 validation ledgers byte-for-byte.

## Resume and boundary

A complete validation-only resume accepted all 150 groups as
`reused_validated`. Pre/post snapshots of all 1,650 private files are
byte-identical, including relative path, SHA-256, size, and modification time.

No labels were opened and no ranks, outcomes, clustering, within-training
pairs, nested-256 calculation, gene/fusion/new-data/Rust/default activation, or
manuscript claim was produced. The result establishes a valid label-closed
nested-192 distance calculation only. A separate prospective gate is required
before any prediction-locked outcome evaluation.
