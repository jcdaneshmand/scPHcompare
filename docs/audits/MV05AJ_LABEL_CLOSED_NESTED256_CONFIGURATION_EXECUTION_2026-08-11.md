# MV5-AJ label-closed nested-256 configuration execution

Date: 2026-08-11

Authorization: MV5-AI commit `23640a48b31b569dd52878c6f7358e6082b297ba`

Prospective engine: `3d553cf086af39a99db4cf4f7c09f2b3ef4bb26c`

Runtime/source binding: `7162d0355dbbd97a4d81caf3edb2c35fc3b8abfe`

## Result

The exact authorized `nested_cells_256_pc30_euclidean_v1` calculation is
complete and accepted. All 150 prospectively frozen groups published
atomically under one worker. The private result tree contains 13,500
256-by-30 views, 70,700 directed heldout-training biological pairs, 27,000
separate H0/H1 landscape summaries, 141,400 exact all-active-consecutive-level
landscape rows, 70,700 matched-energy rows, and 282,800 four-method rows.

The calculation used 7,104.807 monitored worker-seconds (1.974 hours). The
maximum group took 95.101 seconds, peak process-tree RSS was 507,351,040
bytes, and final private storage was 1,664,639,370 bytes. These values are
below the frozen 600-second/4-GiB unit, six-hour configuration, and 4-GiB
storage caps.

## Exact nesting definition

For every sample, seed, and representation, the engine reconstructed the
deterministic `sha256_sample_seed_cell_nested_v1` ordering over sample ID,
accepted seed, and cell ID. The selected 256 point IDs were required to be
exactly positions 1 through 256 of the frozen 384-cell order. Positions 1
through 192 were separately required to be byte-identical to the previously
accepted nested-192 point-ID sequence. All 30 frozen coordinates and Euclidean
geometry were retained without refitting, normalization, or feature access.

## Independent numerical validation

The validator avoided production scientific helpers and independently
checked:

- all 13,500 selected point-ID and transformed-payload identities;
- exact 192-prefix and 256-of-384 subset identities in all 13,500 views;
- all 3,442,500 finite H0 MST deaths, with maximum absolute error
  `4.97e-14`;
- 60 stratified direct H1 diagrams over three fold strata, five seeds, and
  both representations;
- 30 stratified direct energy distances, with maximum absolute error
  `5.33e-15`;
- 60 stratified exact H0/H1 landscape distances using a separate Python
  implementation, with zero observed error;
- all 1,350 artifact-manifest hashes and all 282,800 four-method rows;
- every aggregate cardinality, exact all-active-consecutive-level policy,
  resource limit, and label/ranking/outcome closure field.

All 16 aggregate validation categories passed. The two prospectively selected
repeat groups reproduced all 16 deterministic artifacts byte-for-byte. A
clean second run of the geometry, landscape, and aggregate validators
reproduced all 11 validation ledgers byte-for-byte.

## Resume and boundary

A complete validation-only resume accepted all 150 groups as
`reused_validated`. Pre/post snapshots of all 1,650 private files were exact,
including relative path, SHA-256, byte size, and modification time.

The bound gate contained exactly 43 public source identities. No labels,
Tissue/Approach fields, or outcome data were opened, and no ranks, outcomes,
clustering, within-training pairs, method selection, alternate configuration,
gene/fusion/new-data/Rust/default activation, or manuscript claim was
produced. The result establishes a valid label-closed nested-256 distance
calculation only. The decision ledger authorizes only a separate,
prospectively frozen nested-256 prediction-locked outcome-evaluation sprint.
