# MV5-AV corrected landscape tree consumer implementation

Date: 2026-08-12

Authorization: MV5-AU `5e81118`; API implementation `1be734e`; smoke validator `d8f83c9`.

## Outcome

MV5-AV passes. Two additive public APIs now verify and load completed corrected
sidecars and build deterministic average-linkage trees independently for H0 and
H1. Explicit iteration and cell/gene view identity are mandatory. The combined
matrix is verified as descriptive but is never consumed.

The result contains no combined tree, partition, selected `k`, labels, or
outcomes. Existing workflow formals, defaults, legacy fields, clustering,
visualization, Betti, modular, and cross-iteration paths remain unchanged.

## Validation

Twenty-six focused expectations cover artifact identity, rejection, separate
matrices, descriptive-combined consistency, deterministic serialization,
permutation equivariance, source immutability, and absent partitions.

Two clean read-only smokes each loaded all eight MV5-AT sidecars and built 16
trees: H0/H1 for four cell and four gene views. Public tree ledgers and both
seven-category independent validation outputs are byte-identical. Every source
artifact path, size, modification time, and SHA-256 remained unchanged.

## Decision

The next eligible sprint is MV5-AW: prospectively define whether and how a
tree may become a partition, including a label-closed `k` policy and stability
assessment. No partition, evaluation, visualization, default change, legacy
change, new data, claim, or optimization is authorized yet.
