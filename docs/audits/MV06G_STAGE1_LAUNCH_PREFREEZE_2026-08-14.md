# MV6-G stage-one launch prefreeze

## Disposition

The exact maximum-group scale/ranking sentinel is ready to execute under the
already frozen stage-one authorization. No calculation was launched while
this contract was written or validated. Full production, metadata access,
outcome evaluation, advanced fusion, clustering, defaults, release, and claims
remain closed.

## Bound identities

| Item | Identity |
|---|---|
| Parent prefreeze commit | `95d7615` |
| Parent contract | `f72326c0411c6c954ffb570fbf3e019adb7b09f82b74cd96984409e762077f8b` |
| MV6-F queue root | `f5471633e21d229eeabecadf12989dece2a3a7ab5b5d09f4584b0c3b6410bb5d` |
| Stage-one implementation root | `6a76a11d1b211fcf89fddcd67b7591161619950023420c0bfccbbdc65b76ce82` |
| Source diagrams | `25e9780dc76b0da2fb9e2d503247eb0c6fb915d2d665b41cebc982e69ef0c832` |
| Source query distances | `776bf94392cb6c37a0fdc9243032f1f007ea580d497252bc4cb172ece1dd6a4c` |
| Rust library | `51d3fca4808c6620c46c7b1100ddfad4753af02d37f87ec3eb278008846b160d` |

The stage-one root binds ten sources: the scientific helpers, launch builder,
independent launch validator and repeat checker, atomic runner, external
process-tree monitor, independent scientific reconstructor, grouped-Persim
oracle, immutable-resume checker, and focused tests.

## Exact execution scope

The runner reads the accepted 180-record `diagrams.rds` and 6,500-row
query-to-training distance table. It does not reconstruct views or rerun PH.
For the 65 training samples, it computes all 2,080 unordered sample pairs for
cell H0, cell H1, gene H0, and gene H1: 8,320 exact all-active-level component
rows. It fits four strictly positive training medians, normalizes the 1,625
accepted query pairs, and emits 14,625 rows for the nine fixed methods.

Rankings are ascending by normalized distance with canonical training sample
ID as the exact-tie rule. No query contributes to scale fitting. No label,
tissue, approach, endpoint, or outcome file is an input.

## Atomic and resource contract

The primary and clean repeat each use one worker, 1,800 seconds, 12 GiB of
process-tree RSS, 5 GiB of private storage, and no automatic retry. Each run
publishes one directory atomically only after validating:

- `training-distances.csv` (8,320 rows);
- `scales.csv` (four rows);
- `rankings.csv` (14,625 rows);
- `metrics.csv`; and
- `status.csv` binding every scientific artifact and parent identity.

Any partial directory, stale source, hash drift, nonzero Rust status,
nonfinite/negative distance, degenerate scale, cardinality error, resource
breach, or one-sided artifact set stops execution.

## Required admission evidence

After primary and clean repeat, the committed validators must require:

1. byte identity for training distances, scales, and rankings;
2. independent reconstruction of all four medians and all 14,625 formulas and
   ranks;
3. 12 R and 12 grouped-Persim checks spanning minimum, median, and maximum
   combined interval depth in each cell/gene H0/H1 component;
4. external elapsed/RSS/storage evidence and a complete-production projection
   no greater than 12 worker-hours; and
5. five-file immutable resume with SHA-256, byte size, and modification time
   unchanged.

## Prefreeze validation

All 11 independent launch categories pass, the focused suite passes 19
expectations, the Python oracle script parses, and a clean rebuild reproduces
both generated launch artifacts byte-for-byte. The package-aware suite passes
1,576 tests with zero failures or warnings and two established skips in 164.9
seconds.

| Evidence | SHA-256 |
|---|---|
| Launch contract | `53a99e7bed397101d602e38b65d9b6689c882411cfc9f37f5f08e3f0a686a26d` |
| Source inventory | `05fda40f7a37a2b39d05fec88cda2da0b271c1b4e342bdb6f1dc4dd7c44551c9` |
| 11-category validation | `b3a67d9320401eda11acecb06f0d8c9e4f5a01a4a1279bfc36395528a36c686d` |
| Two-artifact byte repeat | `1039c3b2f3d1aa309d56b65ea226be173e749344ff2146b0ae131ad326194e8e` |

## Resume-portability correction

The first primary/repeat calculation under root `883bbe32…16e2` passed its
scientific, resource, repeat, scale/ranking reconstruction, and R/Persim gates.
The immutable-resume checker then passed an unquoted Rust-library path through
`system2()`, so the child saw the wrong argument count and exited before
validating reuse. All five before/after hashes, byte sizes, and modification
times remained identical. Those outputs are quarantined as scientifically
valid but root-superseded attempt evidence.

The checker now uses `processx::run()` with an argument vector, and a focused
regression requires that safe route. No scientific runner, formula, source
artifact, landscape definition, cap, or output schema changed. The corrected
ten-source root above was independently rebuilt, passed 11/11 categories and
2/2 byte repeat, passed the full 1,576-test suite, and must be committed before
a clean sentinel reexecution.

The only authorized next action is the monitored primary sentinel followed by
one clean repeat and the frozen admission gates.
