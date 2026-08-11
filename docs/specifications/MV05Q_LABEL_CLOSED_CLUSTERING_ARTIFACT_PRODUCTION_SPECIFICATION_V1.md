# MV5-Q label-closed clustering-artifact production specification v1

Date frozen: 2026-08-10  
Starting accepted commit: `eb5c812`  
Parent scientific contract: `MV05N_LABEL_CLOSED_CLUSTERING_RESOURCE_GATE_SPECIFICATION_V1`  
Outcome-label state: closed

## Purpose

MV5-Q materializes auditable sample-clustering artifacts from the accepted MV5-P training distances and the already accepted MV5-D5/MV5-J held-out-to-training distances. It does not reinterpret persistence landscapes, select scientific winners, or inspect biological or technical labels.

## Immutable analysis axis

The queue contains exactly 150 analysis groups: 15 leave-one-study-out folds, two representations (`sct_whole` and `inductive_integrated`), and five frozen distances. Every analysis group uses exactly seeds 20260805 through 20260809 and one common ordered training-sample axis across those seeds.

The frozen distance panel is:

1. exact all-active-level H0 persistence-landscape L2;
2. exact all-active-level H1 persistence-landscape L2;
3. descriptive raw `sqrt(H0^2 + H1^2)` with no component rescaling;
4. representation-matched cell-distribution energy distance;
5. shared training-standardized pseudobulk Euclidean context baseline.

The explicit alias registry is authoritative. In particular, the accepted SCT query representation `sct_fold` maps to training representation `sct_whole`, and both analysis representations reuse the one exactly verified shared pseudobulk training/query distance source.

## Training-only clustering

For every group, PAM is fit for every integer `k = 2:min(10, n_training - 1)` independently in each of the five seed-specific training matrices. Selection uses no held-out sample and no label:

- calculate adjusted Rand agreement for all ten pairs of seed partitions at each k;
- use their mean as stability;
- calculate leave-one-seed-out jackknife Monte Carlo SE;
- define the one-SE threshold as the maximum mean stability minus the SE at the first maximum in ascending k order;
- select the smallest k whose mean stability reaches that threshold.

PAM is then retained at selected k for each seed. Average linkage is fit only as the frozen sensitivity at the PAM-selected k. Cluster IDs are canonicalized solely by lexicographically sorting each cluster's sorted member-ID signature.

## Held-out assignment

Held-out samples never affect training matrices, candidate fits, k selection, medoids, dendrograms, or cluster names.

- PAM: assign to the nearest frozen training medoid; an exact distance tie is resolved by the lexicographically smallest medoid sample ID.
- Average linkage: assign to the frozen cluster with minimum mean held-out-to-member dissimilarity; an exact tie is resolved by the lexicographically smallest canonical sorted-member signature.

## Execution and artifact policy

The prefreeze binds accepted source hashes, implementation hashes, aliases, queue identities, validation requirements, and abort rules before production. Production uses at most two workers, a 900-second group cap, and a 4-GiB process-tree RSS cap. Each group writes private candidate partitions and public label-closed selected artifacts atomically with hash-bound completion status. A valid existing completion is immutable and is reused only after hash validation.

The prespecified clean repeat is every representation-distance group in the fold with the maximum training-sample count (canonical fold-ID tie break). All 150 groups undergo immutable-resume validation.

## Required public evidence

Public evidence may contain only identifiers, distances or clustering artifacts derived without labels, stability agreement among seed partitions, provenance, runtime/resource values, validation state, and closed-label counters. It must not contain tissue, approach, class, biological/technical label, outcome, ARI/NMI against any label, or method-ranking claims.

## Stop boundary

MV5-Q stops before labels are opened; training-alignment or held-out-generalization ARI/NMI; method, representation, component, fold, or tissue ranking; spectral clustering; robustness execution; gene topology; cell/gene fusion; new data; Rust or other optimization; package-default changes; manuscript claims; PDFs, reviewer material, private-cache tracking, `example_run.r`, pushing, or modification of accepted MV5-D through MV5-P artifacts.
