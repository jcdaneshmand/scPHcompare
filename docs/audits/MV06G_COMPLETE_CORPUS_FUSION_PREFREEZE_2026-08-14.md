# MV6-G complete-corpus fusion prefreeze audit

## Disposition

MV6-G prefreeze completes as
`prefreeze_pass_stage1_training_scale_sentinel_only`. The complete accepted
MV6-F corpus supports a prospective, label-closed fusion calculation, but only
the single maximum-group training-scale/ranking sentinel is authorized next.
Full production, label opening, outcome evaluation, advanced fusion,
clustering, defaults, release, and claims remain closed.

## Immutable parent identities

| Item | Accepted identity |
|---|---|
| Parent revision | `bba0b11` |
| MV6-F queue root | `f5471633e21d229eeabecadf12989dece2a3a7ab5b5d09f4584b0c3b6410bb5d` |
| Complete group inventory | `43c4af75cf00ecf8b38309b27415bace9f762fbbc8573ca1961a8a9302488b41` |
| Complete immutable resume | `2742cca6c6abdb0ef41dcea255824d807e2eb630b59344b469ea4eef91d701bb` |
| Accepted Rust library | `51d3fca4808c6620c46c7b1100ddfad4753af02d37f87ec3eb278008846b160d` |
| MV6-G implementation root | `ab7039f1d6f30d5a821557b970a37d62427af8be958d6fb610fd3b5706a531db` |
| Authoritative metadata | `e2b6f1e7d7dc05ba8e80627aad4ea1677c3b718a844397acd06596d7e22a18d0` |

The builder rehashed every accepted `diagrams.rds` and `distances.csv` against
the 75-row complete inventory. During audit development, an unnecessary
re-sort exposed R's compact-versus-explicit row-name serialization as a source
of a false queue digest. The published implementation resets row names before
hashing, reproduces the accepted MV6-F root above, and includes a regression
expectation for that exact identity.

## Frozen training-only calculation

Raw cell/gene and H0/H1 distances are not directly averaged. For every fold,
seed, and component, MV6-G will compute all unordered training-to-training
exact, all-active-level landscape distances and divide by their strictly
positive training median. Held-out-query distances, labels, clipping,
centering, and outcomes cannot affect these 300 scales.

| Quantity | Frozen total |
|---|---:|
| Fold-seed groups | 75 |
| Training biological pairs | 262,675 |
| Training component rows | 1,050,700 |
| Query-to-training biological pairs | 35,350 |
| Component scales | 300 |
| Ranking methods | 9 |
| Query ranking rows | 318,150 |

The nine methods are four normalized components, equal H0/H1 cell and gene
composites, equal cell/gene fusion as the sole primary fusion, and fixed
gene-weight sensitivities 0.25 and 0.75. Weights are never selected from
outcomes. Rankings use ascending distance and canonical training sample ID for
exact ties.

## Endpoint and label firewall

The primary endpoint is cross-study tissue mean reciprocal rank. Fixed 1-NN
balanced accuracy is supportive. The two-member primary MRR family compares
equal fusion separately with the cell and gene composites; fusion benefit
requires both differences to be positive. The contract freezes 2,000
tissue-stratified held-out-study bootstrap replicates (seed `20260815`), 9,999
paired held-out-study sign flips (seed `20260816`), and Holm adjustment across
exactly those two contrasts.

The existing 500-gene panel is explicitly disclosed as technically
transductive because label-free presence and variance were assessed across the
accepted corpus. Fold transforms and all component scales remain
training-only. This is an internal blocked evaluation, not external
validation.

All scales and nine-method rankings must be complete, repeated, resumed, and
bound into a committed prediction manifest before the authoritative metadata
is read. Tissue and approach fields are prohibited from prediction artifacts;
post-label scaling, reranking, weight selection, or method replacement is
prohibited.

## Validation evidence

The independent validator passes 12/12 categories: parent completion,
75 source-group identities, workload, method panel, endpoints/contrasts,
blocked inference, label firewall, resources, implementation identity,
contract summary, downstream firewall, and stage-one-only admission. The
focused test file passes 16 expectations. A clean independent build reproduces
all ten generated prefreeze artifacts byte-for-byte. The supported
package-aware suite passes 1,557 tests with zero failures or warnings and the
two established optional/public-artifact skips in 162.1 seconds.

| Evidence | SHA-256 |
|---|---|
| Contract | `f72326c0411c6c954ffb570fbf3e019adb7b09f82b74cd96984409e762077f8b` |
| 12-category validation | `e2740db4c72aa1c7839f7769e9470c8a389708f96cec513e8e5a1bb043ed71ec` |
| 10-artifact byte-repeat ledger | `5e55bde01c0e79c8e482240bf551e777be5f615dfe6f2a1716837a452efc66e4` |

No outcome label was opened, no biological endpoint was computed, and no
fusion evaluation or clustering job ran.

## Next authorized sprint

Execute only group
`mv06f_group_v1:fbed6ad04f8243313ed439ecb5f29ddd43326a478d9b60fb21ff84be70b6ebf1`:
2,080 training pairs, 8,320 training component rows, four training scales,
1,625 query pairs, and 14,625 label-closed ranking rows. Use one worker, a
1,800-second elapsed cap, 12-GiB process-tree RSS cap, 5-GiB private-storage
cap, atomic publication, and no retry. Full production may be considered only
after exact R/Persim oracles, clean repeat, immutable resume, and a projection
below 12 worker-hours pass.
