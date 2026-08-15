# MV6-G serial completion prefreeze

## Disposition

This sprint prospectively authorizes the remaining 74 label-closed MV6-G
scale/ranking groups under a separately hash-bound serial execution layer.
The admitted scientific runner is unchanged. No production group existed
while this gate was built and validated.

The execution layer must stop before launching a later group after any child
failure, partial directory, stale identity, artifact-validation failure,
30-minute group runtime breach, 12-GiB process-tree RSS breach, 5-GiB private
storage breach, or 12-worker-hour aggregate breach. It uses one worker and no
automatic retry.

## Frozen workload

| Quantity | Remaining total |
|---|---:|
| Groups | 74 |
| Training biological pairs | 260,595 |
| Training component rows | 1,042,380 |
| Query biological pairs | 33,725 |
| Query ranking rows | 303,525 |

The full corpus after successful completion will contain 75 groups, 262,675
training pairs, 1,050,700 training component rows, 35,350 query pairs, and
318,150 nine-method ranking rows. All component scales remain fold-training
only. Cell H0, cell H1, gene H0, and gene H1 remain separate inputs; the fixed
fusion panel changes no landscape definition.

## Bound identities

| Identity | SHA-256 |
|---|---|
| Accepted scientific implementation | `8b0a1e42d9e46234edb847bcab82a2e770a47f3934d21f2a81c352c74c9cec0c` |
| Serial execution implementation | `deb03fbc3e3bcf49d7c6fcbd356256c94aae0e3dae7672ec4150b3702e36f745` |
| Rust library | `51d3fca4808c6620c46c7b1100ddfad4753af02d37f87ec3eb278008846b160d` |
| Parent queue root | `f5471633e21d229eeabecadf12989dece2a3a7ab5b5d09f4584b0c3b6410bb5d` |
| Accepted rebind evidence | `3bc9ac4aa0f2f1b58d166086c4057b853398d1a15047d51906ded7aa7716e771` |

## Validation

The independent prefreeze validator passes 10/10 categories. A clean rebuild
reproduces policy, queue, and source inventory byte-for-byte (3/3). The
focused completion suite passes 7/7 expectations. The package-aware complete
suite passes 1,601 expectations with zero failures and two established skips
in 159.2 seconds. The current root includes private child stdout/stderr capture
and binds both log hashes into every resource metric.

| Evidence | SHA-256 |
|---|---|
| Completion policy | `34b7db132236eb504db57abcb2daa7d10d803691a26b7323e0c09d04201dba63` |
| Completion queue | `5176f7659be8319e033a5431508bb25329047d0242ff50a1d0c65ad66480cf3f` |
| Completion source inventory | `74fcb09cd35fc947f9b88e2b39bd174070137b987210aa35f43191a63abd934b` |
| Independent validation | `fb75cc859ff44f7181b07247683c9937f726a1180a033ee68ea4560a338034fb` |
| Byte-repeat ledger | `334d9d70b283272aba4e07cf55c2bb1e5bd91d8788fc67b00e0b1784c4f5b6b3` |

## Authorized next action

After this prefreeze is committed, execute the 74 groups in exact queue order.
Then require the independent 75-group complete-corpus validator and, in a
separate committed gate, the 445-file SHA-256/byte/mtime immutable-resume
check. Labels, endpoint evaluation, clustering, inference, and claims remain
closed throughout production and resume validation.
