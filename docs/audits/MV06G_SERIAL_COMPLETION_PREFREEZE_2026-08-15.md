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
| Accepted scientific implementation | `9bf8614d8e2dbdfce43792e74f08620712674c8830770c7c8d70b1fea432a71c` |
| Serial execution implementation | `38440b861a04af17e641635a5cf93e4f44bb465c5d3c72edfa864439ee020091` |
| Rust library | `51d3fca4808c6620c46c7b1100ddfad4753af02d37f87ec3eb278008846b160d` |
| Parent queue root | `f5471633e21d229eeabecadf12989dece2a3a7ab5b5d09f4584b0c3b6410bb5d` |
| Accepted rebind evidence | `3bc9ac4aa0f2f1b58d166086c4057b853398d1a15047d51906ded7aa7716e771` |

## Validation

The independent prefreeze validator passes 10/10 categories. A clean rebuild
reproduces policy, queue, and source inventory byte-for-byte (3/3). The
focused completion suite passes 7/7 expectations. The package-aware complete
suite passes 1,592 expectations with zero failures and two established skips
in 202.2 seconds.

| Evidence | SHA-256 |
|---|---|
| Completion policy | `5061a85916c47be2721407c49d067ab0496d846d23f2fd341d36663760a703cb` |
| Completion queue | `dc74ebbfa0feeb56712229ece3b73bafb622ed105da2de008ef2e0bd6f878486` |
| Completion source inventory | `bb39428b8e2616203333544adb7952ac772731b4708b0c99a8dd6198447ac2c3` |
| Independent validation | `fb75cc859ff44f7181b07247683c9937f726a1180a033ee68ea4560a338034fb` |
| Byte-repeat ledger | `2c4531d1022c8fef1735adbb89336f7c2a8b4a9f3a02db62cef522cbd075eb7a` |

## Authorized next action

After this prefreeze is committed, execute the 74 groups in exact queue order.
Then require the independent 75-group complete-corpus validator and, in a
separate committed gate, the 445-file SHA-256/byte/mtime immutable-resume
check. Labels, endpoint evaluation, clustering, inference, and claims remain
closed throughout production and resume validation.
