# MV6-G complete-production rebind prefreeze

## Disposition

This sprint freezes a general dynamic group runner for MV6-G but authorizes
only a byte-equivalence reexecution of the already accepted maximum group.
The remaining 74 groups, metadata, outcomes, advanced fusion, clustering, and
claims remain closed until that rebind passes and a separate serial execution
policy is committed.

## Frozen scope

| Quantity | Remaining total |
|---|---:|
| Groups | 74 |
| Training biological pairs | 260,595 |
| Training component rows | 1,042,380 |
| Query biological pairs | 33,725 |
| Nine-method ranking rows | 303,525 |
| Training samples per group | 65–89 |
| Held-out samples per group | 1–25 |

The implementation root
`9bf8614d8e2dbdfce43792e74f08620712674c8830770c7c8d70b1fea432a71c`
binds the dynamic scale/ranking module, atomic group runner, prefreeze
builder/validator/repeat checker, externally monitored rebind runner,
equivalence validator, and focused tests. The parent fusion contract, queue,
accepted source groups, Rust library, exact all-active-level landscape
definition, nine methods, tie rule, and label firewall are unchanged.

## Resource and stop contract

Any future group remains serial with one worker, a 1,800-second and 12-GiB
process-tree cap, 5-GiB aggregate private-storage cap, 12-worker-hour total cap,
and no automatic retry. The first failure, partial artifact, identity drift,
nonfinite/negative distance, degenerate scale, or firewall violation stops new
work.

Before any of the 74 groups, the general runner must recompute the accepted
stage-one group and reproduce `training-distances.csv`, `scales.csv`, and
`rankings.csv` byte-for-byte while passing the same live caps. This prefreeze
does not itself grant production authorization.

## Validation

The exact 74-row queue independently reconstructs all workload totals. Nine of
nine categories pass, and a clean rebuild reproduces all three generated
prefreeze artifacts byte-for-byte. The focused suite passes 9/9 expectations;
the complete package suite passes 1,585 tests with zero failures or warnings
and two established skips in 167.7 seconds.

| Evidence | SHA-256 |
|---|---|
| Policy | `33e818cba4d6fbf18abd99f0a861221c0a4a1cea2f8bacf095f0e9280bf07c32` |
| Queue | `dc74ebbfa0feeb56712229ece3b73bafb622ed105da2de008ef2e0bd6f878486` |
| Sources | `8cd75d7c24c38eff852f86d69704bf97f2637a31992383c48d8fab0988dfbe5f` |
| Validation | `1d06ad7e368851a5a91ed61738dcbd8fecfe61d740da427bee94237b9f4d65ec` |
| Repeat ledger | `8ab85ac8391a1736c56ea40ff1cf524971cc0f90f863568ef3c092994cb09302` |

The next action is only the monitored maximum-group general-runner rebind and
three-artifact equivalence validator.
