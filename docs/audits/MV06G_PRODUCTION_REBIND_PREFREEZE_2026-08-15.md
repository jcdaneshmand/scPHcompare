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

The corrected implementation root
`8b0a1e42d9e46234edb847bcab82a2e770a47f3934d21f2a81c352c74c9cec0c`
binds the dynamic scale/ranking module, atomic group runner, prefreeze
builder/validator/repeat checker, externally monitored rebind runner,
equivalence validator, and focused tests. It recomputes ranking blocks after
canonical sorting, correcting only the stale-index invariant. The parent
fusion contract, queue,
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
the complete package suite passes 1,601 tests with zero failures and two
established skips in 178.8 seconds.

| Evidence | SHA-256 |
|---|---|
| Policy | `91ccd518c30280226e6a0dd7ceff93706893b123327f4f8e353900001b5b6507` |
| Queue | `dc74ebbfa0feeb56712229ece3b73bafb622ed105da2de008ef2e0bd6f878486` |
| Sources | `2de1ede178b4ccd064b54604b4f332c3cc972b3864ba95654a9f2678b451e212` |
| Validation | `1d06ad7e368851a5a91ed61738dcbd8fecfe61d740da427bee94237b9f4d65ec` |
| Repeat ledger | `ce8d6406b0912ea57eca9b37638e7e9912c4d247a11bb065b1b5a0f50b254974` |

The next action is only the monitored maximum-group general-runner rebind and
three-artifact equivalence validator.
