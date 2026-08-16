# MV7-H source-axis rebind prefreeze

Date: 2026-08-16

Status: accepted forward-only execution rebind; full source/PH production authorized

## What happened

The first v2 production unit, `source__20260805`, stopped after 2.389 seconds
before publishing a source bundle or any PH record. The private receipt and
stderr are retained. Labels and outcomes remained closed.

The parent PCA model stores its 124 ordered `fit_sample_ids` with each value
also repeated as a vector name. The frozen CSV stores the same 124 ordered
values without names. The original use of `identical()` therefore rejected an
attribute difference even though the type, length, membership, values, and
order were equal. After that first guard was corrected, a second construction
guard exposed the same redundant-name issue when canonicalizing the source
record sample axis.

## Correction and real-data proof

Commit `8ceb0ad` introduced two narrowly scoped helpers: one compares ordered
axes after removing names only, and one sorts and removes names for the
canonical source-record axis. The checks still reject order changes, type
changes, membership changes, and length changes. The focused suite passes 26
expectations.

A real seed-20260805 reconstruction then reopened all 124 frozen normalization
caches, rebuilt 248 cell/gene typed views, passed the full source-record
validator, and reproduced all 12 accepted MV7-G sentinel typed views exactly.
The private source record is 146,569,565 bytes with SHA-256
`5f25228552bbf6481b29de12052ae03aa7462590676009edb7770efd38aebf33`.
No private path or matrix is published.

Commit `ecb392b` also corrected the generic prefreeze wording from an historical
claim about zero jobs to the precise invariant that the prefreeze contains zero
execution outputs. This truthfully accommodates the preserved failed receipt
without placing execution state inside a prospective contract.

## Scientific equivalence and final binding

The v4 prefreeze is bound to exact commit
`ecb392b031c9b2887408fc2e130568346db9690a` and implementation root
`7684dbeeffb01c94bd39f82d0dbeae605cd266fe790cb4962ebf47a9bd8fb43b`.
Thirteen scientific/configuration artifacts are byte-identical to v2,
including the panel, cache and sample axes, source/PH/landscape queues, resource
policy, landscape contract, label firewall, runtime, and decision. The
acceptance artifact differs only in the execution-output wording above.

The 1,240-row PH queue remains
`b6b328e045e85f7041d3ae0dd0dbee07736c53b16347adaf423b60c193ac9eb0`;
the 20-group landscape queue remains
`6a85cfe91fa15989d9fe170ad0be492074a22e5bba030d491468d066a461362e`;
and the accepted Rust binary remains
`51d3fca4808c6620c46c7b1100ddfad4753af02d37f87ec3eb278008846b160d`.
The landscape definition remains finite positive intervals, essential H0
excluded, all active consecutive levels, exact squared-L2 integration, no grid,
no level cap, separate H0/H1 components, and streamed groups.

Both v4 builds pass 14/14 independent categories and 17/17 byte-repeat files.
The complete repository suite passes 1,690 expectations in 368 result groups
with zero failures or warnings and the same four established skips.
Only `docs/audits/mv07h-prefreeze-evidence-v4/` may identify the next execution.
The full source/PH run remains authorized; landscape execution still stops after
the single prospectively selected stress group unless its repeat, R-oracle,
analytic-oracle, resource, and immutable-resume gates pass.
