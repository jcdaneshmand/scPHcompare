# MV7-H full-PH repeat-monitor rebind prefreeze

Date: 2026-08-15

Status: accepted forward-only implementation rebind; production authorized

## Finding and correction

A final prelaunch code-path review found that the full-PH monitor attempted to
row-bind the source-repeat queue row and the twelve PH-repeat rows before
constructing their common receipt ledger. Those queue schemas intentionally
differ, so R would stop at the receipt-comparison step after production rather
than close the repeat gate.

No MV7-H source, PH, landscape, clustering, label, or outcome job had started.
Commit `cbc75a4` changed only the repeat receipt axis to select the two shared
fields, `job_id` and `output_file`, before row-binding. A focused regression
test now covers the heterogeneous source/PH queue normalization and the full
MV7-H focused suite passes 22 expectations. The complete regression suite also
passes 367 result groups with zero failures or warnings and the same four
established skips.

## Scientific equivalence

The v2 prefreeze is bound to commit `cbc75a4` and implementation root
`f7e14cc0c7ff8cb693054db7a7d9b8da13c476d34db1327b217ba44996f282bd`.
The original root
`bd68734d0f28ab15e801a74a6ba2092abe21bdc6766d2d46b3302238ea86d5de`
is retained as audit history.

Fourteen scientific/configuration artifacts are byte-identical between v1 and
v2: panel, cache manifest, sentinel axis, full sample-seed axis, source queue,
PH queue, landscape queue, resource caps/projection, landscape contract, label
firewall, runtime, acceptance, and decision. In particular:

- the 1,240-row PH queue SHA-256 remains
  `b6b328e045e85f7041d3ae0dd0dbee07736c53b16347adaf423b60c193ac9eb0`;
- the 20-group landscape queue SHA-256 remains
  `6a85cfe91fa15989d9fe170ad0be492074a22e5bba030d491468d066a461362e`;
- the Rust SHA-256 remains
  `51d3fca4808c6620c46c7b1100ddfad4753af02d37f87ec3eb278008846b160d`;
  and
- the eight-item all-active-level, H0/H1-separated landscape definition is
  unchanged.

## Validation and decision

The v2 builder repeated 16/16 artifacts byte-for-byte and its independent
validator passed 14/14 categories. The correction cannot alter matrices,
typed views, PH, distance values, pair scope, resources, or scientific stopping
rules. It only permits the already required 13-artifact source/PH repeat ledger
to be assembled after execution.

The v2 evidence under `docs/audits/mv07h-prefreeze-evidence-v2/` supersedes v1
for execution identity. The original evidence remains immutable history. The
exact v2 full source/PH run and then one stress landscape group are authorized;
the other 19 landscape groups, H0/H1 combination, clustering, labels, outcomes,
and claims remain closed.
