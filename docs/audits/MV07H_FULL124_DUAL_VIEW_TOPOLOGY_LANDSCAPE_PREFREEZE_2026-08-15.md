# MV7-H full-124 dual-view topology and landscape prefreeze audit

Date: 2026-08-15

Status: prefreeze complete; full source/PH and one landscape stress group
authorized

## Decision

MV7-H passes its prospective implementation and workload gate. The exact
source/PH corpus may run serially, followed by independent validation and only
the single frozen maximum-burden landscape stress group. The other 19 landscape
groups remain closed until the stress group passes repeat, R-reference,
analytic, resource, and immutable-resume validation.

No labels, outcomes, clustering, dimension combination, biological claim, or
new-data operation is authorized.

## Scientific continuity

The implementation is frozen at `1f2960d`. It reuses, rather than refits, the
five accepted MV7-G global standardization and 30-PC models. A real private
sentinel reconstruction reproduced both cell- and gene-view cache keys and
payload hashes exactly. The six sentinel pairs per seed must reproduce 60/60
accepted MV7-G typed views during production.

The two views remain intentionally distinct:

- cell topology: 384 cells as points in the shared 30-PC Euclidean space; and
- gene topology: 500 panel genes as points under explicit Pearson-correlation
  chord distance over the same 384 cells.

Complete Vietoris-Rips H0/H1 persistence uses field 2 and no filtration
truncation. Every full H0 finite-death multiset must match an independently
calculated MST; H1 must be finite, positive, and nonempty.

## Exact workload

The independently reconstructed queues contain:

- 124 samples by five seeds, or 620 sample-seed states;
- five reused-transform source bundles and 1,240 typed views;
- 1,240 PH records, balanced 620 cell and 620 gene;
- 7,626 unordered sample pairs per seed;
- 76,260 view-specific pairs across seeds and views; and
- 152,520 separately stored H0/H1 component-distance rows.

Landscape production is split into 20 atomic groups: five seeds by two views
by two homology dimensions. Each group contains exactly 7,626 rows. The frozen
stress group is seed `20260807`, `gene_topology_v1`, H1; its six-sentinel
interval burden of 10,435 was the greatest of the 20 prospectively ordered
groups. Tissue, approach, study, and outcomes did not enter this ordering.

## Landscape definition

The accepted dissertation-aligned definition is unchanged:

- only finite positive-persistence intervals enter;
- the essential H0 interval is excluded;
- every consecutive active level is used, with zero padding across unequal
  depths;
- H0 and H1 are calculated and published separately;
- integration is exact squared L2 over critical piecewise-linear segments;
- no universal grid or landscape-level cap is permitted; and
- output is streamed by bounded group, never as a dense landscape-function
  matrix.

The accepted Rust library SHA-256 is
`51d3fca4808c6620c46c7b1100ddfad4753af02d37f87ec3eb278008846b160d`.
Rust is the bounded execution engine; the accepted R engine remains the
canonical numerical oracle.

## Resource and recovery contract

Execution is one worker, zero retries, and atomic publication. Source groups
have 1,800-second/8-GiB caps; cell PH 600-second/4-GiB caps; gene PH
1,800-second/8-GiB caps; and landscape groups 3,600-second/12-GiB caps. The
aggregate boundary is 48 worker-hours, 12-GiB peak process-tree RSS, and 4-GiB
retained private state.

The conservative estimate is 6.187 worker-hours for source/PH plus 6.453
worker-hours for landscapes, or 12.640 worker-hours total. The landscape term
scales the complete prior MV6-F combined production time by component-row
count, so it is an admission estimate, not a promised runtime. The prospective
stress group remains mandatory because H1 interval depth can affect cost.

Every source, PH, and landscape output has a hash-and-size receipt. Existing
output without a matching receipt fails closed. One full seed source plus its
12 sentinel PH records must repeat byte-for-byte; the stress distance CSV must
repeat byte-for-byte; and accepted artifacts must retain hash, size, and mtime
across whole-stage resume.

## Validation result

The builder reproduced all 16 public prefreeze artifacts byte-for-byte. The
independent validator passed 14/14 categories, including 32 source hashes,
complete axis/queue reconstruction, balanced component counts, the exact Rust
identity, resource limits, the eight-item landscape contract, and the label
firewall. The focused MV7-H contract suite passed 20 expectations. The complete
source-loaded package suite passed 366 reported test result groups with zero
failures, zero warnings, and the same four established environment/CRAN skips.

## Public/private boundary

Public evidence contains queue identities, relative artifact names, source and
implementation hashes, counts, caps, contract text, and decisions. It contains
no expression matrix, typed-view payload, persistence diagram, private local
path, label value, outcome, PDF, Rust binary, or generated binary. All future
scientific payloads and process logs remain under ignored `tmp/` state.
