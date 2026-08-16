# MV7-H landscape-stress runner rebind

Date: 2026-08-16

Status: accepted; clean v2 stress launch authorized

## Failed v1 launch

The first authorized stress launch failed closed after 3.657 seconds at
127,791,104 bytes peak process-tree RSS. The private receipt has disposition
`failed`, exit status 1, and no distance artifact. No landscape group directory
and no public stress evidence were published. Labels, outcomes, clustering,
dimension combination, and the remaining 19 groups stayed closed.

The stderr was `Error in ph_keys[[first_id]] : subscript out of bounds`. The
runner assigned sample-ID names to the interval list but not to the parallel
PH-cache-key vector, then indexed both by sample ID. The failure occurred before
the first pairwise distance calculation.

## Correction and unchanged contract

Commit `b573d831bde5c163fbbeddcebc524e139b88138b` assigns the frozen ordered
124-sample axis as names on `ph_keys`, matching the existing named interval
list. A static contract test requires both name assignments and both
sample-ID-keyed lookups.

The correction does not change the persistence diagrams, interval extraction,
landscape mathematics, pair axis, Rust kernel, output schema, resource caps, or
scientific estimand. The launch remains restricted to group
`mv07h_landscape_group_v1:eb7602ef5b12f26de930bf1af58d988e216d373f2881b593c68e3ce3e44e94aa`:
seed 20260807, gene topology, H1, all 7,626 unordered sample pairs. The accepted
Rust SHA-256 remains
`51d3fca4808c6620c46c7b1100ddfad4753af02d37f87ec3eb278008846b160d`.

Focused MV7-H validation passes 63/63. The canonical repository suite passes
with only its two established optional MV5-BC skips. The failed v1 private root
is retained as evidence. Relaunch is authorized only in a new v2 private root
under an exact Git HEAD containing this correction and this audit.
