# MV5-P label-closed distance-production launch readiness

## Decision

MV5-P is technically ready to launch the complete frozen distance queues.
The production monitor must stop after distance construction and validation;
clustering, held-out assignment, labels, and outcomes remain prohibited.

## Preserved MV5-O contract

No accepted MV5-D through MV5-O artifact or implementation was modified. The
MV5-O source root remains
`541e7d3aa8acce5d512bbb4819c034735eef47387e91a63abccfa259f53d6de1`.
The exact all-active H0/H1 landscape runner and baseline runner retain their
frozen hashes.

The first full-group admission exposed that the frozen stager expects its PH
manifest argument to contain exactly the samples participating in the
training-pair manifest. MV5-P resolves this at orchestration scope: it derives
a deterministic private training-only PH manifest and matching metrics subset
from the accepted label-closed public sources, validates their complete job and
sample identities, and passes those paths to the unchanged frozen stager.
Held-out diagrams are neither changed nor used in training distances.

## Runtime and admission

An ignored project-local Python environment pins Persim 0.3.8. Its sorted
runtime freeze and 47-file installed Persim manifest are SHA-256 bound. The
installed runtime reproduces the accepted 32-row real-diagram fixture
byte-for-byte.

The complete execution-order-1 integrated group then passed end to end:

- 85 training samples and no held-out intervals;
- 7,140 exact H0/H1 landscape distances in 30 chunks;
- 3,570 energy distances in one unit;
- all distances finite and nonnegative;
- all-active levels, no cap, separate H0/H1, and exact integration preserved;
- labels closed, biological outcomes false, and clustering count zero; and
- two immediate reruns rebuilt zero units and preserved every artifact byte,
  size, and timestamp.

## Orchestration guards

The monitor enforces deterministic queue order, at most two concurrent group
workers, a 900-second group bound, 4-GiB process-tree RSS, 10-GiB private
storage, and projected-plus-observed aggregate worker accounting against the
21.6-hour cap. It validates the 18-file source freeze, runtime hashes, helper
hash, and worker hash before launch; the worker hash is rechecked before every
new group. Partial files, stale status, cap failures, or nonzero downstream
counters stop further launches.

The initial monitor stopped after two completed groups because its metrics
replacement guard misread R's successful `unlink()` return code. Completed
distance artifacts were preserved and no partial file existed. The correction
uses a same-filesystem direct atomic rename, verified against an existing file
on the production mount. A clean two-group replay measured process-tree peaks
of 392,876,032 and 393,621,504 bytes and reproduced every training input,
landscape output, and energy output byte-for-byte. Recovery metrics use the
more conservative replay elapsed times and declare
`clean_resource_replay_byte_identical` RSS provenance before deterministic
resume at execution order 3.

## Verification

- focused MV5-P tests: 17/17;
- complete source suite: 793/793, zero failures, warnings, or skips;
- clean Git-archive package: `R CMD check --no-manual`, `Status: OK`;
- clean tarball SHA-256:
  `bcae7a33a024c1d378f982184eefb194a67caafadb698229aa03fe47c9f572e2`;
- check-log SHA-256:
  `f9c00171104e3bae7aecaab651b1baed567fc060e9628975358e38bfcaa273ce`.

The next authorized action is only the complete label-closed distance run,
followed by the frozen oracles, repeats, resume, matrix, and resource checks.
