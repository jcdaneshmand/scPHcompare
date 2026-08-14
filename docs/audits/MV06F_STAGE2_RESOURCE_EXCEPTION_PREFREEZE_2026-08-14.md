# MV6-F stage-two first-group resource exception prefreeze

| Field | Frozen value |
|---|---|
| Date | 2026-08-14 |
| Stage-two admission commit | `d468d76` |
| Group | `mv06f_group_v1:1fdd87ee5e65ad31c7f286b5b15dc9e13314d98d6ca81e353781a52bde5bbac7` |
| Fold / seed | `large_loso_v1:SRA779509` / `20260805` |
| First disposition | `group_rss_cap_exceeded` |
| Observed process-tree RSS | 8,747,204,608 B |
| Frozen production cap | 8,589,934,592 B |
| Diagnostic cap | 12,884,901,888 B |
| Scope | One-group resource diagnosis only |

## Stop and diagnosis

The serial stage-two monitor stopped on execution-order group 2 after 228.451
seconds when complete process-tree RSS reached 8.747 GB, 157.270 MB above the
prospective 8-GiB per-group cap. Exit status was `-9`; the partial directory
contains no published scientific artifact. No later group launched, no retry
occurred, and labels/fusion/clustering/outcomes remained closed.

The accepted maximum sentinel used the same 25-held-out/65-training workload
but seed `20260807` and peaked near 3.1 GB. The breach is therefore not explained
by pair count alone. It may reflect seed-specific PH complexity or transient
Ripser memory; changing or dropping H1 is not authorized.

## Bounded exception

The complete resource contract already reserves 12 GiB aggregate live RSS and
execution is serial. A separate monitor is therefore prefrozen to rerun only
this exact group with the unchanged scientific runner, queue, implementation
root, Rust binary, inputs, 1,800-second cap, and 10-GiB storage cap, while
raising only the diagnostic process-tree RSS ceiling to 12 GiB.

This is not an automatic retry or a revised production default. The prior
partial state must be moved intact to quarantine before launch. The diagnostic
monitor and authorization are bound by SHA-256. If the group exceeds 12 GiB,
fails numerically, or does not publish a fully validated atomic directory, the
diagnosis stops and the project must profile/optimize PH before further
production.

If it completes, its full measured peak and elapsed time will determine whether
the existing serial monitor can resume through validated reuse or whether a
separate resource-policy revision is required. No second pending group is
authorized by this exception.

Machine-readable evidence is in
`docs/audits/mv06f-stage2-production-evidence/`. Private logs and the failed
partial directory remain excluded from Git.
