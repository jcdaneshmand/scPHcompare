# MV6-F stage-1 launch prefreeze

| Field | Frozen value |
|---|---|
| Date | 2026-08-14 |
| Parent prefreeze | `dbf0f9b` |
| Group | `mv06f_group_v1:fbed6ad04f8243313ed439ecb5f29ddd43326a478d9b60fb21ff84be70b6ebf1` |
| Fold / seed | `large_loso_v1:SRA779509` / `20260807` |
| Sample roles | 25 held-out / 65 training |
| Biological / component pairs | 1,625 / 6,500 |
| Queue root | `f5471633e21d229eeabecadf12989dece2a3a7ab5b5d09f4584b0c3b6410bb5d` |
| Implementation root | `599074b3cd078cf27eb4a85148eb1df2ce3f84a5bdfd3160617b80a78f78c05e` |
| Rust SHA-256 | `51d3fca4808c6620c46c7b1100ddfad4753af02d37f87ec3eb278008846b160d` |
| Monitor SHA-256 | `9bcd7a3cd4462683a151e224e1f994c088bc411b04cc8ff2f04100c5fb36f8a0` |
| Production launched | No |

The external stage-1 monitor is committed before execution. It revalidates the
queue, implementation, and Rust roots, then launches exactly one group runner
with all numerical-library thread counts set to one. Every 0.25 seconds it
measures the complete process-tree RSS and current private-root size.

The process is killed and new work stops if elapsed time exceeds 1,200 seconds,
process-tree RSS exceeds 6 GiB, or the private root exceeds 10 GiB. Output and
error logs remain private. The only public execution artifact is an atomic
technical resource metric containing counts, hashes, resource use, disposition,
closed-label state, and zero downstream counters.

The monitor refuses a pre-existing public metric, unexpected partial group,
stale source, mismatched Rust library, or changed queue/implementation root. A
successful process must also pass the complete group-directory validator before
the public metric can report `completed`.

This launch prefreeze authorizes only the one maximum group. It does not
authorize its repeat automatically, the other 74 groups, fusion, clustering,
outcome access, public Rust adoption, manuscript claims, a release, or a push.
