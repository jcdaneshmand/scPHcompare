# MV8-VB bootstrap type-recovery prefreeze

**Date:** 2026-08-23

**Result:** 12/12 checks pass; no root or scientific artifact was created.

The committed MV8-VA bootstrap compared identical byte values with R's storage-type-sensitive identical() predicate. CSV parsing supplied integer values while file.info supplied doubles, so validation stopped before creating either fresh root.

MV8-VB changes only that predicate to numeric value equality. The first logged launch's manually mistyped head and the mechanically Windows-Git-bound repeat are both preserved. Neither copied, recomputed, or retried a record. After commit, one corrected bootstrap is authorized.
