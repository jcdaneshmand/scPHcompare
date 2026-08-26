# MV17-A source-inventory prefreeze attempt-1 failure

The exact-head `206ad0b164d69d6bc2bce153d2c13f072b60ea0e` MV17-A builder
stopped at its first accepted-closure manifest gate. It created no audit file,
opened no scientific result, and ran no null, localization, or H2 computation.

Read-only diagnosis showed that every MV8-W manifest byte count and SHA-256
matched. The failure was caused solely by R's `vapply()` retaining path names on
the observed hash vector while the expected CSV vector was unnamed; elementwise
equality was complete, but `identical()` correctly rejected the attribute
difference.

Recovery is limited to removing irrelevant vector names before exact equality
checks. Source roots, source values, scientific contracts, counts, tolerances,
and downstream firewalls are unchanged. Attempt 1 remains preserved as an empty
local root; it must not be reused. A distinct v2 audit root may be built only
after this failure record and correction are committed.
