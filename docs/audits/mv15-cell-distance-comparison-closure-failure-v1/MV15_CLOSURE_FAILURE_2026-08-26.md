# MV15 independent-closure failure record

The prospectively authorized MV15 production completed successfully and remains
immutable. The first independent-closure invocation stopped before publishing
any closure artifact. Its failure was confined to the closure validator: three
neighborhood-summary difference vectors were supplied as separate arguments to
`abs()` instead of being combined into one vector first.

No production job was retried, no scientific value changed, and no result was
interpreted. Recovery requires a separately committed closure-only correction
and prospective recovery prefreeze. Production rerun and all downstream stages
remain closed.
