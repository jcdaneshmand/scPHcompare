# MV17-A source-inventory prefreeze attempt-3 failure

The exact-head `73062c20e10598a2075d8041f9bab5d99fd3c1fb` attempt completed
the full immutable-source rehash and produced a local v3 audit. Post-build
aggregate review rejected that audit before commit because the public
`gene_landscape` row reported `axes=28`, inherited from the mixed MV8 closure,
instead of the verified 26 gene-only dimension groups.

The gene PH counts (1,264 artifacts; 2,528 dimension records), gene landscape
counts (626 chunks; 152,688 pairs), cell counts, source hashes, schemas, privacy
state, and all scientific firewalls passed. No unit identifier, path, scientific
value, label, or outcome was published. No null, localization, or H2 job ran.

Recovery changes the one aggregate field to the length of the prospectively
selected gene-only group axis and adds an explicit 26-group validation. The
local v3 root is preserved as rejected evidence and must not be committed or
reused. A distinct v4 audit is required after this record and correction commit.
