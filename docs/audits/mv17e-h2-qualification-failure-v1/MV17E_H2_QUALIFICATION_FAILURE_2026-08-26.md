# MV17-E H2 qualification failure

Both prospectively frozen seven-fixture runs completed successfully under their
resource caps. Ripserr and TDA/GUDHI agree exactly on every H0/H1/H2 interval,
and all seven scientific RDS payloads repeat byte-for-byte.

The closure stopped because complete Vietoris-Rips H2 is not negative on the
finite circle controls. The 42-point circle has 13 H2 intervals with maximum
persistence 0.0698869; its shuffled counterpart has two H2 intervals with
maximum persistence 0.0650376. Both exceed the frozen `1e-8` maximum. Sphere
and torus positive controls pass, and shuffled sphere/torus attenuate.

No threshold may be changed and production must not be rerun. A committed
failure-closure prefreeze must bind all artifacts and authorize only a negative
closure that blocks real-data H2 and MV17-F while allowing H0/H1 work to proceed.
