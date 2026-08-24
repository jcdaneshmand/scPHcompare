# MV8-ZH launcher-interruption recovery prefreeze

- Validation: 25/25 pass.
- Durable public prefix: 163 chunks; one exact completed orphan at order 164.
- Recovery: rehash and adopt the orphan without recomputation, retry, deletion, or overwrite.
- Missing parent telemetry is represented by the frozen 3,600-s/4-GiB caps as conservative upper bounds, not measurements.
- After adoption, strict-prefix production resumes at order 165; every downstream stage remains closed.
