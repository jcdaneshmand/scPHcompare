# MV8-ZT engine-v2 interruption-recovery prefreeze

- Validation: 28/28 pass.
- Durable public prefix: 325 chunks and 80,010 pairs.
- One exact engine-v2 child completed at order 326 before receipt promotion.
- Recovery rehashes and adopts order 326 without recomputation, retry, deletion, or overwrite.
- Missing parent telemetry is represented by the frozen 3,600-s/4-GiB caps as conservative upper bounds, not measurements.
- After adoption, strict-prefix production resumes at order 327; all downstream stages remain closed.
