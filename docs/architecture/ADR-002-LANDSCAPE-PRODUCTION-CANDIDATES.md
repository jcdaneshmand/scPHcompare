# ADR-002: Landscape production candidates and diagram eligibility

## Status

Accepted on 2026-08-05 as a profiling and invalidation decision. No production landscape engine or corrected PH pipeline is activated.

## Context

ADR-001 established `landscape_reference_v1` as a non-default exact/adaptive R oracle and required independent validation, diagram-eligibility review, and production-scale profiling. This sprint completed those gates using an isolated Python environment and the actual current and historical PH code.

The dissertation intended cells to be the observations. Both code lineages pass feature-by-cell assay matrices directly to ripserr, whose point-cloud contract interprets rows as points. Consequently the historical diagrams encode features as points and cannot seed a corrected biological analysis.

## Decision

1. Classify all nine audited historical diagram artifacts as scientifically ineligible due to observation-orientation conflict. Permit them only for clearly labeled stress testing.
2. Retain `landscape_reference_v1` as the package correctness oracle.
3. Reject Persim 0.3.8's built-in exact `p_norm()` for project distances because it fails analytical and independent-quadrature sign-crossing cases.
4. Retain Persim's exact critical-pair construction plus a project-controlled exact linear-segment integral as the leading production prototype, subject to eligible-diagram, worst-case, packaging, and batch-interface validation.
5. Do not use the Python GUDHI grid vectorizer for the primary distance. It is a display/sensitivity candidate.
6. Keep Rust gated. The current evidence identifies an algorithmically promising mature-library-derived route before a new language implementation is justified.

## Activation requirements

- a tested cell-as-observation PH input contract and explicit coordinate space;
- a shared or otherwise defensible filtration/censoring policy;
- freshly generated eligible H0 and H1 diagrams with stable IDs and complete provenance;
- exact agreement with the frozen analytical corpus and the R oracle;
- batch pairwise profiling on eligible diagrams, including worst-case critical-pair growth, runtime, peak memory, determinism, failures, and cross-platform packaging;
- author-team approval and a clean, separately identified scientific rerun.

## Consequences

- Legacy landscape, bottleneck, spectral, clustering, statistics, and figures derived from the invalid diagrams cannot be reused as corrected results.
- No automatic transpose is introduced in this sprint; doing so without fixing the coordinate and computational contract could make PH infeasible or scientifically ambiguous.
- Existing data remain the priority. New data are not needed to resolve this implementation defect.
- The next sprint moves upstream from landscape optimization to corrected diagram generation, then returns to production-engine profiling.

Detailed evidence is in `docs/audits/LANDSCAPE_ORACLE_AND_DIAGRAM_ELIGIBILITY_2026-08-05.md`.
