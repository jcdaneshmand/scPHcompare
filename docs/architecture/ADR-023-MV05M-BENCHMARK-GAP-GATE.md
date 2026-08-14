# ADR-023: Gate label-free clustering before robustness or method expansion

- Status: Accepted
- Date: 2026-08-10
- Decision scope: MV5-M prioritization and MV5-N authorization only

## Context

The biological retrieval axis is complete through MV5-L with null/negative
topology-versus-energy evidence. The project still needs technical-mixing,
clustering, robustness, broader integration, gene, fusion, and external
validation work, but those axes differ substantially in identifiability and
artifact readiness.

Technical mixing is currently not identifiable as a confirmatory held-out
sample endpoint because eligible studies are nested within tissues and LOSO
query pairs exclude same-study comparators. Full clustering is also not ready:
accepted distance artifacts are query-to-training only. Unlike technical
mixing, however, clustering has a coherent leakage-safe path using
training-only partitions and held-out assignment, and its required cell views,
diagrams, query distances, and stable-k helper already exist.

## Decision

1. Use the frozen six-criterion weighted gate, without consulting MV5-L tissue
   heterogeneity, to select one next axis.
2. Select label-free clustering contract/resource gating at score 45; retain
   robustness second at 43.
3. Keep technical mixing validity-blocked until a conditional estimand resolves
   study-within-tissue confounding and same-study comparator availability.
4. Authorize MV5-N only to freeze training-only PAM, label-free stable-k,
   held-out medoid assignment, sensitivity rules, exact pair scope, and bounded
   label-closed resource admission.
5. Require MV5-N to stop before full training-matrix production or any
   biological/technical clustering outcome.

## Consequences

MV5-N can determine whether a fair clustering extension is scientifically and
computationally executable without spending the projected 12.9+ landscape
worker-hours prematurely. It must account for 262,675 missing training pairs
and 525,350 H0/H1 rows per representation.

This decision advances the original sample-clustering aim while preserving the
corrected sample unit and study-held-out design. It does not imply that
clustering will favor PH, integration, any tissue, or any k. Robustness remains
the next eligible axis if clustering fails its contract/resource gate.
