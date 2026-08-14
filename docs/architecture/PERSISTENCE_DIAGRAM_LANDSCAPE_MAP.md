# Persistence Diagram and Landscape Map

## Document control

| Field | Value |
|---|---|
| Status | Conceptual draft; corrected landscape pipeline not yet activated |
| Created | 2026-08-05 |
| Intended uses | Repository documentation, methods development, and eventual manuscript figure source |
| Governing specification | `docs/specifications/PERSISTENCE_LANDSCAPE_SPECIFICATION_V1.md` |
| Required revision gate | Exact/error-controlled landscape engine plus observation-unit and filtration eligibility audits |

This Mermaid diagram is the editable source of truth for the project’s persistence-diagram and persistence-landscape workflow. It distinguishes preprocessing representations, homology dimensions, diagram-derived summaries, distance matrices, and downstream clustering. It describes the intended conceptual architecture; it does not assert that the corrected scientific pipeline has already been activated or rerun.

```mermaid
flowchart TD
    A["One biological sample"] --> B1["Raw point cloud"]
    A --> B2["SCT Individual point cloud"]
    A --> B3["SCT Whole point cloud"]
    A --> B4["Integrated point cloud"]

    B1 --> C["Vietoris–Rips persistent homology"]
    B2 --> C
    B3 --> C
    B4 --> C

    C --> D0["H0 persistence diagram"]
    C --> D1["H1 persistence diagram"]

    D0 --> E0["H0 persistence landscape"]
    D1 --> E1["H1 persistence landscape"]

    D0 --> F0["H0 Betti curve"]
    D1 --> F1["H1 Betti curve"]
    F0 --> FE["Euler curve: Betti0 − Betti1"]
    F1 --> FE

    D0 --> G0["H0 bottleneck distance"]
    D1 --> G1["H1 bottleneck distance"]
    E0 --> H0["H0 landscape L2 distance"]
    E1 --> H1["H1 landscape L2 distance"]

    G0 --> I["BDM and spectral transformation"]
    G1 --> I
    H0 --> J["Separate H0 and H1 LDMs"]
    H1 --> J
    J --> K["Secondary combined LDM"]

    I --> L["Sample-level clustering"]
    J --> L
    K --> L
```

## Interpretation constraints

- Raw, SCT Individual, SCT Whole, and Integrated are separate comparison strata; their filtration scales are not pooled implicitly.
- H0 and H1 are homology dimensions. Landscape levels `lambda_1`, `lambda_2`, and so on exist separately within each dimension.
- The proposed corrected primary landscape calculation uses every active level with exact or error-controlled integration.
- `k=1` is retained only as a clearly labeled paper-compatibility sensitivity.
- H0 and H1 distances are primary separate outputs. Their unweighted combination is secondary and must report the H1 energy contribution.
- Betti, Euler, bottleneck, spectral, and landscape representations are related diagram-derived summaries; they are not interchangeable.
- Sample clustering is downstream of the distance definition. A change to levels, filtration evaluation, distance, or H0/H1 combination can change the clustering.

## Final-figure checklist

Before treating this as a manuscript figure:

1. Replace “intended conceptual architecture” with the approved implemented contract.
2. Add the exact/adaptive integration and numerical-error provenance node.
3. Confirm the eligible persistence-diagram input unit and filtration policy.
4. Decide whether BDM/SDM remain in the revised primary or sensitivity analyses.
5. Add versioned configuration identifiers to the rendered caption.
6. Render deterministic SVG and PNG outputs from this Mermaid source.
7. Review terminology, accessibility, colors, dimensions, and caption with the author team.
