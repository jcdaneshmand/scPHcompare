# Phase 3 Subplan — Literature, References, and Figures

## Objective

Create a source-verifiable scholarly and visual record. This phase may run alongside early technical work, but G3 depends on stable claims and figures.

## Task ledger

| ID | Task | Required evidence/output | Acceptance test | Status |
|---|---|---|---|---|
| P3-01 | Define literature search protocol | Databases, queries, dates, inclusion/exclusion criteria, screening fields | Search is repeatable and covers TDA, integration, benchmarking, and distortion | `not_started` |
| P3-02 | Build prior-work evidence table | Persistent identifiers and extracted problem/unit/method/data/validation/contribution fields | Each included work has an authoritative record and relevance disposition | `in_progress` |
| P3-03 | Build novelty matrix | Comparison across problem, unit, topology, integration, validation, interpretation, software | Proposed novelty has direct supporting contrasts and explicit limits | `not_started` |
| P3-04 | Verify bibliography metadata | DOI/PMID/PMCID/title/authors/venue/year/status/corrections ledger | Every reference resolves to authoritative metadata or is flagged | `not_started` |
| P3-05 | Audit sentence-level support | Claim, citation, quoted evidence location, support strength, action | Every cited factual claim is supported by the cited source | `not_started` |
| P3-06 | Inventory figures | Panel-level provenance, underlying data, script, manual/AI involvement, rights/status | Every panel is reproducible, replaceable, or explicitly marked historical | `in_progress` |
| P3-07 | Regenerate quantitative figures | Versioned tidy data, deterministic script, output hash, visual QA record | Values match analysis outputs and captions; no manual numeric placement | `not_started` |
| P3-08 | Redesign schematic/complex figures | Deterministic source and accessibility checks | Visual meaning is accurate, legible, and regeneration is documented | `not_started` |
| P3-09 | Revise contribution statement | Evidence-bounded contribution and non-claims | Novelty language matches the matrix and current results | `not_started` |
| P3-10 | Evaluate G3 | Literature/reference/figure audit packet | All major claims/citations/figures pass or have approved exceptions | `not_started` |

## Evidence rules

- Prefer primary papers and authoritative bibliographic records.
- Record searches that produce no relevant evidence.
- Do not cite a general TDA paper for a single-cell or data-growth claim it does not support.
- Replace legacy LLM-generated prose through source-grounded rewriting, not cosmetic paraphrase.

## Gate G3 checklist

- [ ] Search protocol and screening record are reproducible.
- [ ] Every reference has verified metadata and sentence-level support.
- [ ] Novelty statement is defensible against the comparison matrix.
- [ ] Every quantitative figure regenerates from versioned outputs.
- [ ] Manual/AI-assisted figure provenance and replacement decisions are recorded.
