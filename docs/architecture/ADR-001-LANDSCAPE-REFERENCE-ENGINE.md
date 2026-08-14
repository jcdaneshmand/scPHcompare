# ADR-001: Persistence-landscape reference engine

## Status

Accepted for non-default correctness and feasibility work on 2026-08-05. Not approved as the production scientific engine and not activated in the public pipeline.

## Context

The approved target contract, `full_l2_error_controlled_v1`, requires all active consecutive landscape levels, exact or error-controlled L2 integration, separate H0/H1 primary distances, and an explicitly secondary combined distance. The convergence audit found H0 depths up to 34,839 and showed that fixed grids can preserve distance ranks while missing narrow-feature energy. A dense grid-by-level representation is therefore not an acceptable reference architecture.

The implementation survey found three relevant established paths:

| Candidate | Strength | Constraint for this project | Decision |
|---|---|---|---|
| Bubenik-Dlotko persistence-landscape toolbox | Published exact and grid algorithms for landscapes and distances | Older standalone C++ integration surface; exact representations may grow quadratically | Retain as algorithmic reference, not an immediate package dependency |
| GUDHI `Persistence_landscape` | Maintained C++ exact representation with distances and other operations | Documentation warns exact representation can be quadratic and memory-prohibitive for large diagrams; R integration would be additional infrastructure | Benchmark before any custom compiled production engine |
| Persim `PersLandscapeExact` | Readable exact critical-pair representation, exact arithmetic operations and p-norm interface | Python dependency is absent from the locked WSL environment; exact representation still has large-diagram scaling limits | Preferred optional independent tractable-case oracle after an approved dependency trial |
| Current R `TDA::landscape` | Already locked and used by scPHcompare | Evaluates supplied levels on a supplied grid; it is not an exact all-level distance engine | Use as an independent pointwise-value cross-check only |

Primary references and official documentation:

- Bubenik, *Statistical Topological Data Analysis using Persistence Landscapes*, JMLR 2015: <https://jmlr.csail.mit.edu/papers/volume16/bubenik15a/bubenik15a.pdf>
- Bubenik and Dlotko, *A persistence landscapes toolbox for topological statistics*: <https://arxiv.org/abs/1501.00179>
- GUDHI persistence-representation documentation: <https://gudhi.inria.fr/doc/latest/group___persistence__representations.html>
- Persim exact-landscape implementation documentation, version 0.3.8 at review time: <https://persim.scikit-tda.org/en/latest/_modules/persim/landscapes/exact.html>

## Decision

Implement a transparent, non-exported R reference interface with two paths:

1. `exact_breakpoint_stream_v1` constructs all tent breakpoints and pairwise crossings for each tractable diagram, then integrates the squared difference of consecutive landscape levels exactly segment by segment. It streams the current level vector and does not serialize a dense node-by-level landscape matrix.
2. `adaptive_quadpack_partitioned_v1` partitions at every finite birth, midpoint, and death, integrates the all-level squared difference on each partition, repeats at fourfold tighter tolerance, and records both QUADPACK's absolute-error estimate and the refinement delta.

`auto` uses the exact path when neither diagram exceeds the provisional 200-finite-interval guard in a dimension and otherwise uses adaptive quadrature. The guard controls an implementation with quadratic breakpoint generation; it is not a scientific landscape-level cap. Both paths exclude invalid, zero-persistence, and infinite intervals before evaluation and report the resulting finite counts.

The result contract includes separate H0 and H1 distances, the secondary Euclidean combined distance, H1 squared-distance contribution, per-dimension numerical diagnostics, and versioned provenance with an explicit `activated_as_scientific_default = FALSE` marker.

## Validation decision

The exact path passes analytical energy, translation, scaling, narrow-feature, overlap, depth, H0/H1, empty/infinite, symmetry, and deterministic-repeat fixtures. Its streamed values agree with the established `TDA::landscape` evaluator on a tractable all-level fixture. The adaptive path agrees with the exact path throughout the benchmark corpus; maximum observed absolute combined-distance difference was approximately `1.12e-18`, and every repeated result was identical.

This is sufficient for a package reference engine and correctness corpus. ADR-002 subsequently records the independent Persim/GUDHI/SciPy cross-check, the Persim 0.3.8 sign-crossing defect, production-scale profiling, and the historical-diagram orientation invalidation. Those findings retain this implementation as the R correctness oracle but do not activate adaptive results as primary scientific output.

## Rust gate

Rust is not approved now. The sprint identifies algorithmic scaling, not language overhead, as the first issue. Mature C++ and Python exact implementations exist, the new R reference corpus is fast for tractable inputs, and no end-to-end production benchmark yet demonstrates that a Rust component is the dominant solution. Before reconsideration:

1. generate eligible cell-oriented diagrams under an approved coordinate and filtration contract;
2. benchmark the corrected exact critical-pair and direct-distance strategies on those eligible diagrams;
3. replace per-evaluation sorting and repeated adaptive scans algorithmically;
4. measure the corrected end-to-end distance stage, including peak memory;
5. demonstrate that an isolated compiled component materially improves the pipeline while passing the frozen corpus on all supported operating systems.

## Consequences

- The existing public K1 behavior remains unchanged.
- No biological result or headline figure is regenerated.
- The reference implementation can now serve as an oracle for optimized implementations.
- Large-diagram production feasibility remains open; adaptive error certification and independent exact-distance validation are explicit activation gates.
