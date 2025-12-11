<!-- .github/copilot-instructions.md for scPHcompare -->
# Copilot / AI Agent Instructions for scPHcompare

These instructions are focused, actionable, and specific to this repository so an AI coding agent can be productive quickly.

1. Purpose & Big Picture
- **This repo**: an R package (`scPHcompare`) that orchestrates persistent homology (PH) calculations and downstream cluster/curve comparisons across many single-cell datasets. The primary entrypoint is `run_unified_pipeline()` in `R/unified_pipeline.R`.
- **Major components**:
  - `R/`: package code (PH computation, postprocessing, helpers).
  - `inst/extdata/`: configuration and templates (`iteration_config.csv`, `custom_iteration_inputs_template.R`).
  - `scripts/` and `tests/`: reproducible examples and unit tests.
  - `Result Tables/`: archived analysis tables (do not modify lightly).

2. Key entry points & examples
- Run the full pipeline from R: `run_unified_pipeline(metadata_path, results_dir, num_cores = 4, integration_methods = "seurat")` (see `R/unified_pipeline.R`).
- Recompute only matrices for a specific iteration: call `process_iteration_calculate_matrices()` in `R/PH_PostProcessing_andAnalysis.R`.
- Use toy data generator for quick experiments: `generate_toy_data()` (files written to `inst/extdata/toy/`).

3. Important files & conventions to reference in code changes
- `inst/extdata/iteration_config.csv`: canonical iteration labels and output name prefixes (e.g. `seurat_integration`). Use these labels when naming iteration outputs.
- `inst/extdata/custom_iteration_inputs_template.R`: template for overriding Seurat objects / precomputed PH/matrices per iteration. Agents should prefer using this template when adding custom input paths.
- `R/apply_all_clustering_methods.R` and `R/PH_PostProcessing_andAnalysis.R`: define clustering and downstream reporting. Column naming conventions for Seurat metadata are strict (see README): e.g. `kmeans_cluster_<dataset>_<group>` and `hierarchical_cluster_<matrix>_ph_<dataset>_<group>_<linkage>`.
- `R/PH_Calculation.R` and `R/PH_Functions.R`: core PH calculations and data structures (persistence diagram lists, landscape lists, BDM/SDM matrices).

4. Caching, long-running jobs and resume behavior
- Bottleneck distance matrix computation writes `BDM_progress_*` files under your `results_dir`. Code resumes from these progress files — preserve them between runs if retrying.
- Betti curve computations use per-dataset caches under `plots/betti_plots/betti_cache/` (hashed filenames). When modifying Betti code, update cache invalidation logic carefully.

5. Dependency & environment workflow
- Uses `renv/` to pin R package versions. Typical setup:
  - In R: `install.packages("renv"); renv::restore()`
  - To run tests locally: `devtools::test()` or `Rscript -e "devtools::test()"` from the project root.
- Prefer running R sessions from the repository root so `renv` and relative paths resolve.

6. Tests and quick verification
- Tests live in `tests/testthat/`. To run a focused test file: `Rscript -e "testthat::test_file('tests/testthat/test-postprocessing_modules.R')"`.
- For broader package checks use: `R CMD check --no-manual .` (run in an Rtools-enabled shell on Windows).

7. Code patterns and style to follow
- Functional R style across `R/` files: small exported wrappers in `R/` with helpers in separate files (e.g. `ph_utils.R`, `betti_utils.R`). Keep new helpers in `R/` and export only user-facing functions via `NAMESPACE` when necessary.
- File naming: saved Seurat objects use `<iteration>_seurat_object.rds` in `results/seurat_objects/`. Respect existing filename prefixes from `iteration_config.csv`.

8. Integration points and external systems
- Seurat and other Bioconductor/CRAN packages are required (see `DESCRIPTION`). Agents should not hardcode alternative package versions — use `renv` for changes.
- Large computations expect adequate memory/cores — examples in README use `num_cores` set when calling `run_unified_pipeline()`.

9. When changing behavior, update these places
- `inst/extdata/iteration_config.csv` when adding/removing iterations.
- `inst/extdata/custom_iteration_inputs_template.R` when exposing new override keys.
- README and `man/` Rd files for public API changes.

10. Safety and repository boundaries
- Do not modify files in `Result Tables/` unless explicitly requested — they are archival outputs used in manuscript analysis.

If anything in here is unclear or you want more details (examples of tests to run, exact lines to update for cache invalidation, or a small PR to demonstrate the pattern), tell me which part to expand and I will iterate.
