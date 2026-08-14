#!/usr/bin/env Rscript

options(warn = 2)

if (!requireNamespace("digest", quietly = TRUE)) {
  stop("digest is required for MV-03 identity binding.", call. = FALSE)
}
source(file.path("R", "toy_baseline.R"), local = FALSE)
source(file.path("R", "dual_view_topology.R"), local = FALSE)
source(file.path("R", "mv03_pilot.R"), local = FALSE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8L) {
  stop(
    "Usage: build_mv03_stage_bc_inputs.R <B|C> <large-residual> ",
    "<large-integrated> <bone-residual> <bone-integrated> <manifest-csv> ",
    "<audit-dir> <view-cache-dir>", call. = FALSE
  )
}

stage <- match.arg(args[[1L]], c("B", "C"))
cache_paths <- list(
  large = list(residual = args[[2L]], integrated = args[[3L]]),
  bone = list(residual = args[[4L]], integrated = args[[5L]])
)
manifest_path <- normalizePath(args[[6L]], winslash = "/", mustWork = TRUE)
audit_dir <- args[[7L]]
view_dir <- args[[8L]]
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(view_dir, recursive = TRUE, showWarnings = FALSE)

manifest <- utils::read.csv(
  manifest_path, stringsAsFactors = FALSE, check.names = FALSE
)
fit_scope_id <- "descriptive_all_pilot_samples"
seeds <- if (stage == "B") 20260805L else 20260805:20260809
job_rows <- list()
subsample_rows <- list()
standardization_rows <- list()
pca_rows <- list()

for (cohort in c("large", "bone")) {
  residual_cache <- readRDS(cache_paths[[cohort]]$residual)
  integrated_cache <- readRDS(cache_paths[[cohort]]$integrated)
  if (!identical(residual_cache$cohort, cohort) ||
      !identical(integrated_cache$cohort, cohort)) {
    stop("Cache cohort mismatch for ", cohort, ".", call. = FALSE)
  }
  panel_fit <- fit_mv03_gene_panel(
    residual_cache$samples, integrated_cache$samples,
    cohort = cohort, fit_scope_id = fit_scope_id, panel_size = 500L
  )
  utils::write.csv(
    data.frame(
      order = seq_along(panel_fit$panel), gene = panel_fit$panel,
      cohort = cohort, fit_scope_id = fit_scope_id,
      panel_sha256 = panel_fit$panel_sha256, stringsAsFactors = FALSE
    ),
    file.path(audit_dir, paste0(
      "mv03-gene-panel-", cohort, "-seed-independent-2026-08-05.csv"
    )), row.names = FALSE
  )
  utils::write.csv(
    panel_fit$eligibility,
    file.path(audit_dir, paste0(
      "mv03-gene-eligibility-", cohort, "-2026-08-05.csv"
    )), row.names = FALSE
  )
  utils::write.csv(
    panel_fit$ranks,
    file.path(audit_dir, paste0(
      "mv03-sct-variance-ranks-", cohort, "-2026-08-05.csv"
    )), row.names = FALSE
  )

  targets <- if (stage == "B") {
    if (cohort == "large") {
      manifest$sample_id[manifest$cohort == cohort &
                           manifest$pilot_role == "core_tissue_median"]
    } else {
      manifest$sample_id[manifest$cohort == cohort]
    }
  } else {
    if (cohort == "large") {
      manifest$sample_id[manifest$cohort == cohort &
                           manifest$pilot_role %in%
                           c("stress_cell_minimum", "stress_cell_maximum")]
    } else {
      manifest$sample_id[manifest$cohort == cohort &
                           manifest$pilot_role == "technical_approach_q25"]
    }
  }
  expected_targets <- if (stage == "B") {
    if (cohort == "large") 8L else 4L
  } else {
    2L
  }
  if (length(targets) != expected_targets) {
    stop("Frozen target count mismatch for Stage ", stage, " ", cohort,
         ".", call. = FALSE)
  }

  sample_ids <- sort(names(residual_cache$samples), method = "radix")
  for (seed in seeds) {
    selected_cells <- lapply(sample_ids, function(sample_id) {
      select_matched_cells(
        colnames(residual_cache$samples[[sample_id]]), n = 384L, seed = seed
      )
    })
    names(selected_cells) <- sample_ids
    for (sample_id in sample_ids) {
      ids <- selected_cells[[sample_id]]
      subsample_rows[[length(subsample_rows) + 1L]] <- data.frame(
        stage = stage, cohort = cohort, sample_id = sample_id, seed = seed,
        selected_order = seq_along(ids), cell_id = as.character(ids),
        eligible_cell_count = attr(ids, "eligible_cell_count"),
        selected_cell_sha256 = attr(ids, "selected_cell_sha256"),
        stringsAsFactors = FALSE
      )
    }

    representations <- if (cohort == "large") {
      c("sct_whole", "seurat_integration")
    } else {
      c("sct_whole", "integrated")
    }
    for (representation in representations) {
      values <- if (representation == "sct_whole") {
        panel_fit$residual_samples
      } else {
        panel_fit$integrated_samples
      }
      prepared <- prepare_mv03_sources(
        values, panel_fit$panel, cohort = cohort,
        representation = representation, fit_scope_id = fit_scope_id,
        seed = seed, selected_cells = selected_cells
      )
      pca_model <- fit_cell_topology_pca(
        prepared$sources, n_components = 30L, pca_seed = seed
      )
      standardization_rows[[length(standardization_rows) + 1L]] <- data.frame(
        stage = stage, cohort = cohort, representation = representation,
        seed = seed, gene = panel_fit$panel,
        center = unname(prepared$center), scale = unname(prepared$scale),
        fit_scope_id = fit_scope_id,
        standardization_id = prepared$standardization_id,
        stringsAsFactors = FALSE
      )
      pca_rows[[length(pca_rows) + 1L]] <- data.frame(
        stage = stage, cohort = cohort, representation = representation,
        seed = seed, component = seq_along(pca_model$singular_values),
        singular_value = pca_model$singular_values,
        pca_cache_key = pca_model$cache_key,
        algorithm = pca_model$algorithm,
        fit_sample_count = length(pca_model$fit_sample_ids),
        fit_scope_id = fit_scope_id,
        stringsAsFactors = FALSE
      )
      for (sample_id in sort(targets, method = "radix")) {
        source_object <- prepared$sources[[sample_id]]
        views <- list(
          cell_topology_v1 = construct_cell_topology_view(
            source_object, pca_model
          ),
          gene_topology_v1 = construct_gene_topology_view(source_object)
        )
        for (view_id in names(views)) {
          path <- file.path(
            view_dir,
            paste(stage, cohort, representation, sample_id, seed, view_id,
                  sep = "__")
          )
          path <- paste0(path, ".rds")
          saveRDS(views[[view_id]], path, compress = FALSE)
          job_rows[[length(job_rows) + 1L]] <- data.frame(
            stage = stage, cohort = cohort, representation = representation,
            sample_id = sample_id, seed = seed, view_id = view_id,
            point_count = length(views[[view_id]]$point_ids),
            coordinate_count = length(views[[view_id]]$coordinate_ids),
            view_cache_key = views[[view_id]]$cache_key,
            view_rds = gsub("\\\\", "/", path), stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  rm(residual_cache, integrated_cache, panel_fit)
  invisible(gc())
}

utils::write.csv(
  do.call(rbind, subsample_rows),
  file.path(audit_dir, paste0(
    "mv03-stage-", tolower(stage), "-matched-cells-2026-08-05.csv"
  )), row.names = FALSE
)
utils::write.csv(
  do.call(rbind, standardization_rows),
  file.path(audit_dir, paste0(
    "mv03-stage-", tolower(stage), "-standardization-2026-08-05.csv"
  )), row.names = FALSE
)
utils::write.csv(
  do.call(rbind, pca_rows),
  file.path(audit_dir, paste0(
    "mv03-stage-", tolower(stage), "-pca-summary-2026-08-05.csv"
  )), row.names = FALSE
)
jobs <- do.call(rbind, job_rows)
utils::write.csv(
  jobs,
  file.path(audit_dir, paste0(
    "mv03-stage-", tolower(stage), "-job-manifest-2026-08-05.csv"
  )), row.names = FALSE
)
message("Built ", nrow(jobs), " scientifically eligible Stage ", stage,
        " view jobs.")
