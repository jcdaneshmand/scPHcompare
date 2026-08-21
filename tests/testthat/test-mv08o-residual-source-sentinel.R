test_that("MV8-O closes the bounded residual source sentinel without opening topology", {
  audit <- file.path("..", "..", "docs", "audits", "mv08o-residual-source-sentinel-closure-v1")
  report_path <- file.path(audit, "MV08O_RESIDUAL_SOURCE_SENTINEL_CLOSURE_2026-08-21.md")
  builder_path <- file.path("..", "..", "scripts", "build_mv08o_residual_source_sentinel_closure.R")
  worker_path <- file.path("..", "..", "scripts", "run_mv08o_residual_source_worker.R")
  sentinel_path <- file.path("..", "..", "scripts", "run_mv08o_residual_source_sentinel.R")
  expect_true(all(file.exists(c(report_path, builder_path, worker_path, sentinel_path))))

  read_audit <- function(name) {
    read.csv(file.path(audit, name), check.names = FALSE, stringsAsFactors = FALSE)
  }
  resource <- read_audit("mv08o-source-sentinel-resource.csv")
  views <- read_audit("mv08o-source-sentinel-view-summary.csv")
  determinism <- read_audit("mv08o-source-sentinel-determinism.csv")
  validation <- read_audit("mv08o-source-sentinel-validation.csv")
  manifest <- read_audit("mv08o-artifact-manifest.csv")

  expect_equal(nrow(resource), 5L)
  expect_identical(as.character(resource$unit_id), c(
    "SRA701877_SRS3279688", "SRA742961_SRS3565197", "HCA_BM_002",
    "SRA742961_SRS3565197", "HCA_BM_002"
  ))
  expect_identical(as.character(resource$run_role), c("primary", "primary", "primary", "repeat", "repeat"))
  expect_true(all(resource$disposition == "completed" & resource$exit_status == 0L))
  expect_true(all(resource$elapsed_seconds <= resource$elapsed_cap_seconds))
  expect_true(all(resource$peak_rss_bytes <= resource$rss_cap_bytes))
  expect_true(all(resource$stderr_class == "known_glmGamPoi_native_fallback"))
  expect_true(all(resource$worker_rows == ifelse(resource$dataset_scope == "internal124", 10L, 4L)))
  expect_true(all(resource$all_geometry_valid))
  expect_false(any(resource$persistence_computed | resource$landscapes_computed |
    resource$clustering_computed | resource$fusion_computed | resource$biological_outcomes_computed))

  expect_equal(nrow(views), 38L)
  required <- views$dataset_scope == "internal124" | views$panel_id == "common475" |
    (views$panel_id == "exact500" & views$representation_id == "sct_pearson_residual_all_qc_fit_selected384")
  diagnostic <- views$dataset_scope == "external8" & views$panel_id == "exact500" &
    views$representation_id == "sct_data_all_qc_fit_selected384"
  expect_true(all(views$values_finite[required]))
  expect_true(all(views$zero_variance_gene_count[required] == 0L))
  expect_true(all(views$correlation_chord_valid[required]))
  expect_equal(sum(diagnostic), 2L)
  expect_true(all(!views$correlation_chord_valid[diagnostic]))
  expect_false(any(views$persistence_computed | views$landscapes_computed |
    views$clustering_computed | views$fusion_computed | views$biological_outcomes_computed))

  expect_equal(nrow(determinism), 14L)
  expect_true(all(determinism$passed))
  expect_equal(nrow(validation), 11L)
  expect_true(all(validation$passed))
  expect_equal(nrow(manifest), 5L)
  expect_true(all(file.exists(file.path(audit, manifest$artifact))))
  live_hashes <- vapply(file.path(audit, manifest$artifact), function(path) {
    digest::digest(file = path, algo = "sha256", serialize = FALSE)
  }, character(1L))
  expect_identical(unname(live_hashes), manifest$sha256)

  public_csv <- list.files(audit, pattern = "\\.csv$", full.names = TRUE)
  public_text <- paste(unlist(lapply(public_csv, readLines, warn = FALSE)), collapse = "\n")
  expect_false(grepl("(?:^|[,\"'])tmp[/\\\\]|/mnt/[a-z]/|[A-Za-z]:[/\\\\]", public_text, perl = TRUE))
  report <- paste(readLines(report_path, warn = FALSE), collapse = "\n")
  expect_match(report, "diagnostic-only invalid view", fixed = TRUE)
  expect_match(report, "separately prefrozen full-source production decision", fixed = TRUE)
})
