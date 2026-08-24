test_that("MV8-ZA binds only live resource monitoring around unchanged workers", {
  root <- file.path("..", "..", "docs", "audits",
                    "mv08za-landscape-sentinel-monitor-prefreeze-v1")
  read <- function(name) read.csv(file.path(root, name), check.names = FALSE,
                                  stringsAsFactors = FALSE)
  binding <- read("mv08za-implementation-bindings.csv")
  checks <- read("mv08za-validation.csv")
  decision <- read("mv08za-decision.csv")
  manifest <- read("mv08za-artifact-manifest.csv")
  expect_equal(nrow(binding), 4L)
  expect_true(all(!binding$scientific_change))
  expect_equal(nrow(checks), 12L)
  expect_true(all(checks$passed))
  expect_equal(decision$authorized_child_processes, 3L)
  expect_equal(decision$Rust_chunks, 2L)
  expect_equal(decision$pairs_per_Rust_chunk, 250L)
  expect_equal(decision$R_oracle_pairs, 1L)
  expect_equal(decision$workers, 1L)
  expect_equal(decision$retries, 0L)
  expect_false(decision$scientific_contract_changed)
  expect_true(all(unlist(decision[c(
    "production_pairs_authorized", "comparison_jobs_authorized",
    "clustering_jobs_authorized", "fusion_jobs_authorized",
    "label_jobs_authorized", "outcome_jobs_authorized"
  )], use.names = FALSE) == 0))
  expect_identical(
    unname(vapply(file.path(root, manifest$artifact), function(path) {
      digest::digest(file = path, algo = "sha256", serialize = FALSE)
    }, character(1L))), manifest$sha256
  )
  monitor <- paste(readLines(file.path("..", "..", "scripts",
                                       "run_mv08za_landscape_sentinel.R")),
                   collapse = "\n")
  expect_match(monitor, "ps_children", fixed = TRUE)
  expect_match(monitor, "kill_tree", fixed = TRUE)
  expect_match(monitor, "MV08ZA_GIT_HEAD", fixed = TRUE)
})
