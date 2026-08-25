test_that("MV8-X implementation parses and retains the isolated rebuild boundary", {
  root <- normalizePath(file.path(testthat::test_path(), "..", ".."),
                        winslash = "/", mustWork = TRUE)
  r_scripts <- file.path(root, "scripts", c(
    "build_mv08x_rust_landscape_admission_prefreeze.R",
    "run_mv08x_rust_landscape_oracles.R",
    "build_mv08y_rust_landscape_admission_closure.R",
    "run_mv08xa_rust_invariant_diagnostic.R",
    "build_mv08xa_oracle_observability_recovery.R"
  ))
  for (script in r_scripts) expect_silent(parse(file = script))

  acquire <- paste(readLines(file.path(
    root, "scripts", "acquire_mv08x_rust_toolchain.sh"
  ), warn = FALSE), collapse = "\n")
  rebuild <- paste(readLines(file.path(
    root, "scripts", "run_mv08x_rust_rebuild.sh"
  ), warn = FALSE), collapse = "\n")
  expect_match(acquire, "--no-modify-path", fixed = TRUE)
  expect_match(acquire, "CARGO_HOME", fixed = TRUE)
  expect_match(acquire, "RUSTUP_HOME", fixed = TRUE)
  expect_match(acquire, "x86_64-unknown-linux-gnu", fixed = TRUE)
  expect_match(acquire, "rustc_host=%s", fixed = TRUE)
  expect_false(grepl("sudo|apt-get|dnf|yum", acquire, perl = TRUE))
  expect_match(rebuild, "CARGO_BUILD_JOBS=1", fixed = TRUE)
  expect_match(rebuild, "--locked -j 1", fixed = TRUE)
  expect_match(rebuild, "-D warnings", fixed = TRUE)
  expect_match(rebuild, "clean Rust builds are not byte-identical", fixed = TRUE)
  expect_false(grepl("sudo|apt-get|dnf|yum", rebuild, perl = TRUE))
  oracle <- paste(readLines(r_scripts[[2L]], warn = FALSE), collapse = "\n")
  expect_match(oracle, "landscape_active_depth", fixed = TRUE)
  expect_match(oracle, "oracle-progress.csv", fixed = TRUE)
  expect_match(oracle, "cumsum(births - deaths)", fixed = TRUE)
})

test_that("MV8-X prospectively freezes one private stress pair per group", {
  root <- file.path(
    "..", "..", "docs", "audits",
    "mv08x-rust-landscape-admission-prefreeze-v1"
  )
  skip_if_not(dir.exists(root), "MV8-X prefreeze has not been produced")
  read <- function(name) utils::read.csv(
    file.path(root, name), check.names = FALSE, stringsAsFactors = FALSE
  )
  contract <- read("mv08x-contract.csv")
  selection <- read("mv08x-oracle-selection.csv")
  resources <- read("mv08x-resource-policy.csv")
  acceptance <- read("mv08x-acceptance.csv")
  implementations <- read("mv08x-implementation-bindings.csv")
  checks <- read("mv08x-validation.csv")
  decision <- read("mv08x-decision.csv")
  manifest <- read("mv08x-artifact-manifest.csv")

  expect_equal(nrow(selection), 28L)
  expect_identical(as.integer(selection$oracle_order), 1:28)
  dimension_counts <- table(selection$homology_dimension)
  expect_identical(
    as.integer(dimension_counts[c("H0", "H1")]),
    c(14L, 14L)
  )
  expect_setequal(selection$dataset_scope, c("external8", "internal124"))
  expect_setequal(selection$panel_id, c("common475", "exact500"))
  expect_setequal(selection$view_kind, c("cell_topology_v1", "gene_topology_v1"))
  expect_equal(length(unique(selection$representation_id)), 3L)
  expect_true(all(selection$reference_route %in% c(
    "r_exact_breakpoint", "r_adaptive_certified"
  )))
  expect_false(any(c(
    "job_id", "unit_id", "output_file", "source_role", "accession",
    "private_path"
  ) %in% names(selection)))
  expect_true(all(resources$workers == 1L))
  expect_true(all(resources$retries == 0L))
  expect_equal(nrow(acceptance), 17L)
  expect_equal(nrow(checks), 20L)
  expect_true(all(checks$passed))
  expect_equal(contract$production_landscape_jobs, 0L)
  expect_equal(decision$oracle_pairs, 28L)
  expect_false(decision$landscape_execution_authorized)
  expect_equal(decision$comparison_jobs_authorized, 0L)
  expect_equal(decision$clustering_jobs_authorized, 0L)
  expect_equal(decision$fusion_jobs_authorized, 0L)
  expect_equal(decision$label_jobs_authorized, 0L)
  expect_equal(decision$outcome_jobs_authorized, 0L)
  expect_identical(decision$outcome_label_state, "closed")
  expect_false(decision$biological_outcomes_computed)

  amendment_root <- file.path(
    "..", "..", "docs", "audits",
    "mv08xa-oracle-observability-recovery-v1"
  )
  amended_roles <- character()
  amendment <- data.frame()
  if (dir.exists(amendment_root)) {
    amendment <- utils::read.csv(
      file.path(amendment_root, "mv08xa-amendment-bindings.csv"),
      check.names = FALSE, stringsAsFactors = FALSE
    )
    amended_roles <- amendment$role[
      !is.na(amendment$old_sha256) & nzchar(amendment$old_sha256)
    ]
    expected_amendment <- amendment$new_sha256
    zp_path <- file.path(
      "..", "..", "docs", "audits",
      "mv08zp-landscape-kernel-repair-prefreeze-v1",
      "mv08zp-implementation-bindings.csv"
    )
    if (file.exists(zp_path)) {
      zp <- utils::read.csv(zp_path, check.names = FALSE,
                            stringsAsFactors = FALSE)
      replacement <- match(amendment$file, zp$file)
      changed <- !is.na(replacement)
      expected_amendment[changed] <- zp$sha256[replacement[changed]]
    }
    zs_path <- file.path(
      "..", "..", "docs", "audits",
      "mv08zs-landscape-oracle-harness-recovery-acceptance-v1",
      "mv08zs-harness-amendment.csv"
    )
    if (file.exists(zs_path)) {
      zs <- utils::read.csv(zs_path, check.names = FALSE,
                            stringsAsFactors = FALSE)
      replacement <- match(amendment$file, zs$file)
      changed <- !is.na(replacement)
      expected_amendment[changed] <- zs$new_sha256[replacement[changed]]
    }
    current_amendment <- vapply(
      file.path("..", "..", amendment$file),
      function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE),
      character(1L)
    )
    expect_identical(unname(current_amendment), expected_amendment)
  }
  retained <- !implementations$role %in% amended_roles
  expected_retained <- implementations$sha256[retained]
  zp_path <- file.path(
    "..", "..", "docs", "audits",
    "mv08zp-landscape-kernel-repair-prefreeze-v1",
    "mv08zp-implementation-bindings.csv"
  )
  if (file.exists(zp_path)) {
    zp <- utils::read.csv(zp_path, check.names = FALSE,
                          stringsAsFactors = FALSE)
    replacement <- match(implementations$file[retained], zp$file)
    changed <- !is.na(replacement)
    expected_retained[changed] <- zp$sha256[replacement[changed]]
  }
  current <- vapply(file.path("..", "..", implementations$file[retained]), function(path) {
    digest::digest(file = path, algo = "sha256", serialize = FALSE)
  }, character(1L))
  expect_identical(unname(current), expected_retained)
  observed <- vapply(file.path(root, manifest$artifact), function(path) {
    digest::digest(file = path, algo = "sha256", serialize = FALSE)
  }, character(1L))
  expect_identical(unname(observed), manifest$sha256)

  public_text <- paste(vapply(
    file.path(root, manifest$artifact),
    function(path) paste(readLines(path, warn = FALSE), collapse = "\n"),
    character(1L)
  ), collapse = "\n")
  expect_false(grepl(
    "SRA[0-9]|HCA_BM|/mnt/|[A-Za-z]:\\\\|output_file|unit_id|job_id",
    public_text, perl = TRUE
  ))
})
