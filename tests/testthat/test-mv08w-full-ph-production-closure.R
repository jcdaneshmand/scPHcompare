test_that("MV8-W independently closes the exact 1,280-record PH inventory", {
  root <- file.path(
    "..", "..", "docs", "audits", "mv08w-full-ph-production-closure-v1"
  )
  expect_true(dir.exists(root))
  read <- function(name) utils::read.csv(
    file.path(root, name), check.names = FALSE, stringsAsFactors = FALSE
  )

  checks <- read("mv08w-validation.csv")
  decision <- read("mv08w-decision.csv")
  inventory <- read("mv08w-full-ph-inventory.csv")
  sources <- read("mv08w-source-cache-rehash.csv")
  production <- read("mv08w-production-ph-rehash.csv")
  resources <- read("mv08w-resource-ledger.csv")
  manifest <- read("mv08w-artifact-manifest.csv")

  expect_equal(nrow(checks), 20L)
  expect_true(all(checks$passed))
  expect_equal(decision$mv08t_records, 23L)
  expect_equal(decision$mv08v_records, 1257L)
  expect_equal(decision$full_ph_records, 1280L)
  expect_equal(decision$fallback_records, 0L)
  expect_equal(decision$validations_passed, decision$validations_total)

  expect_equal(nrow(inventory), 1280L)
  expect_true(all(inventory$independently_validated))
  expect_true(all(inventory$h0_mst_passed))
  expect_equal(nrow(sources), 131L)
  expect_true(all(sources$independently_rehashed))
  expect_equal(nrow(production), 1257L)
  expect_identical(as.integer(production$production_order), 1:1257)
  expect_true(all(production$independently_reconstructed))
  expect_true(all(production$h0_mst_passed))
  expect_equal(sum(production$downstream_jobs), 0)

  expect_equal(nrow(resources), 1257L)
  expect_true(all(resources$disposition == "completed"))
  expect_true(all(resources$workers == 1L))
  expect_true(all(resources$retries == 0L))
  expect_true(all(resources$elapsed_seconds <= resources$elapsed_cap_seconds))
  expect_true(all(
    resources$peak_process_tree_rss_bytes <= resources$rss_cap_bytes
  ))
  expect_true(all(resources$stderr_bytes == 0L))

  firewall_columns <- c(
    "landscape_groups_authorized", "comparison_strata_authorized",
    "clustering_jobs_authorized", "fusion_jobs_authorized",
    "label_jobs_authorized", "outcome_jobs_authorized"
  )
  expect_true(all(as.numeric(unlist(
    decision[firewall_columns], use.names = FALSE
  )) == 0))
  expect_identical(decision$outcome_label_state, "closed")
  expect_false(decision$biological_outcomes_computed)
  expect_true(all(inventory$outcome_label_state == "closed"))
  expect_false(any(inventory$biological_outcomes_computed))

  expect_equal(nrow(manifest), 7L)
  observed <- vapply(file.path(root, manifest$artifact), function(path) {
    expect_equal(as.numeric(file.info(path)$size), as.numeric(
      manifest$bytes[manifest$artifact == basename(path)]
    ))
    digest::digest(file = path, algo = "sha256", serialize = FALSE)
  }, character(1L))
  expect_identical(unname(observed), manifest$sha256)

  public_text <- paste(vapply(
    file.path(root, manifest$artifact),
    function(path) paste(readLines(path, warn = FALSE), collapse = "\n"),
    character(1L)
  ), collapse = "\n")
  expect_false(grepl(
    "/mnt/|[A-Za-z]:\\\\|https?://[^[:space:]]*[?&](token|signature)=",
    public_text, perl = TRUE, ignore.case = TRUE
  ))
})
