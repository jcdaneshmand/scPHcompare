test_that("MV8-ZU independently closes the complete engine-v2 production contract", {
  root <- testthat::test_path("..", "..")
  builder <- file.path(
    root, "scripts", "build_mv08zu_engine_v2_full_landscape_closure.R"
  )
  launcher <- file.path(
    root, "scripts", "run_mv08zu_engine_v2_full_landscape_closure.R"
  )
  prefreeze_v3 <- file.path(
    root, "docs", "audits", "mv08zt-engine-v2-full-landscape-prefreeze-v3"
  )
  recovery_v2 <- file.path(
    root, "docs", "audits", "mv08zt-engine-v2-interruption-recovery-prefreeze-v2"
  )

  expect_true(file.exists(builder))
  expect_true(file.exists(launcher))
  expect_silent(parse(builder))
  expect_silent(parse(launcher))
  text <- paste(readLines(builder, warn = FALSE), collapse = "\n")
  expect_match(text, "if (length(args) != 10L)", fixed = TRUE)
  expect_match(text, "required_engine_version", fixed = TRUE)
  expect_match(text, "distances$rust_engine_version) == 2L", fixed = TRUE)
  expect_match(text, "status$scientific_engine_version == 2L", fixed = TRUE)
  expect_match(text, "rust_scph_landscape_kernel_v2", fixed = TRUE)
  expect_match(text, "recovery_provenance_ok", fixed = TRUE)
  expect_match(text, "orphan_production_order == 326L", fixed = TRUE)
  expect_match(text, "resume_at_production_order == 327L", fixed = TRUE)
  expect_match(text, "fresh_engine_v2_from_zero_no_old_output_reuse", fixed = TRUE)
  expect_match(text, "total_rows == 152744L", fixed = TRUE)
  expect_match(text, "all-active-level", ignore.case = TRUE)
  expect_match(text, "exact streamed", ignore.case = TRUE)
  expect_match(text, "mv08zu-implementation-binding.csv", fixed = TRUE)
  expect_match(text, "hash_bound_terminal_closure_launcher", fixed = TRUE)
  expect_match(text, "public_privacy_ok", fixed = TRUE)
  expect_match(text, "queue_firewall_ok", fixed = TRUE)
  expect_match(text, "completion_contract_ok", fixed = TRUE)
  expect_false(grepl("processx::process|process\\$new|system2\\(", text, perl = TRUE))
  expect_false(grepl("run_mv08z_landscape_chunk", text, fixed = TRUE))

  launcher_lines <- readLines(launcher, warn = FALSE)
  launcher_text <- paste(launcher_lines, collapse = "\n")
  expect_match(launcher_text, "if (length(args) != 0L)", fixed = TRUE)
  expect_match(
    launcher_text, "landscape_production_complete_closure_pending", fixed = TRUE
  )
  expect_match(launcher_text, "progress$completed_chunks == 628L", fixed = TRUE)
  expect_match(launcher_text, "progress$completed_pairs == 152744L", fixed = TRUE)
  expect_match(launcher_text, "nrow(ledger) == 628L", fixed = TRUE)
  expect_match(launcher_text, "nrow(completed) == 628L", fixed = TRUE)
  expect_match(launcher_text, "refusing to overwrite MV8-ZU output", fixed = TRUE)
  expect_match(launcher_text, "system2(", fixed = TRUE)
  expect_match(launcher_text, "expected_builder_sha256", fixed = TRUE)
  expect_false(grepl("run_mv08z_landscape_chunk", launcher_text, fixed = TRUE))
  hash_line <- grep(
    "^expected_builder_sha256 <-", launcher_lines, value = TRUE
  )
  expect_length(hash_line, 1L)
  expected_hash <- sub(".*\"([0-9a-f]+)\".*", "\\1", hash_line)
  expect_identical(
    expected_hash,
    digest::digest(file = builder, algo = "sha256", serialize = FALSE)
  )

  read <- function(path) utils::read.csv(
    path, stringsAsFactors = FALSE, check.names = FALSE
  )
  closure <- read(file.path(prefreeze_v3, "mv08zu-prospective-closure.csv"))
  expect_equal(closure$required_chunks, 628L)
  expect_equal(closure$required_pairs, 152744L)
  expect_equal(closure$required_engine_version, 2L)
  expect_true(closure$require_every_input_rehash)
  expect_true(closure$require_every_pair_identity)
  expect_true(closure$require_every_active_depth)
  expect_true(closure$require_h0_h1_separate)
  expect_true(closure$require_zero_partial)
  expect_true(closure$require_one_worker)
  expect_true(closure$require_zero_retries)
  expect_true(closure$require_resource_caps)
  expect_true(closure$require_aggregate_privacy)
  expect_equal(closure$comparison_jobs_authorized, 0L)
  expect_equal(closure$clustering_jobs_authorized, 0L)
  expect_equal(closure$fusion_jobs_authorized, 0L)
  expect_equal(closure$label_jobs_authorized, 0L)
  expect_equal(closure$outcome_jobs_authorized, 0L)
  expect_equal(closure$manuscript_claim_jobs_authorized, 0L)

  recovery_checks <- read(file.path(recovery_v2, "mv08zt-r1-validation.csv"))
  recovery_policy <- read(file.path(recovery_v2, "mv08zt-r1-recovery-policy.csv"))
  expect_true(all(recovery_checks$passed))
  expect_equal(recovery_policy$accepted_prefix, 325L)
  expect_equal(recovery_policy$orphan_order, 326L)
  expect_equal(recovery_policy$resume_at_order, 327L)
  expect_false(recovery_policy$landscape_recomputation)
  expect_false(recovery_policy$automatic_retry)
})
