test_that("publication evidence tooling freezes the intended boundary", {
  builder_path <- "../../scripts/build_publication_evidence_bundle.py"
  planner_path <- "../../scripts/plan_publication_stack.py"
  skip_if_not(file.exists(builder_path) && file.exists(planner_path),
              "repository-only publication scripts are excluded from source packages")

  builder <- paste(readLines(builder_path, warn = FALSE), collapse = "\n")
  planner <- paste(readLines(planner_path, warn = FALSE), collapse = "\n")

  expect_match(builder, 'path.startswith("results/")', fixed = TRUE)
  expect_match(builder, 'path.startswith("docs/audits/")', fixed = TRUE)
  expect_match(builder, "GitBlobReader", fixed = TRUE)
  expect_match(builder, "FIXED_ZIP_TIME", fixed = TRUE)
  expect_match(builder, "provisional_author_team_review_required", fixed = TRUE)
  expect_match(builder, "Refusing to overwrite", fixed = TRUE)
  expect_match(planner, "max-files", fixed = TRUE)
  expect_match(planner, "max-lines", fixed = TRUE)
  expect_match(planner, "max-raw-bytes", fixed = TRUE)
  expect_match(planner, "docs/audits/*.csv", fixed = TRUE)
  expect_match(planner, "docs/audits/**/*.csv", fixed = TRUE)
  expect_match(planner, "results/**", fixed = TRUE)
})
