test_that("MV13-C independently refits every sentinel stage", {
  root <- testthat::test_path("..", "..")
  script <- paste(readLines(file.path(root, "scripts",
    "build_mv13c_allqc_cell_sentinel_closure.R"), warn = FALSE),
    collapse = "\n")
  expect_match(script, "length(sources) == 124L", fixed = TRUE)
  expect_match(script, "fit_cell_topology_pca", fixed = TRUE)
  expect_match(script, "construct_cell_topology_view", fixed = TRUE)
  expect_match(script, "run_topology_view_ph", fixed = TRUE)
  expect_match(script, "mv07g_validate_ph_against_view_v1", fixed = TRUE)
  expect_match(script, "PCA_rotation_exact", fixed = TRUE)
  expect_match(script, "PH_diagram_exact", fixed = TRUE)
  expect_match(script, "full_execution_authorized_by_this_closure = FALSE",
               fixed = TRUE)
  expect_match(script, "labels_authorized = FALSE", fixed = TRUE)
})
