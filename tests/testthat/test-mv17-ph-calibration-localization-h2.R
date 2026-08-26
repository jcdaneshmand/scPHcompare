test_that("MV17-A freezes four separate H0/H1 calibration estimands", {
  estimands <- mv17a_estimand_registry_v1()
  expect_equal(nrow(estimands), 4L)
  expect_setequal(estimands$view, c("cell", "gene"))
  expect_setequal(estimands$homology_dimension, c("H0", "H1"))
  expect_true(all(estimands$primary_object ==
                    "persistence_diagram_dimension_separate"))
  expect_true(all(estimands$essential_h0[estimands$homology_dimension == "H0"] ==
                    "excluded"))
  expect_false(any(estimands$null_family_selected))
  expect_false(any(estimands$summary_selected))
  expect_false(any(estimands$real_calibration_authorized))
})

test_that("MV17-A freezes localization outputs and H2 fixture classes", {
  localization <- mv17a_localization_registry_v1()
  fixtures <- mv17a_h2_fixture_registry_v1()
  expect_equal(nrow(localization), 4L)
  expect_true(all(grepl("aggregate_", localization$public_output)))
  expect_false(any(localization$method_selected))
  expect_false(any(localization$real_localization_authorized))
  expect_setequal(fixtures$fixture_class, c(
    "sphere", "torus", "circle", "gaussian_cloud", "shuffled_sphere",
    "shuffled_torus", "shuffled_circle"
  ))
  expect_equal(fixtures$expected_H2[fixtures$fixture_class == "sphere"],
               "positive")
  expect_equal(fixtures$expected_H2[fixtures$fixture_class == "circle"],
               "negative")
  expect_false(any(fixtures$execution_authorized))
  expect_false(any(fixtures$real_data_H2_authorized))
})

test_that("MV17-A keeps every downstream surface closed", {
  gates <- mv17a_stage_gate_registry_v1()
  firewall <- mv17a_firewall_v1()
  expect_equal(gates$stage, c("MV17-B", "MV17-C", "MV17-D", "MV17-E",
                              "MV17-F", "MV17-G", "MV17-H", "MV17-I"))
  expect_setequal(gates$stage[gates$implementation_eligible_after_MV17A],
                  c("MV17-B", "MV17-D", "MV17-E"))
  expect_false(any(gates$scientific_execution_authorized))
  expect_true(all(gates$exact_head_prefreeze_required))
  expect_equal(nrow(firewall), 12L)
  expect_true(all(firewall$state == "closed"))
  expect_false(any(firewall$authorized))
  expect_true(mv17a_validate_contract_v1(
    mv17a_estimand_registry_v1(), mv17a_localization_registry_v1(),
    mv17a_h2_fixture_registry_v1(), gates, firewall
  ))
})

test_that("MV17-A artifact-set identity is order-independent and strict", {
  values <- c(
    paste(rep("a", 64L), collapse = ""),
    paste(rep("b", 64L), collapse = "")
  )
  expect_identical(.mv17a_set_sha256(values),
                   .mv17a_set_sha256(rev(values)))
  expect_error(.mv17a_set_sha256("not-a-hash"), "invalid")
})

test_that("MV17-A builder parses", {
  root <- testthat::test_path("..", "..")
  expect_silent(parse(file.path(
    root, "scripts", "build_mv17a_source_inventory_prefreeze.R"
  )))
})
