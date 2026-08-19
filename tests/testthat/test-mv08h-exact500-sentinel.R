test_that("MV8-H exact500 sentinel feasibility is complete and firewalled", {
  audit <- file.path("..", "..", "docs", "audits", "mv08h-exact500-sentinel-v1")
  validation <- read.csv(file.path(audit, "mv08h-exact500-sentinel-validation.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  manifest <- read.csv(file.path(audit, "mv08h-exact500-sentinel-artifact-manifest.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  expect_equal(nrow(validation), 10L)
  expect_true(all(validation$passed %in% c(TRUE, "TRUE", "True")))
  expect_true(any(validation$check_id == "exact500_panel_present" & validation$passed %in% c(TRUE, "TRUE", "True")))
  expect_true(any(validation$check_id == "qc_384_eligibility" & validation$passed %in% c(TRUE, "TRUE", "True")))
  expect_true(any(validation$check_id == "downstream_firewall" & validation$passed %in% c(TRUE, "TRUE", "True")))
  expect_equal(nrow(manifest), 3L)
  expect_true(all(vapply(seq_len(nrow(manifest)), function(i) {
    identical(toupper(digest::digest(file = file.path(audit, manifest$artifact[[i]]), algo = "sha256", serialize = FALSE)), toupper(manifest$sha256[[i]]))
  }, logical(1))))
})
