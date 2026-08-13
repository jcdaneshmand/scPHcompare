args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: mv05bf_rust_ci_fixture.R <library> <output.csv>",
       call. = FALSE)
}
library <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
output <- args[[2L]]
source("R/landscape_rust_prototype.R")

empty <- matrix(numeric(), nrow = 0L, ncol = 2L)
fixtures <- list(
  single_tent = list(matrix(c(0, 2), nrow = 1L), empty, 2 / 3),
  sign_changing_tents = list(
    matrix(c(0, 2), nrow = 1L), matrix(c(0.25, 2.25), nrow = 1L), 7 / 64
  ),
  narrow_feature = list(
    matrix(c(0.499, 0.501), nrow = 1L), empty, 0.002 ^ 3 / 12
  )
)

rows <- lapply(names(fixtures), function(case) {
  fixture <- fixtures[[case]]
  result <- landscape_rust_prototype_dimension(
    fixture[[1L]], fixture[[2L]], 0L, library = library
  )
  error <- abs(result$squared_distance - fixture[[3L]])
  data.frame(
    contract_id = "mv05bf_ci_fixture_v1", case = case,
    expected_squared_distance = fixture[[3L]],
    observed_squared_distance = result$squared_distance,
    absolute_error = error, status = result$status,
    engine_version = result$engine_version, rust_used = result$rust_used,
    fallback_used = FALSE,
    passed = isTRUE(result$rust_used) && result$status == 0L &&
      result$engine_version == 1L && error <= 1e-12,
    stringsAsFactors = FALSE
  )
})

missing <- landscape_rust_prototype_with_fallback(
  matrix(c(0, 2), nrow = 1L), empty, 0L,
  reference_squared = function() 2 / 3,
  library = paste0(library, ".missing")
)
rows[[length(rows) + 1L]] <- data.frame(
  contract_id = "mv05bf_ci_fixture_v1", case = "missing_library_fallback",
  expected_squared_distance = 2 / 3,
  observed_squared_distance = missing$squared_distance,
  absolute_error = abs(missing$squared_distance - 2 / 3),
  status = missing$status, engine_version = missing$engine_version,
  rust_used = missing$rust_used, fallback_used = missing$fallback_used,
  passed = isTRUE(missing$fallback_used) && !isTRUE(missing$rust_used) &&
    missing$status == 9001L && identical(missing$squared_distance, 2 / 3),
  stringsAsFactors = FALSE
)

corrupt <- tempfile("mv05bf-corrupt-library-")
writeBin(charToRaw("not a dynamic library"), corrupt)
on.exit(unlink(corrupt), add = TRUE)
corrupt_result <- landscape_rust_prototype_with_fallback(
  matrix(c(0, 2), nrow = 1L), empty, 0L,
  reference_squared = function() 2 / 3, library = corrupt
)
rows[[length(rows) + 1L]] <- data.frame(
  contract_id = "mv05bf_ci_fixture_v1", case = "corrupt_library_fallback",
  expected_squared_distance = 2 / 3,
  observed_squared_distance = corrupt_result$squared_distance,
  absolute_error = abs(corrupt_result$squared_distance - 2 / 3),
  status = corrupt_result$status,
  engine_version = corrupt_result$engine_version,
  rust_used = corrupt_result$rust_used,
  fallback_used = corrupt_result$fallback_used,
  passed = isTRUE(corrupt_result$fallback_used) &&
    !isTRUE(corrupt_result$rust_used) && corrupt_result$status == 9001L &&
    identical(corrupt_result$squared_distance, 2 / 3),
  stringsAsFactors = FALSE
)

rows <- do.call(rbind, rows)
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
temporary <- paste0(output, ".tmp")
utils::write.csv(rows, temporary, row.names = FALSE)
if (!file.rename(temporary, output)) {
  unlink(temporary)
  stop("Could not atomically publish CI fixture evidence.", call. = FALSE)
}
if (!all(rows$passed)) stop("MV5-BF R CI fixtures failed.", call. = FALSE)
cat("MV5-BF R CI fixtures: ", sum(rows$passed), "/", nrow(rows),
    " passed\n", sep = "")
