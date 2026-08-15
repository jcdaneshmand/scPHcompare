# MV7-F independent-validation helpers.

#' Compare manifest content while ignoring incidental vector names.
mv07f_manifest_matches_v1 <- function(expected_sha, actual_sha,
                                      expected_bytes, actual_bytes) {
  identical(as.character(expected_sha), unname(as.character(actual_sha))) &&
    identical(as.numeric(expected_bytes), unname(as.numeric(actual_bytes)))
}

#' Compare one-row provenance records after CSV round-tripping.
mv07f_provenance_record_matches_v1 <- function(expected, actual) {
  is.data.frame(expected) && is.data.frame(actual) &&
    identical(names(expected), names(actual)) && nrow(expected) == 1L &&
    nrow(actual) == 1L && all(vapply(names(expected), function(name) {
      identical(as.character(expected[[name]]), as.character(actual[[name]]))
    }, logical(1L)))
}

#' Compare serialized mtimes within an explicit sub-millisecond tolerance.
mv07f_mtimes_match_v1 <- function(expected, actual,
                                  tolerance_seconds = 1e-4) {
  expected <- as.numeric(expected)
  actual <- as.numeric(actual)
  length(expected) == length(actual) && length(expected) > 0L &&
    all(is.finite(expected)) && all(is.finite(actual)) &&
    is.numeric(tolerance_seconds) && length(tolerance_seconds) == 1L &&
    is.finite(tolerance_seconds) && tolerance_seconds >= 0 &&
    all(abs(actual - expected) <= tolerance_seconds)
}
