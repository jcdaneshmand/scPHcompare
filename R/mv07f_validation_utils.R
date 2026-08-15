# MV7-F independent-validation helpers.

#' Compare manifest content while ignoring incidental vector names.
mv07f_manifest_matches_v1 <- function(expected_sha, actual_sha,
                                      expected_bytes, actual_bytes) {
  identical(as.character(expected_sha), unname(as.character(actual_sha))) &&
    identical(as.numeric(expected_bytes), unname(as.numeric(actual_bytes)))
}
