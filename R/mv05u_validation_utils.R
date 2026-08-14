#' Parse a strict boolean column from cross-language audit artifacts.
#'
#' Python's CSV writer emits `True` and `False`, while R's CSV type inference
#' commonly leaves those title-case values as character data. Validation must
#' recognize their boolean meaning without accepting arbitrary truthy values.
#'
#' @param values A logical or character vector.
#' @param field_name A field name used in validation errors.
#' @return A logical vector with the same length as `values`.
#' @keywords internal
mv05u_parse_strict_boolean_v1 <- function(values, field_name) {
  normalized <- tolower(trimws(as.character(values)))
  invalid <- is.na(values) | is.na(normalized) |
    !normalized %in% c("true", "false")
  if (any(invalid)) {
    stop(
      "MV5-U boolean field `", field_name,
      "` must contain only true/false values.",
      call. = FALSE
    )
  }
  normalized == "true"
}
