#' Logging utilities
#'
#' Provides a simple logging helper that prepends timestamps to messages.
#'
#' @param msg Character string to log.
#' @return Invisibly returns `NULL` after printing the message.
#' @examples
#' log_message("processing")
log_message <- function(msg) {
  message(sprintf("[%s] %s", Sys.time(), msg))
  invisible(NULL)
}

