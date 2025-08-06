#' Logging utilities
#'
#' Provides a simple logging helper that prepends timestamps to messages.
log_message <- function(msg) {
  message(sprintf("[%s] %s", Sys.time(), msg))
}

