#' Logging utilities
#'
#' Provides a simple logging helper that prepends timestamps to messages.
#'
#' @param msg Character string to log.
#' @return Invisibly returns `NULL` after printing the message.
#' @examples
log_env <- new.env(parent = emptyenv())

#' log_message("processing")
log_message <- function(msg) {
  timestamped <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%OS6"), msg)
  message(timestamped)

  log_con <- if (exists("log_con", envir = log_env, inherits = FALSE)) {
    get("log_con", envir = log_env, inherits = FALSE)
  } else {
    NULL
  }
  if (!is.null(log_con) && isOpen(log_con, rw = "write")) {
    writeLines(timestamped, con = log_con)
    flush(log_con)
  }

  invisible(NULL)
}

#' Write debug-level messages with a consistent prefix.
#'
#' @inheritParams log_message
log_debug <- function(msg) {
  log_message(sprintf("DEBUG: %s", msg))
}

#' Configure logging to also append to a file.
#'
#' Sets up a shared connection used by all logging helpers so that messages are
#' both emitted to the console and persisted to disk.
#'
#' @param log_file_path Path to the log file to append to.
#' @return Invisibly returns the normalized log file path.
set_log_file <- function(log_file_path) {
  if (missing(log_file_path) || is.null(log_file_path) || !nzchar(log_file_path)) {
    stop("log_file_path must be a non-empty string.")
  }

  normalized <- normalizePath(log_file_path, mustWork = FALSE)

  existing_con <- if (exists("log_con", envir = log_env, inherits = FALSE)) {
    get("log_con", envir = log_env, inherits = FALSE)
  } else {
    NULL
  }
  if (!is.null(existing_con) && isOpen(existing_con)) {
    close(existing_con)
  }

  assign("log_con", file(normalized, open = "a"), envir = log_env)
  invisible(normalized)
}

#' Close and clear the shared log file connection if present.
close_log_file <- function() {
  existing_con <- if (exists("log_con", envir = log_env, inherits = FALSE)) {
    get("log_con", envir = log_env, inherits = FALSE)
  } else {
    NULL
  }
  if (!is.null(existing_con) && isOpen(existing_con)) {
    close(existing_con)
  }
  assign("log_con", NULL, envir = log_env)
  invisible(NULL)
}

