#!/usr/bin/env Rscript

# Restore the committed dependency baseline into an explicitly selected library.
#
# This script is intentionally runnable with base R plus renv/activate.R. It
# handles the one known retired-Bioconductor bootstrap gap: BiocVersion 3.19.1
# is no longer available from the Bioconductor 3.19 package archive, while its
# exact official Git commit remains available.

parse_args <- function(args) {
  values <- list(
    project = getwd(),
    library = Sys.getenv("RENV_PATHS_LIBRARY", unset = ""),
    report = ""
  )
  for (arg in args) {
    if (startsWith(arg, "--project=")) {
      values$project <- sub("^--project=", "", arg)
    } else if (startsWith(arg, "--library=")) {
      values$library <- sub("^--library=", "", arg)
    } else if (startsWith(arg, "--report=")) {
      values$report <- sub("^--report=", "", arg)
    } else {
      stop("Unknown argument: ", arg, call. = FALSE)
    }
  }
  values
}

fail_if_empty <- function(value, label) {
  if (!nzchar(value)) {
    stop(label, " must be supplied explicitly.", call. = FALSE)
  }
  value
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
project <- normalizePath(args$project, winslash = "/", mustWork = TRUE)
library <- fail_if_empty(args$library, "--library")
dir.create(library, recursive = TRUE, showWarnings = FALSE)
library <- normalizePath(library, winslash = "/", mustWork = TRUE)
lockfile <- file.path(project, "renv.lock")
activate <- file.path(project, "renv", "activate.R")

if (!file.exists(lockfile) || !file.exists(activate)) {
  stop("The project must contain renv.lock and renv/activate.R.", call. = FALSE)
}

Sys.setenv(RENV_PATHS_LIBRARY = library)
source(activate, local = globalenv())

if (as.character(utils::packageVersion("renv")) != "1.1.5") {
  stop("Expected renv 1.1.5 from the committed bootstrap.", call. = FALSE)
}

lock <- renv::lockfile_read(lockfile)
bioc_record <- lock$Packages$BiocVersion
if (is.null(bioc_record)) {
  stop("renv.lock does not contain the required BiocVersion record.", call. = FALSE)
}

expected_bioc_version <- "3.19.1"
expected_bioc_sha <- "99e637d62c373025b9a757047a144a03c32905cd"
bioc_commit_marker <- file.path(library, ".scphcompare-BiocVersion-commit")
if (!identical(bioc_record$Version, expected_bioc_version) ||
    !startsWith(expected_bioc_sha, bioc_record$git_last_commit)) {
  stop("The BiocVersion bootstrap pin disagrees with renv.lock.", call. = FALSE)
}

installed <- utils::installed.packages(lib.loc = library)
recorded_bioc_sha <- if (file.exists(bioc_commit_marker)) {
  readLines(bioc_commit_marker, n = 1L, warn = FALSE)
} else {
  ""
}
bioc_ready <- "BiocVersion" %in% rownames(installed) &&
  identical(unname(installed["BiocVersion", "Version"]), expected_bioc_version) &&
  identical(recorded_bioc_sha, expected_bioc_sha)

if (!bioc_ready) {
  checkout <- tempfile("scphcompare-BiocVersion-")
  clone_status <- system2(
    "git",
    c("clone", "--quiet", "--branch", "RELEASE_3_19",
      "https://git.bioconductor.org/packages/BiocVersion", shQuote(checkout))
  )
  if (!identical(clone_status, 0L)) {
    stop("Could not clone the official BiocVersion repository.", call. = FALSE)
  }
  checkout_status <- system2(
    "git",
    c("-C", shQuote(checkout), "checkout", "--quiet", expected_bioc_sha)
  )
  if (!identical(checkout_status, 0L)) {
    stop("Could not check out the locked BiocVersion commit.", call. = FALSE)
  }
  install_status <- system2(
    file.path(R.home("bin"), "R"),
    c("CMD", "INSTALL", paste0("--library=", shQuote(library)), shQuote(checkout))
  )
  if (!identical(install_status, 0L)) {
    stop("Could not install the locked BiocVersion bootstrap.", call. = FALSE)
  }
  writeLines(expected_bioc_sha, bioc_commit_marker, useBytes = TRUE)
}

restore_started <- Sys.time()
renv::restore(
  project = project,
  library = library,
  lockfile = lockfile,
  clean = TRUE,
  prompt = FALSE
)
restore_finished <- Sys.time()

records <- lock$Packages
renv_library <- dirname(find.package("renv"))
search_libraries <- unique(c(library, renv_library, .Library))
installed_databases <- lapply(
  search_libraries,
  utils::installed.packages
)

find_installed_field <- function(package, field) {
  for (i in seq_along(installed_databases)) {
    database <- installed_databases[[i]]
    if (package %in% rownames(database)) {
      if (field == "Version") {
        return(unname(database[package, "Version"]))
      }
      description <- utils::packageDescription(
        package,
        lib.loc = search_libraries[[i]],
        fields = field
      )
      if (length(description) == 1L && !is.na(description)) {
        return(unname(description))
      }
      return("")
    }
  }
  NA_character_
}

locked_commit <- vapply(records, function(record) {
  if (!is.null(record$RemoteSha)) {
    record$RemoteSha
  } else if (!is.null(record$git_last_commit)) {
    record$git_last_commit
  } else {
    ""
  }
}, character(1))
installed_commit <- vapply(names(records), function(package) {
  remote_sha <- find_installed_field(package, "RemoteSha")
  if (!is.na(remote_sha) && nzchar(remote_sha)) {
    return(remote_sha)
  }
  git_commit <- find_installed_field(package, "git_last_commit")
  if (!is.na(git_commit) && nzchar(git_commit)) {
    return(git_commit)
  }
  if (identical(package, "BiocVersion") && file.exists(bioc_commit_marker)) {
    return(readLines(bioc_commit_marker, n = 1L, warn = FALSE))
  }
  ""
}, character(1))

verification <- data.frame(
  package = names(records),
  locked_version = vapply(records, `[[`, character(1), "Version"),
  installed_version = vapply(
    names(records), find_installed_field, character(1), field = "Version"
  ),
  locked_source = vapply(records, `[[`, character(1), "Source"),
  locked_commit = locked_commit,
  installed_commit = installed_commit,
  stringsAsFactors = FALSE
)
verification$version_matches <- !is.na(verification$installed_version) &
  verification$locked_version == verification$installed_version
verification$commit_matches <- !nzchar(verification$locked_commit) |
  (nzchar(verification$installed_commit) &
     (startsWith(verification$installed_commit, verification$locked_commit) |
        startsWith(verification$locked_commit, verification$installed_commit)))
verification$matches <- verification$version_matches & verification$commit_matches

if (nzchar(args$report)) {
  report_dir <- dirname(args$report)
  dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(verification, args$report, row.names = FALSE)
}

cat("project=", project, "\n", sep = "")
cat("library=", library, "\n", sep = "")
cat("lock_records=", nrow(verification), "\n", sep = "")
cat("exact_records=", sum(verification$matches), "\n", sep = "")
cat("mismatched_records=", sum(!verification$matches), "\n", sep = "")
cat(
  "restore_elapsed_seconds=",
  sprintf("%.1f", as.numeric(difftime(restore_finished, restore_started, units = "secs"))),
  "\n",
  sep = ""
)

if (!all(verification$matches)) {
  print(verification[!verification$matches, ], row.names = FALSE)
  stop("The restored library does not match renv.lock exactly.", call. = FALSE)
}
