args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: mv05bg_verify_locked_library.R <report.csv>", call. = FALSE)
}

lock <- renv::lockfile_read("renv.lock")
records <- lock$Packages
libraries <- .libPaths()
databases <- lapply(libraries, utils::installed.packages)

installed_field <- function(package, field) {
  for (i in seq_along(databases)) {
    database <- databases[[i]]
    if (package %in% rownames(database)) {
      if (field == "Version") return(unname(database[package, "Version"]))
      value <- utils::packageDescription(
        package, lib.loc = libraries[[i]], fields = field
      )
      if (length(value) == 1L && !is.na(value)) return(unname(value))
      return("")
    }
  }
  NA_character_
}

locked_commit <- vapply(records, function(record) {
  if (!is.null(record$RemoteSha)) return(record$RemoteSha)
  if (!is.null(record$git_last_commit)) return(record$git_last_commit)
  ""
}, character(1L))
installed_commit <- vapply(names(records), function(package) {
  remote <- installed_field(package, "RemoteSha")
  if (!is.na(remote) && nzchar(remote)) return(remote)
  git <- installed_field(package, "git_last_commit")
  if (!is.na(git) && nzchar(git)) return(git)
  if (package == "BiocVersion") {
    for (library in libraries) {
      marker <- file.path(library, ".scphcompare-BiocVersion-commit")
      if (file.exists(marker)) return(readLines(marker, n = 1L, warn = FALSE))
    }
  }
  ""
}, character(1L))

result <- data.frame(
  contract_id = "mv05bg_locked_library_verification_v1",
  package = names(records),
  locked_version = vapply(records, `[[`, character(1L), "Version"),
  installed_version = vapply(
    names(records), installed_field, character(1L), field = "Version"
  ),
  locked_commit = locked_commit,
  installed_commit = installed_commit,
  stringsAsFactors = FALSE
)
result$version_matches <- !is.na(result$installed_version) &
  result$locked_version == result$installed_version
result$commit_matches <- !nzchar(result$locked_commit) |
  (nzchar(result$installed_commit) &
     (startsWith(result$installed_commit, result$locked_commit) |
        startsWith(result$locked_commit, result$installed_commit)))
result$matches <- result$version_matches & result$commit_matches

dir.create(dirname(args[[1L]]), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(result, args[[1L]], row.names = FALSE)
cat("MV5-BG locked library verification: ", sum(result$matches), "/",
    nrow(result), " matched\n", sep = "")
if (!all(result$matches)) {
  print(result[!result$matches, ], row.names = FALSE)
  stop("Active R library does not exactly satisfy renv.lock.", call. = FALSE)
}
