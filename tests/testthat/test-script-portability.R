test_that("repository scripts do not pin a local checkout root", {
  scripts_dir <- "../../scripts"
  skip_if_not(dir.exists(scripts_dir),
              "repository-only scripts are excluded from source packages")

  files <- list.files(scripts_dir, pattern = "[.]R$", full.names = TRUE)
  contents <- stats::setNames(lapply(files, readLines, warn = FALSE), files)
  absolute_setwd <- vapply(contents, function(lines) {
    any(grepl(
      "setwd\\(\\s*['\"](?:/|[A-Za-z]:[/\\\\])",
      lines,
      perl = TRUE
    ))
  }, logical(1L))

  expect_false(any(absolute_setwd), info = paste(
    "absolute setwd() remains in",
    paste(names(absolute_setwd)[absolute_setwd], collapse = ", ")
  ))

  bootstrap <- vapply(contents, function(lines) {
    any(grepl("Unable to resolve the repository root.", lines, fixed = TRUE))
  }, logical(1L))
  expect_gt(sum(bootstrap), 0L)
  expect_true(all(vapply(contents[bootstrap], function(lines) {
    text <- paste(lines, collapse = "\n")
    grepl("~+~", text, fixed = TRUE) &&
      grepl("fixed\\s*=\\s*TRUE", text, perl = TRUE)
  }, logical(1L))))
})
