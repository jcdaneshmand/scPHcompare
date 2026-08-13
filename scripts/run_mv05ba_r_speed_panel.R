#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: run_mv05ba_r_speed_panel.R DIAGRAM_RDS REFERENCES SPEED_PANEL")
}
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Unable to resolve the repository root.")
setwd(normalizePath(file.path(
  dirname(gsub("~+~", " ", sub("^--file=", "", script_arg[[1L]]), fixed = TRUE)), ".."
), mustWork = TRUE))
pkgload::load_all(".", quiet = TRUE)
references <- readRDS(args[[2L]])
panel <- utils::read.csv(args[[3L]], stringsAsFactors = FALSE)
diagrams <- readRDS(args[[1L]])
rows <- vector("list", nrow(panel))
output_path <- file.path(dirname(args[[3L]]), "private-r-speed.csv")
for (index in seq_len(nrow(panel))) {
  item <- panel[index, ]
  first <- diagrams[[item$first_diagram_id]]
  second <- diagrams[[item$second_diagram_id]]
  started <- proc.time()[["elapsed"]]
  value <- persistence_landscape_distance(
    first, second, method = "auto", exact_max_intervals = 500L,
    abs_tol = 1e-8, rel_tol = 1e-8, subdivisions = 200L,
    first_id = item$first_diagram_id, second_id = item$second_diagram_id
  )
  elapsed <- proc.time()[["elapsed"]] - started
  accepted <- references[
    references$stratum_id == item$stratum_id &
      references$pair_order == item$pair_order, ]
  if (nrow(accepted) != 2L ||
      !isTRUE(all.equal(value$dimensions$H0$squared_distance,
                        accepted$reference_squared_distance[
                          accepted$dimension == "H0"], tolerance = 0)) ||
      !isTRUE(all.equal(value$dimensions$H1$squared_distance,
                        accepted$reference_squared_distance[
                          accepted$dimension == "H1"], tolerance = 0))) {
    stop("R speed-panel repeat differs from accepted shard: ", item$panel_order)
  }
  rows[[index]] <- data.frame(
    contract_id = "mv05ba_private_r_speed_v1",
    panel_order = item$panel_order, stratum_id = item$stratum_id,
    pair_order = item$pair_order,
    accepted_cache_key = item$pair_cache_key,
    repeated_cache_key = value$cache_key,
    h0_squared_distance = value$dimensions$H0$squared_distance,
    h1_squared_distance = value$dimensions$H1$squared_distance,
    elapsed_seconds = elapsed, identity_repeated =
      identical(value$cache_key, item$pair_cache_key),
    stringsAsFactors = FALSE
  )
  completed <- do.call(rbind, rows[seq_len(index)])
  temporary <- paste0(output_path, ".tmp")
  utils::write.csv(completed, temporary, row.names = FALSE, quote = TRUE)
  if (!file.rename(temporary, output_path)) stop("Atomic R speed write failed.")
  cat("R panel", index, "of", nrow(panel), "seconds", elapsed, "\n")
}
