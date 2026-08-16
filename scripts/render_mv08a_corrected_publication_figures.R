#!/usr/bin/env Rscript

options(warn = 2, digits = 17, scipen = 999)
packages <- c("digest", "ggplot2", "patchwork", "svglite", "ragg", "scales")
for (package in packages) if (!requireNamespace(package, quietly = TRUE))
  stop(package, " required.")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("usage: render_mv08a_corrected_publication_figures.R PH_RDS H1_PRIVATE_CSV OUTPUT_DIR EXPECTED_HEAD",
       call. = FALSE)
}
ph_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
h1_path <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
output_dir <- args[[3L]]
expected_head <- tolower(trimws(args[[4L]]))
observed_head <- tolower(trimws(system2("git", c("rev-parse", "HEAD"), stdout = TRUE)))
if (!identical(observed_head, expected_head)) stop("MV8-A prospective HEAD mismatch.")
if (dir.exists(output_dir) && length(list.files(output_dir, all.files = TRUE,
                                                 no.. = TRUE))) {
  stop("MV8-A output directory must be empty.")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "data"), showWarnings = FALSE)

readc <- function(path) utils::read.csv(path, stringsAsFactors = FALSE,
                                         check.names = FALSE)
sha <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
write_csv <- function(x, path) utils::write.table(
  x, path, sep = ",", row.names = FALSE, col.names = TRUE, quote = TRUE,
  na = "", qmethod = "double")

paths <- c(
  mv07j_decision = "docs/audits/mv07j-synthesis-evidence/mv07j-decision.csv",
  landscape_contract = "docs/audits/mv07j-synthesis-evidence/mv07j-landscape-observation-contract.csv",
  h0_composite = "docs/audits/mv07j-synthesis-evidence/mv07j-h0-composite-concordance.csv",
  stability = "docs/audits/mv07j-synthesis-evidence/mv07j-complete-stability-curve.csv",
  outcomes = "docs/audits/mv07j-synthesis-evidence/mv07j-complete-outcome-summary.csv",
  algorithm = "docs/audits/mv07j-synthesis-evidence/mv07j-algorithm-sensitivity.csv",
  external_gate = "docs/audits/mv07j-synthesis-evidence/mv07j-external-data-gate.csv",
  tissue_summary = "docs/audits/mv07d-prefreeze-evidence/mv07d-tissue-study-summary.csv",
  metadata_counts = "docs/audits/mv07i-outcome-prefreeze-evidence/mv07i-outcome-metadata-counts.csv",
  ph_metrics = "docs/audits/mv07h-full-ph-evidence/mv07h-ph-metrics.csv",
  label_manifest = "docs/audits/mv07i-label-closed-validation/mv07i-label-closed-artifact-manifest.csv",
  specification = "docs/specifications/MV08A_CORRECTED_PUBLICATION_FIGURE_PREFREEZE_V1.md",
  renderer = "scripts/render_mv08a_corrected_publication_figures.R",
  validator = "scripts/validate_mv08a_corrected_publication_figures.R")
if (any(!file.exists(paths))) stop("MV8-A public source missing.")

decision <- readc(paths[["mv07j_decision"]])
if (nrow(decision) != 1L || decision$decision !=
    "authorize_corrected_figure_implementation_and_read_only_external_dataset_audit")
  stop("MV8-A admission decision mismatch.")
label_manifest <- readc(paths[["label_manifest"]])
h1_manifest <- label_manifest[label_manifest$artifact == "h1_summary", , drop = FALSE]
if (nrow(h1_manifest) != 1L || !identical(tolower(sha(h1_path)),
                                           tolower(h1_manifest$production_sha256)))
  stop("MV8-A H1 source hash mismatch.")
ph_metrics <- readc(paths[["ph_metrics"]])
ph_row <- ph_metrics[ph_metrics$seed == 20260807L &
  ph_metrics$sample_id == "SRA701877_SRS3279688" &
  ph_metrics$view_id == "cell_topology_v1", , drop = FALSE]
if (nrow(ph_row) != 1L || !identical(tolower(sha(ph_path)),
                                     tolower(ph_row$output_sha256)))
  stop("MV8-A fixed PH artifact hash mismatch.")

navy <- "#1B3A57"; teal <- "#238A8D"; orange <- "#E69F00"
rust <- "#C94C4C"; purple <- "#6D597A"; blue <- "#3B82B8"
grey <- "#6B7280"; pale <- "#E8EEF2"; dark <- "#17202A"
view_cols <- c(cell = teal, gene = purple)
rep_cols <- c(cell_H0 = navy, cell_H1 = teal,
  cell_H0_H1_secondary = blue, gene_H0 = "#7A4E2D", gene_H1 = purple,
  gene_H0_H1_secondary = rust)
pretty_rep <- c(cell_H0 = "Cell H0", cell_H1 = "Cell H1",
  cell_H0_H1_secondary = "Cell H0+H1 (secondary)", gene_H0 = "Gene H0",
  gene_H1 = "Gene H1", gene_H0_H1_secondary = "Gene H0+H1 (secondary)")
theme_pub <- function(base_size = 10) ggplot2::theme_minimal(base_size = base_size,
  base_family = "sans") + ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = base_size + 3,
      color = dark, margin = ggplot2::margin(b = 8)),
    plot.subtitle = ggplot2::element_text(color = grey, margin = ggplot2::margin(b = 8)),
    plot.caption = ggplot2::element_text(color = grey, hjust = 0, size = base_size),
    axis.title = ggplot2::element_text(face = "bold", color = dark),
    axis.text = ggplot2::element_text(color = dark),
    strip.text = ggplot2::element_text(face = "bold", color = dark),
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_line(color = "#E5E7EB", linewidth = .25),
    legend.position = "bottom", legend.title = ggplot2::element_blank())

figure_records <- list(); record_index <- 0L
save_pair <- function(plot, id, width, height, title) {
  svg <- file.path(output_dir, paste0(id, ".svg"))
  png <- file.path(output_dir, paste0(id, ".png"))
  ggplot2::ggsave(svg, plot = plot, device = svglite::svglite,
    width = width, height = height, units = "in", bg = "white")
  ggplot2::ggsave(png, plot = plot, device = ragg::agg_png,
    width = width, height = height, units = "in", dpi = 300, bg = "white")
  record_index <<- record_index + 1L
  figure_records[[record_index]] <<- data.frame(
    contract_id = "mv08a_figure_manifest_v1", figure_id = id, title = title,
    svg_filename = basename(svg), svg_sha256 = sha(svg), svg_bytes = file.info(svg)$size,
    png_filename = basename(png), png_sha256 = sha(png), png_bytes = file.info(png)$size,
    width_inches = width, height_inches = height,
    png_width_pixels = round(width * 300), png_height_pixels = round(height * 300),
    dpi = 300L, manuscript_claim_authorized = FALSE, stringsAsFactors = FALSE)
}

# F1: corrected architecture.
nodes <- data.frame(
  id = c("sample", "cell", "gene", "cellph", "geneph", "cellland", "geneland",
    "distance", "matrix", "cluster", "outcome", "composite"),
  x = c(1, 2.4, 2.4, 4, 4, 5.7, 5.7, 7.4, 9, 10.6, 12.2, 7.4),
  y = c(1.5, 2.25, .75, 2.25, .75, 2.25, .75, 1.5, 1.5, 1.5, 1.5, .25),
  label = c("Biological sample\ncomparison unit",
    "Cell view\ncells = observations", "Gene view\ngenes = observations",
    "Cell H0 / H1\ndiagrams", "Gene H0 / H1\ndiagrams",
    "All-level cell\nlandscapes", "All-level gene\nlandscapes",
    "Separate H0 / H1\nexact or controlled L2", "Seed-specific\nsample matrices",
    "Label-free PAM\n+ average sensitivity", "Descriptive\nmetadata alignment",
    "Secondary unweighted\nH0+H1 composite"),
  role = c("unit", "cell", "gene", "cell", "gene", "cell", "gene", "method",
    "method", "method", "outcome", "secondary"), stringsAsFactors = FALSE)
edges <- data.frame(from = c("sample", "sample", "cell", "gene", "cellph", "geneph",
  "cellland", "geneland", "distance", "matrix", "cluster", "distance"),
  to = c("cell", "gene", "cellph", "geneph", "cellland", "geneland", "distance",
    "distance", "matrix", "cluster", "outcome", "composite"),
  secondary = c(rep(FALSE, 11L), TRUE), stringsAsFactors = FALSE)
edges <- merge(edges, nodes[c("id", "x", "y")], by.x = "from", by.y = "id")
names(edges)[names(edges) %in% c("x", "y")] <- c("x", "y")
edges <- merge(edges, nodes[c("id", "x", "y")], by.x = "to", by.y = "id",
               suffixes = c("", "end"))
role_fill <- c(unit = "#DDEAF2", cell = "#DDF2EE", gene = "#ECE4F1",
  method = "#E8EEF2", outcome = "#F7E5DA", secondary = "#F3F4F6")
f1 <- ggplot2::ggplot() +
  ggplot2::geom_segment(data = edges,
    ggplot2::aes(x = x + .45, y = y, xend = xend - .45, yend = yend,
      linetype = secondary), color = grey, linewidth = .55,
    arrow = grid::arrow(length = grid::unit(2.2, "mm"), type = "closed")) +
  ggplot2::geom_label(data = nodes, ggplot2::aes(x, y, label = label, fill = role),
    color = dark, size = 3.1, label.size = .3,
    label.padding = grid::unit(.16, "lines")) +
  ggplot2::scale_fill_manual(values = role_fill, guide = "none") +
  ggplot2::scale_linetype_manual(values = c(`FALSE` = "solid", `TRUE` = "22"),
                                 guide = "none") +
  ggplot2::annotate("text", x = 5.7, y = 3.0,
    label = "Finite positive intervals • essential H0 excluded • all active levels",
    color = dark, fontface = "bold", size = 3.4) +
  ggplot2::annotate("text", x = 8.2, y = -.25,
    label = "No fixed scientific grid • no level cap • samples remain the independent comparison units",
    color = grey, size = 3.1) +
  ggplot2::coord_cartesian(xlim = c(.25, 13), ylim = c(-.45, 3.25), clip = "off") +
  ggplot2::labs(title = "Corrected scPHcompare sample-comparison architecture") +
  ggplot2::theme_void(base_family = "sans") +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 15,
    color = dark, margin = ggplot2::margin(b = 8)), plot.margin = ggplot2::margin(10, 12, 10, 12))
write_csv(nodes, file.path(output_dir, "data", "f1-nodes.csv"))
write_csv(edges[c("from", "to", "secondary")], file.path(output_dir, "data", "f1-edges.csv"))
save_pair(f1, "figure-1-corrected-architecture", 15, 5.2,
          "Corrected sample-comparison and dual-view pipeline")

# F2: cohort and confounding.
tissue <- readc(paths[["tissue_summary"]])
tissue$added_descriptive <- tissue$retained_samples - tissue$primary_samples
tissue_long <- rbind(
  data.frame(tissue = tissue$tissue, group = "Primary 90",
             samples = tissue$primary_samples),
  data.frame(tissue = tissue$tissue, group = "Added descriptive 34",
             samples = tissue$added_descriptive))
tissue_long$tissue <- factor(tissue_long$tissue,
  levels = tissue$tissue[order(tissue$retained_samples)])
cohort <- data.frame(x = c(1, 3, 5, 5, 7), y = c(1, 1, 1.6, .4, 1),
  label = c("127 discovered", "124 retained", "90 primary\n5 cross-study tissues",
    "34 added\n3 single-study tissues", "3 excluded\n<250 cells"),
  kind = c("all", "retained", "primary", "added", "excluded"))
cohort_edges <- data.frame(x = c(1.6, 3.6, 3.6, 1.6), y = c(1, 1, 1, 1),
  xend = c(2.4, 4.4, 4.4, 6.4), yend = c(1, 1.6, .4, 1),
  style = c("solid", "solid", "solid", "dashed"))
cohort_plot <- ggplot2::ggplot() +
  ggplot2::geom_segment(data = cohort_edges,
    ggplot2::aes(x, y, xend = xend, yend = yend, linetype = style),
    color = grey, arrow = grid::arrow(length = grid::unit(2, "mm"))) +
  ggplot2::geom_label(data = cohort, ggplot2::aes(x, y, label = label, fill = kind),
    size = 3.2, label.size = .3, color = dark) +
  ggplot2::scale_fill_manual(values = c(all = pale, retained = "#DDEAF2",
    primary = "#DDF2EE", added = "#F7E5DA", excluded = "#F3F4F6"), guide = "none") +
  ggplot2::scale_linetype_identity() + ggplot2::coord_cartesian(xlim = c(.3, 7.7),
    ylim = c(0, 2.2), clip = "off") + ggplot2::theme_void() +
  ggplot2::labs(title = "A  Cohort accountability") +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 12,
    color = dark, hjust = 0))
tissue_plot <- ggplot2::ggplot(tissue_long,
  ggplot2::aes(samples, tissue, fill = group)) +
  ggplot2::geom_col(width = .72) +
  ggplot2::geom_text(data = tissue, ggplot2::aes(x = retained_samples + .5,
    y = tissue, label = paste0(retained_samples, " / ", retained_studies, " studies")),
    inherit.aes = FALSE, hjust = 0, size = 3, color = dark) +
  ggplot2::scale_fill_manual(values = c("Primary 90" = teal,
    "Added descriptive 34" = orange)) +
  ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, .28))) +
  ggplot2::labs(title = "B  Tissue composition", x = "Retained samples", y = NULL) +
  theme_pub(9)
nesting <- ggplot2::ggplot() +
  ggplot2::annotate("rect", xmin = .5, xmax = 4.5, ymin = .55, ymax = 1.45,
    fill = "#DDEAF2", color = navy) +
  ggplot2::annotate("text", x = 2.5, y = 1.12,
    label = "118 scRNA-seq samples", fontface = "bold", size = 4, color = dark) +
  ggplot2::annotate("rect", xmin = 2.9, xmax = 4.3, ymin = .68, ymax = 1.0,
    fill = "#ECE4F1", color = purple) +
  ggplot2::annotate("text", x = 3.6, y = .84,
    label = "6 snRNA-seq", fontface = "bold", size = 3.2, color = dark) +
  ggplot2::annotate("text", x = 2.5, y = .3,
    label = "All six snRNA-seq = substantia nigra = SRA850958\nApproach effect is not identifiable",
    size = 3.4, color = rust, fontface = "bold") +
  ggplot2::coord_cartesian(xlim = c(.3, 4.7), ylim = c(.05, 1.7)) +
  ggplot2::labs(title = "C  Perfect approach nesting") + ggplot2::theme_void() +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 12,
    color = dark))
f2 <- (cohort_plot / (tissue_plot | nesting)) +
  patchwork::plot_annotation(title = "Cohort structure and non-identifiability",
    subtitle = "Eight tissues and 18 studies; three added tissues are single-study descriptions",
    theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 15,
      color = dark), plot.subtitle = ggplot2::element_text(color = grey)))
write_csv(tissue_long, file.path(output_dir, "data", "f2-tissue-counts.csv"))
write_csv(cohort, file.path(output_dir, "data", "f2-cohort-flow.csv"))
save_pair(f2, "figure-2-cohort-confounding", 12, 8,
          "Cohort accountability and confounding structure")

# F3: one fixed validated diagram and display-only landscape evaluation.
ph <- readRDS(ph_path)
pd <- as.data.frame(ph$topology_result$diagram)
names(pd)[1:3] <- c("dimension", "birth", "death")
pd <- pd[is.finite(pd$birth) & is.finite(pd$death) & pd$death > pd$birth &
         pd$dimension %in% c(0, 1), ]
pd$homology <- factor(paste0("H", pd$dimension), levels = c("H0", "H1"))
limits <- range(c(pd$birth, pd$death))
diag_plot <- ggplot2::ggplot(pd, ggplot2::aes(birth, death, color = homology)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, color = grey, linewidth = .4) +
  ggplot2::geom_point(alpha = .62, size = 1.2) +
  ggplot2::scale_color_manual(values = c(H0 = navy, H1 = rust)) +
  ggplot2::coord_equal(xlim = limits, ylim = limits) +
  ggplot2::labs(title = "A  Finite persistence diagram", x = "Birth",
    y = "Death", caption = "Essential H0 class excluded") + theme_pub(9)
grid_display <- seq(limits[[1]], limits[[2]], length.out = 600L)
land_rows <- list(); cursor <- 0L
for (dimension in 0:1) {
  intervals <- as.matrix(pd[pd$dimension == dimension, c("birth", "death")])
  values <- lapply(grid_display, function(location) {
    v <- pmin(location - intervals[, 1], intervals[, 2] - location)
    sort(v[v > 0], decreasing = TRUE)
  })
  for (level in 1:6) {
    cursor <- cursor + 1L
    land_rows[[cursor]] <- data.frame(homology = paste0("H", dimension),
      level = level, filtration = grid_display,
      value = vapply(values, function(v) if (length(v) >= level) v[[level]] else 0,
                     numeric(1L)))
  }
}
land <- do.call(rbind, land_rows)
land$level_label <- factor(paste0("λ", land$level), levels = paste0("λ", 1:6))
land_plot <- ggplot2::ggplot(land,
  ggplot2::aes(filtration, value, color = level_label, group = level_label)) +
  ggplot2::geom_line(linewidth = .55) +
  ggplot2::facet_wrap(~homology, ncol = 1, scales = "free_y") +
  ggplot2::scale_color_manual(values = c(navy, teal, orange, rust, purple, grey)) +
  ggplot2::labs(title = "B  Consecutive landscape levels",
    subtitle = "λ1–λ6 shown for legibility; every active level enters the scientific distance",
    x = "Filtration value", y = "Landscape height",
    caption = "Dense grid is visualization-only; distance integration is exact or error-controlled") +
  theme_pub(9)
f3 <- (diag_plot | land_plot) + patchwork::plot_layout(widths = c(.82, 1.25))
f3 <- f3 + patchwork::plot_annotation(
  title = "From validated persistence diagrams to corrected all-level landscapes",
  subtitle = "Fixed label-closed cell-view artifact; H0 and H1 remain separate",
  theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 15,
    color = dark), plot.subtitle = ggplot2::element_text(color = grey)))
write_csv(pd[c("dimension", "birth", "death", "homology")],
          file.path(output_dir, "data", "f3-finite-diagram.csv"))
write_csv(land, file.path(output_dir, "data", "f3-display-landscape-levels.csv"))
save_pair(f3, "figure-3-diagram-to-landscape", 13, 7.5,
          "Persistence diagrams to corrected all-level landscapes")

# F4: all-pair H1 contribution and all partition concordances.
h1 <- readc(h1_path)
ecdf_rows <- list(); cursor <- 0L
for (view in c("cell", "gene")) {
  x <- sort(h1$median[h1$view_id == view])
  index <- unique(round(seq(1, length(x), length.out = 300L)))
  cursor <- cursor + 1L
  ecdf_rows[[cursor]] <- data.frame(view = view, contribution = x[index],
    cumulative_fraction = index / length(x))
}
ecdf <- do.call(rbind, ecdf_rows)
med <- aggregate(median ~ view_id, h1, stats::median)
names(med) <- c("view", "median")
med$label_y <- c(cell = .43, gene = .58)[med$view]
med$label <- sprintf("%s median %.3f%%",
  ifelse(med$view == "cell", "Cell", "Gene"), 100 * med$median)
distribution <- ggplot2::ggplot(ecdf,
  ggplot2::aes(contribution, cumulative_fraction, color = view)) +
  ggplot2::geom_line(linewidth = .9) +
  ggplot2::geom_vline(data = med, ggplot2::aes(xintercept = median, color = view),
    linetype = "22", linewidth = .55) +
  ggplot2::geom_label(data = med, ggplot2::aes(x = median, y = label_y,
    label = label, color = view), inherit.aes = FALSE, fill = "white",
    label.size = .25, size = 2.85, show.legend = FALSE) +
  ggplot2::scale_x_log10(labels = function(x) sprintf("%.3g%%", 100 * x)) +
  ggplot2::scale_y_continuous(labels = scales::label_percent()) +
  ggplot2::scale_color_manual(values = view_cols,
    labels = c(cell = "Cell view", gene = "Gene view")) +
  ggplot2::labs(title = "A  H1 share of secondary composite distance",
    subtitle = "All 7,626 median-across-seed sample pairs per view",
    x = "H1 squared-distance contribution (log scale)", y = "Cumulative fraction") +
  theme_pub(9)
hc <- readc(paths[["h0_composite"]])
hc$algorithm <- ifelse(hc$algorithm_id == "pam_stability_k_v1", "PAM", "Average linkage")
hc$view <- factor(hc$view_id, levels = c("cell", "gene"),
  labels = c("Cell view", "Gene view"))
hc$seed_label <- sub("2026", "", hc$seed)
exception <- hc[!hc$exact_partition, , drop = FALSE]
exception$exception_label <- sprintf("Cell seed %s\nARI %.3f",
  exception$seed_label, exception$adjusted_rand_index)
concordance <- ggplot2::ggplot(hc,
  ggplot2::aes(seed_label, adjusted_rand_index, color = view, shape = algorithm)) +
  ggplot2::geom_hline(yintercept = 1, color = grey, linewidth = .35) +
  ggplot2::geom_point(size = 2.6) +
  ggplot2::geom_label(data = exception,
    ggplot2::aes(x = seed_label, y = adjusted_rand_index,
      label = exception_label, color = view), inherit.aes = FALSE,
    nudge_y = .12, fill = "white", label.size = .25, size = 2.85,
    show.legend = FALSE) +
  ggplot2::facet_wrap(~algorithm, ncol = 1) +
  ggplot2::scale_color_manual(values = c("Cell view" = teal, "Gene view" = purple)) +
  ggplot2::scale_shape_manual(values = c(PAM = 16, `Average linkage` = 17)) +
  ggplot2::scale_y_continuous(limits = c(0, 1.04), breaks = c(0, .25, .5, .75, 1)) +
  ggplot2::labs(title = "B  H0 versus H0+H1 partition agreement",
    subtitle = "Secondary composite at the same label-free k=2",
    x = "Technical seed (suffix)", y = "Adjusted Rand index") + theme_pub(9)
f4 <- (distribution | concordance) +
  patchwork::plot_layout(widths = c(1.15, .85))
f4 <- f4 + patchwork::plot_annotation(
  title = "H1 is low-energy in the composite but remains a distinct topology view",
  caption = "No pair or seed was selected; biological cycles are not inferred.",
  theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 15,
    color = dark), plot.caption = ggplot2::element_text(color = grey, hjust = 0)))
write_csv(ecdf, file.path(output_dir, "data", "f4-h1-ecdf.csv"))
write_csv(hc, file.path(output_dir, "data", "f4-h0-composite-concordance.csv"))
save_pair(f4, "figure-4-h1-contribution-concordance", 13, 7.5,
          "H1 contribution and H0/composite concordance")

# F5: complete label-free stability curves.
stab <- readc(paths[["stability"]])
stab$representation <- factor(stab$representation_id, levels = names(pretty_rep),
  labels = unname(pretty_rep))
f5 <- ggplot2::ggplot(stab,
  ggplot2::aes(k, mean_stability, color = representation_id, group = representation_id)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = pmax(0, mean_stability - jackknife_se),
    ymax = pmin(1, mean_stability + jackknife_se), fill = representation_id),
    alpha = .13, color = NA) +
  ggplot2::geom_line(ggplot2::aes(y = one_se_threshold), color = grey,
    linetype = "33", linewidth = .5) +
  ggplot2::geom_line(linewidth = .65) + ggplot2::geom_point(size = 1.6) +
  ggplot2::geom_point(data = stab[stab$k == stab$selected_k, , drop = FALSE],
    shape = 21, fill = "white", stroke = .8, size = 3,
    show.legend = FALSE) +
  ggplot2::geom_vline(xintercept = 2, linetype = "22", color = dark, linewidth = .45) +
  ggplot2::facet_wrap(~representation, ncol = 3) +
  ggplot2::scale_color_manual(values = rep_cols, guide = "none") +
  ggplot2::scale_fill_manual(values = rep_cols, guide = "none") +
  ggplot2::scale_x_continuous(breaks = 2:10) +
  ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .2),
    labels = scales::label_number(accuracy = .1)) +
  ggplot2::labs(title = "Label-free five-seed cluster stability",
    subtitle = "All six representations selected the smallest k within one SE of maximum stability",
    x = "Candidate cluster count k", y = "Mean pairwise seed ARI",
    caption = "Ribbon: delete-one-seed jackknife SE. Horizontal dotted line: one-SE threshold; vertical dashed line and open point: selected k=2. No metadata labels were used.") +
  theme_pub(9)
write_csv(stab, file.path(output_dir, "data", "f5-complete-stability.csv"))
save_pair(f5, "figure-5-label-free-stability", 12, 8,
          "Label-free stability and selected cluster count")

# F6: all 120 complete descriptive units.
out <- readc(paths[["outcomes"]])
out$representation <- factor(out$representation_id, levels = names(pretty_rep),
  labels = unname(pretty_rep))
out$algorithm <- factor(out$algorithm_role, levels = c("primary", "sensitivity"),
  labels = c("PAM (primary)", "Average linkage (sensitivity)"))
out$metric <- factor(out$metric_id,
  levels = c("adjusted_rand_index", "normalized_mutual_information_max"),
  labels = c("ARI", "max-NMI"))
out$population <- factor(out$population_id,
  levels = c("full124_descriptive", "primary90_context_restriction"),
  labels = c("Full 124 descriptive", "Primary-90 context restriction"))
out$axis <- factor(out$label_axis, levels = c("tissue", "study", "canonical_approach"),
  labels = c("Tissue", "Study", "Approach*"))
out$display <- sprintf("%.3f\n±%.3f", out$seed_mean, out$seed_jackknife_se)
f6 <- ggplot2::ggplot(out, ggplot2::aes(representation, axis, fill = seed_mean)) +
  ggplot2::geom_tile(color = "white", linewidth = .7) +
  ggplot2::geom_text(ggplot2::aes(label = display), size = 2.85, color = dark,
                     lineheight = .82) +
  ggplot2::facet_grid(population + algorithm ~ metric, scales = "free_y", space = "free_y") +
  ggplot2::scale_fill_steps2(low = "#C94C4C", mid = "#F7F7F7", high = navy,
    midpoint = .08, limits = c(-.08, .30),
    breaks = c(-.05, 0, .05, .10, .15, .20, .25),
    labels = scales::label_number(accuracy = .01),
    oob = scales::squish, name = "Seed mean") +
  ggplot2::labs(title = "Complete descriptive cluster–metadata alignment",
    subtitle = "All 120 prespecified units; values are five-seed mean ± technical jackknife SE",
    x = NULL, y = NULL,
    caption = "*Full-124 approach is perfectly nested in substantia nigra/SRA850958; primary-90 approach is structurally non-estimable. No p-values or rankings.") +
  theme_pub(8) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30,
    hjust = 1), strip.text.y = ggplot2::element_text(angle = 0),
    legend.position = "right")
write_csv(out, file.path(output_dir, "data", "f6-complete-outcomes.csv"))
save_pair(f6, "figure-6-descriptive-alignment", 14, 11,
          "Complete descriptive tissue and study alignment")

# F7: all 30 algorithm sensitivities.
alg <- readc(paths[["algorithm"]])
alg$representation <- factor(alg$representation_id, levels = names(pretty_rep),
  labels = unname(pretty_rep))
alg$view <- ifelse(grepl("^cell", alg$representation_id), "Cell view", "Gene view")
alg$seed_label <- factor(alg$seed)
f7 <- ggplot2::ggplot(alg,
  ggplot2::aes(representation, adjusted_rand_index, group = representation)) +
  ggplot2::geom_hline(yintercept = c(0, 1), color = grey, linewidth = .35) +
  ggplot2::geom_boxplot(outlier.shape = NA, width = .5, color = dark, fill = NA,
    linewidth = .45) +
  ggplot2::geom_point(ggplot2::aes(color = seed_label),
    position = ggplot2::position_jitter(width = .09, height = 0,
      seed = 20260816L), size = 2.1) +
  ggplot2::facet_wrap(~view, scales = "free_x", nrow = 1) +
  ggplot2::scale_color_manual(values = c(navy, teal, orange, rust, purple),
    labels = function(x) sub("2026", "", x)) +
  ggplot2::scale_y_continuous(limits = c(-.06, 1.03), breaks = c(0, .25, .5, .75, 1)) +
  ggplot2::labs(title = "PAM versus average-linkage partition sensitivity",
    subtitle = "All representations and five technical seeds at the same label-free selected k=2",
    x = NULL, y = "Adjusted Rand index", color = "Seed suffix",
    caption = "High agreement is representation-specific; neither algorithm is selected as a favorable winner.") +
  theme_pub(9) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25,
    hjust = 1))
write_csv(alg, file.path(output_dir, "data", "f7-algorithm-sensitivity.csv"))
save_pair(f7, "figure-7-algorithm-sensitivity", 12, 7,
          "PAM versus average-linkage sensitivity")

source_manifest <- data.frame(
  contract_id = "mv08a_source_manifest_v1", source_order = seq_along(c(paths, ph_path, h1_path)),
  source_id = c(names(paths), "fixed_ph_artifact", "private_h1_pair_table"),
  locator = c(unname(paths), "private_validated_artifact_not_published",
              "private_validated_artifact_not_published"),
  sha256 = vapply(c(paths, ph_path, h1_path), sha, character(1L)),
  bytes = as.numeric(file.info(c(paths, ph_path, h1_path))$size),
  access_class = c(rep("public_repository", length(paths)),
                   "private_hash_validated", "private_hash_validated"),
  expected_head = expected_head, new_ph = FALSE, new_data = FALSE,
  stringsAsFactors = FALSE)
write_csv(source_manifest, file.path(output_dir, "mv08a-source-manifest.csv"))

figure_manifest <- do.call(rbind, figure_records)
write_csv(figure_manifest, file.path(output_dir, "mv08a-figure-manifest.csv"))
data_files <- sort(list.files(file.path(output_dir, "data"), full.names = TRUE))
data_manifest <- data.frame(
  contract_id = "mv08a_data_manifest_v1", data_order = seq_along(data_files),
  filename = basename(data_files), sha256 = vapply(data_files, sha, character(1L)),
  bytes = as.numeric(file.info(data_files)$size), sample_identifiers_published = FALSE,
  confidential_review_text = FALSE, stringsAsFactors = FALSE)
write_csv(data_manifest, file.path(output_dir, "mv08a-data-manifest.csv"))
provenance <- data.frame(
  contract_id = "mv08a_renderer_provenance_v1", accepted_head = expected_head,
  figures_rendered = 7L, svg_files = 7L, png_files = 7L,
  png_dpi = 300L, figure8_status = "deferred_cross_stage_estimand",
  fixed_figure3_seed = 20260807L,
  fixed_figure3_source = "private_identity_redacted_publicly",
  scientific_distance_grid = "none_exact_or_error_controlled",
  display_grid_only_figure3 = 600L, new_ph = FALSE, new_data = FALSE,
  p_values = FALSE, rankings = FALSE, confidential_review_published = FALSE,
  stringsAsFactors = FALSE)
write_csv(provenance, file.path(output_dir, "mv08a-renderer-provenance.csv"))
message("MV8-A rendered: seven SVG/PNG pairs with complete auditable data bundles")
