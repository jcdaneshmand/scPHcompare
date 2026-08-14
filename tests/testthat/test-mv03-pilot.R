test_that("MV-03 gene IDs and technical categories are versioned explicitly", {
  ids <- c(
    "A1CF-ENSG00000148584", "MT-CO1-ENSG00000198888",
    "RPL3-ENSG00000100316", "HBA1-ENSG00000206172",
    "HBP1-ENSG00000105856"
  )
  genes <- canonical_mv03_gene_ids(ids)
  expect_identical(genes, c("A1CF", "MT-CO1", "RPL3", "HBA1", "HBP1"))
  expect_identical(
    mv03_feature_category(ids),
    c(
      "retained_candidate", "mitochondrial_^MT-",
      "ribosomal_protein_^RP[SL]", "hemoglobin_^HB(?!P)",
      "retained_candidate"
    )
  )
})

test_that("MV-03 panel fitting is label-free, paired, and deterministic", {
  genes <- c(
    "A-ENSG000001", "B-ENSG000002", "C-ENSG000003",
    "MT-X-ENSG000004"
  )
  cells <- paste0("cell", 1:5)
  make_matrix <- function(offset) {
    value <- rbind(
      A = c(1, 2, 3, 4, 5) + offset,
      B = c(1, 3, 2, 5, 4) + offset,
      C = c(5, 4, 3, 2, 1) + offset,
      `MT-X` = c(2, 4, 6, 8, 10) + offset
    )
    rownames(value) <- genes
    colnames(value) <- cells
    value
  }
  residual <- list(s1 = make_matrix(0), s2 = make_matrix(1))
  integrated <- list(s1 = make_matrix(2), s2 = make_matrix(3))
  first <- fit_mv03_gene_panel(
    residual, integrated, "fixture", "descriptive", panel_size = 2L
  )
  second <- fit_mv03_gene_panel(
    residual, integrated, "fixture", "descriptive", panel_size = 2L
  )
  expect_identical(first$panel, second$panel)
  expect_identical(first$panel_sha256, second$panel_sha256)
  expect_length(first$panel, 2L)
  expect_false("MT-X" %in% first$panel)
  expect_true(all(first$panel %in% c("A", "B", "C")))
})

test_that("MV-03 source preparation uses equal matched contributions", {
  genes <- paste0("gene", seq_len(500L))
  cells <- paste0("cell", seq_len(384L))
  first <- outer(seq_len(500L), seq_len(384L), function(gene, cell) {
    sin(gene / 17) + cos(cell / 13) + gene * cell / 1e6
  })
  second <- first * 1.7 + 2
  dimnames(first) <- list(genes, cells)
  dimnames(second) <- list(genes, cells)
  selected <- list(s1 = cells, s2 = cells)
  prepared <- prepare_mv03_sources(
    list(s1 = first, s2 = second), genes, "fixture", "sct_whole",
    "descriptive", 20260805L, selected_cells = selected
  )
  pooled <- do.call(cbind, lapply(prepared$sources, `[[`, "matrix"))
  expect_equal(unname(rowMeans(pooled)), rep(0, 500L), tolerance = 1e-10)
  expect_equal(unname(apply(pooled, 1L, stats::sd)), rep(1, 500L),
               tolerance = 1e-10)
  expect_true(all(vapply(
    prepared$sources, inherits, logical(1), "scph_dual_view_source_v1"
  )))
})
