consumer_diagram <- function(shift = 0) {
  h0 <- cbind(0, seq(0, 3.82, length.out = 383) + shift,
              seq(0, 3.82, length.out = 383) + shift + seq(0.3, 0.7, length.out = 383))
  h1 <- cbind(1, seq(0, 0.78, length.out = 79) + shift,
              seq(0, 0.78, length.out = 79) + shift + seq(0.1, 0.25, length.out = 79))
  x <- rbind(h0, h1); colnames(x) <- c("dimension", "birth", "death"); x
}

make_consumer_sidecar <- function() {
  root <- tempfile("consumer-sidecar-"); dir.create(root)
  path <- file.path(root, "pd.rds")
  saveRDS(list(c = consumer_diagram(.08), a = consumer_diagram(0),
               b = consumer_diagram(.03)), path)
  sidecar <- produce_corrected_landscape_artifacts_v1(path, "consumer test", root,
    list(contract_id = "scph_corrected_landscape_workflow_control_v1",
      enabled = TRUE, max_wall_seconds = 120, max_pairs = 3L),
    log_message = function(...) NULL)
  list(root = root, sidecar = sidecar)
}

test_that("verified consumer loader preserves separate matrices and identity", {
  fixture <- make_consumer_sidecar()
  before <- list.files(fixture$sidecar$artifact_dir, recursive = TRUE, full.names = TRUE)
  hashes <- vapply(before, digest::digest, character(1), algo = "sha256", file = TRUE)
  bundle <- read_corrected_landscape_bundle(fixture$sidecar, "iter-1", "cell_topology_v1")
  expect_s3_class(bundle, "scph_corrected_landscape_consumer_bundle_v1")
  expect_identical(names(bundle$matrices), c("H0", "H1", "combined"))
  expect_identical(bundle$view_id, "cell_topology_v1")
  expect_match(bundle$cache_key, "^scph_corrected_consumer_bundle_v1:[0-9a-f]{64}$")
  expect_false(identical(bundle$matrices$H0, bundle$matrices$H1))
  expect_equal(bundle$matrices$combined,
               sqrt(bundle$matrices$H0^2 + bundle$matrices$H1^2), tolerance = 1e-14)
  expect_identical(before, list.files(fixture$sidecar$artifact_dir,
    recursive = TRUE, full.names = TRUE))
  expect_identical(hashes, vapply(before, digest::digest, character(1),
    algo = "sha256", file = TRUE))
  expect_identical(bundle, readRDS({p <- tempfile(fileext = ".rds"); saveRDS(bundle,p); p}))
})

test_that("loader rejects implicit identity and sidecar mismatches", {
  fixture <- make_consumer_sidecar()
  expect_error(read_corrected_landscape_bundle(fixture$sidecar, "", "cell_topology_v1"),
               "iteration_id")
  expect_error(read_corrected_landscape_bundle(fixture$sidecar, "iter", "unknown"),
               "view_id")
  bad <- fixture$sidecar; bad$input_set_sha256 <- paste(rep("0",64), collapse="")
  expect_error(read_corrected_landscape_bundle(bad, "iter", "gene_topology_v1"),
               "input-set key")
  bad <- fixture$sidecar; bad$matrix_cache_key <- "wrong"
  expect_error(read_corrected_landscape_bundle(bad, "iter", "gene_topology_v1"),
               "matrix cache key")
  expect_error(read_corrected_landscape_bundle(list(), "iter", "cell_topology_v1"),
               "lacks required")
})

test_that("average trees are label-free, separate, deterministic, and equivariant", {
  fixture <- make_consumer_sidecar()
  bundle <- read_corrected_landscape_bundle(fixture$sidecar, "iter", "gene_topology_v1")
  result <- corrected_landscape_average_trees(bundle)
  repeated <- corrected_landscape_average_trees(bundle)
  expect_s3_class(result, "scph_corrected_landscape_average_trees_v1")
  expect_identical(result, repeated)
  expect_identical(names(result$trees), c("H0", "H1"))
  expect_null(result$combined_tree); expect_null(result$partitions); expect_null(result$selected_k)
  expect_true(result$provenance$label_free); expect_false(result$provenance$combined_consumed)
  expect_true(all(vapply(result$trees, inherits, logical(1), "hclust")))
  order <- rev(bundle$sample_ids)
  permuted <- bundle
  permuted$matrices <- lapply(bundle$matrices, function(x) x[order, order, drop = FALSE])
  permuted$sample_ids <- order
  permuted$cache_key <- paste0(bundle$cache_key, ":permuted-test")
  ptrees <- corrected_landscape_average_trees(permuted)$trees
  for (dimension in c("H0", "H1")) {
    expect_equal(as.matrix(stats::cophenetic(result$trees[[dimension]]))[order, order],
                 as.matrix(stats::cophenetic(ptrees[[dimension]])), tolerance = 1e-14)
  }
  expect_error(corrected_landscape_average_trees(list()), "not a valid")
})
