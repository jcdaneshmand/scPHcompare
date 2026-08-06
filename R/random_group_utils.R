# Assign sample identities to reproducible, approximately balanced null groups.
assignRandomGroup <- function(seurat_obj, k, new_col_name = "Random_Group",
                              seed = 123, num_bootstraps = 1) {
  if (!methods::is(seurat_obj, "Seurat")) {
    stop("seurat_obj must be a Seurat object.", call. = FALSE)
  }
  if (length(k) != 1L || is.na(k) || k < 1L || k != as.integer(k)) {
    stop("k must be a positive integer.", call. = FALSE)
  }
  if (length(num_bootstraps) != 1L || is.na(num_bootstraps) ||
      num_bootstraps < 1L || num_bootstraps != as.integer(num_bootstraps)) {
    stop("num_bootstraps must be a positive integer.", call. = FALSE)
  }
  if (length(seed) != 1L || is.na(seed) || !is.numeric(seed)) {
    stop("seed must be one finite numeric value.", call. = FALSE)
  }
  if (!is.finite(seed)) {
    stop("seed must be one finite numeric value.", call. = FALSE)
  }
  if (!"orig.ident" %in% colnames(seurat_obj@meta.data)) {
    stop("The Seurat object does not contain a column named 'orig.ident'.", call. = FALSE)
  }

  identities <- as.character(seurat_obj@meta.data$orig.ident)
  if (anyNA(identities) || any(!nzchar(identities))) {
    stop("orig.ident contains missing or empty sample identities.", call. = FALSE)
  }
  unique_ids <- unique(identities)

  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) {
    previous_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    if (had_seed) {
      assign(".Random.seed", previous_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  for (i in seq_len(as.integer(num_bootstraps))) {
    set.seed(as.integer(seed) + i)
    random_groups <- sample(rep(seq_len(as.integer(k)), length.out = length(unique_ids)))
    group_mapping <- stats::setNames(random_groups, unique_ids)
    current_col_name <- if (num_bootstraps == 1L) {
      new_col_name
    } else {
      paste0(new_col_name, "_bootstrap_", i)
    }
    seurat_obj@meta.data[[current_col_name]] <- unname(group_mapping[identities])
  }

  seurat_obj
}
