# Matched historical cell-side clustering benchmark contracts for MV11.

.mv11_required_seeds <- 20260805:20260809
.mv11_cell_components <- c(cell_H0 = "H0", cell_H1 = "H1")
.mv11_cell_stack_id <- "historical_selectedfit_cell_exact500"

mv11_validate_matrix_bundle_v1 <- function(bundle) {
  required <- c(
    "contract_id", "sample_ids", "seeds", "seed_specific",
    "source_inventory", "outcome_label_state", "biological_outcomes_computed"
  )
  if (!is.list(bundle) || !all(required %in% names(bundle)) ||
      !identical(bundle$contract_id, "mv07i_matrix_bundle_v1") ||
      length(bundle$sample_ids) != 124L || anyNA(bundle$sample_ids) ||
      any(!nzchar(bundle$sample_ids)) || anyDuplicated(bundle$sample_ids) ||
      !identical(as.integer(bundle$seeds), .mv11_required_seeds) ||
      !identical(bundle$outcome_label_state, "closed") ||
      !identical(bundle$biological_outcomes_computed, FALSE) ||
      !all(names(.mv11_cell_components) %in% names(bundle$seed_specific))) {
    stop("MV11 historical matrix bundle violates its outer contract",
         call. = FALSE)
  }
  sample_ids <- as.character(bundle$sample_ids)
  for (component in names(.mv11_cell_components)) {
    matrices <- bundle$seed_specific[[component]]
    if (!is.list(matrices) ||
        !identical(names(matrices), as.character(.mv11_required_seeds))) {
      stop("MV11 cell matrix seed axis drift", call. = FALSE)
    }
    for (matrix in matrices) {
      if (!is.matrix(matrix) || !identical(dim(matrix), c(124L, 124L)) ||
          !identical(rownames(matrix), sample_ids) ||
          !identical(colnames(matrix), sample_ids) || anyNA(matrix) ||
          any(!is.finite(matrix)) || any(matrix < 0) ||
          !isTRUE(all.equal(matrix, t(matrix), tolerance = 1e-12,
                            check.attributes = FALSE)) ||
          any(abs(diag(matrix)) > 1e-12)) {
        stop("MV11 cell distance matrix violates its exact axis/value contract",
             call. = FALSE)
      }
    }
  }
  inventory <- bundle$source_inventory
  inventory_required <- c(
    "seed", "view_id", "homology_dimension", "distances_sha256"
  )
  if (!is.data.frame(inventory) ||
      !all(inventory_required %in% names(inventory))) {
    stop("MV11 source inventory is incomplete", call. = FALSE)
  }
  cell <- inventory[
    inventory$view_id == "cell_topology_v1" &
      inventory$homology_dimension %in% unname(.mv11_cell_components), ,
    drop = FALSE
  ]
  if (nrow(cell) != 10L ||
      any(table(cell$seed, cell$homology_dimension) != 1L) ||
      !identical(sort(unique(as.integer(cell$seed))), .mv11_required_seeds) ||
      any(!grepl("^[0-9a-f]{64}$", cell$distances_sha256))) {
    stop("MV11 cell source inventory violates the ten-matrix contract",
         call. = FALSE)
  }
  invisible(bundle)
}

mv11_cell_catalog_v1 <- function(bundle) {
  mv11_validate_matrix_bundle_v1(bundle)
  inventory <- bundle$source_inventory
  inventory <- inventory[
    inventory$view_id == "cell_topology_v1" &
      inventory$homology_dimension %in% unname(.mv11_cell_components), ,
    drop = FALSE
  ]
  rows <- expand.grid(
    homology_dimension = c("H0", "H1"), seed = .mv11_required_seeds,
    stringsAsFactors = FALSE
  )
  rows <- rows[order(rows$homology_dimension, rows$seed,
                     method = "radix"), , drop = FALSE]
  source_key <- paste(inventory$homology_dimension, inventory$seed, sep = "\r")
  row_key <- paste(rows$homology_dimension, rows$seed, sep = "\r")
  source_index <- match(row_key, source_key)
  if (anyNA(source_index)) stop("MV11 catalog source join failed", call. = FALSE)
  data.frame(
    catalog_order = seq_len(nrow(rows)),
    stack_id = .mv11_cell_stack_id,
    dataset_scope = "internal124",
    representation_id = "existing_selectedfit_data_exact500",
    panel_id = "exact500",
    seed = as.integer(rows$seed),
    view_kind = "cell_topology_v1",
    homology_dimension = rows$homology_dimension,
    units = 124L,
    unordered_pairs = choose(124L, 2L),
    source_distances_sha256 = inventory$distances_sha256[source_index],
    k_grid = "2:10",
    methods = 5L,
    outcome_label_state = "closed",
    biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
}

mv11_cell_matrix_v1 <- function(bundle, seed, homology_dimension) {
  mv11_validate_matrix_bundle_v1(bundle)
  seed <- as.integer(seed)
  homology_dimension <- as.character(homology_dimension)
  if (length(seed) != 1L || is.na(seed) ||
      !seed %in% .mv11_required_seeds ||
      length(homology_dimension) != 1L ||
      !homology_dimension %in% c("H0", "H1")) {
    stop("MV11 cell matrix request violates seed/dimension contract",
         call. = FALSE)
  }
  component <- names(.mv11_cell_components)[
    .mv11_cell_components == homology_dimension
  ]
  matrix <- bundle$seed_specific[[component]][[as.character(seed)]]
  .mv05n_validate_distance_matrix(matrix)
}

mv11_select_primary_k_v1 <- function(assignments) {
  required <- c("stack_id", "homology_dimension", "seed", "method_id", "k",
                "sample_id", "cluster")
  if (!is.data.frame(assignments) || !all(required %in% names(assignments))) {
    stop("MV11 primary-K assignments are incomplete", call. = FALSE)
  }
  .mv10_assert_label_closed(assignments, "MV11 primary-K assignments")
  primary <- assignments[
    assignments$stack_id == .mv11_cell_stack_id &
      assignments$method_id == "pam_dissimilarity_v1", , drop = FALSE
  ]
  if (!setequal(primary$homology_dimension, c("H0", "H1")) ||
      !identical(sort(unique(as.integer(primary$seed))),
                 .mv11_required_seeds)) {
    stop("MV11 primary-K input requires five-seed separate H0 and H1",
         call. = FALSE)
  }
  rows <- lapply(c("H0", "H1"), function(dimension) {
    value <- primary[primary$homology_dimension == dimension, , drop = FALSE]
    if (!identical(sort(unique(as.integer(value$k))), .mv10_k_grid)) {
      stop("MV11 primary-K input lacks the complete K=2:10 grid",
           call. = FALSE)
    }
    selected <- mv05_select_stable_k_v1(value)
    if (!identical(selected$status, "selected") || is.na(selected$selected_k)) {
      stop("MV11 primary-K rule did not select a K", call. = FALSE)
    }
    data.frame(
      contract_id = "mv11_primary_k_selection_v1",
      stack_id = .mv11_cell_stack_id,
      homology_dimension = dimension,
      method_id = "pam_dissimilarity_v1",
      selected_k = selected$selected_k,
      threshold = selected$threshold,
      selection_rule =
        "smallest_k_within_one_SE_of_maximum_five_seed_mean_ARI",
      outcome_label_state = "closed",
      biological_outcomes_computed = FALSE,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}
