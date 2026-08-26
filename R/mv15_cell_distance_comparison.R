# Helpers for prospective label-closed MV15 cell-distance comparisons.

mv15_read_bound_stack_v1 <- function(binding, cell_root, gene_root) {
  required <- c(
    "source_kind", "source_group_order", "dataset_scope", "panel_id",
    "seed", "homology_dimension", "units", "unordered_pairs"
  )
  if (!is.data.frame(binding) || nrow(binding) != 1L ||
      !all(required %in% names(binding))) {
    stop("MV15 stack binding must contain one complete row", call. = FALSE)
  }
  source_kind <- as.character(binding$source_kind)
  order <- as.integer(binding$source_group_order)
  root <- if (source_kind == "cell_topology_v2") cell_root else gene_root
  if (!source_kind %in% c("cell_topology_v2", "gene_topology_v1") ||
      is.na(order) || order < 1L || !dir.exists(root)) {
    stop("MV15 stack source binding is invalid", call. = FALSE)
  }
  group <- sprintf("group_%02d", order)
  files <- list.files(
    file.path(root, "production", group), pattern = "^distances[.]csv$",
    recursive = TRUE, full.names = TRUE
  )
  files <- sort(files, method = "radix")
  if (!length(files) || !all(file.exists(files))) {
    stop("MV15 bound distance payload is absent for ", source_kind,
         " group ", order, call. = FALSE)
  }
  payloads <- lapply(files, utils::read.csv, stringsAsFactors = FALSE,
                     check.names = FALSE)
  value <- do.call(rbind, payloads)
  if (!all(c("first_unit_id", "second_unit_id", "distance",
             "homology_dimension") %in% names(value)) ||
      any(as.character(value$homology_dimension) !=
            as.character(binding$homology_dimension))) {
    stop("MV15 distance payload schema or dimension drift", call. = FALSE)
  }
  canonical <- .mv08zy_validate_pairs(value)
  if (nrow(canonical) != as.integer(binding$unordered_pairs) ||
      length(unique(c(canonical$first_unit_id,
                      canonical$second_unit_id))) != as.integer(binding$units)) {
    stop("MV15 distance payload cardinality drift", call. = FALSE)
  }
  file_manifest <- data.frame(
    file_order = seq_along(files),
    bytes = as.numeric(file.info(files)$size),
    sha256 = vapply(files, function(path) digest::digest(
      file = path, algo = "sha256", serialize = FALSE
    ), character(1L)),
    stringsAsFactors = FALSE
  )
  list(
    pairs = canonical,
    file_manifest = file_manifest,
    payload_set_sha256 = digest::digest(
      paste(file_manifest$sha256, collapse = "\n"), algo = "sha256",
      serialize = FALSE
    ),
    pair_axis_sha256 = digest::digest(
      paste(canonical$pair_key, collapse = "\n"), algo = "sha256",
      serialize = FALSE
    )
  )
}

mv15_compare_distance_pairs_v1 <- function(left, right, comparison_id,
                                             neighbor_k) {
  left <- .mv08zy_validate_pairs(left)
  right <- .mv08zy_validate_pairs(right)
  if (!identical(left[c("first_unit_id", "second_unit_id", "pair_key")],
                 right[c("first_unit_id", "second_unit_id", "pair_key")])) {
    stop("MV15 comparison pair axes differ", call. = FALSE)
  }
  comparison_id <- as.character(comparison_id)
  units <- sort(unique(c(left$first_unit_id, left$second_unit_id)),
                method = "radix")
  neighbor_k <- sort(unique(as.integer(neighbor_k)))
  if (length(comparison_id) != 1L || is.na(comparison_id) ||
      !nzchar(comparison_id) || !length(neighbor_k) || anyNA(neighbor_k) ||
      any(neighbor_k < 1L) || any(neighbor_k >= length(units))) {
    stop("MV15 comparison identifiers or neighborhood sizes are invalid",
         call. = FALSE)
  }
  x <- left$distance
  y <- right$distance
  left_median <- stats::median(x)
  right_median <- stats::median(y)
  if (!is.finite(left_median) || !is.finite(right_median) ||
      left_median <= sqrt(.Machine$double.eps) ||
      right_median <= sqrt(.Machine$double.eps) ||
      sum(x ^ 2) <= .Machine$double.eps ||
      sum(y ^ 2) <= .Machine$double.eps ||
      stats::sd(x) <= sqrt(.Machine$double.eps) ||
      stats::sd(y) <= sqrt(.Machine$double.eps)) {
    stop("MV15 comparison distance stack is degenerate", call. = FALSE)
  }
  scale <- max(0, sum(x * y) / sum(x ^ 2))
  stress <- sqrt(sum((y - scale * x) ^ 2) / sum(y ^ 2))
  scaled_change <- abs(y / right_median - x / left_median)
  summary <- data.frame(
    contract_id = "mv15_distance_comparison_summary_v1",
    comparison_id = comparison_id, units = length(units),
    unordered_pairs = nrow(left),
    distance_transform = "sqrt_exact_streamed_squared_L2",
    pearson = stats::cor(x, y, method = "pearson"),
    spearman = stats::cor(x, y, method = "spearman"),
    nonnegative_left_to_right_scale = scale,
    relative_left_to_right_stress = stress,
    left_median_distance = left_median,
    right_median_distance = right_median,
    median_abs_median_scaled_change = stats::median(scaled_change),
    p95_abs_median_scaled_change = as.numeric(stats::quantile(
      scaled_change, 0.95, names = FALSE, type = 7
    )),
    interpretation = "descriptive_no_equivalence_ranking_or_biological_claim",
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    clustering_jobs = 0L, fusion_jobs = 0L, label_jobs = 0L,
    outcome_jobs = 0L, manuscript_claim_jobs = 0L,
    stringsAsFactors = FALSE
  )
  neighbor_rows <- list()
  neighbor_summary <- list()
  for (i in seq_along(neighbor_k)) {
    k <- neighbor_k[[i]]
    left_sets <- .mv08zy_neighbor_sets(left, k)
    right_sets <- .mv08zy_neighbor_sets(right, k)
    overlap <- vapply(units, function(unit) {
      shared <- intersect(left_sets[[unit]], right_sets[[unit]])
      combined <- union(left_sets[[unit]], right_sets[[unit]])
      length(shared) / length(combined)
    }, numeric(1L))
    neighbor_rows[[i]] <- data.frame(
      contract_id = "mv15_unit_neighbor_overlap_v1",
      comparison_id = comparison_id, unit_id = units, k = k,
      neighbor_jaccard = overlap, stringsAsFactors = FALSE
    )
    neighbor_summary[[i]] <- data.frame(
      contract_id = "mv15_neighbor_summary_v1",
      comparison_id = comparison_id, units = length(units), k = k,
      mean_neighbor_jaccard = mean(overlap),
      median_neighbor_jaccard = stats::median(overlap),
      p10_neighbor_jaccard = as.numeric(stats::quantile(
        overlap, 0.10, names = FALSE, type = 7
      )), stringsAsFactors = FALSE
    )
  }
  list(
    summary = summary,
    neighbor_summary = do.call(rbind, neighbor_summary),
    neighbor = do.call(rbind, neighbor_rows),
    pair_axis = left[c("first_unit_id", "second_unit_id", "pair_key")]
  )
}
