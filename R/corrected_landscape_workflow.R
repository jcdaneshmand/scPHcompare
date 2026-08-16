# Additive corrected-landscape workflow artifacts (v1).

.corrected_landscape_control_id_v1 <-
  "scph_corrected_landscape_workflow_control_v1"

.validate_corrected_landscape_control_v1 <- function(control) {
  if (is.null(control)) return(NULL)
  if (!is.list(control) || is.null(names(control)) || any(!nzchar(names(control))) ||
      anyDuplicated(names(control))) {
    stop("corrected_landscape_control must be NULL or a uniquely named list.",
         call. = FALSE)
  }
  allowed <- c("contract_id", "enabled", "max_wall_seconds", "max_pairs",
               "max_rss_bytes", "workers", "existing", "downstream_use",
               "method", "exact_max_intervals", "abs_tol", "rel_tol",
               "subdivisions")
  unknown <- setdiff(names(control), allowed)
  if (length(unknown)) {
    stop("Unknown corrected_landscape_control field(s): ",
         paste(unknown, collapse = ", "), call. = FALSE)
  }
  required <- c("contract_id", "enabled", "max_wall_seconds", "max_pairs")
  if (!all(required %in% names(control))) {
    stop("corrected_landscape_control requires contract_id, enabled, max_wall_seconds, and max_pairs.",
         call. = FALSE)
  }
  scalar_number <- function(value, label, integer = FALSE) {
    if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
        !is.finite(value) || value <= 0 || (integer && value != as.integer(value))) {
      stop(label, " must be one positive ", if (integer) "integer" else "number",
           ".", call. = FALSE)
    }
    if (integer) as.integer(value) else as.numeric(value)
  }
  if (!identical(control$contract_id, .corrected_landscape_control_id_v1) ||
      !identical(control$enabled, TRUE)) {
    stop("Corrected landscape v1 requires its exact contract_id and enabled = TRUE.",
         call. = FALSE)
  }
  fixed <- list(method = "auto", exact_max_intervals = 500L,
                abs_tol = 1e-8, rel_tol = 1e-8, subdivisions = 200L)
  supplied_fixed <- intersect(names(fixed), names(control))
  if (any(!vapply(supplied_fixed, function(name) {
    identical(control[[name]], fixed[[name]])
  }, logical(1L)))) {
    stop("Corrected landscape scientific settings are fixed at auto/500/1e-8/1e-8/200 and cannot be overridden.",
         call. = FALSE)
  }
  result <- list(
    contract_id = .corrected_landscape_control_id_v1,
    enabled = TRUE,
    max_wall_seconds = scalar_number(control$max_wall_seconds,
                                     "max_wall_seconds"),
    max_pairs = scalar_number(control$max_pairs, "max_pairs", integer = TRUE),
    max_rss_bytes = scalar_number(if (is.null(control$max_rss_bytes)) {
      2 * 1024^3
    } else control$max_rss_bytes,
                                  "max_rss_bytes"),
    workers = if (is.null(control$workers)) 1L else control$workers,
    existing = if (is.null(control$existing)) {
      "resume_or_fail"
    } else control$existing,
    downstream_use = if (is.null(control$downstream_use)) {
      "artifacts_only"
    } else control$downstream_use,
    method = fixed$method, exact_max_intervals = fixed$exact_max_intervals,
    abs_tol = fixed$abs_tol, rel_tol = fixed$rel_tol,
    subdivisions = fixed$subdivisions
  )
  if (!identical(as.integer(result$workers), 1L) ||
      !identical(result$existing, "resume_or_fail") ||
      !identical(result$downstream_use, "artifacts_only")) {
    stop("Corrected landscape v1 fixes workers=1, existing='resume_or_fail', and downstream_use='artifacts_only'.",
         call. = FALSE)
  }
  result$workers <- 1L
  if (result$max_rss_bytes < 1.5 * 1024^3) {
    stop("max_rss_bytes must be at least 1.5 GiB for corrected landscape v1.",
         call. = FALSE)
  }
  result
}

.corrected_landscape_resource_plan_v1 <- function(diagrams, control,
                                                   iteration_name) {
  diagrams <- .validate_public_diagram_list_v1(diagrams)
  summaries <- do.call(rbind, lapply(names(diagrams), function(id) {
    .public_diagram_summary_v1(diagrams[[id]], id)
  }))
  rownames(summaries) <- NULL
  if (any(summaries$finite_h0_intervals < 383L) ||
      any(summaries$finite_h0_intervals > 499L) ||
      any(summaries$finite_h1_intervals < 79L) ||
      any(summaries$finite_h1_intervals > 2802L)) {
    stop("Corrected landscape interval counts are outside the profiled H0 383-499 / H1 79-2802 envelope; profiling is required.",
         call. = FALSE)
  }
  pairs <- if (length(diagrams) > 1L) utils::combn(names(diagrams), 2L) else {
    matrix(character(), nrow = 2L, ncol = 0L)
  }
  pair_rows <- if (ncol(pairs)) do.call(rbind, lapply(seq_len(ncol(pairs)), function(i) {
    left <- match(pairs[1L, i], summaries$source_id)
    right <- match(pairs[2L, i], summaries$source_id)
    adaptive <- max(summaries$finite_h0_intervals[c(left, right)]) > 500L ||
      max(summaries$finite_h1_intervals[c(left, right)]) > 500L
    data.frame(
      pair_order = i, first_source_id = pairs[1L, i],
      second_source_id = pairs[2L, i],
      planned_route = if (adaptive) "contains_adaptive_dimension" else "exact_exact",
      planned_wall_seconds = if (adaptive) 240 else 30,
      stringsAsFactors = FALSE
    )
  })) else data.frame(
    pair_order = integer(), first_source_id = character(),
    second_source_id = character(), planned_route = character(),
    planned_wall_seconds = numeric(), stringsAsFactors = FALSE
  )
  planned <- 30 + sum(pair_rows$planned_wall_seconds)
  if (nrow(pair_rows) > control$max_pairs) {
    stop("Corrected landscape pair count exceeds max_pairs.", call. = FALSE)
  }
  if (planned > control$max_wall_seconds) {
    stop("Corrected landscape planned wall time exceeds max_wall_seconds.",
         call. = FALSE)
  }
  list(diagrams = diagrams, summaries = summaries, pairs = pair_rows,
       planned_wall_seconds = planned, iteration_name = iteration_name)
}

.workflow_safe_token_v1 <- function(value) {
  token <- tolower(gsub("[^A-Za-z0-9._-]+", "-", as.character(value)))
  token <- gsub("^-+|-+$", "", token)
  if (!nzchar(token)) stop("Iteration name has no safe artifact token.", call. = FALSE)
  token
}

.atomic_csv_create_or_verify_v1 <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = paste0(".", basename(path), "."),
                        tmpdir = dirname(path))
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  utils::write.csv(value, temporary, row.names = FALSE, na = "", quote = TRUE)
  if (file.exists(path)) {
    if (!identical(digest::digest(path, algo = "sha256", file = TRUE),
                   digest::digest(temporary, algo = "sha256", file = TRUE))) {
      stop("Existing corrected landscape artifact conflicts: ", path,
           call. = FALSE)
    }
    return(invisible(FALSE))
  }
  if (!file.rename(temporary, path)) {
    stop("Atomic corrected landscape CSV rename failed: ", path, call. = FALSE)
  }
  invisible(TRUE)
}

.atomic_rds_create_v1 <- function(value, path) {
  if (file.exists(path)) stop("Refusing to overwrite corrected artifact: ", path,
                              call. = FALSE)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = paste0(".", basename(path), "."),
                        tmpdir = dirname(path))
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  saveRDS(value, temporary, version = 3)
  if (!identical(value, readRDS(temporary))) {
    stop("Corrected landscape temporary RDS failed reload validation.",
         call. = FALSE)
  }
  if (!file.rename(temporary, path)) {
    stop("Atomic corrected landscape RDS rename failed: ", path, call. = FALSE)
  }
  invisible(TRUE)
}

.verify_corrected_pair_shard_v1 <- function(value, first_id, second_id,
                                            first_hash = NULL,
                                            second_hash = NULL) {
  if (!inherits(value, .scph_landscape_pair_class_v1) ||
      !identical(value$mode, "scientific") ||
      !identical(value$specification, .scph_landscape_contract_id_v1) ||
      !identical(value$provenance$first_source_id, first_id) ||
      !identical(value$provenance$second_source_id, second_id) ||
      !all(vapply(value$dimensions, `[[`, logical(1L),
                  "within_requested_tolerance"))) {
    stop("Existing corrected landscape pair shard failed contract validation.",
         call. = FALSE)
  }
  expected <- paste0("scph_landscape_distance_v1:", digest::digest(
    .landscape_result_identity_v1(value), algo = "sha256", serialize = TRUE
  ))
  if (!identical(value$cache_key, expected)) {
    stop("Existing corrected landscape pair shard has a cache-key conflict.",
         call. = FALSE)
  }
  if (!is.null(first_hash) &&
      !identical(value$provenance$first_diagram_sha256, first_hash)) {
    stop("Existing corrected landscape first-diagram hash conflicts.",
         call. = FALSE)
  }
  if (!is.null(second_hash) &&
      !identical(value$provenance$second_diagram_sha256, second_hash)) {
    stop("Existing corrected landscape second-diagram hash conflicts.",
         call. = FALSE)
  }
  invisible(TRUE)
}

.assemble_corrected_landscape_matrix_v1 <- function(diagrams, pair_objects) {
  ids <- names(diagrams)
  matrices <- setNames(lapply(c("H0", "H1", "combined"), function(unused) {
    matrix(0, nrow = length(ids), ncol = length(ids), dimnames = list(ids, ids))
  }), c("H0", "H1", "combined"))
  pair_rows <- vector("list", length(pair_objects))
  pair_cache_keys <- character(length(pair_objects))
  for (index in seq_along(pair_objects)) {
    pair <- pair_objects[[index]]
    first <- pair$provenance$first_source_id
    second <- pair$provenance$second_source_id
    for (dimension in names(matrices)) {
      matrices[[dimension]][first, second] <- pair$distances[[dimension]]
      matrices[[dimension]][second, first] <- pair$distances[[dimension]]
    }
    pair_rows[[index]] <- data.frame(
      first_source_id = first, second_source_id = second,
      h0_distance = unname(pair$distances[["H0"]]),
      h1_distance = unname(pair$distances[["H1"]]),
      combined_distance = unname(pair$distances[["combined"]]),
      h1_squared_distance_fraction = pair$h1_squared_distance_fraction,
      h0_method = pair$dimensions$H0$method,
      h1_method = pair$dimensions$H1$method,
      h0_error_estimate = pair$dimensions$H0$achieved_absolute_error_estimate,
      h1_error_estimate = pair$dimensions$H1$achieved_absolute_error_estimate,
      cache_key = pair$cache_key, stringsAsFactors = FALSE
    )
    pair_cache_keys[[index]] <- pair$cache_key
  }
  pair_diagnostics <- if (length(pair_rows)) do.call(rbind, pair_rows) else {
    data.frame(first_source_id = character(), second_source_id = character(),
      h0_distance = numeric(), h1_distance = numeric(), combined_distance = numeric(),
      h1_squared_distance_fraction = numeric(), h0_method = character(),
      h1_method = character(), h0_error_estimate = numeric(),
      h1_error_estimate = numeric(), cache_key = character(), stringsAsFactors = FALSE)
  }
  diagram_provenance <- do.call(rbind, lapply(ids, function(id) {
    .public_diagram_summary_v1(diagrams[[id]], id)
  }))
  rownames(diagram_provenance) <- NULL
  identity <- list(
    contract_id = "scph_public_landscape_matrix_v1",
    specification = .scph_landscape_contract_id_v1, mode = "scientific",
    sample_ids = ids, diagram_sha256 = diagram_provenance$diagram_sha256,
    pair_cache_keys = pair_cache_keys, matrices = matrices
  )
  structure(list(
    contract_id = identity$contract_id, result_version = 1L,
    specification = identity$specification, mode = "scientific",
    sample_ids = ids, matrices = matrices, pair_diagnostics = pair_diagnostics,
    diagram_provenance = diagram_provenance,
    provenance = list(
      api = "persistence_landscape_distance_matrix", api_version = 1L,
      scientific_contract = .scph_landscape_contract_id_v1,
      method_requested = "auto", exact_max_intervals = 500L,
      absolute_tolerance = 1e-8, relative_tolerance = 1e-8,
      subdivisions = 200L, sample_order_policy = "lexical radix order",
      pair_count = nrow(pair_diagnostics),
      dimension_policy = "H0 and H1 separate; Euclidean combination secondary",
      legacy_reproduction = FALSE, existing_workflow_default_changed = FALSE,
      workflow_artifact_contract = "scph_corrected_landscape_artifact_v1"
    ),
    runtime = list(elapsed_seconds = NA_real_),
    cache_key = paste0("scph_landscape_distance_matrix_v1:", digest::digest(
      identity, algo = "sha256", serialize = TRUE
    ))
  ), class = c(.scph_landscape_matrix_class_v1, "list"))
}

.verify_completion_v1 <- function(artifact_dir, completion) {
  required <- c("resource-plan-v1.csv", "input-manifest-v1.csv",
                "pair-index-v1.csv", "distance-matrix-v1.rds",
                "provenance-v1.csv")
  if (nrow(completion) != length(required) ||
      !setequal(completion$artifact, required)) {
    stop("Corrected landscape completion manifest is incomplete.", call. = FALSE)
  }
  ok <- vapply(seq_len(nrow(completion)), function(index) {
    path <- file.path(artifact_dir, completion$artifact[[index]])
    file.exists(path) && identical(
      digest::digest(path, algo = "sha256", file = TRUE),
      completion$sha256[[index]]
    ) && identical(unname(as.numeric(file.info(path)$size)),
                   as.numeric(completion$bytes[[index]]))
  }, logical(1L))
  if (!all(ok)) stop("Corrected landscape completion hash validation failed.",
                     call. = FALSE)
  pair_index <- utils::read.csv(file.path(artifact_dir, "pair-index-v1.csv"),
                                stringsAsFactors = FALSE, check.names = FALSE)
  pair_ok <- vapply(seq_len(nrow(pair_index)), function(index) {
    path <- file.path(artifact_dir, pair_index$pair_artifact[[index]])
    file.exists(path) && identical(
      digest::digest(path, algo = "sha256", file = TRUE),
      pair_index$pair_sha256[[index]]
    ) && identical(unname(as.numeric(file.info(path)$size)),
                   as.numeric(pair_index$pair_bytes[[index]]))
  }, logical(1L))
  if (!all(pair_ok)) stop("Corrected landscape pair-index hash validation failed.",
                          call. = FALSE)
  invisible(TRUE)
}

produce_corrected_landscape_artifacts_v1 <- function(
    pd_list_path, iteration_name, results_dir, control, log_message = message,
    stop_after_pairs = Inf) {
  control <- .validate_corrected_landscape_control_v1(control)
  if (is.null(control)) return(NULL)
  if (!is.character(pd_list_path) || length(pd_list_path) != 1L ||
      !file.exists(pd_list_path)) {
    stop("Corrected landscape v1 requires an existing pd_list_path.", call. = FALSE)
  }
  diagrams <- readRDS(pd_list_path)
  plan <- .corrected_landscape_resource_plan_v1(diagrams, control, iteration_name)
  diagrams <- plan$diagrams
  input_identity <- list(
    contract_id = "scph_corrected_landscape_input_set_v1",
    control = control[c("contract_id", "method", "exact_max_intervals",
                        "abs_tol", "rel_tol", "subdivisions")],
    sample_ids = names(diagrams),
    diagram_sha256 = plan$summaries$diagram_sha256,
    scientific_contract = .scph_landscape_contract_id_v1,
    engine_version = "landscape_reference_v3"
  )
  input_sha <- digest::digest(input_identity, algo = "sha256", serialize = TRUE)
  artifact_dir <- file.path(results_dir, "corrected_landscape_v1", paste0(
    .workflow_safe_token_v1(iteration_name), "--", input_sha
  ))
  dir.create(file.path(artifact_dir, "pairs"), recursive = TRUE,
             showWarnings = FALSE)
  resource_plan <- data.frame(
    contract_id = "scph_corrected_landscape_resource_plan_v1",
    input_set_sha256 = input_sha, iteration_name = iteration_name,
    samples = length(diagrams), pairs = nrow(plan$pairs),
    exact_exact_pairs = sum(plan$pairs$planned_route == "exact_exact"),
    adaptive_pairs = sum(plan$pairs$planned_route == "contains_adaptive_dimension"),
    planned_wall_seconds = plan$planned_wall_seconds,
    max_wall_seconds = control$max_wall_seconds,
    max_pairs = control$max_pairs, max_rss_bytes = control$max_rss_bytes,
    workers = 1L, admitted = TRUE, stringsAsFactors = FALSE
  )
  input_manifest <- cbind(
    data.frame(contract_id = "scph_corrected_landscape_input_manifest_v1",
               input_set_sha256 = input_sha, sample_order = seq_len(nrow(plan$summaries)),
               stringsAsFactors = FALSE), plan$summaries
  )
  .atomic_csv_create_or_verify_v1(resource_plan,
                                  file.path(artifact_dir, "resource-plan-v1.csv"))
  .atomic_csv_create_or_verify_v1(input_manifest,
                                  file.path(artifact_dir, "input-manifest-v1.csv"))

  completion_path <- file.path(artifact_dir, "completion-v1.csv")
  if (file.exists(completion_path)) {
    completion <- utils::read.csv(completion_path, stringsAsFactors = FALSE,
                                  check.names = FALSE)
    .verify_completion_v1(artifact_dir, completion)
    matrix_object <- readRDS(file.path(artifact_dir, "distance-matrix-v1.rds"))
    if (!inherits(matrix_object, .scph_landscape_matrix_class_v1)) {
      stop("Completed corrected landscape matrix has the wrong class.", call. = FALSE)
    }
    return(list(artifact_dir = artifact_dir, input_set_sha256 = input_sha,
                matrix_cache_key = matrix_object$cache_key,
                pair_count = nrow(plan$pairs), resumed = TRUE,
                downstream_use = "artifacts_only"))
  }

  pair_objects <- vector("list", nrow(plan$pairs))
  pair_index <- vector("list", nrow(plan$pairs))
  for (index in seq_len(nrow(plan$pairs))) {
    row <- plan$pairs[index, ]
    first <- row$first_source_id
    second <- row$second_source_id
    first_hash <- plan$summaries$diagram_sha256[
      match(first, plan$summaries$source_id)]
    second_hash <- plan$summaries$diagram_sha256[
      match(second, plan$summaries$source_id)]
    pattern <- sprintf("^pair-%06d--[0-9a-f]{64}\\.rds$", index)
    existing_paths <- list.files(file.path(artifact_dir, "pairs"),
                                 pattern = pattern, full.names = TRUE)
    if (length(existing_paths) > 1L) {
      stop("Multiple corrected landscape shards claim one pair order.",
           call. = FALSE)
    }
    if (length(existing_paths) == 1L) {
      pair_path <- existing_paths[[1L]]
      value <- readRDS(pair_path)
      .verify_corrected_pair_shard_v1(value, first, second,
                                      first_hash, second_hash)
    } else {
      value <- persistence_landscape_distance(
        diagrams[[first]], diagrams[[second]], method = "auto",
        exact_max_intervals = 500L, abs_tol = 1e-8, rel_tol = 1e-8,
        subdivisions = 200L, mode = "scientific",
        first_id = first, second_id = second
      )
      value$runtime$elapsed_seconds <- NA_real_
      cache_sha <- sub("^scph_landscape_distance_v1:", "", value$cache_key)
      pair_path <- file.path(artifact_dir, "pairs", sprintf(
        "pair-%06d--%s.rds", index, cache_sha
      ))
      .atomic_rds_create_v1(value, pair_path)
      value <- readRDS(pair_path)
      .verify_corrected_pair_shard_v1(value, first, second,
                                      first_hash, second_hash)
    }
    relative_path <- file.path("pairs", basename(pair_path))
    pair_objects[[index]] <- value
    pair_index[[index]] <- data.frame(
      contract_id = "scph_corrected_landscape_pair_index_v1",
      input_set_sha256 = input_sha, pair_order = index,
      first_source_id = first, second_source_id = second,
      pair_artifact = gsub("\\\\", "/", relative_path),
      pair_cache_key = value$cache_key,
      pair_sha256 = digest::digest(pair_path, algo = "sha256", file = TRUE),
      pair_bytes = as.numeric(file.info(pair_path)$size),
      h0_method = value$dimensions$H0$method,
      h1_method = value$dimensions$H1$method,
      h0_error_estimate = value$dimensions$H0$achieved_absolute_error_estimate,
      h1_error_estimate = value$dimensions$H1$achieved_absolute_error_estimate,
      h0_certified = value$dimensions$H0$within_requested_tolerance,
      h1_certified = value$dimensions$H1$within_requested_tolerance,
      stringsAsFactors = FALSE
    )
    if (is.finite(stop_after_pairs) && index >= stop_after_pairs) {
      stop("MV5-AS intentional interruption after verified pair shard.",
           call. = FALSE)
    }
  }
  pair_index <- if (length(pair_index)) do.call(rbind, pair_index) else {
    data.frame(contract_id = character(), input_set_sha256 = character(),
      pair_order = integer(), first_source_id = character(),
      second_source_id = character(), pair_artifact = character(),
      pair_cache_key = character(), pair_sha256 = character(), pair_bytes = numeric(),
      h0_method = character(), h1_method = character(), h0_error_estimate = numeric(),
      h1_error_estimate = numeric(), h0_certified = logical(), h1_certified = logical(),
      stringsAsFactors = FALSE)
  }
  .atomic_csv_create_or_verify_v1(pair_index,
                                  file.path(artifact_dir, "pair-index-v1.csv"))
  matrix_object <- .assemble_corrected_landscape_matrix_v1(diagrams, pair_objects)
  matrix_path <- file.path(artifact_dir, "distance-matrix-v1.rds")
  if (file.exists(matrix_path)) {
    if (!identical(matrix_object, readRDS(matrix_path))) {
      stop("Existing corrected landscape matrix conflicts with verified shards.",
           call. = FALSE)
    }
  } else {
    .atomic_rds_create_v1(matrix_object, matrix_path)
  }
  if (!identical(matrix_object, readRDS(matrix_path))) {
    stop("Corrected landscape matrix failed serialization validation.",
         call. = FALSE)
  }
  provenance <- data.frame(
    contract_id = "scph_corrected_landscape_provenance_v1",
    input_set_sha256 = input_sha, matrix_cache_key = matrix_object$cache_key,
    scientific_contract = .scph_landscape_contract_id_v1,
    engine_version = "landscape_reference_v3", method = "auto",
    exact_max_intervals = 500L, abs_tol = 1e-8, rel_tol = 1e-8,
    subdivisions = 200L, level_policy = "all_consecutive_active_levels",
    dimension_policy = "H0_H1_separate_combined_descriptive",
    downstream_use = "artifacts_only", legacy_artifacts_changed = FALSE,
    stringsAsFactors = FALSE
  )
  .atomic_csv_create_or_verify_v1(provenance,
                                  file.path(artifact_dir, "provenance-v1.csv"))
  completion_artifacts <- c("resource-plan-v1.csv", "input-manifest-v1.csv",
                            "pair-index-v1.csv", "distance-matrix-v1.rds",
                            "provenance-v1.csv")
  completion <- data.frame(
    contract_id = "scph_corrected_landscape_completion_v1",
    input_set_sha256 = input_sha, artifact = completion_artifacts,
    sha256 = vapply(completion_artifacts, function(path) digest::digest(
      file.path(artifact_dir, path), algo = "sha256", file = TRUE
    ), character(1)),
    bytes = vapply(completion_artifacts, function(path) as.numeric(file.info(
      file.path(artifact_dir, path))$size), numeric(1)),
    stringsAsFactors = FALSE
  )
  .atomic_csv_create_or_verify_v1(completion, completion_path)
  .verify_completion_v1(artifact_dir, completion)
  log_message(sprintf("Completed additive corrected landscape artifacts for %s at %s",
                      iteration_name, artifact_dir))
  list(artifact_dir = artifact_dir, input_set_sha256 = input_sha,
       matrix_cache_key = matrix_object$cache_key, pair_count = nrow(plan$pairs),
       resumed = FALSE, downstream_use = "artifacts_only")
}
