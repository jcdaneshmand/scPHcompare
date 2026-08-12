# MV5-AP realistic public-landscape compatibility gate helpers.

mv05ap_select_depth_triplets_v1 <- function(manifest) {
  required <- c(
    "diagram_id", "stratum_id", "cohort", "representation", "sample_id",
    "view_id", "h0_finite_intervals", "h1_finite_intervals",
    "diagram_sha256", "result_file_sha256", "result_file"
  )
  if (!is.data.frame(manifest) || !all(required %in% names(manifest)) ||
      anyNA(manifest[, required]) || anyDuplicated(manifest$diagram_id)) {
    stop("MV5-AP requires one complete row per unique accepted diagram.",
         call. = FALSE)
  }
  groups <- split(manifest, manifest$stratum_id)
  group_sizes <- vapply(groups, nrow, integer(1L))
  if (length(groups) < 2L || any(group_sizes < 3L)) {
    stop("Every MV5-AP stratum must contain at least three diagrams.",
         call. = FALSE)
  }
  selected <- lapply(groups, function(group) {
    group <- group[order(
      group$h1_finite_intervals, group$sample_id, group$diagram_id,
      method = "radix"
    ), , drop = FALSE]
    indices <- c(1L, ceiling(nrow(group) / 2), nrow(group))
    roles <- c("minimum_h1_depth", "middle_order_h1_depth", "maximum_h1_depth")
    result <- group[indices, , drop = FALSE]
    result$selection_role <- roles
    result$selection_rule <- paste0(
      "within_stratum_order_h1_then_sample_then_diagram_take_1_ceiling_n2_n"
    )
    result
  })
  selected <- do.call(rbind, selected)
  rownames(selected) <- NULL
  role_order <- match(
    selected$selection_role,
    c("minimum_h1_depth", "middle_order_h1_depth", "maximum_h1_depth")
  )
  selected[order(selected$stratum_id, role_order, method = "radix"),
           , drop = FALSE]
}

mv05ap_decide_v1 <- function(provenance_verified, default_exact_status,
                             raised_guard_exact_status,
                             adaptive_default_status,
                             serialization_identical) {
  scalar_flag <- function(value, label) {
    if (!is.logical(value) || length(value) != 1L || is.na(value)) {
      stop(label, " must be one non-missing logical.", call. = FALSE)
    }
    value
  }
  provenance_verified <- scalar_flag(provenance_verified,
                                     "provenance_verified")
  serialization_identical <- scalar_flag(serialization_identical,
                                         "serialization_identical")
  statuses <- c(default_exact_status, raised_guard_exact_status,
                adaptive_default_status)
  if (anyNA(statuses) || any(!statuses %in% c("success", "error"))) {
    stop("MV5-AP method statuses must be success or error.", call. = FALSE)
  }
  safe <- provenance_verified && serialization_identical &&
    identical(default_exact_status, "success") &&
    identical(adaptive_default_status, "success")
  blocking_issues <- c(
    if (!provenance_verified) "diagram_provenance_failure",
    if (!identical(default_exact_status, "success")) {
      "default_exact_guard_rejects_realistic_diagrams"
    },
    if (!identical(adaptive_default_status, "success")) {
      "adaptive_default_not_error_controlled_on_realistic_diagrams"
    },
    if (!serialization_identical) "serialization_drift"
  )
  if (!length(blocking_issues)) blocking_issues <- "none"
  data.frame(
    contract_id = "mv05ap_realistic_compatibility_gate_v1",
    decision = if (safe) {
      "authorize_later_opt_in_workflow_integration_prefreeze"
    } else "abort_before_full_subset_and_require_numeric_engine_remediation",
    representative_subset_complete = FALSE,
    opt_in_integration_authorized = safe,
    workflow_default_change_authorized = FALSE,
    artifact_rewrite_authorized = FALSE,
    blocking_issue = if (!provenance_verified) {
      "diagram_provenance_failure"
    } else if (!identical(default_exact_status, "success")) {
      "default_exact_guard_rejects_realistic_diagrams"
    } else if (!identical(adaptive_default_status, "success")) {
      "adaptive_default_not_error_controlled_on_realistic_diagrams"
    } else if (!serialization_identical) {
      "serialization_drift"
    } else "none",
    blocking_issues = paste(blocking_issues, collapse = ";"),
    raised_guard_exact_status = raised_guard_exact_status,
    stringsAsFactors = FALSE
  )
}
