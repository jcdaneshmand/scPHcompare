# MV5-AP-R1 realistic compatibility rerun helpers.

mv05apr1_pair_plan_v1 <- function(subset) {
  required <- c(
    "diagram_id", "stratum_id", "sample_id", "selection_role",
    "h0_finite_intervals", "h1_finite_intervals", "diagram_sha256",
    "result_file_sha256", "result_file"
  )
  if (!is.data.frame(subset) || !all(required %in% names(subset)) ||
      anyNA(subset[, required]) || anyDuplicated(subset$diagram_id)) {
    stop("MV5-AP-R1 requires the complete unique frozen subset.", call. = FALSE)
  }
  groups <- split(subset, subset$stratum_id)
  if (length(groups) != 8L || any(vapply(groups, nrow, integer(1L)) != 3L)) {
    stop("MV5-AP-R1 requires exactly three diagrams in each of eight strata.",
         call. = FALSE)
  }
  rows <- lapply(groups, function(group) {
    group <- group[order(group$diagram_id, method = "radix"), , drop = FALSE]
    pairs <- utils::combn(seq_len(nrow(group)), 2L)
    do.call(rbind, lapply(seq_len(ncol(pairs)), function(index) {
      left <- group[pairs[1L, index], ]
      right <- group[pairs[2L, index], ]
      data.frame(
        stratum_id = left$stratum_id,
        first_diagram_id = left$diagram_id,
        second_diagram_id = right$diagram_id,
        first_sample_id = left$sample_id,
        second_sample_id = right$sample_id,
        pair_id = paste(left$stratum_id, left$diagram_id, right$diagram_id,
                        sep = "::"),
        stringsAsFactors = FALSE
      )
    }))
  })
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result[order(result$stratum_id, result$first_diagram_id,
               result$second_diagram_id, method = "radix"), , drop = FALSE]
}

mv05apr1_decide_v1 <- function(provenance_verified, pairs_complete,
                               certificates_pass, strict_agreement_pass,
                               repeat_pass, serialization_pass,
                               resources_pass, immutability_pass) {
  values <- c(
    provenance_verified = provenance_verified,
    pairs_complete = pairs_complete,
    certificates_pass = certificates_pass,
    strict_agreement_pass = strict_agreement_pass,
    repeat_pass = repeat_pass,
    serialization_pass = serialization_pass,
    resources_pass = resources_pass,
    immutability_pass = immutability_pass
  )
  if (!is.logical(values) || anyNA(values)) {
    stop("All MV5-AP-R1 decision inputs must be non-missing logical values.",
         call. = FALSE)
  }
  passed <- all(values)
  failures <- names(values)[!values]
  if (!length(failures)) failures <- "none"
  data.frame(
    contract_id = "mv05apr1_continuation_decision_v1",
    decision = if (passed) {
      "authorize_opt_in_workflow_integration_prefreeze_only"
    } else "stop_and_remediate_realistic_landscape_gate",
    realistic_gate_passed = passed,
    opt_in_integration_prefreeze_authorized = passed,
    workflow_integration_authorized = FALSE,
    workflow_default_change_authorized = FALSE,
    legacy_artifact_rewrite_authorized = FALSE,
    blocking_issues = paste(failures, collapse = ";"),
    next_sprint = if (passed) "MV5-AR" else "none",
    stringsAsFactors = FALSE
  )
}
