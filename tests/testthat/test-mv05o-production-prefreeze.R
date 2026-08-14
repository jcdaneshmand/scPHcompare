test_that("MV5-O source freeze rejects prohibited paths", {
  expect_error(mv05o_source_freeze_v1(
    "docs/private/reviewer.pdf", "bad", strrep("a", 40L)),
    "incomplete|prohibited"
  )
})

test_that("MV5-O inventory and queue contracts are exact", {
  make_groups <- function() {
    one <- expand.grid(fold_index = 1:15, seed = 20260805:20260809,
                       stringsAsFactors = FALSE)
    rows <- do.call(rbind, lapply(
      c("sct_whole", "inductive_integrated"), function(representation) {
        transform(one, representation = representation)
      }
    ))
    pair_counts <- rep(
      c(2080L, 3160L, 3240L, 3321L, 3486L, 3570L, 3655L, 3828L, 3916L),
      c(5L, 5L, 5L, 5L, 5L, 10L, 20L, 10L, 10L)
    )
    rows$source_group_id <- paste(rows$representation, rows$fold_index,
                                  rows$seed, sep = ":")
    rows$group_order <- ave(seq_len(nrow(rows)), rows$representation, FUN = seq_along)
    rows$fold_id <- paste0("fold:", rows$fold_index)
    rows$unordered_training_pairs <- rep(pair_counts, 2L)
    rows$training_samples <- as.integer(
      (1 + sqrt(1 + 8 * rows$unordered_training_pairs)) / 2
    )
    rows$h0_h1_request_rows <- 2L * rows$unordered_training_pairs
    rows$chunk_count <- rows$h0_h1_request_rows
    rows$request_identity_set_sha256 <- vapply(rows$source_group_id,
                                               .mv05o_digest, character(1L))
    rows$pair_scope <- "training_training_unordered"
    rows$landscape_definition_id <- "all_active_exact_critical_pairs_v1"
    rows$outcome_label_state <- "closed"
    rows$biological_outcomes_computed <- FALSE
    rows
  }
  groups <- make_groups()
  chunks <- do.call(rbind, lapply(seq_len(nrow(groups)), function(i) {
    do.call(rbind, lapply(c("H0", "H1"), function(dimension) {
      count <- groups$unordered_training_pairs[[i]]
      data.frame(
        chunk_id = paste0("chunk:", i, ":", dimension),
        source_group_id = groups$source_group_id[[i]],
        group_order = groups$group_order[[i]], fold_id = groups$fold_id[[i]],
        seed = groups$seed[[i]], representation = groups$representation[[i]],
        homology_dimension = dimension, request_rows = count,
        request_identity_set_sha256 = .mv05o_digest(paste(i, dimension)),
        first_pair_request_id = paste0("pair:", i, ":", dimension, ":1"),
        last_pair_request_id = paste0("pair:", i, ":", dimension, ":", count),
        outcome_label_state = "closed", biological_outcomes_computed = FALSE,
        stringsAsFactors = FALSE
      )
    }))
  }))
  expanded <- lapply(seq_len(nrow(chunks)), function(i) {
    n <- chunks$request_rows[[i]]
    sizes <- c(rep(250L, n %/% 250L), if (n %% 250L) n %% 250L)
    sizes <- sizes[sizes > 0L]
    rows <- chunks[rep(i, length(sizes)), ]
    rows$request_rows <- sizes
    rows$chunk_id <- paste0(rows$chunk_id, ":", seq_along(sizes))
    rows$request_identity_set_sha256 <- vapply(rows$chunk_id, .mv05o_digest,
                                               character(1L))
    rows
  })
  chunks <- do.call(rbind, expanded)
  groups$chunk_count <- as.integer(table(paste(chunks$source_group_id,
                                                chunks$representation, sep = "\r"))[
    paste(groups$source_group_id, groups$representation, sep = "\r")])
  expect_invisible(mv05o_validate_mv05n_inventories_v1(groups, chunks))
  queues <- mv05o_build_queues_v1(
    groups, chunks, strrep("a", 64L),
    c(prefreeze = strrep("b", 64L), stager = strrep("c", 64L),
      landscape = strrep("d", 64L), baseline = strrep("e", 64L))
  )
  expect_equal(nrow(queues$groups), 150L)
  expect_equal(nrow(queues$landscape), 4340L)
  expect_equal(nrow(queues$baseline), 225L)
  expect_equal(sum(queues$landscape$request_rows), 1050700L)
  expect_equal(sum(queues$baseline$pair_rows), 788025L)
  expect_false(any(queues$groups$production_executed))
  plan <- mv05o_build_validation_plan_v1(
    queues$groups, queues$landscape, queues$baseline
  )
  expect_equal(nrow(plan), 15L)
  expect_equal(sum(grepl("exact_r_oracle", plan$validation_id)), 12L)
  expect_error(.mv05o_assert_label_closed(data.frame(tissue = "x")), "prohibited")
})

test_that("MV5-O abort registry freezes ten non-retrying guards", {
  rules <- mv05o_abort_rules_v1()
  expect_equal(nrow(rules), 10L)
  expect_false(any(rules$automatic_retry))
  expect_true(all(rules$outcome_label_state == "closed"))
  expect_true(any(grepl("21.6", rules$trigger, fixed = TRUE)))
  expect_true(any(grepl("4294967296", rules$trigger, fixed = TRUE)))
  expect_true(any(grepl("10737418240", rules$trigger, fixed = TRUE)))
})
