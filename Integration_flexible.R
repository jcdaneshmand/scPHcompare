perform_integration <- function(
  seurat_list,
  metadata_column = "", # default is empty so the whole list is one group
  integration_features_n = 3000,
  dims = 20,
  num_cores = 16,
  min_cells_threshold = 200,
  output_path = "final_integrated_seurat.rds",
  checkpoint_dir = "./integration_checkpoints",
  anchor_batches_dir = "./anchor_batches",
  clear_existing = FALSE,
  gene_subset_method = "intersection",
  k_anchor = 5,
  k_filter = 50,
  k_weight = 50,
  min_final_anchors = 1000
) {
  ##############################
  ## 0) LOGGING & UTILS
  ##############################
  log_message <- function(msg) {
    message(sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%OS6"), msg))
  }
  log_debug <- function(msg) {
    message(sprintf("[%s] DEBUG: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%OS6"), msg))
  }
  clean_memory <- function() invisible(gc(full = TRUE))
  debug_error_handler <- function(e) {
    cat("An error occurred in a parallel worker. sys.calls():\n")
    print(sys.calls())
    stop(e)
  }

  ##############################
  ## 1) INIT CHECKPOINT SYSTEM
  ##############################
  ensure_dir <- function(d) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  }
  init_checkpoint_system <- function() {
    ensure_dir(checkpoint_dir)
    ensure_dir(anchor_batches_dir)
    if (clear_existing) {
      log_debug("Clearing existing checkpoints and caches.")
      try(unlink(file.path(checkpoint_dir, "*"), force = TRUE))
      try(unlink(file.path(anchor_batches_dir, "*"), force = TRUE))
    }
  }

  ##############################
  ## 2) CHECKPOINT & ANCHOR CACHING
  ##############################
  get_checkpoint_file <- function(prefix, name_id) {
    # Change the extension to .h5seurat for SeuratDisk compatibility.
    file.path(checkpoint_dir, sprintf("%s_checkpoint_%s.h5seurat", prefix, name_id))
  }

  save_checkpoint <- function(obj, prefix, name_id) {
    f <- get_checkpoint_file(prefix, name_id)
    log_debug(sprintf("Saving checkpoint (H5Seurat) => %s", basename(f)))
    SeuratDisk::SaveH5Seurat(obj, f, overwrite = TRUE)
  }

  load_checkpoint <- function(prefix, name_id) {
    f <- get_checkpoint_file(prefix, name_id)
    if (!file.exists(f)) {
      return(NULL)
    }
    log_debug(sprintf("Loading checkpoint (H5Seurat) => %s", basename(f)))
    tryCatch({
      SeuratDisk::LoadH5Seurat(f)
    },
      error = function(e) {
        log_debug(sprintf("Checkpoint load error => %s", e$message))
        NULL
      }
    )
  }

  get_anchor_file <- function(prefix, name_id) {
    file.path(anchor_batches_dir, sprintf("%s_anchors_%s.rds", prefix, name_id))
  }
  save_anchors <- function(anchors, prefix, name_id) {
    anchor_file <- get_anchor_file(prefix, name_id)
    log_debug(sprintf("Saving anchors => %s", basename(anchor_file)))
    saveRDS(anchors, anchor_file)
  }
  load_anchors <- function(prefix, name_id) {
    anchor_file <- get_anchor_file(prefix, name_id)
    if (!file.exists(anchor_file)) return(NULL)
    log_debug(sprintf("Loading anchors => %s", basename(anchor_file)))
    tryCatch({
      readRDS(anchor_file)
    }, error = function(e) {
      log_debug(sprintf("Error loading anchors => %s", e$message))
      NULL
    })
  }

  ##############################
  ## 3) MATRIX COMPLETION
  ##############################
  complete_matrix <- function(mat, target_genes) {
    missing <- setdiff(target_genes, rownames(mat))
    if (length(missing) > 0) {
      zeros <- Matrix::Matrix(0, nrow = length(missing), ncol = ncol(mat), sparse = TRUE)
      rownames(zeros) <- missing
      mat <- rbind(mat, zeros)
    }
    mat <- mat[target_genes,, drop = FALSE]
    mat
  }
  safe_complete_matrix <- function(mat, target_genes) {
    if (!inherits(mat, "dgCMatrix")) stop("Input matrix must be a dgCMatrix.")
    result <- complete_matrix(mat, target_genes)
    if (!identical(sort(rownames(result)), sort(target_genes))) {
      stop("Feature mismatch after matrix completion.")
    }
    result
  }

  ##############################
  ## 4) LOG PROBLEMATIC OBJECTS/GROUPS
  ##############################
  log_problematic_object <- function(obj) {
    nm <- obj@project.name
    n_cells <- ncol(obj)
    sct_data <- GetAssayData(obj, "SCT", "data")
    n_genes_sct <- nrow(sct_data)
    n_var_feats <- length(VariableFeatures(obj))
    log_message(sprintf("* Problematic Object '%s': %d cells, %d genes (SCT[data]), %d var.features",
                        nm, n_cells, n_genes_sct, n_var_feats))
    if ("pca" %in% names(obj@reductions)) {
      pca_dims <- ncol(obj@reductions$pca@cell.embeddings)
      log_message(sprintf("  PCA dimension => %d", pca_dims))
    } else {
      log_message("  No PCA reduction found.")
    }
  }
  log_problematic_group <- function(lst, group_id) {
    log_message(sprintf("Problematic group '%s' => %d objects:", group_id, length(lst)))
    for (o in lst)
      log_problematic_object(o)
  }

  ##############################
  ## 5) REMOVE EXTRANEOUS LAYERS
  ##############################
  remove_extraneous_layers <- function(sct_assay) {
    if (!requireNamespace("SeuratObject", quietly = TRUE)) {
      log_debug("SeuratObject not available; skipping layer removal.")
      return(sct_assay)
    }
    keep <- c("data", "counts", "scale.data")
    layer_names <- names(SeuratObject::Layers(sct_assay))
    for (ly in setdiff(layer_names, keep)) {
      sct_assay <- SeuratObject::RemoveLayer(sct_assay, ly)
    }
    sct_assay
  }

  ##############################
  ## 6) KITCHEN SINK CLEANUP (WITH SCT REBUILD)
  ##############################
  kitchen_sink_cleanup <- function(obj, unify_feats, dims_use = 20, verbose = FALSE) {
    # Remove old PCA and extraneous layers
    obj@reductions$pca <- NULL
    if ("SCT" %in% names(obj@assays)) {
      sct_assay <- obj[["SCT"]]
      sct_assay <- remove_extraneous_layers(sct_assay)
      obj[["SCT"]] <- sct_assay
    }
    # Restrict unify_feats to model features if available
    if ("SCT" %in% names(obj@assays) && length(obj@assays$SCT@SCTModel.list) > 0) {
      model_feats <- rownames(obj@assays$SCT@SCTModel.list[[1]]@feature.attributes)
      unify_feats <- intersect(unify_feats, model_feats)
    }
    # Rebuild the SCT assay from the raw RNA assay (mirroring main flow)
    obj <- SCTransform(
      obj,
      assay = "RNA",
      new.assay.name = "SCT",
      residual.features = unify_feats,
      verbose = verbose
    )
    if (length(obj@assays$SCT@SCTModel.list) > 1) {
      obj@assays$SCT@SCTModel.list <- tail(obj@assays$SCT@SCTModel.list, 1)
    }
    DefaultAssay(obj) <- "SCT"
    all_sct_genes <- rownames(GetAssayData(obj, "SCT", "data"))
    obj <- ScaleData(obj, assay = "SCT", features = all_sct_genes, verbose = verbose)
    if (length(VariableFeatures(obj)) == 0) {
      VariableFeatures(obj) <- all_sct_genes
    }
    obj@reductions$pca <- NULL
    obj <- RunPCA(obj, assay = "SCT", features = VariableFeatures(obj), npcs = dims_use, verbose = verbose)
    obj
  }

  ##############################
  ## 7) VERIFY FEATURES (ENSURE PCA BEFORE)
  ##############################
  ensure_pca <- function(obj, dims_use = 20) {
    if (!("pca" %in% names(obj@reductions))) {
      var_feats <- VariableFeatures(obj)
      if (length(var_feats) == 0) {
        mat <- GetAssayData(obj, "SCT", "data")
        var_feats <- rownames(mat)
        VariableFeatures(obj) <- var_feats
      }
      obj <- RunPCA(obj, assay = "SCT", features = var_feats, npcs = dims_use, verbose = FALSE)
    }
    obj
  }
  verify_features_on_subscript_error <- function(obj_list, required_feats, dims_use = 20) {
    log_message(sprintf("Verifying features => required_feats => %d", length(required_feats)))
    new_list <- lapply(obj_list, function(o) {
      o <- ensure_pca(o, dims_use = dims_use)
      mat <- GetAssayData(o, "SCT", "data")
      missing <- setdiff(required_feats, rownames(mat))
      if (length(missing) > 0) {
        new_mat <- safe_complete_matrix(mat, union(rownames(mat), required_feats))
        o <- SetAssayData(o, "SCT", "data", new.data = new_mat)
        all_genes <- rownames(new_mat)
        o <- ScaleData(o, assay = "SCT", features = all_genes, verbose = FALSE)
      }
      o
    })
    new_list
  }

  ##############################
  ## 8) PREP FOR ANCHORS
  ##############################
  prep_sct_for_anchors <- function(obj_list, anchor_feats) {
    obj_list <- PrepSCTIntegration(
      object.list = obj_list,
      anchor.features = anchor_feats,
      assay = "SCT"
    )
    obj_list
  }

  ##############################
  ## 9) FIND ANCHORS (RPCA -> CCA FALLBACK)
  ##############################
  find_anchors_with_fallback <- function(obj_list, anchor_feats, dims_use, prefix, name_id) {
    loadedA <- load_anchors(prefix, name_id)
    if (!is.null(loadedA)) {
      log_debug(sprintf("Using cached anchors => prefix=%s, name_id=%s", prefix, name_id))
      return(loadedA)
    }
    obj_list <- prep_sct_for_anchors(obj_list, anchor_feats)
    obj_list <- lapply(obj_list, ensure_pca, dims_use = dims_use)
    obj_list <- lapply(obj_list, function(o) {
      all_genes <- rownames(GetAssayData(o, "SCT", "data"))
      o <- ScaleData(o, assay = "SCT", features = all_genes, verbose = FALSE)
      o
    })
    scale_rows <- Reduce(intersect, lapply(obj_list, function(o) {
      rownames(o[["SCT"]]@scale.data)
    }))
    anchor_feats <- intersect(anchor_feats, scale_rows)
    if (length(anchor_feats) < min_final_anchors) {
      stop(sprintf("Too few anchor.features remain => %d left (need at least %d).",
                   length(anchor_feats), min_final_anchors))
    }
    anchors <- tryCatch({
      FindIntegrationAnchors(
        object.list = obj_list,
        normalization.method = "SCT",
        anchor.features = anchor_feats,
        dims = 1:dims_use,
        reduction = "rpca",
        k.anchor = k_anchor,
        k.filter = k_filter
      )
    }, error = function(e) {
      message(sprintf("RPCA error => fallback to CCA => prefix=%s, name_id=%s => %s", prefix, name_id, e$message))
      obj_list <- prep_sct_for_anchors(obj_list, anchor_feats)
      obj_list <- lapply(obj_list, ensure_pca, dims_use = dims_use)
      scale_rows <- Reduce(intersect, lapply(obj_list, function(o) {
        rownames(o[["SCT"]]@scale.data)
      }))
      anchor_feats <- intersect(anchor_feats, scale_rows)
      if (length(anchor_feats) < min_final_anchors) {
        stop(sprintf("After fallback to CCA => anchor.features => %d => below %d",
                     length(anchor_feats), min_final_anchors))
      }
      FindIntegrationAnchors(
        object.list = obj_list,
        normalization.method = "SCT",
        anchor.features = anchor_feats,
        dims = 1:dims_use,
        reduction = "cca",
        k.anchor = k_anchor,
        k.filter = k_filter
      )
    })
    save_anchors(anchors, prefix, name_id)
    anchors
  }

  ##############################
  ## 10) SINGLE GROUP INTEGRATION
  ##############################
  split_group_in_half <- function(objs) {
    n <- length(objs)
    if (n < 2) stop("Cannot split a group with <2 objects.")
    half <- as.integer(n / 2)
    list(objs[1:half], objs[(half + 1):n])
  }

  integrate_single_group <- function(obj_list, group_id, dims_use, integration_feats_n,
                                      global_feats = NULL, min_feats = 500) {
    cached <- load_checkpoint("group", group_id)
    if (!is.null(cached) && inherits(cached, "Seurat")) {
      log_message(sprintf("Loaded cached integrated single group => %s", group_id))
      return(cached)
    }
    if (length(obj_list) == 1) {
      log_message(sprintf("Group '%s' => only 1 object => returning as-is.", group_id))
      return(obj_list[[1]])
    }

    # Use Seurat's feature selection if global_feats is not provided.
    if (is.null(global_feats)) {
      global_feats <- SelectIntegrationFeatures(object.list = obj_list, nfeatures = integration_feats_n)
    }

    integrated_res <- tryCatch({
      anchors <- find_anchors_with_fallback(obj_list, global_feats, dims_use, prefix = "group", name_id = group_id)
      IntegrateData(
        anchorset = anchors,
        dims = 1:dims_use,
        new.assay.name = "integrated",
        k.weight = k_weight
      )
    }, error = function(e) {
      if (grepl("sparseMatrix", e$message, ignore.case = TRUE) ||
          grepl("p\\[length\\(p\\)\\] cannot exceed", e$message, ignore.case = TRUE)) {
        log_message(sprintf("Matrix size error in group '%s': %s", group_id, e$message))
        iteration <- 1
        repeat {
          new_nfeats <- floor(integration_feats_n * (0.95) ^ iteration)
          if (new_nfeats < min_feats) {
            stop("Unable to integrate even with the minimum feature set.")
          }
          log_message(sprintf("Retrying integration for group '%s' with %d features (iteration %d)",
                              group_id, new_nfeats, iteration))
          new_global_feats <- SelectIntegrationFeatures(object.list = obj_list, nfeatures = new_nfeats)
          result <- tryCatch({
            anchors_retry <- find_anchors_with_fallback(obj_list, new_global_feats, dims_use,
                                                        prefix = "group",
                                                        name_id = paste0(group_id, "_retry_", iteration))
            IntegrateData(
              anchorset = anchors_retry,
              dims = 1:dims_use,
              new.assay.name = "integrated",
              k.weight = k_weight
            )
          }, error = function(e2) e2)
          if (!inherits(result, "error")) {
            return(result)
          }
          iteration <- iteration + 1
        }
      } else if (grepl("subscript out of bounds", e$message, ignore.case = TRUE)) {
        log_message(sprintf("Subscript error in group '%s' => verifying features.", group_id))
        log_problematic_group(obj_list, group_id)
        verified_objs <- verify_features_on_subscript_error(obj_list, global_feats, dims_use = dims_use)
        verified_objs <- lapply(verified_objs, function(o) {
          kitchen_sink_cleanup(o, unify_feats = global_feats, dims_use = dims_use, verbose = FALSE)
        })
        verified_objs <- prep_sct_for_anchors(verified_objs, global_feats)
        anchors2 <- tryCatch({
          loaded2 <- load_anchors("group", paste0(group_id, "_retry"))
          if (!is.null(loaded2)) {
            loaded2
          } else {
            FindIntegrationAnchors(
              object.list = verified_objs,
              normalization.method = "SCT",
              anchor.features = global_feats,
              dims = 1:dims_use,
              reduction = "rpca",
              k.anchor = k_anchor,
              k.filter = k_filter
            )
          }
        }, error = function(e2) {
          message(sprintf("RPCA fallback => CCA in group '%s'_retry => %s", group_id, e2$message))
          verified_objs <- lapply(verified_objs, function(o) {
            kitchen_sink_cleanup(o, unify_feats = global_feats, dims_use = dims_use, verbose = FALSE)
          })
          verified_objs <- prep_sct_for_anchors(verified_objs, global_feats)
          FindIntegrationAnchors(
            object.list = verified_objs,
            normalization.method = "SCT",
            anchor.features = global_feats,
            dims = 1:dims_use,
            reduction = "cca",
            k.anchor = k_anchor,
            k.filter = k_filter
          )
        })
        save_anchors(anchors2, "group", paste0(group_id, "_retry"))
        integrated2 <- IntegrateData(
          anchorset = anchors2,
          dims = 1:dims_use,
          new.assay.name = "integrated",
          k.weight = k_weight
        )
        return(integrated2)
      } else {
        stop(e)
      }
    })

    if (inherits(integrated_res, "Seurat")) {
      save_checkpoint(integrated_res, "group", group_id)
    }
    integrated_res
  }

  ##############################
  ## 11) PAIRWISE MERGE
  ##############################
  merge_pairwise_integration <- function(objA, objB, dims_use, prefix = "pairwise",
                                       name_id = "mergeX", universal_feats = NULL,
                                       integration_feats_n = 6000, min_feats = 500) {
    cached_merge <- load_checkpoint(prefix, name_id)
    if (!is.null(cached_merge) && inherits(cached_merge, "Seurat")) {
      log_message(sprintf("Loaded cached pairwise merge => %s", name_id))
      return(cached_merge)
    }

    # If universal_feats is not provided, compute it using SelectIntegrationFeatures.
    if (is.null(universal_feats)) {
      universal_feats <- SelectIntegrationFeatures(object.list = list(objA, objB),
                                                 nfeatures = integration_feats_n)
    }

    # Helper function to try the integration with a given feature set.
    try_integration <- function(nfeats, feats) {
      # Uses variables dims_use, prefix, name_id, and k_weight from the outer scope.
      anchors <- find_anchors_with_fallback(list(objA, objB), feats, dims_use, prefix, name_id)
      IntegrateData(
      anchorset = anchors,
      dims = 1:dims_use,
      new.assay.name = "merged_pairwise",
      k.weight = k_weight
    )
    }

    merged <- tryCatch({
      try_integration(integration_feats_n, universal_feats)
    }, error = function(e) {
      if (grepl("sparseMatrix", e$message, ignore.case = TRUE) ||
        grepl("p\\[length\\(p\\)\\] cannot exceed", e$message, ignore.case = TRUE)) {
        log_message(sprintf("Matrix size error in merge: %s", e$message))
        # Iteratively reduce features by 10% until the merge succeeds or we hit the minimum threshold.
        new_nfeats <- integration_feats_n
        repeat {
          new_nfeats <- floor(new_nfeats * 0.9)
          if (new_nfeats < min_feats) {
            stop("Unable to merge even with the minimum feature set.")
          }
          log_message(sprintf("Retrying merge with %d features", new_nfeats))
          new_universal_feats <- SelectIntegrationFeatures(object.list = list(objA, objB),
                                                         nfeatures = new_nfeats)
          result <- tryCatch({
            try_integration(new_nfeats, new_universal_feats)
          }, error = function(e2) e2)
          if (!inherits(result, "error")) {
            return(result)
          }
        }
      } else if (grepl("subscript out of bounds", e$message, ignore.case = TRUE)) {
        log_message(sprintf("Subscript error in merge => verifying features => prefix=%s, name_id=%s", prefix, name_id))
        log_problematic_group(list(objA, objB), name_id)
        featsA <- rownames(GetAssayData(objA, "SCT", "data"))
        featsB <- rownames(GetAssayData(objB, "SCT", "data"))
        if (is.null(universal_feats)) {
          unify_feats <- union(featsA, featsB)
        } else {
          unify_feats <- universal_feats
        }
        newA <- verify_features_on_subscript_error(list(objA), unify_feats, dims_use = dims_use)[[1]]
        newB <- verify_features_on_subscript_error(list(objB), unify_feats, dims_use = dims_use)[[1]]
        newA <- kitchen_sink_cleanup(newA, unify_feats = unify_feats, dims_use = dims_use, verbose = FALSE)
        newB <- kitchen_sink_cleanup(newB, unify_feats = unify_feats, dims_use = dims_use, verbose = FALSE)
        newAB <- prep_sct_for_anchors(list(newA, newB), unify_feats)
        anchors2 <- tryCatch({
          loaded2 <- load_anchors(prefix, paste0(name_id, "_retry"))
          if (!is.null(loaded2)) {
            loaded2
          } else {
            FindIntegrationAnchors(
            object.list = newAB,
            normalization.method = "SCT",
            anchor.features = unify_feats,
            dims = 1:dims_use,
            reduction = "rpca",
            k.anchor = k_anchor,
            k.filter = k_filter
          )
          }
        }, error = function(e2) {
          message(sprintf("RPCA fallback => CCA in merge '%s'_retry => %s", name_id, e2$message))
          newAB <- lapply(newAB, function(o) {
            kitchen_sink_cleanup(o, unify_feats = unify_feats, dims_use = dims_use, verbose = FALSE)
          })
          newAB <- prep_sct_for_anchors(newAB, unify_feats)
          FindIntegrationAnchors(
          object.list = newAB,
          normalization.method = "SCT",
          anchor.features = unify_feats,
          dims = 1:dims_use,
          reduction = "cca",
          k.anchor = k_anchor,
          k.filter = k_filter
        )
        })
        save_anchors(anchors2, prefix, paste0(name_id, "_retry"))
        integrated2 <- IntegrateData(
        anchorset = anchors2,
        dims = 1:dims_use,
        new.assay.name = "merged_pairwise",
        k.weight = k_weight
      )
        return(integrated2)
      } else {
        stop(e)
      }
    })

    if (inherits(merged, "Seurat")) {
      save_checkpoint(merged, prefix, name_id)
    }
    rm(objA, objB);
    gc()
    merged
  }

  ##############################
  ## 12) HANDLE GENE SETS
  ##############################
  handle_gene_sets <- function(valid_objects, nfeatures = 6000) {
    # Simply select the top nfeatures across the objects.
    SelectIntegrationFeatures(object.list = valid_objects, nfeatures = nfeatures)
  }

  ##############################
  ## 13) MAIN EXECUTION BLOCK
  ##############################
  if (file.exists(output_path)) {
    log_message(sprintf("Loading existing final result => %s", output_path))
    return(readRDS(output_path))
  }
  init_checkpoint_system()
  log_message(sprintf("Filtering & converting %d Seurat objects, min_cells_threshold=%d",
                      length(seurat_list), min_cells_threshold))
  valid_objects <- parallel::mclapply(seq_along(seurat_list), function(i) {
    tryCatch({
      obj <- seurat_list[[i]]
      if (!inherits(obj, "Seurat")) return(NULL)
      if (ncol(obj) < min_cells_threshold) {
        log_debug(sprintf("Skipping object %d => fewer than %d cells", i, min_cells_threshold))
        return(NULL)
      }
      new_names <- paste0("Sample", i, "_", colnames(obj))
      obj <- RenameCells(obj, new.names = new_names)
      obj
    }, error = debug_error_handler)
  }, mc.cores = num_cores, mc.preschedule = FALSE)
  valid_objects <- purrr::compact(valid_objects)
  rm(seurat_list)
  clean_memory()
  if (length(valid_objects) == 0) stop("No valid Seurat objects remain after filtering.")

  log_message("Ensuring each valid object has PCA if not already.")
  valid_objects <- parallel::mclapply(valid_objects, function(obj) {
    tryCatch({
      dat <- GetAssayData(obj, "SCT", "data")
      if (nrow(dat) == 0) stop(sprintf("Object %s => empty SCT data layer", obj@project.name))
      if (length(VariableFeatures(obj)) == 0) {
        log_debug(sprintf("Object %s => no var.features; using all genes.", obj@project.name))
        VariableFeatures(obj) <- rownames(dat)
      }
      if (!("pca" %in% names(obj@reductions))) {
        log_debug(sprintf("Object %s => no PCA. Running PCA.", obj@project.name))
        var_feats <- if (length(VariableFeatures(obj)) == 0) rownames(dat) else VariableFeatures(obj)
        obj <- RunPCA(obj, assay = "SCT", features = var_feats, verbose = TRUE)
      }
      obj
    }, error = debug_error_handler)
  }, mc.cores = num_cores)
  clean_memory()

  log_message("Handling gene sets (intersection or union).")
  target_genes <- handle_gene_sets(valid_objects)

  log_message("Rebuilding SCT assay for each object to unify gene sets.")
  valid_objects <- parallel::mclapply(valid_objects, function(obj) {
    tryCatch({
      if (!("SCT" %in% names(obj@assays))) {
        log_debug(sprintf("Object %s => no SCT assay => skipping rebuild", obj@project.name))
        return(obj)
      }
      old_data <- GetAssayData(obj, "SCT", "data")
      old_counts <- tryCatch({ GetAssayData(obj, "SCT", "counts") }, error = function(e) NULL)
      if (is.null(old_counts) || nrow(old_counts) == 0) {
        old_counts <- expm1(old_data)
        log_debug(sprintf("Object %s => 'counts' missing => using expm1(data).", obj@project.name))
      }
      if (gene_subset_method == "intersection") {
        common_genes <- intersect(rownames(old_data), target_genes)
        old_data <- old_data[common_genes,, drop = FALSE]
        old_counts <- old_counts[common_genes,, drop = FALSE]
      } else {
        completed_data <- safe_complete_matrix(old_data, target_genes)
        old_data <- completed_data
        if (!all(rownames(old_data) %in% rownames(old_counts))) {
          old_counts <- safe_complete_matrix(old_counts, rownames(old_data))
        } else {
          old_counts <- old_counts[rownames(old_data),, drop = FALSE]
        }
      }
      common_genes <- rownames(old_data)
      meta_feats <- obj@assays$SCT@meta.features
      if (!is.null(meta_feats)) {
        missing_in_meta <- setdiff(common_genes, rownames(meta_feats))
        if (length(missing_in_meta) > 0) {
          extra_meta <- data.frame(matrix(NA, nrow = length(missing_in_meta), ncol = ncol(meta_feats)))
          rownames(extra_meta) <- missing_in_meta
          colnames(extra_meta) <- colnames(meta_feats)
          meta_feats <- rbind(meta_feats, extra_meta)
        }
        meta_feats <- meta_feats[common_genes,, drop = FALSE]
      }
      new_var_features <- intersect(VariableFeatures(obj), common_genes)
      if (length(new_var_features) < 1000) {
        new_var_features <- union(new_var_features, common_genes)
      }
      new_sct_assay <- new("SCTAssay",
                           data = old_data,
                           counts = old_counts,
                           scale.data = matrix(nrow = 0, ncol = ncol(old_data)),
                           meta.features = meta_feats,
                           var.features = new_var_features,
                           SCTModel.list = obj[["SCT"]]@SCTModel.list,
                           misc = obj[["SCT"]]@misc,
                           key = obj[["SCT"]]@key)
      obj[["SCT"]] <- new_sct_assay
      if (!is.null(obj@assays$SCT@SCTModel.list)) {
        model_feats <- rownames(obj@assays$SCT@SCTModel.list[[1]]@feature.attributes)
        data_feats <- rownames(GetAssayData(obj, "SCT", "data"))
        if (!all(model_feats %in% data_feats)) {
          obj <- SCTransform(obj, assay = "RNA", new.assay.name = "SCT", residual.features = common_genes)
          if (length(obj@assays$SCT@SCTModel.list) > 1) {
            obj@assays$SCT@SCTModel.list <- tail(obj@assays$SCT@SCTModel.list, 1)
          }
          DefaultAssay(obj) <- "SCT"
        }
      }
      all_sct_genes <- rownames(GetAssayData(obj, "SCT", "data"))
      obj <- ScaleData(obj, assay = "SCT", features = all_sct_genes, verbose = FALSE)
      obj
    }, error = debug_error_handler)
  }, mc.cores = num_cores)
  clean_memory()

  ##############################
  ## 14) GROUPING OF OBJECTS FOR INTEGRATION
  ##############################
  # New option: if no metadata column is provided (empty string), treat the whole list as one group.
  if (metadata_column == "") {
    log_message("No metadata column provided. Integrating the entire Seurat list as one group.")
    group_list <- list("all" = valid_objects)
  } else {
    group_ids <- sapply(valid_objects, function(o) unique(o@meta.data[[metadata_column]]))
    group_list <- split(valid_objects, f = group_ids)
    log_message(sprintf("Generated %d groups based on '%s'.", length(group_list), metadata_column))
  }

  ##############################
  ## 15) INTEGRATION OF GROUPS
  ##############################
  # Step 1: Integrate each group separately.
  integrated_groups <- lapply(names(group_list), function(grp) {
    log_message(sprintf("Integrating group => %s", grp))
    integrate_single_group(
      obj_list = group_list[[grp]],
      group_id = grp,
      dims_use = dims,
      integration_feats_n = integration_features_n
    )
  })
  names(integrated_groups) <- names(group_list)

  # If only one group exists, no further merging is needed.
  if (length(integrated_groups) == 1) {
    log_message("Only one integrated group present. Skipping merge steps.")
    final_res <- integrated_groups[[1]]
  } else {
    # Step 2: Select universal features using Seurat's integration feature selection.
    universal_feats <- SelectIntegrationFeatures(object.list = integrated_groups, nfeatures = integration_features_n)

    # Step 3: Iteratively merge integrated groups one by one.
    final_res <- integrated_groups[[1]]
    for (i in 2:length(integrated_groups)) {
      cur_grp_id <- names(integrated_groups)[i]
      merge_id <- paste0("merge_iter_", i)
      log_message(sprintf("Merging integrated group '%s' into current merged object.", cur_grp_id))
      merge_attempt <- tryCatch({
        merge_pairwise_integration(
          final_res,
          integrated_groups[[i]],
          dims_use = dims,
          prefix = "pairwise",
          name_id = merge_id,
          universal_feats = universal_feats
        )
      }, error = function(e) {
        log_message(sprintf("Pairwise merge error for merge_id '%s': %s", merge_id, e$message))
        # Attempt one more round of verification & cleanup before giving up on merge.
        verified_objs <- verify_features_on_subscript_error(
          list(final_res, integrated_groups[[i]]),
          universal_feats,
          dims_use = dims
        )
        verified_objs <- lapply(verified_objs, function(o) {
          kitchen_sink_cleanup(o, unify_feats = universal_feats, dims_use = dims, verbose = FALSE)
        })
        merged_retry <- tryCatch({
          merge_pairwise_integration(
            verified_objs[[1]],
            verified_objs[[2]],
            dims_use = dims,
            prefix = "pairwise",
            name_id = paste0(merge_id, "_retry"),
            universal_feats = universal_feats
          )
        }, error = function(e2) {
          log_message(sprintf("Retry merge failed for merge_id '%s': %s", merge_id, e2$message))
          NULL
        })
        if (is.null(merged_retry)) {
          log_message("Merging halted for this pair; returning detailed list of integrated groups.")
          return(list(
            error_message = e$message,
            growing_object = final_res,
            problematic_object = integrated_groups[[i]],
            integrated_groups = integrated_groups
          ))
        } else {
          merged_retry
        }
      })
      if (is.list(merge_attempt)) {
        # merge failed; return list with details
        final_res <- merge_attempt
        break
      } else {
        final_res <- merge_attempt
        integrated_groups[[i]] <- merge_attempt
      }
      clean_memory()
    }
  }

  ##############################
  ## 16) SAVE & RETURN RESULT
  ##############################
  saveRDS(final_res, output_path)
  log_message(sprintf("Integration complete. Saved => %s", output_path))
  final_res}
