# --- Load necessary libraries ---
load_integration_libs <- function() {
  packages <- c("Seurat", "SeuratObject", "SeuratDisk", "Matrix", "parallel", "purrr", "plyr")
  loaded <- sapply(packages, requireNamespace, quietly = TRUE)
  if (!all(loaded)) {
    missing_pkgs <- packages[!loaded]
    stop(paste("Required packages for Seurat integration missing:", paste(missing_pkgs, collapse = ", ")))
  }
  library(Seurat)
  library(purrr)
  library(plyr)
}

# Define %||% operator
`%||%` <- function(a, b) if (!is.null(a)) a else b

  perform_integration <- function(
    seurat_list, # Input MUST be a list for Seurat integration
    integration_method = c("seurat", "none"), # Should always be "seurat" if called correctly now
    integration_params = list(), # For Seurat: dims, k.anchor, k.filter, k.weight, nfeatures etc.
    metadata_column = "", # Optional grouping column name
    num_cores = parallel::detectCores() %/% 2,
    min_cells_threshold = 200, # Filtering should happen before calling this
    output_path = "final_integrated_seurat.rds", # Specific name for Seurat output
# --- Seurat specific parameters (can be overridden by integration_params) ---
    integration_features_n = 3000,
    dims = 30,
    k_anchor = 5,
    k_filter = 200,
    k_weight = 100,
    min_final_anchors = 1000,
# --- Checkpointing ---
    checkpoint_dir = "./integration_checkpoints",
    anchor_batches_dir = "./anchor_batches",
    clear_existing = FALSE
) {

    # Ensure valid method (should only be seurat or none now)
    integration_method <- match.arg(integration_method)

    # Load required libraries
    load_integration_libs()

    # --- Internal Helper Functions (Keep relevant ones for Seurat) ---
    log_message <- function(msg, type = "INFO") message(sprintf("[%s] [%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%OS6"), type, msg))
    log_debug <- function(msg) log_message(msg, type = "DEBUG")
    clean_memory <- function() invisible(gc(full = TRUE))
    ensure_dir <- function(d) if (!dir.exists(d)) dir.create(d, recursive = TRUE)
    get_checkpoint_file <- function(prefix, name_id) file.path(checkpoint_dir, sprintf("%s_checkpoint_%s.h5seurat", prefix, name_id))
    save_checkpoint <- function(obj, prefix, name_id) { f <- get_checkpoint_file(prefix, name_id); log_debug(sprintf("Saving checkpoint => %s", basename(f))); tryCatch(SeuratDisk::SaveH5Seurat(obj, f, overwrite = TRUE), error = function(e) log_message(paste("Failed checkpoint save:", e$message), type = "WARN")) }
    load_checkpoint <- function(prefix, name_id) { f <- get_checkpoint_file(prefix, name_id); if (!file.exists(f)) return(NULL); log_debug(sprintf("Loading checkpoint => %s", basename(f))); tryCatch(SeuratDisk::LoadH5Seurat(f), error = function(e) { log_debug(sprintf("Checkpoint load error: %s", e$message)); NULL }) }
    get_anchor_file <- function(prefix, name_id) file.path(anchor_batches_dir, sprintf("%s_anchors_%s.rds", prefix, name_id))
    save_anchors <- function(anchors, prefix, name_id) { anchor_file <- get_anchor_file(prefix, name_id); log_debug(sprintf("Saving anchors => %s", basename(anchor_file))); tryCatch(saveRDS(anchors, anchor_file), error = function(e) log_message(paste("Failed anchor save:", e$message), type = "WARN")) }
    load_anchors <- function(prefix, name_id) { anchor_file <- get_anchor_file(prefix, name_id); if (!file.exists(anchor_file)) return(NULL); log_debug(sprintf("Loading anchors => %s", basename(anchor_file))); tryCatch(readRDS(anchor_file), error = function(e) { log_debug(sprintf("Anchor load error: %s", e$message)); NULL }) }
    # Seurat specific helpers needed:
    ensure_pca <- function(obj, dims_use = 30) {
      if (!"SCT" %in% Assays(obj)) return(obj)
      if (!("pca" %in% names(obj@reductions))) {
        log_debug(sprintf("Object %s => no PCA. Running PCA.", obj@project.name))
        var_feats <- VariableFeatures(obj)
        if (length(var_feats) == 0) { var_feats <- rownames(GetAssayData(obj, "SCT")) }
        if (length(var_feats) > 0) {
          npcs_calc <- min(dims_use, length(var_feats) - 1, ncol(obj) - 1)
          if (npcs_calc < 2) { log_message(paste("Cannot run PCA on", Project(obj)), type = "WARN"); return(obj) }
          obj <- tryCatch(RunPCA(obj, assay = "SCT", features = var_feats, npcs = npcs_calc, verbose = FALSE), error = function(e) { log_message(paste("RunPCA failed for", Project(obj), ":", e$message), type = "WARN"); return(obj) })
        } else { log_message(paste("No variable features found for", Project(obj), ", skipping PCA."), type = "WARN") }
      }
      obj
    }
    prep_sct_for_anchors <- function(obj_list, anchor_feats) { tryCatch(PrepSCTIntegration(object.list = obj_list, anchor.features = anchor_feats, assay = "SCT", verbose = TRUE), error = function(e) { log_message(paste("Error PrepSCTIntegration:", e$message), type = "ERROR"); return(obj_list) }) }
    find_anchors_with_fallback <- function(obj_list, anchor_feats, dims_use, prefix, name_id, k_anchor, k_filter, min_final_anchors) {
      # Pass params explicitly
      loadedA <- load_anchors(prefix, name_id);
      if (!is.null(loadedA)) { log_debug(sprintf("Cached anchors: %s", name_id)); return(loadedA) }
      obj_list <- prep_sct_for_anchors(obj_list, anchor_feats)
      obj_list <- lapply(obj_list, ensure_pca, dims_use = dims_use)
      anchor_feats_used <- Reduce(intersect, lapply(obj_list, function(o) VariableFeatures(o, assay = "SCT")))
      if (length(anchor_feats_used) < min_final_anchors) anchor_feats_used <- anchor_feats
      if (length(anchor_feats_used) < min_final_anchors) stop(sprintf("Need %d anchor feats, only %d remain.", min_final_anchors, length(anchor_feats_used)))

      anchors <- tryCatch({
        FindIntegrationAnchors(object.list = obj_list, normalization.method = "SCT", anchor.features = anchor_feats_used, dims = 1:dims_use, reduction = "rpca", k.anchor = k_anchor, k.filter = k_filter, verbose = TRUE)
      }, error = function(e) { log_message(sprintf("RPCA->CCA fallback: %s => %s", name_id, e$message), type = "WARN"); FindIntegrationAnchors(object.list = obj_list, normalization.method = "SCT", anchor.features = anchor_feats_used, dims = 1:dims_use, reduction = "cca", k.anchor = k_anchor, k.filter = k_filter, verbose = TRUE) })
      save_anchors(anchors, prefix, name_id);
      anchors
    }
    integrate_single_group <- function(obj_list, group_id, dims_use, integration_feats_n, k_anchor, k_filter, k_weight, min_final_anchors, global_feats = NULL, min_feats = 500) {
      cached <- load_checkpoint("group", group_id);
      if (!is.null(cached)) { log_message(sprintf("Cached group: %s", group_id)); return(cached) }
      if (length(obj_list) == 1) { log_message(sprintf("Group '%s' single object.", group_id)); return(obj_list[[1]]) }
      if (is.null(global_feats)) { global_feats <- SelectIntegrationFeatures(object.list = obj_list, nfeatures = integration_feats_n) }
      integrated_res <- tryCatch({
        anchors <- find_anchors_with_fallback(obj_list, global_feats, dims_use, prefix = "group", name_id = group_id, k_anchor, k_filter, min_final_anchors);
        if (nrow(anchors@anchors) < 10) log_message(paste("WARN: Low anchors for group", group_id), type = "WARN");
        IntegrateData(anchorset = anchors, dims = 1:dims_use, new.assay.name = "integrated", k.weight = k_weight, verbose = TRUE)
      }, error = function(e) { log_message(paste("Error integrating group", group_id, ":", e$message), type = "ERROR"); stop(e) }) # Add more sophisticated error handling if needed
      if (inherits(integrated_res, "Seurat")) { DefaultAssay(integrated_res) <- "integrated"; save_checkpoint(integrated_res, "group", group_id) };
      integrated_res
    }
    merge_pairwise_integration <- function(objA, objB, dims_use, prefix = "pairwise", name_id = "mergeX", universal_feats = NULL, integration_feats_n = 3000, k_anchor, k_filter, k_weight, min_final_anchors, min_feats = 500) {
      cached_merge <- load_checkpoint(prefix, name_id);
      if (!is.null(cached_merge)) { log_message(sprintf("Cached merge: %s", name_id)); return(cached_merge) }
      if (is.null(universal_feats)) { universal_feats <- SelectIntegrationFeatures(object.list = list(objA, objB), nfeatures = integration_feats_n) }
      merged_res <- tryCatch({
        anchors <- find_anchors_with_fallback(list(objA, objB), universal_feats, dims_use, prefix, name_id, k_anchor, k_filter, min_final_anchors);
        if (nrow(anchors@anchors) < 10) log_message(paste("WARN: Low anchors for merge", name_id), type = "WARN");
        IntegrateData(anchorset = anchors, dims = 1:dims_use, new.assay.name = "integrated", k.weight = k_weight, verbose = TRUE)
      }, error = function(e) { log_message(paste("Error merging pair", name_id, ":", e$message), type = "ERROR"); stop(e) })
      if (inherits(merged_res, "Seurat")) { DefaultAssay(merged_res) <- "integrated"; save_checkpoint(merged_res, prefix, name_id) };
      rm(objA, objB);
      gc();
      merged_res
    }
    handle_gene_sets <- function(valid_objects, nfeatures = 3000) { SelectIntegrationFeatures(object.list = valid_objects, nfeatures = nfeatures) }


    ##############################
    ## MAIN EXECUTION BLOCK
    ##############################
    final_res <- NULL # Initialize result

    if (integration_method == "none") {
      log_message("Integration method is 'none'. No integration performed.")
      return(NULL)
    }

    # --- Only Seurat Workflow Remains ---
    if (integration_method == "seurat") {
      log_message("Starting Seurat Integration Workflow...")
      if (!is.list(seurat_list) || length(seurat_list) < 2) {
        log_message("Seurat integration requires a list of at least two Seurat objects.", type = "ERROR")
        return(NULL)
      }

      # --- Use parameters passed or defaults ---
      # Use %||% safely by defining it or using rlang::`%||%`
      integration_features_n <- integration_params$integration_features_n %||% integration_features_n
      dims <- integration_params$dims %||% dims
      k_anchor <- integration_params$k_anchor %||% k_anchor
      k_filter <- integration_params$k_filter %||% k_filter
      k_weight <- integration_params$k_weight %||% k_weight
      min_final_anchors <- integration_params$min_final_anchors %||% min_final_anchors


      # --- 1) Initialize & Filter (Seurat) ---
      init_checkpoint_system() # Setup checkpoint dirs
      log_message(sprintf("Filtering & preparing %d input Seurat objects (min_cells=%d)", length(seurat_list), min_cells_threshold))
      valid_objects <- parallel::mclapply(seq_along(seurat_list), function(i) {
        obj <- seurat_list[[i]]
        if (!inherits(obj, "Seurat") || ncol(obj) < min_cells_threshold) return(NULL)
        if (!"SCT" %in% Assays(obj)) { log_message(paste("WARN: Object", Project(obj), "missing SCT assay"), type = "WARN"); return(NULL) }
        return(obj)
      }, mc.cores = num_cores)
      valid_objects <- purrr::compact(valid_objects)
      rm(seurat_list);
      clean_memory()
      if (length(valid_objects) < 2) { log_message("Seurat integration requires >= 2 valid objects.", type = "ERROR"); return(NULL) }
      log_message(sprintf("%d objects available for Seurat integration.", length(valid_objects)))


      # --- 2) Grouping Logic (Seurat) ---
      if (metadata_column == "" || !metadata_column %in% unlist(lapply(valid_objects, function(o) names(o@meta.data)))) {
        log_message(paste("No valid metadata column ('", metadata_column, "') provided/found. Integrating all as one group.", sep = ""), type = "INFO")
        group_list <- list("all_samples" = valid_objects)
      } else {
        group_ids <- sapply(valid_objects, function(o) tryCatch(as.character(unique(o@meta.data[[metadata_column]])[1]), error = function(e) NA_character_))
        if (any(is.na(group_ids))) log_message("WARN: NA group IDs found.", type = "WARN")
        group_list <- split(valid_objects, f = group_ids)
        log_message(sprintf("Generated %d groups based on '%s'.", length(group_list), metadata_column))
      }
      rm(valid_objects);
      gc()

      # --- 3) Seurat: Integrate within groups ---
      log_message("Step 1: Integrating within each group...")
      integrated_groups <- lapply(names(group_list), function(grp) {
        log_message(sprintf("--> Processing Group: %s (%d objects)", grp, length(group_list[[grp]])))
        integrate_single_group(
              obj_list = group_list[[grp]], group_id = grp,
              dims_use = dims, integration_feats_n = integration_features_n,
        # Pass Seurat specific params needed by helpers called within
              k_anchor = k_anchor, k_filter = k_filter, k_weight = k_weight, min_final_anchors = min_final_anchors
          )
      })
      names(integrated_groups) <- names(group_list)
      rm(group_list);
      gc()
      failed_groups <- sapply(integrated_groups, function(x)!inherits(x, "Seurat"))
      if (any(failed_groups)) { log_message(paste("WARN: Integration failed for groups:", paste(names(integrated_groups)[failed_groups], collapse = ", ")), type = "WARN"); integrated_groups <- integrated_groups[!failed_groups] }
      if (length(integrated_groups) == 0) { log_message("Integration failed for all groups.", type = "ERROR"); return(NULL) }

      # --- 4) Seurat: Merge across integrated groups ---
      if (length(integrated_groups) == 1) {
        log_message("Only one integrated group. Assigning as final result.")
        final_res <- integrated_groups[[1]]
      } else {
        log_message(paste("Step 2: Iteratively merging", length(integrated_groups), "integrated groups..."))
        universal_feats <- SelectIntegrationFeatures(object.list = integrated_groups, nfeatures = integration_features_n)
        final_res <- integrated_groups[[1]]
        for (i in 2:length(integrated_groups)) {
          cur_grp_id <- names(integrated_groups)[i]
          merge_id <- paste0("across_groups_iter_", i)
          log_message(sprintf("--> Merging group '%s' (iter %d)", cur_grp_id, i))
          merged_pair <- merge_pairwise_integration(
                                  final_res, integrated_groups[[i]], dims_use = dims, prefix = "pairwise",
                                  name_id = merge_id, universal_feats = universal_feats,
                                  integration_feats_n = integration_features_n,
          # Pass Seurat specific params needed by helpers
                                  k_anchor = k_anchor, k_filter = k_filter, k_weight = k_weight, min_final_anchors = min_final_anchors
                              )
          if (!inherits(merged_pair, "Seurat")) { log_message(paste("Pairwise merging failed at iter", i), type = "ERROR"); break }
          final_res <- merged_pair;
          rm(merged_pair);
          gc()
        }
        log_message("Finished merging integrated groups.")
      }

    } else {
      # integration_method is not "seurat" or "none"
      log_message(paste("Integration method '", integration_method, "' is not supported."), type = "ERROR")
      return(NULL)
    }

    # --- Save & Return ---
    if (!is.null(final_res) && inherits(final_res, "Seurat")) {
      if ("integrated" %in% Assays(final_res)) DefaultAssay(final_res) <- "integrated"
      else if (length(Assays(final_res)) > 0) DefaultAssay(final_res) <- Assays(final_res)[1]

      tryCatch({
        log_message(sprintf("Saving final integrated Seurat object => %s", output_path))
        saveRDS(final_res, output_path)
        # if(integration_method == "seurat") save_checkpoint(final_res, "final", "merged_all") # Optional
      }, error = function(e) {
        log_message(paste("Error saving final Seurat result:", e$message), type = "ERROR")
      })
    } else if (integration_method != "none") {
      log_message(paste("Seurat integration resulted in NULL or non-Seurat object."), type = "WARN")
    }

    log_message(paste("perform_integration finished for method:", integration_method))
    return(final_res)

  }
# --- End of perform_integration function ---
