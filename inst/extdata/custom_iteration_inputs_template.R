# Template for supplying custom iteration inputs to `run_unified_pipeline()` and `run_postprocessing_pipeline()`.
# Replace the placeholder file paths with the Seurat object and PH list RDS files
# you want to use for each iteration. Leave any entries as NULL if you do not want
# to override that component.

custom_iteration_inputs_template <- list(
  "Seurat Integration" = list(
    seurat_object_path = "path/to/seurat_integration_seurat_object.rds",
    ph_list_path = "path/to/seurat_integration_ph_list.rds",
    bdm_matrix_path = "path/to/seurat_integration_bdm.rds",
    sdm_matrix_path = "path/to/seurat_integration_sdm.rds",
    landscape_list_path = "path/to/seurat_integration_landscape_list.rds",
    landscape_l2_distance_matrix_path = "path/to/seurat_integration_landscape_l2.rds"
  ),
  "Harmony Integration" = list(
    seurat_object_path = "path/to/harmony_integration_seurat_object.rds",
    ph_list_path = "path/to/harmony_integration_ph_list.rds",
    bdm_matrix_path = "path/to/harmony_integration_bdm.rds",
    sdm_matrix_path = "path/to/harmony_integration_sdm.rds",
    landscape_list_path = "path/to/harmony_integration_landscape_list.rds",
    landscape_l2_distance_matrix_path = "path/to/harmony_integration_landscape_l2.rds"
  ),
  Raw = list(
    seurat_object_path = "path/to/raw_seurat_object.rds",
    ph_list_path = "path/to/raw_ph_list.rds",
    bdm_matrix_path = "path/to/raw_bdm.rds",
    sdm_matrix_path = "path/to/raw_sdm.rds",
    landscape_list_path = "path/to/raw_landscape_list.rds",
    landscape_l2_distance_matrix_path = "path/to/raw_landscape_l2.rds"
  ),
  SCT_Individual = list(
    seurat_object_path = "path/to/sct_individual_seurat_object.rds",
    ph_list_path = "path/to/sct_individual_ph_list.rds",
    bdm_matrix_path = "path/to/sct_individual_bdm.rds",
    sdm_matrix_path = "path/to/sct_individual_sdm.rds",
    landscape_list_path = "path/to/sct_individual_landscape_list.rds",
    landscape_l2_distance_matrix_path = "path/to/sct_individual_landscape_l2.rds"
  ),
  SCT_Whole = list(
    seurat_object_path = "path/to/sct_whole_seurat_object.rds",
    ph_list_path = "path/to/sct_whole_ph_list.rds",
    bdm_matrix_path = "path/to/sct_whole_bdm.rds",
    sdm_matrix_path = "path/to/sct_whole_sdm.rds",
    landscape_list_path = "path/to/sct_whole_landscape_list.rds",
    landscape_l2_distance_matrix_path = "path/to/sct_whole_landscape_l2.rds"
  )
)
