# load libraries

library(Seurat)
library(SeuratObject)


load("~/Cardiac_organoids/code/entropy_functions/clean_nodatasets_060720.RData")

min_cells <- 20
celltype_col <- "celltype"


# compute entropy scores for single cell datasets:

hca.sc <- readRDS("~/Cardiac_organoids/processed_data/HCA_v1_sc_processed_filtered.rds")
sc.organoids <- readRDS("~/Cardiac_organoids/processed_data/sc_organoids.rds")
dataset_list <- c(sc.organoids, list(hca.sc = hca.sc))

dataset_names <- names(dataset_list)

entropy_df_list_sc <- lapply(
  dataset_names,
  function(dataset_name) {
    compute_entropy_for_dataset(
      seu = dataset_list[[dataset_name]],
      dataset_name = dataset_name,
      celltype_col = celltype_col,
      min_cells = min_cells
    )
  }
)

names(entropy_df_list_sc) <- dataset_names
combined_entropy_df_sc <- do.call(rbind, entropy_df_list_sc)



# compute entropy scores for single nuclei datasets:

hca.sn <- readRDS("~/Cardiac_organoids/processed_data/HCA_v1_sn_processed_filtered.rds")
sn.organoids <- readRDS("~/Cardiac_organoids/processed_data/sn_organoids.rds")
dataset_list <- c(sn.organoids, list(hca.sn = hca.sn))

dataset_names <- names(dataset_list)

entropy_df_list_sn <- lapply(
  dataset_names,
  function(dataset_name) {
    compute_entropy_for_dataset(
      seu = dataset_list[[dataset_name]],
      dataset_name = dataset_name,
      celltype_col = celltype_col,
      min_cells = min_cells
    )
  }
)

names(entropy_df_list_sn) <- dataset_names
combined_entropy_df_sn <- do.call(rbind, entropy_df_list_sn)