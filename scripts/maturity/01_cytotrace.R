# load libraries

library(CytoTRACE2)


# run CytoTRACE2 on single cell datasets

hca.sc <- readRDS("~/Cardiac_organoids/processed_data/HCA_v1_sc_processed_filtered.rds")
sc.organoids <- readRDS("~/Cardiac_organoids/processed_data/sc_organoids.rds")
dataset_list <- c(sc.organoids, list(hca.sc = hca.sc))


# find gene intersection
gene_lists <- lapply(dataset_list, function(seu) {
  rownames(seu@assays$RNA$counts)
})
common_genes <- Reduce(intersect, gene_lists)
length(common_genes)

count_matrices <- list()
metadata_list <- list()

for (i in seq_along(dataset_list)) {
  
  seu <- dataset_list[[i]]
  dataset <- names(dataset_list)[i]
  
  new_cellnames <- paste(dataset, colnames(seu), sep = "_") # to avoid duplicate cell names
  
  colnames(seu) <- new_cellnames
  
  counts <- seu@assays$RNA$counts[common_genes, ]
  
  count_matrices[[i]] <- counts
  
  meta <- data.frame(
    cell = new_cellnames,
    celltype = seu$celltype,
    dataset = dataset
  )
  
  metadata_list[[i]] <- meta
}

# merge counts and metadata
combined_counts_sc <- do.call(cbind, count_matrices)
combined_metadata_sc <- do.call(rbind, metadata_list)

# ensure ordering matches
combined_metadata_sc <- combined_metadata_sc[match(colnames(combined_counts_sc), combined_metadata_sc$cell), ]
identical(colnames(combined_counts_sc), combined_metadata_sc$cell)

cytotrace2_result_sc <- cytotrace2(combined_counts_sc, species = "human")

idx <- match(rownames(cytotrace2_result_sc), combined_metadata_sc$cell)
stopifnot(!anyDuplicated(combined_metadata_sc$cell), !anyNA(idx))

cytotrace2_all_sc <- data.frame(
  combined_metadata_sc[idx, c("cell", "dataset", "celltype"), drop = FALSE],
  cytotrace2_result_sc,
  row.names = NULL
)

cell_counts <- as.data.frame(with(cytotrace2_all_sc, table(dataset, celltype)))
valid_groups <- subset(cell_counts, Freq >= 20, select = c(dataset, celltype))

cytotrace2_filtered_sc <- merge(
  cytotrace2_all_sc,
  valid_groups,
  by = c("dataset", "celltype"),
  sort = FALSE
)




# run CytoTRACE2 on single nuclei datasets

hca.sn <- readRDS("~/Cardiac_organoids/processed_data/HCA_v1_sn_processed_filtered.rds")
sn.organoids <- readRDS("~/Cardiac_organoids/processed_data/sn_organoids.rds")
dataset_list <- c(sn.organoids, list(hca.sn = hca.sn))

# find gene intersection
gene_lists <- lapply(dataset_list, function(seu) {
  rownames(seu@assays$RNA$counts)
})
common_genes <- Reduce(intersect, gene_lists)
length(common_genes)

count_matrices <- list()
metadata_list <- list()

for (i in seq_along(dataset_list)) {
  
  seu <- dataset_list[[i]]
  dataset <- names(dataset_list)[i]
  
  new_cellnames <- paste(dataset, colnames(seu), sep = "_") # to avoid duplicate cell names
  
  colnames(seu) <- new_cellnames
  
  counts <- seu@assays$RNA$counts[common_genes, ]
  
  count_matrices[[i]] <- counts
  
  meta <- data.frame(
    cell = new_cellnames,
    celltype = seu$celltype,
    dataset = dataset
  )
  
  metadata_list[[i]] <- meta
}

# merge counts and metadata
combined_counts_sn <- do.call(cbind, count_matrices)
combined_metadata_sn <- do.call(rbind, metadata_list)

# ensure ordering matches
combined_metadata_sn <- combined_metadata_sn[match(colnames(combined_counts_sn), combined_metadata_sn$cell), ]
identical(colnames(combined_counts_sn), combined_metadata_sn$cell)

cytotrace2_result_sn <- cytotrace2(combined_counts_sn, species = "human")

idx <- match(rownames(cytotrace2_result_sn), combined_metadata_sn$cell)
stopifnot(!anyDuplicated(combined_metadata_sn$cell), !anyNA(idx))

cytotrace2_all_sn <- data.frame(
  combined_metadata_sn[idx, c("cell", "dataset", "celltype"), drop = FALSE],
  cytotrace2_result_sn,
  row.names = NULL
)

cell_counts <- as.data.frame(with(cytotrace2_all_sn, table(dataset, celltype)))
valid_groups <- subset(cell_counts, Freq >= 20, select = c(dataset, celltype))

cytotrace2_filtered_sn <- merge(
  cytotrace2_all_sn,
  valid_groups,
  by = c("dataset", "celltype"),
  sort = FALSE
)