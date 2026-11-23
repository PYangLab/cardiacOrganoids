# load libraries

library(Cepo)
library(SeuratObject)
library(dplyr)
library(tibble)
library(tidyr)
library(matrixStats)
library(purrr)


############ run once for sc datasets and once for sn datasets #####################

# single cell:

hca.sc <- readRDS("~/Cardiac_organoids/processed_data/HCA_v1_sc_processed_filtered.rds")
hca.sc$batch_final <- hca.sc$sangerID
hca.sc$group_label <- paste(hca.sc$batch_final, hca.sc$region, hca.sc$donor, hca.sc$gender, sep = "_") # these are the groups Cepo will run over separately
hca.sc$average_across1 <- paste(hca.sc$batch_final, hca.sc$region, hca.sc$donor, sep = "_") # averages across genders for each batch_donor_chamber
hca.sc$average_across2 <- paste(hca.sc$batch_final, hca.sc$region, sep = "_") # averages across donors for each batch_region

sc.organoids <- readRDS("~/Cardiac_organoids/processed_data/sc_organoids.rds")

seu <- hca.sc
reference_name <- "hca.sc"
organoids <- sc.organoids


# single nuclei:

hca.sn <- readRDS("~/Cardiac_organoids/processed_data/HCA_v1_sn_processed_filtered.rds")

hca.sn$batch_final <- hca.sn$sangerID
hca.sn$group_label <- paste(hca.sn$batch_final, hca.sn$region, hca.sn$donor, hca.sn$gender, sep = "_") # these are the groups Cepo will run over separately
hca.sn$average_across1 <- paste(hca.sn$batch_final, hca.sn$region, hca.sn$donor, sep = "_") # averages across genders for each batch_donor_chamber
hca.sn$average_across2 <- paste(hca.sn$batch_final, hca.sn$region, sep = "_") # averages across donors for each batch_region

sn.organoids <- readRDS("~/Cardiac_organoids/processed_data/sn_organoids.rds")

seu <- hca.sn
reference_name <- "hca.sn"
organoids <- sn.organoids

###############################################################



# subset for cell types with at least 20 cells

organoids <- lapply(organoids, function(seu_obj) {
  
  cell_counts <- table(seu_obj$celltype)
  
  keep_celltypes <- names(cell_counts)[cell_counts >= 20]
  
  seu_obj <- subset(seu_obj, subset = celltype %in% keep_celltypes)
  
  return(seu_obj)
})



# subset for common genes

# set threshold
exprsPct <- 0.05

# define function for computing "keep" genes for a Seurat object
get_keep_genes <- function(seu, celltype_col = "celltype") {
  
  cts <- unique(seu[[celltype_col]][, 1])
  cellTypes <- seu[[celltype_col]][, 1]
  
  exprsMat <- as.matrix(GetAssayData(seu, assay = "RNA", layer = "data"))
  
  meanPct.list <- lapply(cts, function(ct) {
    idx <- which(cellTypes == ct)
    (matrixStats::rowSums2(exprsMat[, idx, drop = FALSE] > 0) / length(idx)) > exprsPct
  })
  names(meanPct.list) <- cts
  
  keep <- matrixStats::rowSums2(do.call(cbind, meanPct.list)) > 0 # keep if expressed in at least one cell type
  names(keep) <- rownames(exprsMat)
  
  return(names(keep)[keep])
}

# apply function to organoid datasets
keep_list <- lapply(organoids, function(obj) {
  get_keep_genes(obj, "celltype")
})

# apply function to HCA 
keep_seu <- get_keep_genes(seu, "celltype")

# intersect common genes across organoids
common_genes_organoids <- Reduce(intersect, keep_list)

# intersect genes with HCA
common_genes <- intersect(common_genes_organoids, keep_seu)

# subset organoid and HCA genes
organoids <- lapply(organoids, function(obj) obj[common_genes, ])
seu <- seu[common_genes, ]



# run Cepo on HCA

cepo_seu <- multiCepo(
  seu = seu,
  batch_var = "group_label",
  celltype_var = "celltype",
  min_cells = 20,
  ncores = 10,
  exprs_pct = NULL
)
# remove NULLs
cepo_seu <- cepo_seu[!sapply(cepo_seu, is.null)]

averaged_by_x <- lapply(unique(seu$average_across1), function(x) {
  gather_list <- cepo_seu[grep(paste0(x, "_"), names(cepo_seu))] 
  averageCepoStats(gather_list)
})
names(averaged_by_x) <- unique(seu$average_across1) 

averaged_by_y <- lapply(unique(seu$average_across2), function(x) {
  gather_list <- averaged_by_x[grep(paste0(x, "_"), names(averaged_by_x))]
  averageCepoStats(gather_list)
})
names(averaged_by_y) <- unique(seu$average_across2)

averaged_by_z <- lapply(unique(seu$batch_final), function(x) {
  gather_list <- averaged_by_y[grep(paste0(x, "_"), names(averaged_by_y))]
  averageCepoStats(gather_list)
})
names(averaged_by_z) <- unique(seu$batch_final)

ref_cepo_average <- averageCepoStats(averaged_by_z)



# run Cepo on organoids

set.seed(2026)

organoid_cepo_bootstrap <- bootstrap_organoid_cepo(
  organoids = organoids,
  reps = 10,
  frac = 0.8,
  celltype_var = "celltype",
  min_cells = 15,
  seed = 2026
)




# correlate Cepo scores

# correlate each subsampled organoid Cepo matrix to the fixed reference Cepo matrix

bootstrap_correlations <- list()

for (organoid_name in names(organoid_cepo_bootstrap)) {
  for (rep_name in names(organoid_cepo_bootstrap[[organoid_name]])) {
    
    organoid_cepo <- organoid_cepo_bootstrap[[organoid_name]][[rep_name]]
    
    common_genes <- intersect(
      rownames(organoid_cepo),
      rownames(ref_cepo_average)
    )
    
    cor_mat <- cor(
      as.matrix(organoid_cepo[common_genes, ]),
      as.matrix(ref_cepo_average[common_genes, ]),
      method = "pearson",
      use = "pairwise.complete.obs"
    )
    
    bootstrap_correlations[[paste(organoid_name, rep_name, sep = ".")]] <-
      as.data.frame(cor_mat) %>%
      rownames_to_column("organoid_celltype") %>%
      pivot_longer(
        cols = -organoid_celltype,
        names_to = "reference_celltype",
        values_to = "cor"
      ) %>%
      mutate(
        organoid = organoid_name,
        rep = rep_name
      )
  }
}

bootstrap_correlations_df <- bind_rows(bootstrap_correlations)