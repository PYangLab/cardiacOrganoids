# load libraries

library(SingleCellExperiment)
library(Seurat)
library(zellkonverter)
library(Matrix)
library(scReClassify)


# load the log-normalised counts

h5 <- "~/Cardiac_organoids/raw_data/HCA/Global_lognormalised.hdf"

# Read sparse X matrix
x_data <- as.numeric(rhdf5::h5read(h5, "/X/data"))
x_indices <- as.integer(rhdf5::h5read(h5, "/X/indices"))
x_indptr <- as.integer(rhdf5::h5read(h5, "/X/indptr"))

lognorm_mat <- sparseMatrix(
  i = x_indices + 1, # AnnData stores zero-based indices, so add 1 for R
  p = x_indptr,
  x = x_data,
  dims = c(32732, 704296)
)


# add gene IDs / names

gene_ids <- h5read(h5, "/var/gene_id")
gene_names <- h5read(h5, "/var/gene_name-new")
rownames(lognorm_mat) <- gene_ids


# load metadata

read_cat <- function(h5, path) {
  cats <- as.vector(h5read(h5, paste0(path, "/categories")))
  codes <- as.integer(h5read(h5, paste0(path, "/codes")))
  
  out <- rep(NA_character_, length(codes))
  keep <- codes >= 0
  out[keep] <- cats[codes[keep] + 1]
  out
}

read_obs_field <- function(h5, name, otype) {
  path <- paste0("/obs/", name)
  
  if (otype == "H5I_GROUP") {
    read_cat(h5, path)
  } else {
    as.vector(h5read(h5, path))
  }
}

obs_info <- h5ls(h5)
obs_fields <- obs_info[obs_info$group == "/obs", c("name", "otype")]

obs_list <- setNames(
  lapply(seq_len(nrow(obs_fields)), function(i) {
    read_obs_field(h5, obs_fields$name[i], obs_fields$otype[i])
  }),
  obs_fields$name
)

obs <- as.data.frame(obs_list, stringsAsFactors = FALSE, check.names = FALSE)

rownames(obs) <- obs$barcode
colnames(lognorm_mat) <- obs$barcode


# load the raw counts

sce <- zellkonverter::readH5AD("~/Cardiac_organoids/raw_data/HCA/HCA_v2_global_raw.h5ad", reader = "R")

raw_counts <- assay(sce, "X")
raw_counts <- as(raw_counts, "dgCMatrix")


# set the same gene order

lognorm_mat <- lognorm_mat[rownames(raw_counts), ]

print(identical(rownames(lognorm_mat), rownames(raw_counts)))
identical(
  unname(as.character(rownames(lognorm_mat))),
  unname(as.character(rownames(raw_counts)))
)
stopifnot(all(rownames(lognorm_mat) == rownames(raw_counts)))


# set the same cell order

lognorm_mat <- lognorm_mat[, colnames(raw_counts)]
obs <- obs[colnames(raw_counts), , drop = FALSE]

stopifnot(all(colnames(lognorm_mat) == colnames(raw_counts)))
stopifnot(all(rownames(obs) == colnames(raw_counts)))


# swap Ensembl IDs for gene symbols

gene_lookup <- data.frame(
  gene_id = as.character(gene_ids),
  gene_symbol = as.character(gene_names),
  stringsAsFactors = FALSE
)

new_symbols <- gene_lookup$gene_symbol[match(rownames(raw_counts), gene_lookup$gene_id)]

keep <- !is.na(new_symbols) & new_symbols != "" # TRUE for all
raw_counts <- raw_counts[keep, ]
lognorm_mat <- lognorm_mat[keep, ]
new_symbols <- new_symbols[keep]

# conservative duplicate handling
keep_unique <- !duplicated(new_symbols) # TRUE for all
raw_counts <- raw_counts[keep_unique, ]
lognorm_mat <- lognorm_mat[keep_unique, ]
new_symbols <- new_symbols[keep_unique]

rownames(raw_counts) <- new_symbols
rownames(lognorm_mat) <- new_symbols

stopifnot(all(rownames(raw_counts) == rownames(lognorm_mat)))


# create Seurat object

seu <- CreateSeuratObject(
  counts = raw_counts,
  meta.data = obs,
  assay = "RNA"
)

seu <- SetAssayData(
  seu,
  assay = "RNA",
  layer = "data",
  new.data = lognorm_mat
)


# subset for HCA v1 (sn and sc; multiome excluded)

hca.sn <- seu[,seu$modality == "snRNA"]
hca.sc <- seu[,seu$modality == "scRNA"]


# filtering according to Litvinukova et al. 

filter_hca_object <- function(seu, modality = c("sc", "sn")) {
  modality <- match.arg(modality)
  
  meta <- seu@meta.data
  
  required_cols <- c(
    "nCount_RNA",
    "nFeature_RNA",
    "pct_counts_mt",
    "pct_counts_ribo",
    "scrublet_score"
  )
  
  missing_cols <- setdiff(required_cols, colnames(meta))
  if (length(missing_cols) > 0) {
    stop("Missing metadata columns: ", paste(missing_cols, collapse = ", "))
  }
  
  if (modality == "sc") {
    keep <- with(meta,
                 nCount_RNA > 500 &
                   nCount_RNA < 15000 &
                   nFeature_RNA > 200 &
                   pct_counts_mt < 20 &
                   pct_counts_ribo < 20 &
                   scrublet_score < 0.3
    )
  }
  
  if (modality == "sn") {
    keep <- with(meta,
                 nCount_RNA > 500 &
                   nCount_RNA < 15000 &
                   nFeature_RNA > 300 &
                   nFeature_RNA < 6000 &
                   pct_counts_mt < 5 &
                   pct_counts_ribo < 5 &
                   scrublet_score < 0.3
    )
  }
  
  keep[is.na(keep)] <- FALSE
  
  message(modality, ": keeping ", sum(keep), " / ", length(keep), " cells")
  
  subset(seu, cells = rownames(meta)[keep])
}

hca.sc <- filter_hca_object(hca.sc, modality = "sc")
hca.sn <- filter_hca_object(hca.sn, modality = "sn")

# drop cell types with only 1 cell
drop_types <- c("Adipocyte", "Mast cell")
keep <- !(hca.sc$cell_type %in% drop_types)
hca.sc <- subset(hca.sc, cells = colnames(hca.sc)[keep])
hca.sc$cell_type <- droplevels(factor(hca.sc$cell_type))


# data processing

hca.sn <- FindVariableFeatures(hca.sn)
hca.sn <- ScaleData(hca.sn)
hca.sn <- RunPCA(hca.sn)
hca.sn <- RunUMAP(hca.sn, reduction = "pca", dims = 1:30)

hca.sc <- FindVariableFeatures(hca.sc)
hca.sc <- ScaleData(hca.sc)
hca.sc <- RunPCA(hca.sc)
hca.sc <- RunUMAP(hca.sc, reduction = "pca", dims = 1:30)


# run scReClassify on single nuclei:

# determine number of PCs to cover 70% variance
pcs <- Embeddings(hca.sn, reduction = "pca") 
stdev <- hca.sn@reductions$pca@stdev
var_per_pc <- stdev^2
cum_var <- cumsum(var_per_pc) / sum(var_per_pc)
d <- which(cum_var >= 0.7)[1]

# clamp to 10–20 PCs
M_prime <- max(10, min(d, 20))
matPCs <- pcs[, 1:M_prime]
dat <- as.SingleCellExperiment(hca.sn)
reducedDim(dat, "matPCs") <- matPCs

# run scReClassify
set.seed(1)
cellTypes.reclassify <- multiAdaSampling(
  dat,
  hca.sn$cell_type,
  reducedDimName = "matPCs",
  classifier = "svm",
  percent = 1,
  L = 10
)

# percentage of cells that change:
paste(mean(cellTypes.reclassify$final != hca.sn$cell_type) * 100, "%") # 0.316740544488035 %

hca.sn$scReClassify_label <- cellTypes.reclassify$final


# run scReClassify on single cell:

# remove cell types with only one cell
drop_types <- c("Adipocyte", "Mast cell")
keep <- !(hca.sc$cell_type %in% drop_types)
hca.sc <- subset(hca.sc, cells = colnames(hca.sc)[keep])
hca.sc$cell_type <- droplevels(factor(hca.sc$cell_type))

# determine number of PCs to cover 70% variance
pcs <- Embeddings(hca.sc, reduction = "pca") 
stdev <- hca.sc@reductions$pca@stdev
var_per_pc <- stdev^2
cum_var <- cumsum(var_per_pc) / sum(var_per_pc)
d <- which(cum_var >= 0.7)[1]

# clamp to 10–20 PCs
M_prime <- max(10, min(d, 20))
matPCs <- pcs[, 1:M_prime]
dat <- as.SingleCellExperiment(hca.sc)
reducedDim(dat, "matPCs") <- matPCs

# run scReClassify
set.seed(1)
cellTypes.reclassify <- multiAdaSampling(
  dat,
  hca.sc$cell_type,
  reducedDimName = "matPCs",
  classifier = "svm",
  percent = 1,
  L = 10
)

# percentage of cells that change:
paste(mean(cellTypes.reclassify$final != hca.sc$cell_type) * 100, "%") # 0.748459383753501 %

hca.sc$scReClassify_label <- cellTypes.reclassify$final


# harmonise cell type labels

hca.sc$celltype <- as.character(hca.sc$scReClassify_label)
hca.sc$celltype[grep("Neural cell", hca.sc$celltype)] <- "Neural.cell"
hca.sc$celltype[grep("Mural cell", hca.sc$celltype)] <- "Mural.cell"
hca.sc$celltype[grep("Endothelial cell", hca.sc$celltype)] <- "Endothelial.cell"
hca.sc$celltype[grep("Mesothelial cell", hca.sc$celltype)] <- "Mesothelial.cell"
hca.sc$celltype[grep("Lymphatic Endothelial cell", hca.sc$celltype)] <- "Lymphatic.endothelial.cell"
hca.sc$celltype[grep("Ventricular Cardiomyocyte", hca.sc$celltype)] <- "Ventricular.cardiomyocyte"
hca.sc$celltype[grep("Atrial Cardiomyocyte", hca.sc$celltype)] <- "Atrial.cardiomyocyte"

hca.sc$celltype_region <- paste(hca.sc$celltype, hca.sc$region, sep = "_")

saveRDS(hca.sc, file = "~filepath/HCA_v1_sc_processed_filtered.rds")

hca.sn$celltype <- as.character(hca.sn$scReClassify_label)
hca.sn$celltype[grep("Mast cell", hca.sn$celltype)] <- "Mast.cell"
hca.sn$celltype[grep("Neural cell", hca.sn$celltype)] <- "Neural.cell"
hca.sn$celltype[grep("Mural cell", hca.sn$celltype)] <- "Mural.cell"
hca.sn$celltype[grep("Endothelial cell", hca.sn$celltype)] <- "Endothelial.cell"
hca.sn$celltype[grep("Mesothelial cell", hca.sn$celltype)] <- "Mesothelial.cell"
hca.sn$celltype[grep("Lymphatic Endothelial cell", hca.sn$celltype)] <- "Lymphatic.endothelial.cell"
hca.sn$celltype[grep("Ventricular Cardiomyocyte", hca.sn$celltype)] <- "Ventricular.cardiomyocyte"
hca.sn$celltype[grep("Atrial Cardiomyocyte", hca.sn$celltype)] <- "Atrial.cardiomyocyte"

hca.sn$celltype_region <- paste(hca.sn$celltype, hca.sn$region, sep = "_")

saveRDS(hca.sn, file = "~filepath/HCA_v1_sn_processed_filtered.rds")