# load libraries

library(Matrix)
library(purrr)
library(Seurat)
library(readr)
library(dplyr)


# set seed

set.seed(100)


# load the organoid datasets and convert to Seurat objects:


#### Giacomelli 2020 ######################################

# load replicate 1

expr_matrix1 <- readMM("~/Cardiac_organoids/raw_data/Giacomelli/GSM4443384_MTC_S2_L007_matrix.mtx.gz")

barcodes1 <- read.table("~/Cardiac_organoids/raw_data/Giacomelli/GSM4443384_MTC_S2_L007_barcodes.tsv.gz", header = FALSE, sep = "\t")
colnames(expr_matrix1) <- barcodes1$V1

features1 <- read.delim("~/Cardiac_organoids/raw_data/Giacomelli/GSM4443384_MTC_S2_L007_genes.tsv.gz", header = FALSE, stringsAsFactors = FALSE)
rownames(expr_matrix1) <- make.unique(features1$V2)

giacomelli_cmec_cluster1 <- read.csv("~/Cardiac_organoids/raw_data/Giacomelli/GSM4443384_MTC_S2_L007_cluster.tsv.gz", sep="")
metadata1 <- giacomelli_cmec_cluster1 %>%
  dplyr::select(cell.id, celltype) %>%
  dplyr::distinct()
rownames(metadata1) <- sub("^[^_]*_", "", metadata1$cell.id)
metadata1$cell.id <- NULL
metadata1$batch <- 1

expr_matrix1 <- expr_matrix1[, colnames(expr_matrix1) %in% rownames(metadata1)]
identical(rownames(metadata1), colnames(expr_matrix1))
rownames(metadata1) <- paste0("Batch1_", sub("^[^_]*_", "", rownames(metadata1))) 
colnames(expr_matrix1) <- paste0("Batch1_", sub("^[^_]*_", "", colnames(expr_matrix1))) # add the batch number back in to avoid duplicate barcode names between replicates


# load replicate 2

expr_matrix2 <- readMM("~/Cardiac_organoids/raw_data/Giacomelli/GSM4443385_CMEC_S2_L002_matrix.mtx.gz")

barcodes2 <- read.table("~/Cardiac_organoids/raw_data/Giacomelli/GSM4443385_CMEC_S2_L002_barcodes.tsv.gz", header = FALSE, sep = "\t")
colnames(expr_matrix2) <- barcodes2$V1

features2 <- read.delim("~/Cardiac_organoids/raw_data/Giacomelli/GSM4443385_CMEC_S2_L002_genes.tsv.gz", header = FALSE, stringsAsFactors = FALSE)
rownames(expr_matrix2) <- make.unique(features2$V2)

giacomelli_cmec_cluster2 <- read.csv("~/Cardiac_organoids/raw_data/Giacomelli/GSM4443385_CMEC_S2_L002_cluster.tsv.gz", sep="")
metadata2 <- giacomelli_cmec_cluster2 %>%
  dplyr::select(cell.id, celltype) %>%
  dplyr::distinct()
rownames(metadata2) <- sub("^[^_]*_", "", metadata2$cell.id)
metadata2$cell.id <- NULL
metadata2$batch <- 2

expr_matrix2 <- expr_matrix2[, colnames(expr_matrix2) %in% rownames(metadata2)]
identical(rownames(metadata2), colnames(expr_matrix2))
rownames(metadata2) <- paste0("Batch2_", sub("^[^_]*_", "", rownames(metadata2))) 
colnames(expr_matrix2) <- paste0("Batch2_", sub("^[^_]*_", "", colnames(expr_matrix2))) # add the batch number back in to avoid duplicate barcode names between replicates


# combine both replicates

identical(rownames(expr_matrix1), rownames(expr_matrix2))
expr_matrix <- cbind(expr_matrix1, expr_matrix2)
metadata <- rbind(metadata1, metadata2)
identical(rownames(metadata), colnames(expr_matrix))

giacomelli <- CreateSeuratObject(counts = expr_matrix,
                                 assay = "RNA",
                                 meta.data = metadata)



#### Zhang 2022 ######################################

expr_matrix <- readMM("~/Cardiac_organoids/raw_data/Zhang/GSM5071427_MCO_Mix.matrix.mtx.gz")

barcodes <- read.table("~/Cardiac_organoids/raw_data/Zhang/GSM5071427_MCO_Mix.barcodes.tsv.gz", header = FALSE, sep = "\t")
colnames(expr_matrix) <- barcodes$V1

features <- read.delim("~/Cardiac_organoids/raw_data/Zhang/GSM5071427_MCO_Mix.genes.tsv.gz", header = FALSE, stringsAsFactors = FALSE)
rownames(expr_matrix) <- make.unique(features$V2)

zhang <- CreateSeuratObject(counts = expr_matrix, assay = "RNA")



#### Abilez 2025 ######################################

expr_matrix1 <- readMM("~/Cardiac_organoids/raw_data/Abilez/matrix-3.mtx")
expr_matrix2 <- readMM("~/Cardiac_organoids/raw_data/Abilez/matrix-2.mtx")
expr_matrix3 <- readMM("~/Cardiac_organoids/raw_data/Abilez/matrix.mtx")

barcodes1 <- read.table("~/Cardiac_organoids/raw_data/Abilez/barcodes.tsv", header = FALSE, sep = "\t")
barcodes2 <- read.table("~/Cardiac_organoids/raw_data/Abilez/barcodes-2.tsv", header = FALSE, sep = "\t")
barcodes3 <- read.table("~/Cardiac_organoids/raw_data/Abilez/barcodes-3.tsv", header = FALSE, sep = "\t")

colnames(expr_matrix1) <- paste0("Batch1_", barcodes1$V1)
colnames(expr_matrix2) <- paste0("Batch2_", barcodes2$V1)
colnames(expr_matrix3) <- paste0("Batch3_", barcodes3$V1)

features <- read.delim("~/Cardiac_organoids/raw_data/Abilez/features.tsv", header = FALSE, stringsAsFactors = FALSE)
rownames(expr_matrix1) <- make.unique(features$V2)
rownames(expr_matrix2) <- make.unique(features$V2)
rownames(expr_matrix3) <- make.unique(features$V2)

# combine replicates
identical(rownames(expr_matrix1), rownames(expr_matrix2))
identical(rownames(expr_matrix1), rownames(expr_matrix3))
expr_matrix <- cbind(expr_matrix1, expr_matrix2)
expr_matrix <- cbind(expr_matrix, expr_matrix3)

metadata <- data.frame(
  batch = sub("_.*$", "", colnames(expr_matrix)),
  row.names = colnames(expr_matrix)
)

abilez <- CreateSeuratObject(counts = expr_matrix,
                             assay = "RNA",
                             meta.data = metadata)



#### Ho 2022 ######################################

expr_matrix <- read.delim("~/Cardiac_organoids/raw_data/Ho/GSE168464_CCO_Project_merged_UMI_counts.txt", row.names = 1)

day <- sub("^.*(d[0-9]+)$", "\\1", colnames(expr_matrix))

metadata <- data.frame(
  day = day,
  row.names = colnames(expr_matrix)
)

ho <- CreateSeuratObject(
  counts = expr_matrix,
  assay = "RNA",
  meta.data = metadata
)



#### Schmidt 2023 ######################################

schmidt <- readRDS("~/Cardiac_organoids/raw_data/Schmidt/GSE239890_seurat_object.rds")



#### Zou 2026 ######################################

zou_dir <- path.expand("~/Cardiac_organoids/raw_data/Zou")

matrix_files <- list.files(
  zou_dir,
  pattern = "_matrix\\.mtx\\.gz$",
  full.names = TRUE
)

sample_ids <- sub("_matrix\\.mtx\\.gz$", "", basename(matrix_files))
names(matrix_files) <- sample_ids

read_one_10x_flat <- function(sample_id) {
  mtx_file <- file.path(zou_dir, paste0(sample_id, "_matrix.mtx.gz"))
  bc_file  <- file.path(zou_dir, paste0(sample_id, "_barcodes.tsv.gz"))
  ft_file  <- file.path(zou_dir, paste0(sample_id, "_features.tsv.gz"))
  
  stopifnot(file.exists(mtx_file), file.exists(bc_file), file.exists(ft_file))
  
  mat <- readMM(mtx_file)
  
  barcodes <- read.delim(
    bc_file,
    header = FALSE,
    stringsAsFactors = FALSE
  )
  
  features <- read.delim(
    ft_file,
    header = FALSE,
    stringsAsFactors = FALSE
  )
  
  gene_names <- if (ncol(features) >= 2) features$V2 else features$V1
  
  rownames(mat) <- make.unique(gene_names)
  colnames(mat) <- paste0(sample_id, "_", barcodes$V1)
  
  mat
}

mat_list <- lapply(sample_ids, read_one_10x_flat)
names(mat_list) <- sample_ids

expr_matrix <- do.call(cbind, mat_list)

metadata <- data.frame(
  sample_id = sub("_[^_]+$", "", colnames(expr_matrix)),
  row.names = colnames(expr_matrix)
)

metadata$gsm <- sub("_.*$", "", metadata$sample_id)
metadata$sample_name <- sub("^GSM[0-9]+_", "", metadata$sample_id)

metadata$condition <- ifelse(
  grepl("D25\\.5", metadata$sample_name),
  "SCPD_D25.5",
  ifelse(
    grepl("^SCPD", metadata$sample_name),
    "SCPD_D42.5",
    ifelse(
      grepl("^SCWOPD", metadata$sample_name),
      "SCWOPD_D42.5",
      NA_character_
    )
  )
)

metadata$replicate <- ifelse(
  grepl("-[0-9]+$", metadata$sample_name),
  sub(".*-([0-9]+)$", "\\1", metadata$sample_name),
  NA
)

identical(rownames(metadata), colnames(expr_matrix))

zou <- CreateSeuratObject(
  counts = expr_matrix,
  assay = "RNA",
  meta.data = metadata,
  project = "Zou"
)



#### Feng 2022 ######################################

expr_matrix <- readMM("~/Cardiac_organoids/raw_data/Feng/GSE163619_GEO_unsorted_counts.mtx")

barcodes <-read.table("~/Cardiac_organoids/raw_data/Feng/GSE163619_GEO_unsorted_barcodes.txt", header = FALSE, sep = "\t")

colnames(expr_matrix) <- make.unique(barcodes$V1)

features <- read.delim("~/Cardiac_organoids/raw_data/Feng/GSE163619_GEO_unsorted_genes.txt", header = FALSE, stringsAsFactors = FALSE)

rownames(expr_matrix) <- features$V1

metadata <- read_csv("~/Cardiac_organoids/raw_data/Feng/GSE163619_GEO_unsorted_metadata.csv")
metadata <- as.data.frame(metadata)
rownames(metadata) <- colnames(expr_matrix)

feng <- CreateSeuratObject(counts = expr_matrix,
                           assay = "RNA",
                           meta.data = metadata)

# remove monolayer samples, disease cell lines, and cell line WTC
feng <- feng[,feng$cell.line == "SCVI114"]
feng <- feng[,feng$differentiation.method == "Org"]



#### O'Hern 2025 ######################################

raw_dir <- path.expand("~/Cardiac_organoids/raw_data/OHern")
extract_base <- file.path(raw_dir, "extracted_cellranger")

tar_files <- list.files(
  raw_dir,
  pattern = "^GSM.*_cellranger\\.tar\\.gz$",
  full.names = TRUE
)

tar_sample_ids <- basename(tar_files) |>
  sub("_cellranger\\.tar\\.gz$", "", x = _) |>
  sub("^GSM[0-9]+_", "", x = _)

names(tar_files) <- tar_sample_ids

extracted_dirs <- list.dirs(
  extract_base,
  recursive = FALSE,
  full.names = TRUE
)

extracted_sample_ids <- basename(extracted_dirs)

sample_ids <- union(tar_sample_ids, extracted_sample_ids)

find_10x_matrix <- function(extract_dir) {
  all_files <- list.files(extract_dir, recursive = TRUE, full.names = TRUE)
  
  h5 <- all_files[grepl("filtered_feature_bc_matrix\\.h5$|raw_feature_bc_matrix\\.h5$", all_files)]
  if (length(h5) > 0) {
    return(list(type = "h5", path = h5[1]))
  }
  
  dirs <- unique(dirname(all_files))
  
  has_mtx <- file.exists(file.path(dirs, "matrix.mtx")) |
    file.exists(file.path(dirs, "matrix.mtx.gz"))
  
  has_barcodes <- file.exists(file.path(dirs, "barcodes.tsv")) |
    file.exists(file.path(dirs, "barcodes.tsv.gz"))
  
  has_features <- file.exists(file.path(dirs, "features.tsv")) |
    file.exists(file.path(dirs, "features.tsv.gz")) |
    file.exists(file.path(dirs, "genes.tsv")) |
    file.exists(file.path(dirs, "genes.tsv.gz"))
  
  hits <- dirs[has_mtx & has_barcodes & has_features]
  
  if (length(hits) == 0) {
    stop("No 10x matrix found inside: ", extract_dir)
  }
  
  filtered <- hits[grepl("filtered_feature_bc_matrix", hits)]
  if (length(filtered) > 0) hits <- filtered
  
  list(type = "dir", path = hits[1])
}

load_cellranger_sample <- function(sample_id) {
  extract_dir <- file.path(extract_base, sample_id)
  
  tar_file <- if (sample_id %in% names(tar_files)) {
    tar_files[[sample_id]]
  } else {
    NA_character_
  }
  
  dir.create(extract_dir, recursive = TRUE, showWarnings = FALSE)
  
  if (length(list.files(extract_dir, recursive = TRUE)) == 0) {
    if (is.na(tar_file) || !file.exists(tar_file)) {
      stop("No tar file and no extracted data found for sample: ", sample_id)
    }
    
    untar(tar_file, exdir = extract_dir)
  }
  
  mat_info <- find_10x_matrix(extract_dir)
  
  counts <- if (mat_info$type == "h5") {
    Read10X_h5(mat_info$path)
  } else {
    Read10X(data.dir = mat_info$path)
  }
  
  if (is.list(counts)) {
    counts <- counts[["Gene Expression"]] %||% counts[["RNA"]] %||% counts[[1]]
  }
  
  seu <- CreateSeuratObject(
    counts = counts,
    project = sample_id
  )
  
  seu$sample_id <- sample_id
  seu$condition <- sample_id
  
  seu
}

seu_list <- lapply(sample_ids, load_cellranger_sample)
names(seu_list) <- sample_ids

ohern <- merge(
  x = seu_list[[1]],
  y = seu_list[-1],
  add.cell.ids = names(seu_list),
  project = "OHern"
)

ohern$condition <- ohern$sample_id
ohern <- JoinLayers(ohern, assay = "RNA")



#### Pocock 2025 ######################################

# load replicate 1

expr_matrix1 <- readMM("~/Cardiac_organoids/raw_data/Pocock/GSM8739294_DM_1_matrix.mtx.gz")

barcodes1 <- read.table("~/Cardiac_organoids/raw_data/Pocock/GSM8739294_DM_1_barcodes.tsv.gz", header = FALSE, sep = "\t")
colnames(expr_matrix1) <- paste0("Batch1_", barcodes1$V1)

features1 <- read.delim("~/Cardiac_organoids/raw_data/Pocock/GSM8739294_DM_1_features.tsv.gz", header = FALSE, stringsAsFactors = FALSE)
rownames(expr_matrix1) <- make.unique(features1$V2)

metadata1 <- data.frame(batch = 1, V1 = barcodes1)
rownames(metadata1) <- paste0("Batch1_", metadata1$V1)
metadata1$V1 <- NULL

# load replicate 2

expr_matrix2 <- readMM("~/Cardiac_organoids/raw_data/Pocock/GSM8739295_DM_2_matrix.mtx.gz")

barcodes2 <- read.table("~/Cardiac_organoids/raw_data/Pocock/GSM8739295_DM_2_barcodes.tsv.gz", header = FALSE, sep = "\t")
colnames(expr_matrix2) <- paste0("Batch2_", barcodes2$V1)

features2 <- read.delim("~/Cardiac_organoids/raw_data/Pocock/GSM8739295_DM_2_features.tsv.gz", header = FALSE, stringsAsFactors = FALSE)
rownames(expr_matrix2) <- make.unique(features2$V2)

metadata2 <- data.frame(batch = 2, V1 = barcodes2)
rownames(metadata2) <- paste0("Batch2_", metadata2$V1)
metadata2$V1 <- NULL


# combine both replicates

identical(rownames(expr_matrix1), rownames(expr_matrix2))
expr_matrix <- cbind(expr_matrix1, expr_matrix2)
metadata <- rbind(metadata1, metadata2)
identical(rownames(metadata), colnames(expr_matrix))

pocock <- CreateSeuratObject(counts = expr_matrix,
                             assay = "RNA",
                             meta.data = metadata)



#### Voges 2023 ######################################

# load replicate 1

expr_matrix1 <- readMM("~/Cardiac_organoids/raw_data/Pocock/GSM8739292_SF_1_matrix.mtx.gz")

barcodes1 <- read.table("~/Cardiac_organoids/raw_data/Pocock/GSM8739292_SF_1_barcodes.tsv.gz", header = FALSE, sep = "\t")
colnames(expr_matrix1) <- paste0("Batch1_", barcodes1$V1)

features1 <- read.delim("~/Cardiac_organoids/raw_data/Pocock/GSM8739292_SF_1_features.tsv.gz", header = FALSE, stringsAsFactors = FALSE)
rownames(expr_matrix1) <- make.unique(features1$V2)

metadata1 <- data.frame(batch = 1, V1 = barcodes1)
rownames(metadata1) <- paste0("Batch1_", metadata1$V1)
metadata1$V1 <- NULL

# load replicate 2

expr_matrix2 <- readMM("~/Cardiac_organoids/raw_data/Pocock/GSM8739293_SF_2_matrix.mtx.gz")

barcodes2 <- read.table("~/Cardiac_organoids/raw_data/Pocock/GSM8739293_SF_2_barcodes.tsv.gz", header = FALSE, sep = "\t")
colnames(expr_matrix2) <- paste0("Batch2_", barcodes2$V1)

features2 <- read.delim("~/Cardiac_organoids/raw_data/Pocock/GSM8739293_SF_2_features.tsv.gz", header = FALSE, stringsAsFactors = FALSE)
rownames(expr_matrix2) <- make.unique(features2$V2)

metadata2 <- data.frame(batch = 2, V1 = barcodes2)
rownames(metadata2) <- paste0("Batch2_", metadata2$V1)
metadata2$V1 <- NULL


# combine both replicates

identical(rownames(expr_matrix1), rownames(expr_matrix2))
expr_matrix <- cbind(expr_matrix1, expr_matrix2)
metadata <- rbind(metadata1, metadata2)
identical(rownames(metadata), colnames(expr_matrix))

voges <- CreateSeuratObject(counts = expr_matrix,
                            assay = "RNA",
                            meta.data = metadata)



#### Landau 2024 ######################################

expr_matrix_mac <- readMM("~/Cardiac_organoids/raw_data/Landau/GSM8196011_filtered_2A_matrix.mtx.gz")
expr_matrix_non <- readMM("~/Cardiac_organoids/raw_data/Landau/GSM8196012_filtered_3A_matrix.mtx.gz")

barcodes_mac <-read.table("~/Cardiac_organoids/raw_data/Landau/GSM8196011_filtered_2A_barcodes.tsv.gz", header = FALSE, sep = "\t")
barcodes_non <-read.table("~/Cardiac_organoids/raw_data/Landau/GSM8196012_filtered_3A_barcodes.tsv.gz", header = FALSE, sep = "\t")

colnames(expr_matrix_mac) <- paste0("mac_", barcodes_mac$V1)
colnames(expr_matrix_non) <- paste0("non_", barcodes_non$V1)

features <- read.delim("~/Cardiac_organoids/raw_data/Landau/GSM8196011_filtered_2A_features.tsv.gz", header = FALSE, stringsAsFactors = FALSE)

rownames(expr_matrix_mac) <- make.unique(features$V2)
rownames(expr_matrix_non) <- make.unique(features$V2)

# combine replicates
identical(rownames(expr_matrix_mac), rownames(expr_matrix_non))
expr_matrix <- cbind(expr_matrix_mac, expr_matrix_non)

metadata <- data.frame(
  condition = sub("_.*$", "", colnames(expr_matrix)),
  row.names = colnames(expr_matrix)
)

landau <- CreateSeuratObject(counts = expr_matrix,
                             assay = "RNA",
                             meta.data = metadata)







# cell filtering and processing for single cell datasets:

seu_list_sc <- list(giacomelli = giacomelli,
                    zhang = zhang,
                    schmidt = schmidt,
                    ho = ho,
                    abilez = abilez,
                    feng = feng,
                    zou = zou,
                    ohern = ohern)

# set thresholds
gene_pct_threshold <- 2
mt_pct_threshold <- 20

# add percent mitochondrial and percent genes quantified per dataset
seu_list_sc <- imap(seu_list_sc, function(seu, nm) {

  assay_use <- if ("RNA" %in% Assays(seu)) "RNA" else DefaultAssay(seu)
  
  counts <- GetAssayData(seu, assay = assay_use, layer = "counts")
  
  n_genes_quantified <- sum(Matrix::rowSums(counts > 0) > 0)
  
  mito_pattern <- if (any(grepl("^MT-", rownames(seu[[assay_use]])))) {
    "^MT-"
  } else {
    "^mt-"
  }
  
  seu[["percent.mt"]] <- PercentageFeatureSet(
    seu,
    assay = assay_use,
    pattern = mito_pattern
  )
  
  seu[["pct.genes.quantified"]] <-
    100 * seu@meta.data[[paste0("nFeature_", assay_use)]] / n_genes_quantified
  
  seu[["n.genes.quantified.dataset"]] <- n_genes_quantified
  
  seu
})

# filter cells within each dataset using that dataset's own gene denominator
seu_list_sc_filtered <- imap(seu_list_sc, function(seu, nm) {
  subset(
    seu,
    subset = pct.genes.quantified >= gene_pct_threshold &
      percent.mt <= mt_pct_threshold
  )
})

seu_list_sc_filtered <- lapply(names(seu_list_sc_filtered), function(sample_id) {
  message("Processing ", sample_id)
  
  seu <- seu_list_sc_filtered[[sample_id]]
  
  seu <- NormalizeData(seu)
  seu <- FindVariableFeatures(seu)
  seu <- ScaleData(seu)
  seu <- RunPCA(seu)
  seu <- RunUMAP(seu, reduction = "pca", dims = 1:50)
  
  seu
})





# cell filtering and processing for single nuclei datasets:

seu_list_sn <- list(pocock = pocock,
                    voges = voges,
                    landau = landau)

# set thresholds
gene_pct_threshold <- 2
mt_pct_threshold <- 5

# add percent mitochondrial and percent genes quantified per dataset
seu_list_sn <- imap(seu_list_sn, function(seu, nm) {
  assay_use <- if ("RNA" %in% Assays(seu)) "RNA" else DefaultAssay(seu)
  
  counts <- GetAssayData(seu, assay = assay_use, layer = "counts")
  
  n_genes_quantified <- sum(Matrix::rowSums(counts > 0) > 0)
  
  mito_pattern <- if (any(grepl("^MT-", rownames(seu[[assay_use]])))) {
    "^MT-"
  } else {
    "^mt-"
  }
  
  seu[["percent.mt"]] <- PercentageFeatureSet(
    seu,
    assay = assay_use,
    pattern = mito_pattern
  )
  
  seu[["pct.genes.quantified"]] <-
    100 * seu@meta.data[[paste0("nFeature_", assay_use)]] / n_genes_quantified
  
  seu[["n.genes.quantified.dataset"]] <- n_genes_quantified
  
  seu
})

# filter cells within each dataset using that dataset's own gene denominator
seu_list_sn_filtered <- imap(seu_list_sn, function(seu, nm) {
  subset(
    seu,
    subset = pct.genes.quantified >= gene_pct_threshold &
      percent.mt <= mt_pct_threshold
  )
})

seu_list_sn_filtered <- lapply(names(seu_list_sn_filtered), function(sample_id) {
  message("Processing ", sample_id)
  
  seu <- seu_list_sn_filtered[[sample_id]]
  
  seu <- NormalizeData(seu)
  seu <- FindVariableFeatures(seu)
  seu <- ScaleData(seu)
  seu <- RunPCA(seu)
  seu <- RunUMAP(seu, reduction = "pca", dims = 1:50)
  
  seu
})