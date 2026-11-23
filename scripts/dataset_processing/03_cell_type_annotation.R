# load libraries

library(scClassify)
library(BiocParallel)
library(SingleR)
library(Seurat)


# train scClassify on HCA sc:

hca.sc <- readRDS("~/Cardiac_organoids/processed_data/HCA_v1_sc_processed_filtered.rds")

# cap the cell type counts
seu <- hca.sc
celltypes <- unique(seu$celltype)
Idents(seu) <- seu$celltype
cells_to_keep <- c()

for (celltype in celltypes) {
  
  cells_of_type <- WhichCells(seu, idents = celltype)
  
  if (length(cells_of_type) > 500) {
    cells_of_type <- sample(cells_of_type, 500)
  }
  cells_to_keep <- c(cells_to_keep, cells_of_type)
}

if (length(cells_to_keep) > 0) {
  seu <- subset(seu, cells = cells_to_keep)
}

# train model
model_hca_sc <- train_scClassify(
  exprsMat_train = seu@assays$RNA$data,
  cellTypes_train = seu$celltype,
  tree = "HOPACH",
  selectFeatures = "limma",
  topN = 50,
  hopach_kmax = 5,
  pSig = 0.05,
  cellType_tree = NULL,
  weightsCal = FALSE,
  parallel = TRUE,
  BPPARAM = MulticoreParam(workers=15),
  verbose = TRUE,
  returnList = FALSE
)


# train scClassify on HCA sn:

hca.sn <- readRDS("~/Cardiac_organoids/processed_data/HCA_v1_sn_processed_filtered.rds")

# cap the cell type counts
seu <- hca.sn
celltypes <- unique(seu$celltype)
Idents(seu) <- seu$celltype
cells_to_keep <- c()

for (celltype in celltypes) {
  
  cells_of_type <- WhichCells(seu, idents = celltype)
  
  if (length(cells_of_type) > 500) {
    cells_of_type <- sample(cells_of_type, 500)
  }
  cells_to_keep <- c(cells_to_keep, cells_of_type)
}

if (length(cells_to_keep) > 0) {
  seu <- subset(seu, cells = cells_to_keep)
}

# train model
model_hca_sn <- train_scClassify(
  exprsMat_train = seu@assays$RNA$data,
  cellTypes_train = seu$celltype,
  tree = "HOPACH",
  selectFeatures = "limma",
  topN = 50,
  hopach_kmax = 5,
  pSig = 0.05,
  cellType_tree = NULL,
  weightsCal = FALSE,
  parallel = TRUE,
  BPPARAM = MulticoreParam(workers=15),
  verbose = TRUE,
  returnList = FALSE
)


# run cell type annotation using scClassify, Seurat, and SingleR

annotate_one_seurat <- function(obj, ref = hca.sc, scclassify_model = model_hca_sc) {
  
  # scClassify
  
  sc_preds <- predict_scClassify(
    obj@assays$RNA$data,
    scclassify_model,
    k = 10,
    prob_threshold = 0.7,
    cor_threshold_static = 0.5,
    cor_threshold_high = 0.7,
    features = "limma",
    algorithm = "WKNN",
    similarity = "pearson",
    cutoff_method = c("dynamic", "static"),
    parallel = TRUE,
    BPPARAM = MulticoreParam(workers = 15),
    verbose = TRUE
  )
  
  obj$scClassify_label <- sc_preds$pearson_WKNN_limma$predLabelMat[, 4]
  
  # Seurat label transfer

  common_features <- intersect(
    rownames(ref@assays$RNA$data),
    rownames(obj@assays$RNA$data)
  )
  
  print(length(common_features))
  
  anchors <- FindTransferAnchors(
    reference = ref,
    query = obj,
    query.assay = "RNA",
    normalization.method = "LogNormalize",
    reference.reduction = "pca",
    dims = 1:30,
    features = common_features
  )
  
  transfer_preds <- TransferData(
    anchorset = anchors,
    refdata = ref$celltype,
    dims = 1:30
  )
  
  colnames(transfer_preds)[1] <- "seurat_transfer_label_raw"
  colnames(transfer_preds)[which(colnames(transfer_preds) == "prediction.score.max")] <- "seurat.score"
  
  obj <- AddMetaData(obj, transfer_preds)
  
  # SingleR 
  
  ref_labels <- ref$celltype
  names(ref_labels) <- colnames(ref)
  
  common_genes <- intersect(rownames(ref), rownames(obj))
  
  ref_expr <- ref@assays$RNA$data[common_genes, , drop = FALSE]
  test_expr <- obj@assays$RNA$data[common_genes, , drop = FALSE]
  
  singleR <- SingleR(
    test = test_expr,
    ref = ref_expr,
    labels = ref_labels,
    de.method = "wilcox",
    fine.tune = TRUE,
    prune = TRUE
  )
  
  singleR_df <- as.data.frame(singleR)
  
  obj$singleR_label_raw <- singleR_df[colnames(obj), "labels"]
  obj$singleR_label_pruned <- singleR_df[colnames(obj), "pruned.labels"]
  obj$singleR_delta_next <- singleR_df[colnames(obj), "delta.next"]
  
  return(obj)
}

seu_list_sc_filtered_annotated <- lapply(
  seu_list_sc_filtered,
  annotate_one_seurat,
  ref = hca.sc,
  scclassify_model = model_hca_sc
)

seu_list_sn_filtered_annotated <- lapply(
  seu_list_sn_filtered,
  annotate_one_seurat,
  ref = hca.sn,
  scclassify_model = model_hca_sn
)


# create consensus labels:

consensus_cols <- c(
  "scClassify_label",
  "seurat_transfer_label_raw",
  "singleR_label_raw"
)

make_consensus_celltype <- function(seu, cols = consensus_cols) {
  missing_cols <- setdiff(cols, colnames(seu@meta.data))
  
  consensus <- apply(seu@meta.data[, cols, drop = FALSE], 1, function(x) {
    x <- as.character(x)
    x <- x[!is.na(x) & x != ""]
    
    tab <- table(x)
    
    if (length(tab) > 0 && max(tab) >= 2) {
      names(tab)[which.max(tab)]
    } else {
      "unassigned"
    }
  })
  
  seu$celltype <- consensus
  seu
}

seu_list_sc_filtered_annotated <- lapply(seu_list_sc_filtered_annotated, make_consensus_celltype)
seu_list_sn_filtered_annotated <- lapply(seu_list_sn_filtered_annotated, make_consensus_celltype)


# merge cardiomyocyte labels

cm_labels <- c("Atrial.cardiomyocyte", "Ventricular.cardiomyocyte")

seu_list_sc_filtered_annotated <- lapply(seu_list_sc_filtered_annotated, function(seu) {
  if (!"celltype" %in% colnames(seu@meta.data)) {
    stop("Missing metadata column: celltype")
  }
  seu$celltype <- as.character(seu$celltype)
  seu$celltype[seu$celltype %in% cm_labels] <- "Cardiomyocyte"
  seu
})

seu_list_sn_filtered_annotated <- lapply(seu_list_sn_filtered_annotated, function(seu) {
  if (!"celltype" %in% colnames(seu@meta.data)) {
    stop("Missing metadata column: celltype")
  }
  seu$celltype <- as.character(seu$celltype)
  seu$celltype[seu$celltype %in% cm_labels] <- "Cardiomyocyte"
  seu
})