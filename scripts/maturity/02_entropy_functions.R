# load libraries

library(Seurat)
library(SeuratObject)


# define functions

master_entropy = function(expressionMatrix, n = 1000, species = "mouse"){
  #This portion splits up large datasets to go easier on memory usage.
  if(ncol(expressionMatrix) > 5000){
    lst = split(colnames(expressionMatrix), (seq(ncol(expressionMatrix))-1) %/% 5000)
    result = unlist(lapply(lst, function(x){master_entropy(expressionMatrix[, x], n = n, species = species)}))
  }
  else{
    if(species == "mouse"){
      expressionMatrix = mito_correct(expressionMatrix) 
    }
    if(species == "human"){
      G_list = G_list_human
    }
    expressionMatrix = expressionMatrix[rownames(expressionMatrix) %in% G_list[G_list$use_for_entropy == TRUE, ]$symbol & !startsWith(rownames(expressionMatrix), "Rps") & !startsWith(rownames(expressionMatrix), "Rpl") & !startsWith(rownames(expressionMatrix), "RPS") & !startsWith(rownames(expressionMatrix), "RPL"), ]
    expressionMatrix = apply(expressionMatrix, 2, function(x) tail(sort(x), n))
    p=sweep(expressionMatrix,MARGIN=2,FUN="/",STATS=colSums(expressionMatrix));
    logp=log(p)
    logp[logp==-Inf]=0
    plogp=p*logp
    result=-1*colSums(plogp)
  }
  return(result)
}

data_qc = function(dataset, study, timepoint_list = NA, scn_calc = TRUE, species = "mouse", sample_type = NA, append_to_alldata = FALSE, append_only = FALSE, isolation = NA, sequencing = NA, mapping = NA, datatype = NA, doi = NA, other_meta = NA){
  data = get(dataset)
  qc = as.data.frame(matrix(nrow = ncol(data), ncol = 27))
  colnames(qc) = c("timepoint", "study", "depth", "genes", "entropy", "top5", "top5_norm", "depth_norm", "mito", "mito_all", "ribo", "cm_score", "max_score", "max_celltype", "good_cell", "data", "cellname", "other_meta", "include_dataset", "reason", "isolation", "sequencing", "mapping", "datatype", "species", "sample_type", "full_label")
  rownames(qc) = colnames(data)
  qc$study = study
  qc$timepoint = timepoint_list
  qc$depth = colSums(data)
  qc$genes = colSums(data > 0)
  if(append_to_alldata == TRUE){
    alldata[study, ] <<- c(isolation, sequencing, mapping, datatype, as.numeric(median(qc$depth)), as.numeric(median(qc$genes)), as.numeric(nrow(qc)), paste(unique(qc$timepoint), collapse = ", "), doi, dataset, species, sample_type)
  }
  if(append_only == FALSE){
    qc$entropy = master_entropy(data, species = species)
    
    qc$top5 = top_n_genes(data, 5)
    qc$top5_norm = 0
    if(!is.na(unique(qc$timepoint))){
      for (j in unique(qc$timepoint)){
        qc[qc$timepoint == j, ]$top5_norm = qc[qc$timepoint == j, ]$top5/median(qc[qc$timepoint == j, ]$top5)
      }
    }
    
    qc$depth_norm = 0
    if (!is.na(unique(qc$timepoint))){
      for (j in unique(qc$timepoint)){
        qc[qc$timepoint == j, ]$depth_norm = log(qc[qc$timepoint == j, ]$depth/median(qc[qc$timepoint == j, ]$depth))
      }
    }
    
    qc$mito = mito(data, species = species)
    qc$mito_all = mito(data, species = species, all = TRUE)
    qc$ribo = ribo(data, species = species)
    
    if(scn_calc == TRUE){
      if(species == "human"){
        scn_data = scn_human(data, qc, study, qc)
      }
      else {
        scn_data = scn(data, qc, study, qc)
      }
      qc$cm_score = scn_data$cm_score
      qc$max_score = scn_data$max_score
      qc$max_celltype = scn_data$max_celltype
    }
    
    #Please note that the gene number isn't hardcoded into the computation of "good cell" - this is because other uses may want to compute other entropy scores with different gene numbers but still relevantly compute a good cell
    qc$good_cell = ((qc$top5_norm < 1.3) & (qc$depth_norm > -0.5) & (qc$max_celltype == "cardiac muscle cell"))
    qc$data = dataset
    qc$cellname = colnames(data)
    qc$include_dataset = TRUE
    qc$isolation = isolation
    qc$sequencing = sequencing
    qc$mapping = mapping
    qc$datatype = datatype
    qc$species = species
    qc$sample_type = sample_type
    qc$full_label = paste(qc$species, qc$sample_type)
    qc$other_meta = other_meta
    return(qc)
  }
}

get_counts_matrix <- function(seu, assay = "RNA") {
  tryCatch(
    GetAssayData(seu, assay = assay, layer = "counts"),
    error = function(e) GetAssayData(seu, assay = assay, slot = "counts")
  )
}

compute_entropy_for_dataset <- function(seu, dataset_name, celltype_col = "celltype", min_cells = 20) {
  if (!celltype_col %in% colnames(seu@meta.data)) {
    stop("Column '", celltype_col, "' not found in metadata for dataset: ", dataset_name)
  }
  
  celltypes <- sort(unique(seu@meta.data[[celltype_col]]))
  celltypes <- celltypes[!is.na(celltypes)]
  
  entropy_df_list <- list()
  
  for (ct in celltypes) {
    message("Processing ", dataset_name, " / ", ct)
    
    cells_use <- colnames(seu)[seu@meta.data[[celltype_col]] == ct]
    
    if (length(cells_use) < min_cells) {
      message(
        "Skipping celltype '", ct, "' in ", dataset_name,
        " (", length(cells_use), " cells < ", min_cells, ")"
      )
      next
    }
    
    subset_seu <- seu[, cells_use]
    expr_mat <- get_counts_matrix(subset_seu, assay = "RNA")
    
    entropy_scores <- master_entropy(
      expr_mat,
      species = "human"
    )
    
    entropy_df_list[[ct]] <- data.frame(
      cell_id = colnames(subset_seu),
      celltype = ct,
      entropy = entropy_scores,
      dataset = dataset_name,
      stringsAsFactors = FALSE
    )
  }
  
  do.call(rbind, entropy_df_list)
}