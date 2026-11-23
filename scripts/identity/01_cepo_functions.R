# load libraries

library(Cepo)
library(SeuratObject)
library(dplyr)
library(matrixStats)
library(purrr)


# define Cepo functions

multiCepo <- function(seu, batch_var, celltype_var, min_cells=20, ncores=5, exprs_pct=NULL) {
  
  # helper function to compute cepo when there is more than one cell type per condition after min cell filtering
  singleCepo <- function(seu, celltype_var=celltype_var, min_cells=min_cells, ncores=ncores) {
    
    celltype_col <- seu@meta.data[[celltype_var]]
    valid_ct <- names(which(table(celltype_col) > min_cells))
    ct_idx <- celltype_col %in% valid_ct
    
    if (sum(ct_idx) == 0) return(NULL) # added line
    
    seu <- seu[, ct_idx]
    
    if (length(unique(seu@meta.data[[celltype_var]])) > 1) {
      mat <- GetAssayData(seu, assay = "RNA", layer = "data")
      cepo <- Cepo(mat, cellType = seu@meta.data[[celltype_var]], minCells = min_cells, workers = ncores,
                   exprsPct = exprs_pct)
      return(as.data.frame(cepo$stats))
    } else {
      return(NULL)
    }
  }
  
  stats_list = lapply(unique(seu@meta.data[[batch_var]]), function(batch) {
    
    batch_idx = seu@meta.data[[batch_var]] == batch
    filt <- seu[, batch_idx]
    
    singleCepo(seu=filt, celltype_var=celltype_var, min_cells=min_cells, ncores=ncores)
  })
  names(stats_list) = unique(seu@meta.data[[batch_var]])
  stats_list
}

averageCepoStats <- function(cepo_list) {
  
  cepo_list <- cepo_list[!sapply(cepo_list, is.null)]
  
  common_genes <- Reduce(intersect, lapply(cepo_list, rownames))
  all_ct <- unique(Reduce(c, lapply(cepo_list, colnames)))
  
  by_ct <- lapply(all_ct, function(ct) {
    ct_vals <- lapply(cepo_list, function(x) if (ct %in% colnames(x)) x[common_genes, ct])
    
    if (length(ct_vals) > 1)
      matrixStats::rowMeans2(do.call(cbind, ct_vals))
    else
      do.call(cbind, ct_vals) # make sure it's numeric
  })
  
  result <- do.call(cbind, by_ct)
  colnames(result) <- all_ct
  rownames(result) <- common_genes
  return(result)
}

run_cepo <- function(seu_obj, min_cells = 15) {
  message("Running single Cepo...")
  
  cepo_res <- Cepo(
    seu_obj@assays$RNA$data,
    cellTypes = seu_obj$celltype,
    minCells = min_cells,
    exprsPct = 0.05
  )
  
  as.matrix(cepo_res$stats)
}

subsample_organoid <- function(
    seu_obj,
    frac = 0.8,
    celltype_var = "celltype"
) {
  meta <- seu_obj[[]]
  meta$cell <- rownames(meta)
  
  sampled_cells <- meta %>%
    dplyr::mutate(celltype_for_sampling = as.character(.data[[celltype_var]])) %>%
    dplyr::group_by(celltype_for_sampling) %>%
    dplyr::group_modify(function(.x, .y) {
      n <- nrow(.x)
      n_keep <- ceiling(frac * n)
      
      .x[base::sample(seq_len(n), size = n_keep, replace = FALSE), , drop = FALSE]
    }) %>%
    dplyr::ungroup() %>%
    dplyr::pull(cell)
  
  subset(seu_obj, cells = sampled_cells)
}

bootstrap_organoid_cepo <- function(
    organoids,
    reps = 10,
    frac = 0.8,
    celltype_var = "celltype",
    min_cells = 15,
    seed = 2026
) {
  organoid_cepo_bootstrap <- purrr::imap(
    organoids,
    function(seu_obj, organoid_name) {
      message("Bootstrapping organoid: ", organoid_name)
      
      org_i <- match(organoid_name, names(organoids))
      
      rep_results <- lapply(seq_len(reps), function(rep_i) {
        message("  Replicate ", rep_i, "/", reps)
        
        set.seed(seed + org_i * 1000 + rep_i)
        
        seu_sub <- subsample_organoid(
          seu_obj = seu_obj,
          frac = frac,
          celltype_var = celltype_var
        )
        
        run_cepo(
          seu_obj = seu_sub,
          min_cells = min_cells
        )
      })
      
      names(rep_results) <- paste0("rep", seq_len(reps))
      rep_results
    }
  )
  
  organoid_cepo_bootstrap
}