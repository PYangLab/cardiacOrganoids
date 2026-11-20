#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(BiocParallel)
  library(Matrix)
  library(scClassify)
  library(Seurat)
  library(SingleR)
})

if (packageVersion("Seurat") < "5.0.0") {
  stop("This script requires Seurat >= 5.0")
}


source(file.path("./pipeline_utils.R"))

config <- parse_pipeline_args(
  defaults = list(
    data_dir = "../data",
    out_dir = "./offtarget_results",
    workers = 15L
  )
)

data_dir <- path.expand(config$data_dir)
out_dir <- path.expand(config$out_dir)
workers <- config$workers

if (workers < 1L) {
  stop("--workers must be at least 1")
}

message("Pipeline configuration:")
message("  Data directory: ", data_dir)
message("  Output directory: ", out_dir)
message("  Workers: ", workers)



set.seed(100)
options(future.globals.maxSize = 16 * 1024^3)

manifest_path <- file.path(out_dir, "organoid_input_summary.csv")
reference_path <- file.path(data_dir, "fetal_hca_sc_processed.rds")



required_inputs <- c(
  manifest_path, reference_path
)

if (!all(file.exists(required_inputs))) {
  stop(
    "Step 2 - Missing: ",
    paste(required_inputs[!file.exists(required_inputs)], collapse = ", ")
  )
}

stage_dir <- file.path(out_dir, "02_fetal_hca_classification")
prediction_dir <- file.path(stage_dir, "predictions")
checkpoint_dir <- file.path(stage_dir, "checkpoints")
dir.create(prediction_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

run_contract <- list(
  refined_reference_path = normalizePath(reference_path),
  query_population = "supplied manuscript celltype == 'unassigned' only",
  reference_label = "author_mid_grain (the atlas authors' original fetal mid_grain label)",
  training_cap = "at most 500 reference cells per original author_mid_grain label",
  training_seed = 100L,
  scClassify = list(
    tree = "HOPACH", selectFeatures = "limma", topN = 50L,
    hopach_kmax = 5L, pSig = 0.05, weightsCal = FALSE,
    k = 10L, prob_threshold = 0.7, cor_threshold_static = 0.5,
    cor_threshold_high = 0.7, algorithm = "WKNN",
    similarity = "pearson", cutoff_method = c("dynamic", "static")
  ),
  Seurat = list(
    normalization.method = "LogNormalize",
    reference.reduction = "pca", dims = 1:30
  ),
  SingleR = list(
    de.method = "wilcox", fine.tune = TRUE, prune = TRUE,
    consensus_uses = "raw labels; pruned labels retained as diagnostics"
  ),
  consensus = "exact agreement of at least 2 of the 3 terminal-label calls"
)
run_contract_path <- file.path(stage_dir, "classification_contract.rds")
if (file.exists(run_contract_path)) {
  if (!identical(readRDS(run_contract_path), run_contract)) {
    stop(
      "Existing Step 2 checkpoints were made under a different contract. ",
      "Use a new output directory or inspect/remove only the stale Step 2 run."
    )
  }
} else {
  saveRDS(run_contract, run_contract_path)
}

message("Loading prepared fetal HCA reference; primary labels are author_mid_grain")
reference <- readRDS(reference_path)

required_reference_metadata <- c(
  "mid_grain"
)
missing_reference_metadata <- setdiff(
  required_reference_metadata, colnames(reference@meta.data)
)
if (length(missing_reference_metadata)) {
  stop(
    "Refined fetal reference metadata missing: ",
    paste(missing_reference_metadata, collapse = ", ")
  )
}
reference@meta.data$author_mid_grain <- as.character(reference@meta.data$mid_grain)
if (anyNA(reference@meta.data$author_mid_grain) ||
    any(!nzchar(as.character(reference@meta.data$author_mid_grain)))) {
  stop("Fetal author_mid_grain contains missing or empty labels")
}

if (!"pca" %in% names(reference@reductions) ||
    ncol(Embeddings(reference, "pca")) < 30L) {
  stop("Refined fetal reference must contain at least 30 PCA dimensions")
}

reference_data <- get_rna_layer(reference, "data")
reference_cells <- colnames(reference_data)

reference_labels <- as.character(reference@meta.data$author_mid_grain)
names(reference_labels) <- reference_cells
label_vocabulary <- sort(unique(reference_labels))


bp <- if (.Platform$OS.type == "unix" && workers > 1L) {
  MulticoreParam(workers = workers)
} else {
  SerialParam()
}

model_path <- file.path(checkpoint_dir, "fetal_author_mid_grain_scClassify_model.rds")
training_cells_path <- file.path(
  stage_dir, "fetal_author_mid_grain_scClassify_training_cells.csv"
)
if (file.exists(model_path)) {
  if (!file.exists(training_cells_path)) {
    stop("scClassify model checkpoint exists without its training-cell ledger")
  }
  message("Loading existing fetal scClassify model checkpoint")
  scclassify_model <- readRDS(model_path)
} else {
  message("Training fetal HCA scClassify model with original author_mid_grain labels")
  set.seed(100)
  training_cells <- unlist(lapply(label_vocabulary, function(label) {
    candidates <- reference_cells[reference_labels == label]
    if (length(candidates) > 500L) sample(candidates, 500L) else candidates
  }), use.names = FALSE)
  training_ledger <- data.frame(
    cell = training_cells,
    fetal_hca_label = unname(reference_labels[training_cells]),
    stringsAsFactors = FALSE
  )
  if (anyNA(training_ledger$fetal_hca_label)) {
    stop("Failed to align fetal scClassify training cells and labels")
  }
  write.csv(training_ledger, training_cells_path, row.names = FALSE)
  
  scclassify_model <- train_scClassify(
    exprsMat_train = reference_data[, training_cells, drop = FALSE],
    cellTypes_train = setNames(
      training_ledger$fetal_hca_label, training_ledger$cell
    ),
    tree = "HOPACH",
    selectFeatures = "limma",
    topN = 50,
    hopach_kmax = 5,
    pSig = 0.05,
    cellType_tree = NULL,
    weightsCal = FALSE,
    parallel = workers > 1L,
    BPPARAM = bp,
    verbose = TRUE,
    returnList = FALSE
  )
  saveRDS(scclassify_model, model_path)
}

message("Reading the 22-dataset Step 1 query manifest")
query_manifest <- read.csv(manifest_path, stringsAsFactors = FALSE)
required_manifest_fields <- c("dataset", "cells", "unassigned", "query_rds")
missing_manifest_fields <- setdiff(required_manifest_fields, colnames(query_manifest))
if (length(missing_manifest_fields)) {
  stop("Query manifest missing: ", paste(missing_manifest_fields, collapse = ", "))
}
if (nrow(query_manifest) != 22L || anyDuplicated(query_manifest$dataset)) {
  stop("Expected 22 unique organoid datasets in the Step 0 query manifest")
}
if (!all(file.exists(query_manifest$query_rds))) {
  stop("One or more split query RDS files from Step 0 are missing")
}
if (sum(query_manifest$unassigned) != 56436L) {
  stop(
    "Expected 56,436 manuscript-unassigned cells; manifest contains ",
    format(sum(query_manifest$unassigned), big.mark = ",")
  )
}


empty_result <- function(dataset) {
  data.frame(
    dataset = character(), cell = character(), original_celltype = character(),
    fetal_hca_scClassify = character(), fetal_hca_seurat = character(),
    fetal_hca_singleR_pruned = character(),
    fetal_hca_consensus = character(), fetal_hca_vote_support = integer(),
    fetal_hca_recovered = logical(),
    stringsAsFactors = FALSE
  )
}

extract_terminal_scclassify <- function(predictions, query_cells, dataset) {
  hierarchy <- predictions$pearson_WKNN_limma$predLabelMat
  if (is.null(hierarchy) || nrow(hierarchy) != length(query_cells)) {
    stop(dataset, ": malformed scClassify prediction hierarchy")
  }
  if (is.null(rownames(hierarchy)) || !all(query_cells %in% rownames(hierarchy))) {
    stop(dataset, ": cannot align scClassify predictions to query cells")
  }
  hierarchy <- hierarchy[query_cells, , drop = FALSE]
  accepted <- c(label_vocabulary, "unassigned")
  valid_terminal_columns <- which(vapply(seq_len(ncol(hierarchy)), function(j) {
    values <- as.character(hierarchy[, j])
    !anyNA(values) && all(nzchar(values)) && all(values %in% accepted)
  }, logical(1)))
  if (!length(valid_terminal_columns)) {
    observed <- unique(unlist(lapply(seq_len(ncol(hierarchy)), function(j) {
      setdiff(unique(as.character(hierarchy[, j])), accepted)
    })))
    stop(
      dataset,
      ": no scClassify hierarchy column was demonstrably terminal. ",
      "Non-terminal/unknown values included: ", paste(observed, collapse = ", ")
    )
  }
  terminal_column <- max(valid_terminal_columns)
  labels <- as.character(hierarchy[, terminal_column])
  list(
    labels = labels,
    terminal_column = terminal_column,
    tree_depth = ncol(hierarchy)
  )
}

run_one_dataset <- function(object, dataset, dataset_index, expected_unassigned) {
  output_path <- file.path(prediction_dir, paste0(dataset, ".csv.gz"))
  if (file.exists(output_path)) {
    message("Skipping completed dataset: ", dataset)
    result <- read.csv(gzfile(output_path), check.names = FALSE,
                       stringsAsFactors = FALSE)
    if (nrow(result) != expected_unassigned) {
      stop(dataset, ": completed prediction row count does not match manifest")
    }
    return(result)
  }
  
  all_cells <- rownames(object@meta.data)
  supplied_labels <- as.character(object@meta.data$celltype)
  if (anyNA(supplied_labels)) stop(dataset, ": supplied celltype contains NA")
  query_cells <- all_cells[supplied_labels == "unassigned"]
  if (length(query_cells) != expected_unassigned) {
    stop(dataset, ": manuscript-unassigned count differs from Step 1 manifest")
  }
  if (!length(query_cells)) {
    result <- empty_result(dataset)
    write_csv_gz(result, output_path)
    return(result)
  }
  
  message(
    "Classifying ", format(length(query_cells), big.mark = ","),
    " adult-unassigned cells from ", dataset,
    " (", dataset_index, "/", nrow(query_manifest), ")"
  )
  query <- subset(object, cells = query_cells)
  query_data <- get_rna_layer(query, "data")
  
  common_features <- intersect(rownames(reference_data), rownames(query_data))
  if (length(common_features) < 1000L) {
    stop(dataset, ": fewer than 1,000 shared genes with fetal HCA")
  }
  
  sc_path <- file.path(checkpoint_dir, paste0(dataset, "_scClassify.rds"))
  if (file.exists(sc_path)) {
    sc_predictions <- readRDS(sc_path)
  } else {
    set.seed(1000L + dataset_index)
    sc_predictions <- predict_scClassify(
      query_data,
      scclassify_model,
      k = 10,
      prob_threshold = 0.7,
      cor_threshold_static = 0.5,
      cor_threshold_high = 0.7,
      features = "limma",
      algorithm = "WKNN",
      similarity = "pearson",
      parallel = workers > 1L,
      BPPARAM = bp,
      verbose = TRUE
    )
    saveRDS(sc_predictions, sc_path)
  }
  sc_result <- extract_terminal_scclassify(sc_predictions, query_cells, dataset)
  
  seurat_path <- file.path(checkpoint_dir, paste0(dataset, "_seurat_transfer.rds"))
  if (file.exists(seurat_path)) {
    transfer <- readRDS(seurat_path)
  } else {
    set.seed(2000L + dataset_index)
    anchors <- FindTransferAnchors(
      reference = reference,
      query = query,
      query.assay = "RNA",
      normalization.method = "LogNormalize",
      reference.reduction = "pca",
      dims = 1:30,
      features = common_features
    )
    transfer <- TransferData(
      anchorset = anchors,
      refdata = reference_labels,
      dims = 1:30
    )
    saveRDS(transfer, seurat_path)
    rm(anchors)
    gc()
  }
  
  singler_path <- file.path(checkpoint_dir, paste0(dataset, "_SingleR.rds"))
  if (file.exists(singler_path)) {
    single_r <- readRDS(singler_path)
  } else {
    set.seed(3000L + dataset_index)
    single_r <- SingleR(
      test = query_data[common_features, , drop = FALSE],
      ref = reference_data[common_features, , drop = FALSE],
      labels = reference_labels,
      de.method = "wilcox",
      fine.tune = TRUE,
      prune = TRUE,
      num.threads = min(4L, workers)
    )
    saveRDS(single_r, singler_path)
  }
  
  if (is.null(rownames(transfer)) || !all(query_cells %in% rownames(transfer))) {
    stop(dataset, ": cannot align Seurat transfer results to query cells")
  }
  transfer <- transfer[query_cells, , drop = FALSE]
  single_r <- as.data.frame(single_r)
  if (is.null(rownames(single_r)) || !all(query_cells %in% rownames(single_r))) {
    stop(dataset, ": cannot align SingleR results to query cells")
  }
  single_r <- single_r[query_cells, , drop = FALSE]
  
  vote_table <- data.frame(
    fetal_hca_scClassify = sc_result$labels,
    fetal_hca_seurat = as.character(transfer$predicted.id),
    fetal_hca_singleR_raw = as.character(single_r$labels),
    stringsAsFactors = FALSE,
    row.names = query_cells
  )

  consensus <- apply(vote_table, 1L, exact_consensus)
  support <- ifelse(
    consensus == "unassigned",
    0L,
    rowSums(vote_table == consensus, na.rm = TRUE)
  )
  
  
  result <- data.frame(
    dataset = dataset,
    cell = query_cells,
    original_celltype = "unassigned",
    fetal_hca_scClassify = vote_table$fetal_hca_scClassify,
    fetal_hca_seurat = vote_table$fetal_hca_seurat,
    fetal_hca_singleR_pruned = as.character(single_r$pruned.labels),
    fetal_hca_consensus = consensus,
    fetal_hca_vote_support = as.integer(support),
    fetal_hca_recovered = consensus != "unassigned",
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  
  assert_unique_keys(result, c("dataset", "cell"), dataset)
  write_csv_gz(result, output_path)
  rm(query, sc_predictions, transfer, single_r)
  gc()
  result
}

results <- vector("list", nrow(query_manifest))
names(results) <- query_manifest$dataset
for (i in seq_len(nrow(query_manifest))) {
  dataset <- query_manifest$dataset[[i]]
  object <- readRDS(query_manifest$query_rds[[i]])
  if (nrow(object@meta.data) != query_manifest$cells[[i]]) {
    stop(dataset, ": split query cell count differs from Step 0 manifest")
  }
  results[[i]] <- run_one_dataset(
    object = object,
    dataset = dataset,
    dataset_index = i,
    expected_unassigned = query_manifest$unassigned[[i]]
  )
  rm(object)
  gc()
}

all_predictions <- do.call(rbind, results)
rownames(all_predictions) <- NULL

assert_unique_keys(
  all_predictions, c("dataset", "cell"), "combined fetal-HCA ledger"
)
write_csv_gz(
  all_predictions,
  file.path(stage_dir, "fetal_hca_all_adult_unassigned_cells.csv.gz")
)

dataset_summary <- do.call(rbind, lapply(names(results), function(dataset) {
  x <- results[[dataset]]
  data.frame(
    dataset = dataset,
    adult_hca_unassigned = nrow(x),
    fetal_hca_recovered = sum(x$fetal_hca_recovered),
    doubly_unassigned = sum(!x$fetal_hca_recovered),
    fetal_hca_recovered_percent = if (nrow(x)) {
      100 * mean(x$fetal_hca_recovered)
    } else {
      NA_real_
    },
    stringsAsFactors = FALSE
  )
}))
write.csv(
  dataset_summary,
  file.path(stage_dir, "fetal_hca_dataset_summary.csv"),
  row.names = FALSE
)

consensus_counts <- as.data.frame(table(
  dataset = factor(all_predictions$dataset, levels = query_manifest$dataset),
  fetal_hca_consensus = all_predictions$fetal_hca_consensus
))
consensus_counts <- consensus_counts[consensus_counts$Freq > 0L, , drop = FALSE]
write.csv(
  consensus_counts,
  file.path(stage_dir, "fetal_hca_consensus_counts_by_dataset.csv"),
  row.names = FALSE
)


doubly_unassigned <- all_predictions[
  !all_predictions$fetal_hca_recovered,
  c("dataset", "cell"),
  drop = FALSE
]
assert_unique_keys(
  doubly_unassigned, c("dataset", "cell"), "doubly-unassigned key table"
)
write_csv_gz(
  doubly_unassigned,
  file.path(stage_dir, "doubly_unassigned_cell_keys.csv.gz")
)


# Join the fetal results back to complete organoid objects for whole-dataset
# visualization. Cells already assigned by the adult-HCA manuscript consensus
# remain in the objects but receive NA for fetal classifier outputs because
# they were deliberately not queried. The original celltype is never changed.
annotated_list_path <- file.path(
  stage_dir, "sc_organoid_list_with_fetal_hca_labels.rds"
)
joined_prediction_fields <- c(
  "fetal_hca_scClassify", "fetal_hca_seurat", "fetal_hca_singleR_pruned",
  "fetal_hca_consensus", "fetal_hca_vote_support",
  "fetal_hca_recovered"
)
joined_status_fields <- c(
  "fetal_hca_queried", "fetal_hca_reference_label_source",
  "offtarget_analysis_status", "sequential_hca_celltype",
  "sequential_hca_label_source"
)

validate_annotated_list <- function(objects) {
  if (!is.list(objects) || length(objects) != 22L ||
      !identical(names(objects), query_manifest$dataset)) {
    stop("Annotated organoid output is not the expected ordered 22-object list")
  }
  total_cells <- sum(vapply(objects, function(x) nrow(x@meta.data), integer(1)))
  total_queried <- sum(vapply(objects, function(x) {
    sum(x@meta.data$fetal_hca_queried)
  }, integer(1)))
  missing_fields <- unique(unlist(lapply(objects, function(x) {
    setdiff(
      c("celltype", joined_prediction_fields, joined_status_fields),
      colnames(x@meta.data)
    )
  })))
  if (total_cells != 211236L || total_queried != 56436L || length(missing_fields)) {
    stop(
      "Annotated organoid list failed validation: cells=", total_cells,
      ", queried=", total_queried,
      if (length(missing_fields)) {
        paste0(", missing metadata: ", paste(missing_fields, collapse = ", "))
      } else ""
    )
  }
  invisible(TRUE)
}

# The classifiers and summary tables are complete; free the large reference
# before reading all 22 full organoid objects into one list.
rm(reference, reference_data, scclassify_model)
gc()

if (file.exists(annotated_list_path)) {
  message("Validating existing annotated organoid-list output")
  annotated_organoids <- readRDS(annotated_list_path)
  validate_annotated_list(annotated_organoids)
  rm(annotated_organoids)
  gc()
} else {
  message("Joining fetal-HCA predictions to all cells in the 22 organoid objects")
  annotated_organoids <- vector("list", nrow(query_manifest))
  names(annotated_organoids) <- query_manifest$dataset
  for (i in seq_len(nrow(query_manifest))) {
    dataset <- query_manifest$dataset[[i]]
    object <- readRDS(query_manifest$query_rds[[i]])
    cells <- rownames(object@meta.data)
    original_manuscript_celltype <- as.character(object@meta.data$celltype)
    predictions <- results[[dataset]]
    if (anyDuplicated(cells) || anyDuplicated(predictions$cell) ||
        !all(predictions$cell %in% cells)) {
      stop(dataset, ": cannot uniquely join fetal predictions to the full object")
    }
    existing_conflicts <- intersect(
      c(joined_prediction_fields, joined_status_fields),
      colnames(object@meta.data)
    )
    if (length(existing_conflicts)) {
      stop(
        dataset, ": input object already contains Step 2 metadata fields: ",
        paste(existing_conflicts, collapse = ", ")
      )
    }
    
    match_index <- match(cells, predictions$cell)
    for (field in joined_prediction_fields) {
      object@meta.data[[field]] <- predictions[[field]][match_index]
    }
    queried <- original_manuscript_celltype == "unassigned"
    if (!identical(queried, !is.na(match_index))) {
      stop(dataset, ": joined query membership differs from manuscript celltype")
    }
    object@meta.data$fetal_hca_queried <- queried
    object@meta.data$fetal_hca_reference_label_source <- ifelse(
      queried, "author_mid_grain", NA_character_
    )
    
    status <- rep("adult_HCA_assigned_not_queried", length(cells))
    status[queried & object@meta.data$fetal_hca_recovered] <-
      "fetal_HCA_recovered"
    status[queried & !object@meta.data$fetal_hca_recovered] <-
      "doubly_unassigned"
    object@meta.data$offtarget_analysis_status <- factor(
      status,
      levels = c(
        "adult_HCA_assigned_not_queried", "fetal_HCA_recovered",
        "doubly_unassigned"
      )
    )
    
    sequential_label <- as.character(object@meta.data$celltype)
    recovered <- queried & object@meta.data$fetal_hca_recovered
    sequential_label[recovered] <- object@meta.data$fetal_hca_consensus[recovered]
    object@meta.data$sequential_hca_celltype <- sequential_label
    sequential_source <- rep("adult_HCA_manuscript_consensus", length(cells))
    sequential_source[recovered] <- "fetal_HCA_author_label_consensus"
    sequential_source[queried & !object@meta.data$fetal_hca_recovered] <-
      "unresolved_after_adult_and_fetal_HCA"
    object@meta.data$sequential_hca_label_source <- sequential_source
    
    
    annotated_organoids[[i]] <- object
    rm(object)
    gc()
  }
  validate_annotated_list(annotated_organoids)
  temporary_output <- tempfile(
    pattern = "annotated_organoid_list_", tmpdir = stage_dir, fileext = ".rds"
  )
  saveRDS(annotated_organoids, temporary_output, compress = "gzip")
  if (!file.rename(temporary_output, annotated_list_path)) {
    stop("Failed to move the completed annotated organoid list into place")
  }
  rm(annotated_organoids)
  gc()
}

