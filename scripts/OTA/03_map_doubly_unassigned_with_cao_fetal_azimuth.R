#!/usr/bin/env Rscript

# Map only cells that remain unassigned after both adult- and fetal-HCA
# consensus classification to the complete Cao human fetal-development
# Azimuth reference. 

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
})

if (!requireNamespace("Azimuth", quietly = TRUE)) {
  stop(
    "Package 'Azimuth' is required. Install the current official package in ",
    "a fresh R session with: remotes::install_github('satijalab/azimuth', ",
    "ref = 'master')"
  )
}
if (packageVersion("Azimuth") < "0.5.1") {
  stop(
    "Azimuth >= 0.5.1 is required; installed version is ",
    as.character(packageVersion("Azimuth"))
  )
}
if (packageVersion("Seurat") < "5.0.0") {
  stop("Seurat >= 5.0 is required")
}

source(file.path("./pipeline_utils.R"))

config <- parse_pipeline_args(
  defaults = list(
    out_dir = "./offtarget_results",
    cao_reference_dir = "../data/cao_fetal_azimuth"
  )
)

out_dir <- path.expand(config$out_dir)
reference_dir <- path.expand(config$cao_reference_dir)

message("Pipeline configuration:")
message("  Output directory: ", out_dir)
message("  Cao fetal Azimuth reference: ", reference_dir)

join_full_list <- TRUE
set.seed(20092026)
options(future.globals.maxSize = 16 * 1024^3)
if (requireNamespace("future", quietly = TRUE)) {
  future::plan("sequential")
}

manifest_path <- file.path(out_dir, "organoid_input_summary.csv")
step2_dir <- file.path(out_dir, "02_fetal_hca_classification")

key_path <- file.path(step2_dir, "doubly_unassigned_cell_keys.csv.gz")
step2_summary_path <- file.path(step2_dir, "fetal_hca_dataset_summary.csv")
annotated_step2_path <- file.path(
  step2_dir, "sc_organoid_list_with_fetal_hca_labels.rds"
)
required_inputs <- c(
  manifest_path, key_path, step2_summary_path
)
if (!all(file.exists(required_inputs))) {
  stop(
    "Missing: ",
    paste(required_inputs[!file.exists(required_inputs)], collapse = ", ")
  )
}
if (join_full_list && !file.exists(annotated_step2_path)) {
  stop(
    "join_full_list=TRUE requires the Step 2 annotated organoid list: ",
    annotated_step2_path
  )
}

reference_files <- file.path(reference_dir, c("ref.Rds", "idx.annoy"))
if (!all(file.exists(reference_files))) {
  stop(
    "The Cao reference directory must contain ref.Rds and idx.annoy. Missing: ",
    paste(basename(reference_files[!file.exists(reference_files)]),
          collapse = ", ")
  )
}

# Published MD5 values for the official Human Fetus Azimuth v1.0.0 reference,
# Zenodo record 4738021. This prevents accidental use of a different reference.
expected_reference_md5 <- c(
  `ref.Rds` = "3a7f09a2bbff76f940c6c33209875665",
  `idx.annoy` = "966421e3e4dd28b96f0a75ba96569dd9"
)
observed_reference_md5 <- unname(tools::md5sum(reference_files))
names(observed_reference_md5) <- basename(reference_files)
if (!identical(observed_reference_md5[names(expected_reference_md5)],
               expected_reference_md5)) {
  stop(
    "Cao fetal reference checksum mismatch. Observed: ",
    paste(
      names(observed_reference_md5), observed_reference_md5,
      sep = "=", collapse = "; "
    )
  )
}

stage_dir <- file.path(out_dir, "03_cao_fetal_azimuth")
prediction_dir <- file.path(stage_dir, "predictions")
score_dir <- file.path(stage_dir, "score_matrices")
embedding_dir <- file.path(stage_dir, "embeddings")
dir.create(prediction_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(score_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(embedding_dir, recursive = TRUE, showWarnings = FALSE)

manifest <- read.csv(manifest_path, stringsAsFactors = FALSE)
required_manifest_fields <- c("dataset", "cells", "query_rds")
missing_manifest_fields <- setdiff(required_manifest_fields, colnames(manifest))
if (length(missing_manifest_fields)) {
  stop("Input manifest is missing: ", paste(missing_manifest_fields, collapse = ", "))
}
if (nrow(manifest) != 22L || anyDuplicated(manifest$dataset)) {
  stop("Expected 22 unique organoid datasets in the input manifest")
}
if (!all(file.exists(manifest$query_rds))) {
  stop("One or more split organoid query RDS files are unavailable")
}

doubly_unassigned <- read.csv(
  gzfile(key_path), stringsAsFactors = FALSE, check.names = FALSE
)
if (!all(c("dataset", "cell") %in% colnames(doubly_unassigned))) {
  stop("Doubly-unassigned key table must contain dataset and cell")
}
assert_unique_keys(
  doubly_unassigned, c("dataset", "cell"), "Step 2 doubly-unassigned keys"
)


step2_summary <- read.csv(step2_summary_path, stringsAsFactors = FALSE)
if (!all(c("dataset", "doubly_unassigned") %in% colnames(step2_summary))) {
  stop("Step 2 summary lacks dataset or doubly_unassigned")
}
expected_by_dataset <- setNames(
  step2_summary$doubly_unassigned, step2_summary$dataset
)
observed_by_dataset <- table(factor(
  doubly_unassigned$dataset, levels = manifest$dataset
))
if (!identical(
  as.integer(observed_by_dataset),
  as.integer(expected_by_dataset[manifest$dataset])
)) {
  stop("Doubly-unassigned key counts do not match the Step 2 summary")
}

message("Checking the official Cao fetal Azimuth reference")
reference <- readRDS(reference_files[[1L]])
required_reference_levels <- c("annotation.l1", "annotation.l2", "organ")
missing_reference_levels <- setdiff(
  required_reference_levels, colnames(reference@meta.data)
)
if (length(missing_reference_levels)) {
  stop(
    "Cao reference metadata is missing: ",
    paste(missing_reference_levels, collapse = ", ")
  )
}

reference_hierarchy <- unique(data.frame(
  annotation_l1 = as.character(reference@meta.data$annotation.l1),
  annotation_l2 = as.character(reference@meta.data$annotation.l2),
  organ = as.character(reference@meta.data$organ),
  stringsAsFactors = FALSE
))
if (anyNA(reference_hierarchy) || any(!nzchar(unlist(reference_hierarchy)))) {
  stop("Cao reference hierarchy contains missing or empty labels")
}
write.csv(
  reference_hierarchy,
  file.path(stage_dir, "cao_fetal_reference_label_hierarchy.csv"),
  row.names = FALSE
)
reference_hierarchy_key <- do.call(
  paste, c(reference_hierarchy, sep = "\r")
)
reference_features <- rownames(reference)
reference_summary <- data.frame(
  reference = "Azimuth Human Fetus",
  source_study = "Cao et al., Science 2020",
  zenodo_record = "4738021",
  reference_cells = ncol(reference),
  reference_features = nrow(reference),
  annotation_l1_labels = length(unique(reference_hierarchy$annotation_l1)),
  annotation_l2_labels = length(unique(reference_hierarchy$annotation_l2)),
  organ_labels = length(unique(reference_hierarchy$organ)),
  ref_Rds_md5 = observed_reference_md5[["ref.Rds"]],
  idx_annoy_md5 = observed_reference_md5[["idx.annoy"]],
  stringsAsFactors = FALSE
)
write.csv(
  reference_summary,
  file.path(stage_dir, "cao_fetal_reference_audit.csv"),
  row.names = FALSE
)
rm(reference)
gc()

run_contract <- list(
  reference_directory = reference_dir,
  reference_md5 = observed_reference_md5,
  query_population = paste(
    "Only cells with offtarget_analysis_status == 'doubly_unassigned'",
    "after adult- and fetal-HCA consensus"
  ),
  annotation_levels = required_reference_levels,
  input_layer = "original integer RNA counts",
  dataset_handling = "each of 22 datasets mapped independently",
  seed = 20092026L,
  azimuth_version = as.character(packageVersion("Azimuth")),
  seurat_version = as.character(packageVersion("Seurat")),
  interpretation = paste(
    "Raw reference-mapping only"
  )
)
contract_path <- file.path(stage_dir, "mapping_contract.rds")
if (file.exists(contract_path)) {
  if (!identical(readRDS(contract_path), run_contract)) {
    stop(
      "Existing Step 3 checkpoints were generated under a different mapping ",
      "contract. Use a new results directory or inspect the stale Step 3 run."
    )
  }
} else {
  saveRDS(run_contract, contract_path)
}

get_assay_matrix <- function(object, assay) {
  value <- tryCatch(
    SeuratObject::LayerData(object = object, assay = assay, layer = "data"),
    error = function(e) NULL
  )
  if (is.null(value)) {
    value <- Seurat::GetAssayData(object = object, assay = assay, slot = "data")
  }
  value
}

validate_prediction_result <- function(result, dataset, expected_cells) {
  required <- c(
    "dataset", "cell",
    "cao_fetal_annotation_l1", "cao_fetal_annotation_l1_score",
    "cao_fetal_annotation_l2", "cao_fetal_annotation_l2_score",
    "cao_fetal_organ", "cao_fetal_organ_score",
    "cao_fetal_mapping_score", "cao_fetal_hierarchy_consistent"
  )
  missing <- setdiff(required, colnames(result))
  if (length(missing)) {
    stop(dataset, ": prediction checkpoint is missing: ",
         paste(missing, collapse = ", "))
  }
  if (nrow(result) != length(expected_cells) ||
      anyDuplicated(result$cell) ||
      !setequal(result$cell, expected_cells) ||
      any(result$dataset != dataset)) {
    stop(dataset, ": prediction checkpoint does not match expected query cells")
  }
  if (anyNA(result[, c(
    "cao_fetal_annotation_l1", "cao_fetal_annotation_l1_score",
    "cao_fetal_annotation_l2", "cao_fetal_annotation_l2_score",
    "cao_fetal_organ", "cao_fetal_organ_score",
    "cao_fetal_hierarchy_consistent"
  )])) {
    stop(dataset, ": core Cao predictions contain missing values")
  }
  result
}

azimuth_k_weight <- 20L
azimuth_mapping_score_k <- 20L

run_azimuth_with_consistent_k <- function(query) {
  embedding_method <- "IntegrateEmbeddings.TransferAnchorSet"
  mapping_method <- "MappingScore.AnchorSet"
  method_namespace <- asNamespace("Seurat")
  
  required_methods <- c(embedding_method, mapping_method)
  missing_methods <- required_methods[!vapply(
    required_methods,
    exists,
    logical(1),
    envir = method_namespace,
    inherits = FALSE
  )]
  
  if (length(missing_methods)) {
    stop(
      "Cannot locate Seurat method(s): ",
      paste(missing_methods, collapse = ", ")
    )
  }
  
  trace(
    what = embedding_method,
    tracer = substitute({
      k.weight <- K_WEIGHT
      message(
        "Applying k.weight = ", k.weight,
        " to Azimuth's IntegrateEmbeddings step"
      )
    }, list(K_WEIGHT = azimuth_k_weight)),
    print = FALSE,
    where = method_namespace
  )
  
  on.exit(
    untrace(
      what = embedding_method,
      where = method_namespace
    ),
    add = TRUE
  )
  
  trace(
    what = mapping_method,
    tracer = substitute({
      kanchors <- K_WEIGHT
      ksmooth <- K_SMOOTH
      
      message(
        "Applying kanchors = ", kanchors,
        " and ksmooth = ", ksmooth,
        " to Azimuth's MappingScore step"
      )
    }, list(
      K_WEIGHT = azimuth_k_weight,
      K_SMOOTH = azimuth_mapping_score_k
    )),
    print = FALSE,
    where = method_namespace
  )
  
  on.exit(
    untrace(
      what = mapping_method,
      where = method_namespace
    ),
    add = TRUE
  )
  
  Azimuth::RunAzimuth(
    query = query,
    reference = reference_dir,
    annotation.levels = required_reference_levels,
    k.weight = azimuth_k_weight,
    mapping.score.k = azimuth_mapping_score_k,
    verbose = TRUE
  )
}

run_one_dataset <- function(dataset, query_path, expected_cells,
                            dataset_index, total_datasets) {
  prediction_path <- file.path(prediction_dir, paste0(dataset, ".csv.gz"))
  score_path <- file.path(score_dir, paste0(dataset, ".rds"))
  embedding_path <- file.path(embedding_dir, paste0(dataset, ".rds"))
  umap_path <- file.path(embedding_dir, paste0(dataset, "_ref_umap.csv.gz"))
  feature_path <- file.path(prediction_dir, paste0(dataset, "_feature_audit.csv"))
  error_path <- file.path(prediction_dir, paste0(dataset, ".error.txt"))
  
  if (file.exists(prediction_path) && file.exists(score_path) &&
      file.exists(embedding_path)) {
    message("Validating completed Cao mapping: ", dataset)
    result <- read.csv(
      gzfile(prediction_path), stringsAsFactors = FALSE, check.names = FALSE
    )
    return(validate_prediction_result(result, dataset, expected_cells))
  }
  
  message(
    "Mapping ", format(length(expected_cells), big.mark = ","),
    " doubly-unassigned cells from ", dataset,
    " (", dataset_index, "/", total_datasets, ")"
  )
  object <- readRDS(query_path)
  object_cells <- get_object_cells(object)
  if (anyDuplicated(object_cells) || !all(expected_cells %in% object_cells)) {
    stop(dataset, ": doubly-unassigned cells cannot be aligned to query object")
  }
  counts <- get_rna_layer(object, "counts")[, expected_cells, drop = FALSE]
  if (!identical(colnames(counts), expected_cells)) {
    stop(dataset, ": count columns changed order during subsetting")
  }
  nonzero_counts <- if (inherits(counts, "sparseMatrix")) {
    counts@x
  } else {
    as.vector(counts[counts != 0])
  }
  if (length(nonzero_counts) &&
      any(abs(nonzero_counts - round(nonzero_counts)) > 1e-6)) {
    stop(dataset, ": Azimuth input counts are not integer-valued")
  }
  write.csv(
    data.frame(
      dataset = dataset,
      query_cells = ncol(counts),
      query_features = nrow(counts),
      reference_features = length(reference_features),
      overlapping_features = sum(rownames(counts) %in% reference_features),
      stringsAsFactors = FALSE
    ),
    feature_path, row.names = FALSE
  )
  
  query <- Seurat::CreateSeuratObject(
    counts = counts, project = dataset, min.cells = 0, min.features = 0
  )
  if (!identical(colnames(query), expected_cells)) {
    stop(dataset, ": CreateSeuratObject changed query cell names or order")
  }
  rm(object, counts, nonzero_counts)
  gc()
  
  mapped <- tryCatch(
    run_azimuth_with_consistent_k(query),
    error = function(e) {
      writeLines(
        c(as.character(Sys.time()), conditionMessage(e)),
        error_path
      )
      stop(
        dataset,
        ": Cao fetal Azimuth mapping failed: ",
        conditionMessage(e)
      )
    }
  )
  if (!identical(colnames(mapped), expected_cells)) {
    stop(dataset, ": Azimuth output cells differ from the submitted query")
  }
  
  native_fields <- c(
    l1 = "predicted.annotation.l1",
    l1_score = "predicted.annotation.l1.score",
    l2 = "predicted.annotation.l2",
    l2_score = "predicted.annotation.l2.score",
    organ = "predicted.organ",
    organ_score = "predicted.organ.score"
  )
  missing_native <- setdiff(unname(native_fields), colnames(mapped@meta.data))
  if (length(missing_native)) {
    stop(dataset, ": Azimuth output is missing: ",
         paste(missing_native, collapse = ", "))
  }
  md <- mapped@meta.data[expected_cells, , drop = FALSE]
  hierarchy_key <- paste(
    as.character(md[[native_fields[["l1"]]]]),
    as.character(md[[native_fields[["l2"]]]]),
    as.character(md[[native_fields[["organ"]]]]),
    sep = "\r"
  )
  mapping_score <- if ("mapping.score" %in% colnames(md)) {
    as.numeric(md$mapping.score)
  } else {
    rep(NA_real_, nrow(md))
  }
  result <- data.frame(
    dataset = dataset,
    cell = expected_cells,
    cao_fetal_annotation_l1 = as.character(md[[native_fields[["l1"]]]]),
    cao_fetal_annotation_l1_score = as.numeric(
      md[[native_fields[["l1_score"]]]]
    ),
    cao_fetal_annotation_l2 = as.character(md[[native_fields[["l2"]]]]),
    cao_fetal_annotation_l2_score = as.numeric(
      md[[native_fields[["l2_score"]]]]
    ),
    cao_fetal_organ = as.character(md[[native_fields[["organ"]]]]),
    cao_fetal_organ_score = as.numeric(md[[native_fields[["organ_score"]]]]),
    cao_fetal_mapping_score = mapping_score,
    cao_fetal_hierarchy_consistent = hierarchy_key %in%
      reference_hierarchy_key,
    stringsAsFactors = FALSE
  )
  result <- validate_prediction_result(result, dataset, expected_cells)
  
  score_assays <- grep(
    "^prediction.score", names(mapped@assays), value = TRUE
  )
  expected_score_assays <- paste0(
    "prediction.score.", required_reference_levels
  )
  if (!all(expected_score_assays %in% score_assays)) {
    stop(dataset, ": not all per-class prediction-score assays were retained")
  }
  score_matrices <- lapply(expected_score_assays, function(assay) {
    value <- get_assay_matrix(mapped, assay)
    if (!identical(colnames(value), expected_cells)) {
      stop(dataset, ": score matrix cells differ for assay ", assay)
    }
    value
  })
  names(score_matrices) <- expected_score_assays
  

  embeddings <- lapply(mapped@reductions, Seurat::Embeddings)
  if (!length(embeddings)) {
    stop(dataset, ": Azimuth output did not retain projected embeddings")
  }
  # Verify that every embedding contains exactly the submitted cells,
  # then restore the original query-cell order.
  for (reduction_name in names(embeddings)) {
    reduction_cells <- rownames(embeddings[[reduction_name]])
    
    if (
      anyDuplicated(reduction_cells) ||
      !setequal(reduction_cells, expected_cells)
    ) {
      stop(
        dataset, ": ", reduction_name,
        " embedding cells differ from the submitted query"
      )
    }
    
    embeddings[[reduction_name]] <- embeddings[[reduction_name]][
      expected_cells,
      ,
      drop = FALSE
    ]
  }
  
  if ("ref.umap" %in% names(embeddings)) {
    ref_umap <- embeddings[["ref.umap"]]
    write_csv_gz(
      data.frame(
        dataset = dataset, cell = rownames(ref_umap), ref_umap,
        check.names = FALSE, stringsAsFactors = FALSE
      ),
      umap_path
    )
  }
  
  temporary_prediction <- paste0(prediction_path, ".tmp")
  temporary_scores <- paste0(score_path, ".tmp")
  temporary_embeddings <- paste0(embedding_path, ".tmp")
  write_csv_gz(result, temporary_prediction)
  saveRDS(score_matrices, temporary_scores, compress = "gzip")
  saveRDS(embeddings, temporary_embeddings, compress = "gzip")
  if (!file.rename(temporary_prediction, prediction_path) ||
      !file.rename(temporary_scores, score_path) ||
      !file.rename(temporary_embeddings, embedding_path)) {
    stop(dataset, ": failed to move completed checkpoint files into place")
  }
  if (file.exists(error_path)) unlink(error_path)
  rm(query, mapped, md, score_matrices, embeddings)
  gc()
  result
}

mapping_order <- manifest$dataset[order(
  as.integer(observed_by_dataset[manifest$dataset]),
  manifest$dataset
)]
results <- vector("list", length(mapping_order))
names(results) <- mapping_order
for (i in seq_along(mapping_order)) {
  dataset <- mapping_order[[i]]
  expected_cells <- doubly_unassigned$cell[
    doubly_unassigned$dataset == dataset
  ]
  if (!length(expected_cells)) {
    stop(dataset, ": no doubly-unassigned cells were supplied")
  }
  query_path <- manifest$query_rds[match(dataset, manifest$dataset)]
  results[[i]] <- run_one_dataset(
    dataset = dataset,
    query_path = query_path,
    expected_cells = expected_cells,
    dataset_index = i,
    total_datasets = length(mapping_order)
  )
}

all_predictions <- do.call(rbind, results)
rownames(all_predictions) <- NULL
all_predictions$dataset <- factor(
  all_predictions$dataset, levels = manifest$dataset
)
all_predictions <- all_predictions[
  order(all_predictions$dataset, all_predictions$cell), , drop = FALSE
]
all_predictions$dataset <- as.character(all_predictions$dataset)


expected_key <- paste(doubly_unassigned$dataset, doubly_unassigned$cell, sep = "\r")
observed_key <- paste(all_predictions$dataset, all_predictions$cell, sep = "\r")

write_csv_gz(
  all_predictions,
  file.path(stage_dir, "cao_fetal_all_doubly_unassigned_cells.csv.gz")
)

dataset_summary <- do.call(rbind, lapply(manifest$dataset, function(dataset) {
  x <- all_predictions[all_predictions$dataset == dataset, , drop = FALSE]
  data.frame(
    dataset = dataset,
    doubly_unassigned_cells = nrow(x),
    hierarchy_consistent_cells = sum(x$cao_fetal_hierarchy_consistent),
    hierarchy_consistent_percent = 100 * mean(x$cao_fetal_hierarchy_consistent),
    median_l1_score = median(x$cao_fetal_annotation_l1_score),
    median_l2_score = median(x$cao_fetal_annotation_l2_score),
    median_organ_score = median(x$cao_fetal_organ_score),
    cells_with_mapping_score = sum(is.finite(x$cao_fetal_mapping_score)),
    stringsAsFactors = FALSE
  )
}))
write.csv(
  dataset_summary,
  file.path(stage_dir, "cao_fetal_dataset_summary.csv"),
  row.names = FALSE
)

for (field in c(
  "cao_fetal_annotation_l1", "cao_fetal_annotation_l2", "cao_fetal_organ"
)) {
  counts <- as.data.frame(table(
    dataset = factor(all_predictions$dataset, levels = manifest$dataset),
    label = all_predictions[[field]]
  ))
  counts <- counts[counts$Freq > 0L, , drop = FALSE]
  write.csv(
    counts,
    file.path(stage_dir, paste0(field, "_counts_by_dataset.csv")),
    row.names = FALSE
  )
}

joined_output_path <- file.path(
  stage_dir, "sc_organoid_list_with_fetal_hca_and_cao_labels.rds"
)
joined_fields <- setdiff(
  colnames(all_predictions), c("dataset", "cell")
)
validate_joined_list <- function(objects) {
  if (length(objects) != 22L || !identical(names(objects), manifest$dataset)) {
    stop("Joined organoid list does not contain the expected 22 datasets")
  }
  for (dataset in manifest$dataset) {
    object <- objects[[dataset]]

    if (!all(c("celltype", "cao_fetal_queried", joined_fields) %in%
             colnames(object@meta.data))) {
      stop(dataset, ": joined object lacks Cao metadata fields")
    }
    expected <- doubly_unassigned$cell[doubly_unassigned$dataset == dataset]
    observed <- rownames(object@meta.data)[object@meta.data$cao_fetal_queried]
    if (!setequal(expected, observed)) {
      stop(dataset, ": joined Cao query membership is incorrect")
    }
  }
  invisible(TRUE)
}

if (join_full_list) {
  if (file.exists(joined_output_path)) {
    message("Validating existing Cao-annotated full organoid list")
    joined_organoids <- readRDS(joined_output_path)
    validate_joined_list(joined_organoids)
    rm(joined_organoids)
    gc()
  } else {
    message(
      "Joining Cao predictions to all cells in the 22 Step 2 annotated objects"
    )
    joined_organoids <- readRDS(annotated_step2_path)
    if (length(joined_organoids) != 22L ||
        !identical(names(joined_organoids), manifest$dataset)) {
      stop("Step 2 annotated list does not match the 22-dataset manifest")
    }
    for (dataset in manifest$dataset) {
      object <- joined_organoids[[dataset]]
      cells <- rownames(object@meta.data)
      predictions <- all_predictions[
        all_predictions$dataset == dataset, , drop = FALSE
      ]
      if (anyDuplicated(cells) || anyDuplicated(predictions$cell) ||
          !all(predictions$cell %in% cells)) {
        stop(dataset, ": Cao predictions cannot be uniquely joined")
      }
      conflicts <- intersect(
        c("cao_fetal_queried", joined_fields),
        colnames(object@meta.data)
      )
      if (length(conflicts)) {
        stop(
          dataset, ": input object already contains Step 3 fields: ",
          paste(conflicts, collapse = ", ")
        )
      }
      match_index <- match(cells, predictions$cell)
      object@meta.data$cao_fetal_queried <- !is.na(match_index)
      for (field in joined_fields) {
        object@meta.data[[field]] <- predictions[[field]][match_index]
      }
      if (!identical(
        object@meta.data$cao_fetal_queried,
        as.character(object@meta.data$offtarget_analysis_status) ==
        "doubly_unassigned"
      )) {
        stop(dataset, ": Cao query membership differs from Step 2 status")
      }
      joined_organoids[[dataset]] <- object
      rm(object)
      gc()
    }
    validate_joined_list(joined_organoids)
    temporary_output <- tempfile(
      pattern = "cao_annotated_organoid_list_",
      tmpdir = stage_dir, fileext = ".rds"
    )
    saveRDS(joined_organoids, temporary_output, compress = "gzip")
    if (!file.rename(temporary_output, joined_output_path)) {
      stop("Failed to move completed Cao-annotated organoid list into place")
    }
    rm(joined_organoids)
    gc()
  }
}


