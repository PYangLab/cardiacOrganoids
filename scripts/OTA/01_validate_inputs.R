
suppressPackageStartupMessages({
  library(Seurat)
})


source(file.path("./pipeline_utils.R"))

config <- parse_pipeline_args(
  defaults = list(
    data_dir = "../data",
    out_dir = "./offtarget_results"
  )
)

data_dir <- path.expand(config$data_dir)
out_dir <- path.expand(config$out_dir)

message("Pipeline configuration:")
message("  Data directory: ", data_dir)
message("  Output directory: ", out_dir)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
query_dir <- file.path(out_dir, "00_queries")
dir.create(query_dir, recursive = TRUE, showWarnings = FALSE)

paths <- c(
  organoids = file.path(data_dir, "sc_organoid_list.rds"), #list of sc organoids
  adult_hca = file.path(data_dir, "hca_sc_processed.rds"), #HCA adult reference (processed)
  fetal_hca = file.path(data_dir, "fetal_hca_sc_processed.rds") #HCA fetal reference (processed)
)

missing_files <- paths[!file.exists(paths)]
if (length(missing_files)) {
  stop("Missing required input files: ", paste(missing_files, collapse = ", "))
}

file_manifest <- data.frame(
  input = names(paths),
  path = unname(paths),
  bytes = as.numeric(file.info(paths)$size),
  stringsAsFactors = FALSE
)
write.csv(file_manifest, file.path(out_dir, "input_file_manifest.csv"), row.names = FALSE)


message("Reading organoid list")
organoids <- readRDS(paths[["organoids"]])
if (!is.list(organoids) || length(organoids) != 22L) {
  stop("Expected a list of 22 organoid datasets; observed ", length(organoids))
}
if (is.null(names(organoids)) || anyDuplicated(names(organoids))) {
  stop("Organoid list must have 22 unique dataset names")
}

required_query_metadata <- c(
  "celltype"
)

query_rows <- lapply(names(organoids), function(dataset) {
  object <- organoids[[dataset]]
  cells <- get_object_cells(object)
  if (length(cells) != nrow(object@meta.data)) {
    stop(dataset, ": cell names and metadata rows differ")
  }
  if (anyDuplicated(cells)) {
    stop(dataset, ": duplicated cell names")
  }
  missing_metadata <- setdiff(required_query_metadata, colnames(object@meta.data))
  if (length(missing_metadata)) {
    stop(dataset, ": missing metadata: ", paste(missing_metadata, collapse = ", "))
  }
  data_matrix <- get_rna_layer(object, "data")

  labels <- as.character(object@meta.data[["celltype"]])
  query_path <- file.path(query_dir, paste0(dataset, ".rds"))
  if (!file.exists(query_path)) {
    # Splitting the multi-gigabyte list once lets later classification scripts load one query at a time.
    saveRDS(object, query_path, compress = FALSE)
  }
  data.frame(
    dataset = dataset,
    cells = length(cells),
    genes = nrow(data_matrix),
    assigned = sum(labels != "unassigned", na.rm = TRUE),
    unassigned = sum(labels == "unassigned", na.rm = TRUE),
    unassigned_percent = 100 * mean(labels == "unassigned", na.rm = TRUE),
    assay_class = class(object@assays$RNA)[1],
    has_pca = "pca" %in% names(object@reductions),
    has_umap = "umap" %in% names(object@reductions),
    query_rds = normalizePath(query_path),
    stringsAsFactors = FALSE
  )
})
query_summary <- do.call(rbind, query_rows)
write.csv(query_summary, file.path(out_dir, "organoid_input_summary.csv"), row.names = FALSE)

if (sum(query_summary$cells) != 211236L) {
  stop("Expected 211,236 total cells; observed ", sum(query_summary$cells))
}

rm(organoids, query_rows)
gc()

message("Reading adult HCA reference")
adult <- readRDS(paths[["adult_hca"]])
adult_required <- c("scReClassify_label", "age", "donor")
adult_missing <- setdiff(adult_required, colnames(adult@meta.data))
if (length(adult_missing)) {
  stop("Adult HCA metadata missing: ", paste(adult_missing, collapse = ", "))
}
adult_data <- get_rna_layer(adult, "data")
adult_labels <- standardize_adult_hca_labels(
  adult@meta.data[["scReClassify_label"]]
)
adult_summary <- data.frame(
  cells = ncol(adult_data),
  genes = nrow(adult_data),
  labels = length(unique(adult_labels)),
  ages = length(unique(adult@meta.data[["age"]])),
  donors = length(unique(adult@meta.data[["donor"]])),
  stringsAsFactors = FALSE
)
write.csv(adult_summary, file.path(out_dir, "adult_hca_input_summary.csv"), row.names = FALSE)
write.csv(
  as.data.frame(sort(table(adult_labels), decreasing = TRUE)),
  file.path(out_dir, "adult_hca_label_counts.csv"),
  row.names = FALSE
)
rm(adult, adult_data, adult_labels)
gc()

message("Reading fetal HCA reference")
fetal <- readRDS(paths[["fetal_hca"]])
fetal_required <- c(
  "coarse_grain", "mid_grain", "fine_grain", "age", "week", "donor",
  "cell_or_nuclei"
)
fetal_missing <- setdiff(fetal_required, colnames(fetal@meta.data))
if (length(fetal_missing)) {
  stop("Fetal HCA metadata missing: ", paste(fetal_missing, collapse = ", "))
}
if (!all(as.character(fetal@meta.data[["cell_or_nuclei"]]) == "cell")) {
  stop("Fetal reference contains entries not marked as single cells")
}
fetal_data <- get_rna_layer(fetal, "data")
fetal_summary <- data.frame(
  cells = ncol(fetal_data),
  genes = nrow(fetal_data),
  coarse_labels = length(unique(fetal@meta.data[["coarse_grain"]])),
  mid_labels = length(unique(fetal@meta.data[["mid_grain"]])),
  fine_labels = length(unique(fetal@meta.data[["fine_grain"]])),
  ages = length(unique(fetal@meta.data[["age"]])),
  donors = length(unique(fetal@meta.data[["donor"]])),
  stringsAsFactors = FALSE
)
write.csv(fetal_summary, file.path(out_dir, "fetal_hca_input_summary.csv"), row.names = FALSE)

write.csv(
  as.data.frame(sort(
    table(fetal@meta.data[["mid_grain"]]),
    decreasing = TRUE
  )),
  file.path(out_dir, "fetal_hca_mid_grain_counts.csv"),
  row.names = FALSE
)

