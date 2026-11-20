#!/usr/bin/env Rscript


# Run section 1 once; run sections 2-3 for SC and then repeat them for SN.

# Section 1. Setup ----
code_dir <- '.'
source(file.path(code_dir, 'AV_mapping_args.R'), local=TRUE)
# E.g. Rscript 04_score_organoids.R --modality SN.
config <- parse_pipeline_args(defaults=list(
  data_dir='../data',
  out_dir='./results',
  gene_file='../data/Cardiomyocyte_atrial_and_ventricular_module_gene_sets.xlsx',
  modality='SC'
))
print_pipeline_config(config)
source(file.path(code_dir, 'config_AV_mapping.R'), local=TRUE)
source(file.path(code_dir, 'AV_mapping_functions.R'), local=TRUE)
modules <- readRDS(file.path(cfg$cache_dir, 'module_gene_sets.rds'))

# Section 2. Choose a modality and load its original organoid list ----
# Interactive: set config$modality <- 'SN' before rerunning sections 2-3.
modality <- config$modality 
sources <- cfg$sources[[modality]]
organoids <- readRDS(cfg$organoid_files[[modality]])
datasets <- sort(names(organoids))
writeLines(datasets, file.path(cfg$cache_dir, paste0(modality, '_organoid_datasets.txt')))
score_dir <- file.path(cfg$cache_dir, 'organoid_scores', modality)
dir.create(score_dir, recursive=TRUE, showWarnings=FALSE)


# 3. Select CMs, normalize and score each organoid separately ----
for (dataset in datasets) {
  obj <- organoids[[dataset]]
  md <- obj[[]]
  keep <- as.character(md[[cfg$organoid_cm_column]]) %in% cfg$organoid_cm_labels
  md <- md[keep, , drop=FALSE]
  counts <- rna_counts(obj)[, rownames(md), drop=FALSE]
  # Release the full organoid object before normalization/scoring.
  organoids[dataset] <- list(NULL)
  rm(obj)
  gc()
  seu <- normalize_cm_counts(counts, md)
  rm(counts)
  gc()
  scores <- score_av_modules(seu, modules, sources)
  scores$dataset <- dataset
  scores$modality <- modality
  saveRDS(scores, file.path(score_dir, paste0(dataset, '_scores.rds')))
  cat('Scored:', dataset, '|', nrow(md), 'cardiomyocytes\n')
  rm(seu, md, scores)
  gc()
}
capture.output(sessionInfo(), file=file.path(cfg$log_dir, paste0(modality, '_organoid_sessionInfo.txt')))
rm(organoids)
gc()

